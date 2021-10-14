function [  ] = DTI_cluster_perm_2groups_LBPD( S )

% It calculates clusters (neighbouring voxels) of significant voxels after
% statistics computed comparing the DTI of two groups of subjects.
% It permutes the significant voxels and compute the clusters also for the
% permuted data to get the significant clusters within the original data.
% Parallel computing is done for original clusters and for blocks of
% permutations.

%  INPUT:   -S.V_p:             data in 3D (p-values), obtained by contrasting the DTI of 2 groups of subjects.
%                               if empty [], it tries to load it from S.outdir where it should have been automatically saved it the first time.
%           -S.V_t:             same but for t-values
%                               if empty [], it tries to load it from S.outdir where it should have been automatically saved it the first time.
%           -S.p_thresh:        threshold for binarising the p-values (e.g. 0.05)
%           -S.thresh_mc:       threshold for Monte-Carlo simulations (e.g. 0.001)
%           -S.perm_numb:       number of permutations
%           -S.outdir:          output directory
%           -S.analysis_name:   name for the analysis
%           -S.maskFA:          path to mask FA created by FSL TBSS (or equivalent)
%           -S.sizemass:        specifying if your want to do analysis considering cluster size (number of voxels forming the cluster)
%                               or mass (sum of statistics of the voxels forming the cluster)
%                                   'size'; having cluster size
%                                   'mass'; having cluster mass
%           -S.clsopt:          cluster specifications:
%                                   i.   .scheduler: 'cluster' (parallel computing)
%                                                    'none' (local)
%                                   ii.  .slot: number of slots (1 slot = 1 GB of memory)
%                                   iii. .queue: 1 = all; 0 = short; 2 = long
%           -S.fun_run:         1x4 vector to indicate which steps you want to run (e.g. [1 1 0 0])
%                               The steps are:
%                                   original clusters (i) (parallel computing)
%                                   permuted clusters (ii) (parallel computing)
%                                   combining original and permuted clusters (iii) (local)
%                                   plotting (iv)


%  OUTPUT:  -saved statistics:  'significant clusters.mat' is a 3x2 cell matrix file with the significant clusters stored in (considering both cluster size and mass)
%                               Clusters in the matrix: 1,1: cond1>cond2 size; 2,1: cond1>cond2 mass; 1,2: cond2>cond1 size; 2,2: cond2>cond1 mass; 3x2: legend





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 20/12/2020


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






if S.fun_run(1) == 1 || S.fun_run(2) == 1
    %setting cluster (parallel computing)
    addpath('/projects/MINDLAB2020_MEG-AuditoryPatternRecognition/scripts/leonardo')
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/DTI_Portis')
    clusterconfig('slot', S.clsopt.slot); % set manually the job cluster slots
    clusterconfig('scheduler', S.clsopt.scheduler);
    clusterconfig('long_running', S.clsopt.queue); % 0 = short; 1 = all; 2 = long
    clusterconfig('wait',0)
end


%getting some inputs
if ~isempty(S.V_p) %checking if user provided PVAL (p-values)
    V_p = S.V_p;
    ap = dir([S.outdir '/PVAL.mat*']);
    if isempty(ap) %checking if PVAL is already in the outdir folder
        disp('saving p-values..')
        if ~exist([S.outdir],'dir')
            mkdir([S.outdir])
        end
        save([S.outdir '/PVAL.mat'],'V_p'); %if not, saving it for later usages..
    end
else %otherwise trying to load it
    ap = dir([S.outdir '/PVAL.mat*']);
    if ~isempty(ap)
        load([ap(1).folder '/' ap(1).name]);
        V_p = V_p;
    else
        disp('p-values are not provided nor in the outdir folder..')
    end
end
if ~isempty(S.V_t) %same concept but for t-values
    V_t = S.V_t;
    at = dir([S.outdir '/TVAL.mat*']);
    if isempty(at)
        disp('saving t-values..')
        if ~exist([S.outdir],'dir')
            mkdir([S.outdir])
        end
        save([S.outdir '/TVAL.mat'],'V_t');
    end
else
    at = dir([S.outdir '/TVAL.mat*']);
    if ~isempty(at)
        load([at(1).folder '/' at(1).name]);
        V_t = V_t;
    else
        disp('t-values are not provided nor in the outdir folder..')
    end
end
%getting size of the 3D data
SS = size(V_p); 

%cond 1 > cond 2
%binarizing matrix according to p-values and t-values
V_bin_1 = zeros(SS(1),SS(2),SS(3));
V_p(V_p == 0) = NaN;
V_bin_1(V_p < S.p_thresh) = 1; %assigning 1 to significant values
V_bin_1(V_t < 0) = 0; %deleting p-values for negative contrasts (where cond2 > cond1)
Vbinp_1 = V_bin_1;
V_bin_1(V_bin_1==0) = NaN; %assigning 'NaN' to 0 values (that corresponds to the non-significant voxels)

%cond 2 > cond 1
%binarizing matrix according to p-values and t-values
V_bin_2 = zeros(SS(1),SS(2),SS(3));
V_bin_2(V_p < S.p_thresh) = 1; %assigning 1 to significant values
V_bin_2(V_t > 0) = 0; %deleting p-values for positive contrasts (where cond1 > cond2)
VTOT = V_bin_2 + Vbinp_1; %summing significant voxels for cond1 > cond2 and for cond2 > cond1 for later permutation test
VTOT(VTOT==0) = NaN; %same for VTOT
V_bin_2(V_bin_2==0) = NaN; %assigning 'NaN' to 0 values (that corresponds to the non-significant voxels)

%preparing structure with inputs for clustering function (original clusters)
S_ori = [];
S_ori.V_bin_1 = V_bin_1;
S_ori.V_bin_2 = V_bin_2;
S_ori.outdir = S.outdir;
S_ori.analysis_name = S.analysis_name;
S_ori.V_t = S.V_t;

    
if S.fun_run(1) == 1
    jobid = job2cluster(@cluster_DTI_perm_2groups_1,S_ori); %running with parallel computing
end


%permutations
if S.fun_run(2) == 1
    %preparing structure with inputs for clustering function (permuted clusters)
    S_perm = [];
    S_perm.outdir = S.outdir;
    S_perm.analysis_name = S.analysis_name;
    S_perm.VTOT = VTOT; %significant voxels (binary) for permutations
    S_perm.V_t = abs(S.V_t); %here we have absolute values since we want to take into account both directions of the contrast (both cond1>cond2 and cond2>cond1 which has originally negative t-values) when doing the permutations
    if S.perm_numb < 100 %if you ask less than 100 permutations (e.g. n permutations with n < 100), I split the permutations in n jobs)
        S_perm.perm = 1;
        ff = S.perm_numb;
    else %otherwise I split the permutations up to 100 different jobs (this could be changed in the future, but I guess it is more convenient this way for the current settings of our parallel computing system)
        S_perm.perm = round(S.perm_numb/100);
        ff = 100;
    end
    %submitting jobs
    for kk = 1:ff
        S_perm.kk = kk;
%         if kk == ff
%             clusterconfig('wait',1)
%         end
        jobid = job2cluster(@cluster_DTI_perm_2groups_2,S_perm); %running with parallel computing
    end
end

disp(['check your results in: ' S.outdir])

if S.fun_run(3) == 1
    %combining original clusters and permuted clusters
    S_fin.outdir = S.outdir;
    S_fin.analysis_name = S.analysis_name;
    S_fin.thresh_mc = S.thresh_mc;
    S_fin.perm = S.perm_numb;
    clust = cluster_DTI_perm_2groups_3(S_fin);
    save([S.outdir '/' S.analysis_name '_significant_clusters.mat'],'clust')
end


%plotting results
if S.fun_run(4) == 1
    disp('reshaping cluster information')
    if S.fun_run(3) == 0
        load([S.outdir '/' S.analysis_name '_significant_clusters.mat'])
    end
    if strcmp(S.sizemass,'size')
        cl1 = clust{1,1}; %cond1>cond2
        cl2 = clust{1,2}; %cond2>cond1
    elseif strcmp(S.sizemass,'mass')
        cl1 = clust{2,1}; %cond1>cond2
        cl2 = clust{2,2}; %cond2>cond1
    end
    %loading mask for nifti image (%%%
    fname = S.maskFA;
%     fname = '/scratch5/MINDLAB2017_MEG-LearningBach/DTI_Portis/tbss_bis/stats/mean_FA_skeleton_mask.nii.gz';
    aalnii = load_nii(fname);
    %getting MNI coordinates (you may get a warning since the skeletonized image is not exactly in MNI space, but close enough; I would not worry too much about that at the moment)
    [ mni_coords, xform ] = osl_mnimask2mnicoords(fname);
    %getting AAL information to provide some help to understand the results 
    parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_1mm.nii.gz'; %load this from the provided codes folder
%     parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_2mm_try.nii.gz'; %load this from the provided codes folder    
    %loading AAL labels
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_label_cerebellum.mat'); %load this from the provided codes folder
    %extracting AAL coordinates information
    K = nii.load(parcelfile);
    %getting non-zero indices of original data (useful for later calculation of clusters)
    ind = find(V_t(:,:,:,1) ~= 0); %getting indices of non-0 voxels
    [i1,i2,i3] = ind2sub(SS,ind); %reshaping the vector into a matrix shaped as the original data
    if ~exist([S.outdir '/signclust_images'],'dir')
        mkdir([S.outdir '/signclust_images'])
    end
    %creating image with only significant clusters (cond1>cond2)
    load([S_ori.outdir '/' S_ori.analysis_name '_origclust_cond1cond2.mat'])
    disp(['saving cluster ' S.sizemass ' cond1>cond2 images and excel files..'])
    for ii = 1:(length(cl1)-1) %over significant clusters; cl1 corresponds to significant clusters (-1 because the last one would actually be non-significant.. long but for now pointless issue..)
        Vsign = zeros(SS(1),SS(2),SS(3)); %initializing matrix for each significant cluster
        %create excel file with cluster information
        PDn = cell((cl1(ii)+2),6); %creating cell (all significant voxels plus two lines of legend)
        %legend
        PDn{1,1} = 'Brain region'; PDn{1,2} = 'Hemisphere'; PDn{1,3} = 'T-stat'; PDn{1,4} = 'MNI coordinates'; %1st row
        PDn{2,4} = 'X'; PDn{2,5} = 'Y'; PDn{2,6} = 'Z'; %2nd row
        for jj = 1:length(Idx1t{ii}) %over voxels of the significant clusters
            Vsign(Idx1t{ii}(jj),Idx2t{ii}(jj),Idx3t{ii}(jj)) = V_t(Idx1t{ii}(jj),Idx2t{ii}(jj),Idx3t{ii}(jj));
            PDn{jj+2,3} = V_t(Idx1t{ii}(jj),Idx2t{ii}(jj),Idx3t{ii}(jj)); %storing t-values
            %finding correspondance between non-zero indices of original data and the specific significant voxel jj
            az = double(Idx1t{ii}(jj) == i1);
            bz = double(Idx2t{ii}(jj) == i2);
            cz = double(Idx3t{ii}(jj) == i3);
            indl = find((az + bz + cz) == 3); %such value will be used to get the MNI coordinates which are stored in 1D vector correspondent to the non-zero vector of 1D cooridnates of the original data
            %using such index to get the correspondent MNI coordinates (and storing such coordinates)
            PDn{jj+2,4} = (mni_coords(indl,1)*(-1))-1; %x (trick to get the proper coordinate.. for some reasons we otherwise get a very sligthly different coordinate..)
            PDn{jj+2,5} = mni_coords(indl,2); %y
            PDn{jj+2,6} = mni_coords(indl,3); %z
            hjk = K(Idx1t{ii}(jj),Idx2t{ii}(jj),Idx3t{ii}(jj));
            if hjk ~= 0
                PDn{jj+2,1} = lab(hjk,3:end); %getting AAL label (CONSIDER TO CHANGE THIS TO ANOTHER ATLAS..)
            end
            if (mni_coords(indl,1)*(-1))-1 > 0 %storing hemisphere
                PDn{jj+2,2} = 'R';
            else
                PDn{jj+2,2} = 'L';
            end
        end
        %elaborated way to get voxels sorted in descendent order..
        [lk,lki] = sort(cell2mat(PDn(3:size(PDn,1),3)),'descend');
        caz = zeros(size(PDn,1),1); caz(1) = 1; caz(2) = 2; caz(3:size(PDn,1)) = lki + 2;
        PDn = PDn(caz,:);
        %saving image (nifti)
        aalnii.img = Vsign;
        save_nii(aalnii,[S.outdir '/signclust_images/' S.analysis_name '_cond1cond2_clust_' num2str(ii) '.nii.gz']);
        disp([num2str(ii) ' / ' num2str(length(cl1)-1)])
        %saving excel file
%         PDnn = cell2table(PDn(~any(cellfun('isempty',PDn),2),:)); %remove the possible empty cell
        PDnn = cell2table(PDn); %remove the possible empty cell
        writetable(PDnn,[S.outdir '/signclust_images/' S.analysis_name '_cond1_cond2.xlsx'],'Sheet',ii) %printing excel file
    end

    %creating image with only significant clusters (cond2>cond1)
    load([S_ori.outdir '/' S_ori.analysis_name '_origclust_cond2cond1.mat'])
    disp(['saving cluster ' S.sizemass ' cond2>cond1 images and excel files..'])
    for ii = 1:(length(cl2)-1) %over significant clusters; cl1 corresponds to significant clusters (-1 because the last one would actually be non-significant.. long but for now pointless issue..)
        Vsign = zeros(SS(1),SS(2),SS(3)); %initializing matrix for each significant cluster
        %create excel file with cluster information
        PDn = cell((cl2(ii)+2),6); %creating cell (all significant voxels plus two lines of legend)
        %legend
        PDn{1,1} = 'Brain region'; PDn{1,2} = 'Hemisphere'; PDn{1,3} = 'T-stat'; PDn{1,4} = 'MNI coordinates'; %1st row
        PDn{2,4} = 'X'; PDn{2,5} = 'Y'; PDn{2,6} = 'Z'; %2nd row
        for jj = 1:length(Idx1t{ii}) %over voxels of the significant clusters
            Vsign(Idx1t{ii}(jj),Idx2t{ii}(jj),Idx3t{ii}(jj)) = V_t(Idx1t{ii}(jj),Idx2t{ii}(jj),Idx3t{ii}(jj));
            PDn{jj+2,3} = V_t(Idx1t{ii}(jj),Idx2t{ii}(jj),Idx3t{ii}(jj)); %storing t-values
            %finding correspondance between non-zero indices of original data and the specific significant voxel jj
            az = double(Idx1t{ii}(jj) == i1);
            bz = double(Idx2t{ii}(jj) == i2);
            cz = double(Idx3t{ii}(jj) == i3);
            %%%SISTEMARE MNI COORDINATES FLIPPED RIGHT AND LEFT!!!
            %%%ANCHE LEFT AND RIGHT QUELLO A FIANCO 6 INVECE CHE 7
            indl = find((az + bz + cz) == 3); %such value will be used to get the MNI coordinates which are stored in 1D vector correspondent to the non-zero vector of 1D cooridnates of the original data
            %using such index to get the correspondent MNI coordinates (and storing such coordinates)
            PDn{jj+2,4} = (mni_coords(indl,1)*(-1))-1; %x
            PDn{jj+2,5} = mni_coords(indl,2); %y
            PDn{jj+2,6} = mni_coords(indl,3); %z
            hjk = K(Idx1t{ii}(jj),Idx2t{ii}(jj),Idx3t{ii}(jj));
            if hjk ~= 0
                PDn{jj+2,1} = lab(hjk,3:end); %getting AAL label (CONSIDER TO CHANGE THIS TO ANOTHER ATLAS..)
            end
            if (mni_coords(indl,1)*(-1))-1 > 0 %storing hemisphere
                PDn{jj+2,2} = 'R';
            else
                PDn{jj+2,2} = 'L';
            end
        end
        %elaborated way to get voxels sorted in descendent order..
        [lk,lki] = sort(cell2mat(PDn(3:size(PDn,1),3)));
        caz = zeros(size(PDn,1),1); caz(1) = 1; caz(2) = 2; caz(3:size(PDn,1)) = lki + 2;
        PDn = PDn(caz,:);
        %saving image (nifti)
        aalnii.img = Vsign;
        save_nii(aalnii,[S.outdir '/signclust_images/' S.analysis_name '_cond2cond1_clust_' num2str(ii) '.nii.gz']);
        disp([num2str(ii) ' / ' num2str(length(cl2)-1)])
        %saving excel file
%         PDnn = cell2table(PDn(~any(cellfun('isempty',PDn),2),:)); %remove the possible empty cell
        PDnn = cell2table(PDn); %remove the possible empty cell
        writetable(PDnn,[S.outdir '/signclust_images/' S.analysis_name '_cond2_cond1.xlsx'],'Sheet',ii) %printing excel file
    end
end




end


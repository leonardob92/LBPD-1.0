function [ PP ] = BrainSources_MonteCarlosim_3D_LBPD_D( S )

% It identifies 3D spatial clusters (by using islands3D_LBPD_D.m) on a binarized
% 3D matrix and then computes Monte Carlo simulation for testing the significance of the
% clusters in the original data. The clusters always refer to only one direction of the
% contrast (e.g. cond1 > cond2 or viceversa, but not both at the same
% time).
% This function does not support data evolving in time (the best use of it
% is with data that has been compressed within one time-window to only one time-point).
% It uses functions for reading and writing nifti images that are connected
% to OSL/FSL.




%  INPUT:   -S.data:        actual data (binarized p-values) (double).
%                           1st x 2nd x3rd spatial dimension
%           -S.T:           path to image with t-values or 3D matrix with actual t-values stored in the 3D brain layout
%           -S.permut:      number of permutations for Monte Carlo simulation
%           -S.clustmax:    set 1 for only max cluster size of each permutation MCS (more strict).
%                           set 0 for every size of each cluster detected for each permutation MCS (less strict).
%           -S.permthresh:  threshold for considering significant the size of the original clusters
%                           (expressed between 0 and 1; e.g. 5% = 0.05)
%           -S.parcelfile:  path to nifti images of MNI brain (e.g. AAL in 8-mm brain T1 MNI152)
%           -S.labels:      path to file with labels of parcellation (e.g. AAL)
%           -S.MNIcoords:   MNI coordinates of MNI brain correspondent to the one used for the parcellation and
%                           provided in S.parcelfile(e.g. 8-mm brain T1 MNI152).
%           -S.mask:        mask of the brain layout you have your results in
%           -S.outdir:      directory where you want to save the results image in
%           -S.anal_name:   name for the analysis (used to identify and save image and results)


%  OUTPUT:  -PP:            matrix with significant clusters information:
%                               -1st col: cluster #
%                               -2nd col: cluster size
%                               -3th col: MCS p-value for the cluster
%                               -4th col: cell with information about the voxels forming the cluster (e.g. stats, label and MNI coordinates)
%           -PP is also saved on disk together with the structure S
%           -Image with t-values for the voxels forming the significant clusters are also saved on disk (now all clusters combined together in only one image)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 27/01/2021
% Leonardo Bonetti, Aarhus, DK, 27/02/2023


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    








%input variables (think if you want to check them later..) and setting some parameters
dummyt = S.data;
permut = S.permut;
clustmax = S.clustmax;
permthresh = S.permthresh;
mask = S.mask;
SS = size(mask);
if ischar(S.T) %if you provide path
    T = load_nii(S.T); %loading image
    T2 = T.img; %extracting data
else
    T2 = S.T; %otherwise assigning the data that was provided by the user
%     T = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %loading a rnadom mask just to have the structure for saving the nifti image later..
end

%actual computation
%getting data
dummyt(dummyt==0) = NaN; %assigning 'NaN' to 0 values (that corresponds to the non-significant voxels)
[~,k_info_t,Idx1t,Idx2t,Idx3t] = islands3D_LBPD_D(dummyt); %extracting binary patterns (islands)
if k_info_t ~= 0
    disp(['clusters found in tvalues are n = ' num2str(length(k_info_t(:,1)))])
    %permutations
    ind = find(mask > 0); %getting indices of non-0s channels in brain layout
    clustperm = [];
    maxclust = zeros(permut,1);
    for pp = 1:permut
        disp(['computing permutation ' num2str(pp) ' - ' num2str(permut)])
        idx_dummy = randperm(length(ind)); %create a permuted array from 1 to length(ind)
        idx_2 = ind(idx_dummy); %shuffle the indexes in ind according to the order in idx_dummy
        r_dummy = zeros(SS(1)*SS(2)*SS(3),1); %initialise vector
        r_dummy(ind) = dummyt(idx_2); %taking element in matrix data with index idx_2
        r_resh = reshape(r_dummy,[SS(1),SS(2),SS(3)]); %reshaping the vector into a matrix shaped as the original data
        r_resh(r_resh==0) = NaN; %removing 0 values as done above for the actual data
        [~,k_info_r,~,~,~] = islands3D_LBPD_D(r_resh); %calculating clusters (islands) for permuted data
        if k_info_r == 0 %checking if some clusters have been found
            maxclust(pp,1) = 0;
            k_info_r = [0, 0, 0];
        else
            maxclust(pp,1) = max(k_info_r(:,2)); %storing maximum size of permuted clusters
        end
        clustperm = cat(1,clustperm,k_info_r); %storing all of the sizes of the detected permuted clusters
    end
    %sorting parameters of permuted clusters
    if clustmax == 1 %maximum sizes of permuted clusters (this usually tends to a decente normal distribution a bit shifted towards left..); to my understanding it is fine
        dummysort = sort(maxclust,'descend');
    else %or all of the sizes of the permuted clusters (this usually tends to only the right side of a normal distribution with a very peaky mean..); to my understanding it is fine
        dummysort = sort(clustperm(:,2),'descend');
    end
    threshfinal = dummysort(floor((permthresh*100*length(dummysort)/100) + 1)); %getting the final threshold for significance
    %storing the output
    d = find(k_info_t(:,2) > threshfinal); %looking for significant clusters (according to their sizes)
    cc2 = zeros(length(d),1);
    for kk = 1:length(d) %over significant clusters
        cc = find(k_info_t(d(kk),2) > dummysort); %getting a specific p-value by.. finding how many times original significant cluster was larger than permuted ones
        cc2(kk) = (length(dummysort) - length(cc))/length(dummysort); %then getting false positive (amount of simulated clusters - times when original significant cluster was larger than permuted ones) and dividing them by amount of simulated clusters
    end
    %getting AAL information to provide some help to understand the results
    load(S.labels);
    %extracting AAL coordinates information
    K = nii.load(S.parcelfile);
    mkdir([S.outdir '/Temp'])
    %storing significant clusters and preparing output
    PP = cell(length(d)+1,4);
    PP{1,1} = 'Cluster #'; PP{1,2} = 'Cluster size'; PP{1,3} = 'MCS p-value'; PP{1,4} = 'Voxels information';
    for hh = 1:length(d) %over significant clusters
        PP(hh+1,1) = {hh}; %new ID of the clusters
        PP(hh+1,2) = {k_info_t(d(hh),2)}; %sizes
        PP(hh+1,3) = {cc2(hh,1)}; %p-values of each cluster
        %preparing image with results to be saved for cluster dd
        fname = [S.outdir '/Temp/' S.anal_name '_SignClust_' num2str(hh) '_Tvals.nii.gz']; %path and name
        V = zeros(SS(1),SS(2),SS(3));
        for jj = 1:length(Idx1t{d(hh)}) %over voxels forming the cluster dd
            V(Idx1t{d(hh)}(jj),Idx2t{d(hh)}(jj),Idx3t{d(hh)}(jj)) = T2(Idx1t{d(hh)}(jj),Idx2t{d(hh)}(jj),Idx3t{d(hh)}(jj)); %for the voxels forming the cluster dd, storing the correspondent t-values
        end
        %saving results as nifti image
        if ischar(S.T) %if you provide path
            T.img = V;
            save_nii(T,fname);
        else
            maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
            nii = make_nii(V,[8 8 8]);
            nii.img = V; %storing matrix within image structure
            nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
            disp(['saving nifti image'])
            save_nii(nii,fname); %printing image
        end
        %getting MNI coordinates of significant voxels within the provided image
        [ mni_coords, ~ ] = osl_mnimask2mnicoords(fname);
        %extracting statistics
        VV = V(V~=0);
        %indices of non-zero values of nifti image
        VI = find(V~=0);
        %sorting results in order to have strongest voxels at the top (positive t-values) or at the bottom (negative t-values)
        [VV2, II] = sort(VV,'descend');
        VI = VI(II);
        mni_coords = mni_coords(II,:);
        %getting AAL indices
        ROI = zeros(length(VI),1);
        cnt = 2;
        %final cell
        PDn = cell(length(Idx1t{d(hh)})+2,6);
        %setting legend
        PDn{1,1} = 'Brain Region'; PDn{1,2} = 'Hemisphere'; PDn{1,3} = 't-value'; PDn{1,4} = 'MNI Coordinates';
        PDn{2,4} = 'x'; PDn{2,5} = 'y'; PDn{2,6} = 'z';
        for ii = 1:length(VI)
            ROI(ii) = K(VI(ii));
            cnt = cnt + 1;
            PDn(cnt,4) = {mni_coords(ii,1)}; %storing MNI coordinates
            PDn(cnt,5) = {mni_coords(ii,2)}; %storing MNI coordinates
            PDn(cnt,6) = {mni_coords(ii,3)}; %storing MNI coordinates
            if mni_coords(ii,1) > 0 %storing hemisphere
                PDn(cnt,2) = {'R'};
            else
                PDn(cnt,2) = {'L'};
            end
            PDn(cnt,3) = {round(VV2(ii),2)}; %storing t-statistics
            if ROI(ii) > 0 && ROI(ii) < 91
                PDn(cnt,1) = {lab(ROI(ii),3:end)}; %storing ROI label (if it belongs to AAL parcellation; consider to change atlas in the future..)
            end
        end
        PP(hh+1,4) = {PDn}; %storing coordinates and label names
    end
    %Combining images with significant clusters
    if ~isempty(d) %if there are significant clusters
        path = [S.outdir '/Temp']; %general path
        listt = dir([S.outdir '/Temp/*gz']); %list of images (one image for significant MEG source cluster)
        cmd = ['fslmaths /']; %starting the command to use 'flsmaths'
        for ll = 1:length(listt) %over images
            if ll ~= length(listt)
                dumo = [path(2:end) '/' listt(ll).name ' -add /']; %getting path and name of each image
            else
                dumo = [path(2:end) '/' listt(ll).name ' /']; %getting path and name of each image
            end
            cmd = strcat(cmd,dumo); %combining path of different images
        end
        output = [S.outdir(2:end) '/' S.anal_name '.nii.gz']; %creating output name
        cmd = strcat(cmd,output); %final command line
        system(cmd) %submitting command line to terminal
        rmdir([S.outdir '/Temp'],'s') %removing temporary folder with independent images for each figure
    end
    %saving cluster results
    save([S.outdir '/' S.anal_name '_SignClust.mat'],'PP','S');
end

end


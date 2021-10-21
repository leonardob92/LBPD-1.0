function [ idk3, kkk2, KJK ] = FunctionalSpatialClustering_voxels2ROIs_LBPD_D( S )

% It computes kmeans clustering to define a functionally-based
% parcellation. It usually reduces 3559 voxels into a set of approximately
% 10-30 parcels.
% The algorithm works combining the temporal and spatial profiles of each
% voxel. The temporal profile refers to the peak value or the time-index of
% the peak value of each voxel.
% I suggest to use the function as follows:
% 1) Run the temporal clustering (based on the maximum value or the
% correspondent time-index) for each voxel for a number of different
% clustering solutions (e.g. 2:20).
% 2) Evaluate the best clustering solution (evaluation can be made by using
% the elbow method and the silhouette method that are automatically
% calculated by the current function. Please note that there is no perfect
% way to clusterize data and therefore your choices should be made on the
% basis of the suggestions of the current function as well as on your
% critical thinking and specific necessities.
% 3) Set the temporal clustering to the number of clusters that you decided
% and then run again the kmeans clustering algorithm, but this time on the
% spatial profile of the voxels. You may run here cluster solutions k =
% 2:10 or 2:20.
% 4) Evaluate for each of the temporal cluster which is the best clustering
% solution on a spatial level.
% 5) Choose the best spatial clustering solution for each of the temporal
% clusters.
% 6) Run the function one last time by inserting the proper number of
% clusters for temporal and spatial profiles.
% The algorithm requires functions from nifti toolbox, if you want to plot
% the nifti images. You can donwload it here: https://github.com/mcnablab/NIFTI_toolbox





%   INPUT:

%%% computation settings %%%
%           -S.data:                        matrix with data (voxels x time-points
%           -S.timep:                       time-samples to do the computation on (useful if the data matrix has too large time-series)
%           -S.time:                        time in seconds (its length can be bigger than S.timep, but not bigger!)
%           -S.randsim:                     simulation with random numbers:
%                                               -0: no simulations
%                                               -1: simulation on random numbers (rand function in Matlab)
%                                               -2: simulation by randomising the order of the actual data
%           -S.ROItimeseries_rand:          1 for plottinh time-series of randomly simulated data
%           -S.rperm:                       number of permutation for randomly simulated data
%           -S.clustsol:                    number of temporal clusters for different kmeans clustering solutions
%                                           (e.g. 3:5 means that you want to test cluster solutions with 3, 4 and 5 clusters)
%                                           YOU SHOULD START BY TESTING AND CHOOSING THE TEMPORAL CLUSTER AS A FIRST STEP
%           -S.clustspatial:                spatial clusters solutions to be tested (e.g. 2:10 computes spatial clustering from k = 2 to k = 10 cluster)
%                                           Based on the 3D coordinates of each voxel the algorithm clusterizes them spatially.
%                                           This operation is done for each of the
%                                           previously defined temporal clusters.
%                                           Set to 0 not to run it (useful when you evaluate the temporal clustering).
%                                           YOU SHOULD USE THIS AS SECOND STEP TO EVALUATE THE BEST SPATIAL SOLUTIONS OF THE CLUSTERS AFTER THAT YOU COMPUTED THE TEMPORAL ONES
%           -S.clspatvect:                  %specify different spatial cluster solutions for different temporal clusters (NOTE that this makes sense if clustsol has only 1 cluster
%                                           solution and if that solution is equal to the length of clspatvect!!).
%                                           YOU SHOULD USE THIS AS THIRD STEP TO DEFINE THE PRECISE CLUSTERING SOLUTIONS THAT YOU WANT FOR THE SPATIAL DOMAIN (WHILE KEEPING FIXED
%                                           THE PREVIOUSLY CHOSEN FOR THE TEMPORAL DOMAIN).
%           -S.PCA (ONLY EXPLORATORY!!):    1 to use Monte-Carlo based PCA to obtain parcels timeseries (reducing the timeseries of their voxels)
%                                           0 to use simple arithmetic mean
%           -S.PCA_method:                  method for normalizing the sign of the eigenvectors if you use PCA for "averaging" the voxels of each spatio-temporal parcel
%                                           You can choose either: 'occurrences', 'max_abs', 'average'.
%                                           SUGGESTED: either 'max_abs' or 'average'
%           -S.max_value:                   max value or time-index of the max value for clustering algorithms..
%                                               -1 for temporal clusterization on maximum value;
%                                               -0 for temporal clusterization on time index of maximum value
%           -S.percentmax:                  value (percentage) of maximum values to be considered for clustering purposes (e.g 100 uses all voxels; 5 only the top 5% of the voxels)


%%% plotting settings %%%
%           -S.ROItimescol:                 0 for different colors (time-series).. 1 for different shades (hues) of red
%           -S.plot.scatterl:               1 for scatterplots (of maximum values or indices of voxels).. 0 otherwise
%           -S.plot.viol:                   if you set 1 for "S.plot.scatterl', then you can choose to have violins and scatter plots together (S.plot.viol = 1) or scatter plots only (S.plot.viol = 0)
%           -S.plot.brainv:                 1 for plotting in the brain (nifti file).. 0 otherwise.
%                                           Please, note that the final clustering solution (obtained when "S.clspatvect" is defined) automatically generates this particular plot.
%                                           Set this to 0 can be useful when you test several different temporal and spatial clustering solutions.. in that case you may want to
%                                           produce only few particular plotting within the brain..
%           -S.plot.ROItimeseries:          1 to get timeseries of clusters (for now by computing the mean over the voxels belonging to each different cluster..)
%           -S.plot.ROItimeseries_rand:     same but for randomised data..
%           -S.plot.brain_inv_each_parcel:  1 to get plotting within the brain template, separately for each parcel.
%                                           This works only for spatio-temporal clustering solutions (if you want temporal clustering only, you should write as in the following examples:
%                                           (e.g. S.clustsol = 6; S.clspatvect = [1 1 1 1 1 1])).



%   OUTPUT: -idk3:                          voxels stored into different parcels (kmeans clustering on temporal indices/maximum values)
%           -kkk2:                          voxels stored into different parcels (kmeans clustering on temporal indices/maximum values (i) and on spatial coordinates (ii))
%           -KJK:                           timeseries for each parcel (works only for spatio-temporal clustering solutions).
%                                           (if you want temporal clustering only, you should write as in the following examples: (e.g. S.clustsol = 6; S.clspatvect = [1 1 1 1 1 1])).
%                                           3th dimension: 1st: mean over voxels; 2nd: standard error over voxels 







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 17/06/2020


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%preparing mask for later production of nifti image(s)
maskk = load_nii(S.options.mask_fname); %getting the mask for creating the figure
% idm = find(maskk.img ~= 0); %getting indices of non-0 elements of the matrix
%preparing data and options..
KJK = 0;
idk3 = 0;
kkk2 = 0;
P = S.data(:,S.timep)'; %getting and transposing actual data matrix
timesel = S.time(S.timep); %selecting proper time-points in seconds
if S.randsim == 1 %randomised numbers (rand function in Matlab)
    R = rand(size(P,1),size(P,2))'; %generating random data matrix
elseif S.randsim == 2 %randomised order of actual data
    %randomizing data..
    idx_dummy = randperm(size(P,1)*size(P,2)); %create a permuted array from 1 to size of original data vector
    r_dummy = zeros(size(P,1)*size(P,2),1); %initialise new vector
    r_dummy(1:size(P,1)*size(P,2)) = P(idx_dummy); %taking element in matrix data with index idx_dummy
    r_resh = reshape(r_dummy,[size(P,1),size(P,2)]); %reshaping the vector into a matrix shaped as the original data
    R = r_resh; %this is because I do later operations on P..
end
%maximum values or temporal indices of maximum values
if S.max_value == 1
    [im,mc] = max(abs(P)); %getting maximum values and indices of peaks of matrix P (absolute values since it makes sense only to work on positive absolute value data..)
    thck = round((length(im)*S.percentmax)/100); %finding indices o percentmax % of maximum values
    thck2 = sort(im,'descend');
    im(im<thck2(thck))=0; %thresholding values (assigning 0s to low-value voxels in order to help the kmeans clustering algorithm to work more properly)
else
    [mc,im] = max(P); %getting maximum values and indices of peaks of matrix P
    thck = round((length(im)*S.percentmax)/100); %finding indices o percentmax % of maximum values
    thck2 = sort(im,'descend');
    im(im<thck2(thck))=0; %thresholding values (assigning 0s to low-value voxels in order to help the kmeans clustering algorithm to work more properly)
    im = im + S.timep(1) - 1; %this is a trik to visualize better the results.. if you select data from time-point 1 there is no problem, but if you select data from time-point x (x~=1) then in the scatter plots you would not visualize the proper time-point of maxium value (with regards wth the original data indices).. therefore I have to sum the first time-point of the selected time-window (and subtract 1..)
end

%some default color specifications for later plotting..
cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
cl(1, :) = cb(4, :); %same..
cnt = 0; %initialising counter
if length(S.clustspatial) == 1
    cntkkk = -S.clustspatial;
    KKK = zeros(size(S.data,1),1);
end
Ak = cell(length(S.clustsol),4);
SFkk = cell(length(S.clustsol)-1,1);
ckv = 0;
%checking if user request about spatial clustering makes sense
if length(S.clustsol) == 1 && ~isempty(S.clspatvect)
    if length(S.clspatvect) ~= S.clustsol(1)
        error('if you ask for different spatial cluster solutions you need to have them in the same number as the temporal clusters.. so length(clspatvect) must be equal to clustsol..')
    end
end

%actual computation
EVA = cell(length(S.clustsol),2); %spatial cluster evaluation matrix
for ol = S.clustsol
    cnt = cnt + 1;
    clustnum = ol;
    %preallocating space    
    CO = zeros(clustnum,1);
    A = zeros(1,1);
    %%%%%%%%%%%%%%%% 1) kmeans on temporal domain.. %%%%%%%%%%%%%%%%
%     opts = statset('Display','final');
%     [idk2,C,smd] = kmeans(im',clustnum,'Replicates',100,'Options',opts); %kmeans
    [idk2,C,smd] = kmeans(im',clustnum,'Replicates',100); %kmeans
    A(1,1) = sum(smd);
    %storing cluster sizes for later calculations..
    poi = zeros(clustnum,1);
    imclustI = cell(clustnum,1);
    for jj = 1:clustnum
        imclustI{jj} = im(idk2'==jj); %selecting voxels of cluster jj
        poi(jj,1) = length(imclustI{jj}); %storing cluster size
    end
    [poiP,poiI] = sort(poi,'descend'); %sorting clusters
    CO(:,1) = C(poiI);
    %since different kmeans runs can easily return the same clusters
    %with different indices, here I assign the same indices to the
    %clusters sorted in descended order.. this is very important for
    %later calculations as well as for clarity (e.g. now cluster
    %numbers correspond between scatter plot, information on the
    %clusters (sorted in descendent order) and plots in the brain (if
    %run..)
    idk3 = zeros(length(idk2),1);
    for hh = 1:length(poiI)
        idk3(idk2==poiI(hh)) = hh;
    end
    %plotting within brain template
    if S.plot.brainv == 1
        if ~exist(S.outdir,'dir')
            mkdir(S.outdir)
        end
        fname2 = [S.outdir '/tempclust_k_' num2str(clustnum) '_time_' num2str(S.timep(1)) '_' num2str(S.timep(end)) '_actualdata'];
        %using nifti toolbox
        dumimg = maskk.img; %getting the image matrix
        dumimg(find(maskk.img ~= 0)) = idk3; %assigning values of the new variable to the non-0 indices of the image matrix
        nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
        nii.img = dumimg; %storing matrix within image structure
        nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
        disp('saving nifti image')
        save_nii(nii,fname2); %printing image
        %previous solution using functions from OSL
%         fname_out = nii.quicksave(idk3,fname2,S.options) %plotting as nifti..
    end
    %scatter plots (temporal clustering)
    if S.plot.scatterl == 1
        %preparing colormap..
        if S.ROItimescol ~= 1
            cmp = colormap(lines(S.clustsol(end))); %extracting some colours from a colormap
        else
            cmp = zeros(clustnum,3); %producing a colormap with proressive darker shades of red
            q = 1/clustnum;
            for uu = 1:clustnum
                cmp(uu,1) = 1-(q*uu);
            end
        end
        %producing scatter plots (with clusters sorted by their cluster sizes (2 separate loops since we want to have also the clusters in the scatter plots to be sorted by their sizes)
        figure
        d2 = cell(1,clustnum);
        for jj = 1:clustnum %here there is probably a better way to do the scatter plot.. without doing the loop.. but for now this is fine..
            imclustI2 = imclustI{poiI(jj)}; %plotting clusters sorted by their size (descending order)
            if S.plot.viol ~= 1
                %codes for scatter plot
                scatter(ones(length(imclustI2),1)*jj,imclustI2,[],cmp(jj,:)) %scatter plotting..
                hold on
                scatter(jj,CO(jj),'kx','LineWidth',1);
                hold on
            else
                %codes for violins and scatter plots together
                d2(1,jj) = {imclustI2}; %within each cell a px1 double of observations where p is the number of observations
            end
        end
        set(gcf,'Color','w')
        grid minor
        if S.plot.viol == 1
            d2 = d2'; %transposing it..
            %making figure
            h = rm_raincloud2(d2,cl); %running the function. This function was develped by Micah and colleagues. I called it "2" simply because I think I made some changes and therefore I did not want to make any confusion.. however the total credit for this function belongs to Micah and colleagues
            grid minor
            set(gcf,'Color','w')
            grid minor
        end
    end
    %preparing colormap..
    if S.ROItimescol ~= 1
        cmp = colormap(lines(S.clustsol(end))); %extracting some colours from a colormap
    else
        cmp = zeros(clustnum,3); %producing a colormap with proressive darker shades of red
        q = 1/clustnum;
        for uu = 1:clustnum
            cmp(uu,1) = 1-(q*uu);
        end
    end
    %getting timeseries of brain areas
    if S.plot.ROItimeseries == 1
        P2 = P'; %transposing data matrix
        jg = zeros(clustnum,size(P2,2)); %preallocating space
        figure
        for iii = 1:clustnum
            %getting voxels of cluster ii and combining them into a single ROI
            jg(iii,:) = mean(P2(idk3==iii,:),1); %for now trying with mean..
            plot(timesel,jg(iii,:),'Color',cmp(iii,:),'DisplayName',[num2str(iii)])
            hold on
        end
        grid minor
        legend('show')
        xlim([timesel(1) timesel(end)])
        set(gcf,'Color','w')
    end
    
    %%%%%%%%%%%%%%%% 2) kmeans clustering on spatial domain.. %%%%%%%%%%%%%%%%
    if sum(double(S.clustspatial ~= 0)) ~= 0 || sum(double(S.clspatvect ~= 0)) ~= 0
        disp('computing spatial clustering on different temporal cluster solutions')
        SFk = cell(clustnum,1);
        eva2 = zeros(1,clustnum);
        for zz = 1:clustnum %over kmeans temporal clusters
            MNI = S.coords_template(idk3==zz,:); %getting MNI coordinates of cluster zz
            SF = zeros(1,length(S.clustspatial));
            cnf = 0;
            if isempty(S.clspatvect) %if user does not want to have different spatial cluster solutions for different temporal clusters..
                for ff = S.clustspatial %over spatial cluster solutions
                    cnf = cnf + 1;
                    if size(MNI,1) >= ff
                        disp(['temporal cluster ' num2str(zz) ' / ' num2str(clustnum) ' - spatial cluster ' num2str(ff) ' / ' num2str(S.clustspatial)])
%                         [idk2f,Cf,smdf] = kmeans(MNI,ff,'Replicates',100,'Options',opts); %kmeans
                        [idk2f,Cf,smdf] = kmeans(MNI,ff,'Replicates',100); %kmeans
                        SF(1,cnf) = sum(smdf); %saving sum of square distances
                    end
                end
                %trying out some cluster evaluation algorithms (then I chose silhouette..)
                eva = evalclusters(MNI,'kmeans','silhouette','KList',S.clustspatial);
                eva2(1,zz) = eva.OptimalK;
                figure
                plot(eva)
                grid minor
                title(['silhouette - spatial clustering k = ' num2str(zz) ' of temporal clustering k = ' num2str(clustnum)])
                %                     eva2 = evalclusters(MNI,'kmeans','CalinskiHarabasz','KList',clustspatial)
                %                     eva3 = evalclusters(MNI,'kmeans','DaviesBouldin','KList',clustspatial)
                %                     eva4 = evalclusters(MNI,'kmeans','gap','KList',clustspatial)
                set(gcf,'Color','w')
                
                if length(S.clustspatial) == 1 %working only if the requested spatial cluster solutions are = 1 (otherwise it would be messy with many calculation and plots.. here the idea is that first user decides how many temporal and spatial clusters and then it asks for plotting..)
                    cntkkk = cntkkk + S.clustspatial;
                    KKK(idk3==zz,1) = idk2f + cntkkk; %assigning a new progressive ID to the definitive clusters
                end
            else %otherwise user provides a specific request of different spatial cluster solutions for different temporal clusters
%                 [idk2fv] = kmeans(MNI,S.clspatvect(zz),'Replicates',100,'Options',opts); %kmeans
                [idk2fv] = kmeans(MNI,S.clspatvect(zz),'Replicates',100); %kmeans
                KKK(idk3==zz,1) = idk2fv + ckv;
                ckv = ckv + S.clspatvect(zz); %assigning a new progressive ID to the definitive clusters
            end
            SFk(zz,1) = {SF}; %storing sum of square distances of spatial clusters for each of the clusters emerged for kmeans temporal..
        end
        if length(S.clustspatial) == 1 || ~isempty(S.clspatvect) %for now it makes sense only if we have 1 cluster spatial configuration.. (otherwise it would be quite messy..)
            lop = zeros(max(KKK),1); %working with maximum value of KKK that would be the maximum ID number of the spatial clusters
            for hhk = 1:max(KKK)
                lop(hhk,1) = length(find(KKK==hhk));
            end
            [kkkn,kkki] = sort(lop,'descend'); %sorting cluster sizes (temporal and spatial clusters together)
            kkk2 = zeros(size(S.data,1),1);
            for hhk = 1:max(KKK)
                kkk2(KKK==kkki(hhk)) = hhk; %reassinging cluster ID in order to have them sorted (descending order)
            end
            if isempty(S.clspatvect) %if empty, each temporal cluster is divided into an equal number of spatial clusters
                if ~exist(S.outdir,'dir')
                    mkdir(S.outdir)
                end
                fname2 = [S.outdir '/tempclust_k_' num2str(clustnum) '_spatialclust_k_' num2str(S.clustspatial) '_time_' num2str(S.timep(1)) '_' num2str(S.timep(end)) '_actualdata'];
            else %otherwise the number of spatial clusters varies depending on the user request
                if ~exist(S.outdir,'dir')
                    mkdir(S.outdir)
                end
                fname2 = [S.outdir '/tempclust_k_' num2str(clustnum) '_spatialclust_clspatvect_time_' num2str(S.timep(1)) '_' num2str(S.timep(end)) '_actualdata'];
            end
            %using nifti toolbox
            dumimg = maskk.img; %getting the image matrix
            dumimg(find(maskk.img ~= 0)) = kkk2; %assigning values of the new variable to the non-0 indices of the image matrix
            nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
            nii.img = dumimg; %storing matrix within image structure
            nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
            disp('saving nifti image')
            save_nii(nii,fname2); %printing image
            %previous solution using functions from OSL
%             fname_out = nii.quicksave(kkk2,fname2,S.options) %plotting as nifti..
            %timeseries plotting
            P2 = P'; %transposing data matrix
            if S.ROItimescol ~= 1
                cmpk = colormap(lines(max(KKK))); %extracting some colours from a colormap
            else
                cmpk = zeros(max(KKK),3); %producing a colormap with proressive darker shades of red
                q = 1/(max(KKK));
                for uu = 1:max(KKK)
                    cmpk(uu,1) = 1-(q*uu);
                end
            end
            KJK = zeros(max(KKK),length(timesel),2); %3D = 1: mean; 2: standard error between voxels
            figure
            for iiik = 1:max(KKK) %over spatio-temporal clusters
                %getting voxels of cluster (parcel) iiik and combining them into a single ROI
                H = P2(kkk2==iiik,:); %extracting timeseries of voxels forming parcel iiik
                %PCA is not really suggested, in practice.. even though in principle should be a very good approach; more work here is needed
                if S.PCA == 1
                    %last input is the method for normalizing the sign of the eigenvectors
                    TM = PCA_LBPD(H,10,0,S.PCA_method); %running MCS-based PCA; 10 refers to number of permutations (10 looks small but should be good enough; feel free to change it to a bigger number; 0 refers to not having plotting solutions)
                    jgj = TM.PCA_timeseries; %extracting PCA timeseries
                else
                    jgj = mean(H,1); %arithmetic mean
                end
%                 jgj = mean(P2(kkk2==iiik,:),1); %for now trying with mean.. THEN THINK ABOUT HAVING PCA MAYBE.. OR MEDIAN OR SOMETHING LIKE THAT..
                KJK(iiik,:,1) = jgj; %storing timeseries of parcel iiik..
                KJK(iiik,:,2) = (std(H,0,1))./sqrt(length(find(kkk2==iiik))); %storing timeseries of parcel iiik..
                plot(timesel,jgj,'Color',cmpk(iiik,:),'DisplayName',[num2str(iiik)])
                hold on
                if S.plot.brain_inv_each_parcel == 1
                    kkk3 = zeros(size(S.data,1),1); %zero vector
                    if ~exist(S.outdir,'dir')
                        mkdir(S.outdir)
                    end
                    kkk3(kkk2==iiik) = max(abs(P2(kkk2==iiik,:)),[],2); %assigning amplitude values (maximum value over time) of the voxels composing each parcel
                    % kkk3(kkk2==iiik) = mean(P2(kkk2==iiik,:),2); %assigning amplitude values (averaged over time) of the voxels composing each parcel
                    fname2 = [S.outdir '/tempclust_k_' num2str(clustnum) '_spatialclust_clspatvect_time_' num2str(S.timep(1)) '_' num2str(S.timep(end)) '_parcel_' num2str(iiik) '_of_' num2str(max(KKK)) '_parcels']; %creating name for each parcel of the current spatio-temporal parcellation
                    %using nifti toolbox
                    dumimg = maskk.img; %getting the image matrix
                    dumimg(find(maskk.img ~= 0)) = kkk3; %assigning values of the new variable to the non-0 indices of the image matrix
                    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
                    nii.img = dumimg; %storing matrix within image structure
                    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
                    disp('saving nifti image')
                    save_nii(nii,fname2); %printing image
                    %previous solution using functions from OSL
%                     fname_out = nii.quicksave(kkk3,fname2,S.options) %plotting as nifti..
                end
            end
            grid minor
            legend('show')
            xlim([timesel(1) timesel(end)])
            set(gcf,'Color','w')
        end 
        SFkk(ol,1) = {SFk}; %..in each of the different clustering solution of kmeans temporal
    end
    Ak(cnt,1) = {ol};
    Ak(cnt,2) = {A};
    if sum(double(S.clustspatial ~= 0)) ~= 0
        EVA{cnt,1} = clustnum; %storing number of k temporal clusters requested
        EVA{cnt,2} = eva2; %storing correspondent best spatial cluster solution within each cnt temporal solution
    end
    disp(['temporal cluster ' num2str(cnt) ' / ' num2str(length(S.clustsol))])
end

%randomised data (possibly useful to make some simulations..)
if S.randsim ~= 0
    cnt = 0; %initialising counter
    Akr = cell(length(S.clustsol),2);
    [mcr,imr] = max(R); %getting maximum values and indices of peaks of matrix R
    for ol = S.clustsol
        cnt = cnt + 1;
        clustnum = ol;
        %preallocating space
        Arr = zeros(1,S.rperm);
        for ppp = 1:S.rperm %over number of permutations
            %             opts = statset('Display','final'); %if you want to show the progrssion of the replications of kmeans clustering..
            %             [idk2,C,sumdist] = kmeans(imr',clustnum,'Replicates',repl,'Options',opts); kmeans (same as above)
            [idk2,C,sumdist] = kmeans(imr',clustnum,'Replicates',100); %kmeans (otherwise run this)
            Arr(1,ppp) = sum(sumdist); %sum of square distance from ech voxel time-point of maximum amplitude and the centroid of its cluster
            disp(['cluster solution ' num2str(ol) ' / ' num2str(S.clustsol(end)) ' - permutation ' num2str(ppp) ' / ' num2str(S.rperm)])
            %getting timeseries of randomised brain areas
            if S.ROItimeseries_rand == 1
                R2 = R'; %transposing data matrix
                jgr = zeros(clustnum,size(R2,2)); %preallocating space
                figure
                for iii = 1:clustnum
                    %getting voxels of cluster ii and combining them into a single ROI
                    jgr(iii,:) = mean(R2(idk2==iii,:),1); %for now trying with mean..
                    plot(timesel,jgr(iii,:),'Color',cmp(iii,:),'DisplayName',[num2str(iii)])
                    hold on
                end
                grid minor
                legend('show')
                xlim([timesel(1) timesel(end)])
                set(gcf,'Color','w')
            end
        end
        %storing SSD for randomised data
        Akr(cnt,1) = {ol};
        Akr(cnt,2) = {Arr};
    end
end


%plotting..
%plotting cluster evaluations..
if length(S.clustsol) > 1
    %cluster evaluation algorithms silhouette
    eva = evalclusters(im','kmeans','silhouette','KList',S.clustsol);
    %     evat = eva.OptimalK;
    figure
    plot(eva)
    grid minor
    title(['silhouette - temporal clustering - k = ' num2str(S.clustsol(1)) ' : ' num2str(S.clustsol(end))])
    set(gcf,'Color','w')
end

if length(S.clustsol) > 1 %doing this plotting only if you request more than one cluster solution.. otherwise it is meaningless..
    %plotting the mean of different cluster solution SSD
    figure
    Am2 = cell2mat(Ak(:,2)); %actual data
    %     Am2 = mean(cell2mat(Ak(:,4))); %actual data
    plot(Am2,'DisplayName','data')
    if S.randsim ~= 0
        hold on
        %         if randmean == 1
        %             Am2r = mean(cell2mat(Akr(:,2)),2); %mean SSD for random data (mean over permutations..)
        %             %more strict option.. just getting the minimum value..
        %         else
        Am2r = min(cell2mat(Akr(:,2)),[],2);
        %         end
        plot(Am2r,'DisplayName','rand')
    end
    grid minor
    legend('show');
    set(gcf,'Color','w')
    xticks(1:length(Am2)); %specifying x-axis ticks (where to put the labels..)
    a = cell(1,length(S.clustsol));
    for ii = 1:length(S.clustsol)
        a(1,ii) = {num2str(S.clustsol(ii))}; %assigning labels corresponding to the above-specified ticks (here we want to show the cluster sizes on the x-axis)
    end
    xticklabels([a]) %plotting the labels
    %codes to plot difference between SSD actual data and SSD randomised data and weighting them by the number of clusters..)
    if S.randsim ~= 0
        figure
        plot(abs((Am2-Am2r)'.*(S.clustsol)))
        grid minor
        set(gcf,'Color','w')
    end
    xticks(1:length(Am2)); %specifying x-axis ticks (where to put the labels..)
    xticklabels([a]) %plotting the labels
end
%plotting SSD for spatial clustering solutions (1 plot for each of the different temporal cluster solutions)
if length(S.clustspatial) > 1 && isempty(S.clspatvect)
    for jj = S.clustsol(1):S.clustsol(end)
        figure
        for ii = 1:length(SFkk{jj,1})
            plot(SFkk{jj,1}{ii,1},'DisplayName',num2str(ii))
            hold on
        end
        grid minor
        legend('show')
        set(gcf,'Color','w')
        xticks(1:length(S.clustspatial)); %specifying x-axis ticks (where to put the labels..)
        a = cell(1,length(S.clustspatial));
        for ii = 1:length(S.clustspatial)
            a(1,ii) = {num2str(S.clustspatial(ii))}; %assigning labels corresponding to the above-specified ticks (here we want to show the cluster sizes on the x-axis)
        end
        xticklabels([a]) %plotting the labels
        title(num2str(jj))
    end
end



end


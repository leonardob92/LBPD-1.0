function [ O ] = MEG_SR_Stats_twogroups_LBPD( S )
O = [];

% It computes statistics after beamforming (computed using MEG_SR_Beam_LBPD.m).
% It tests the difference between two experimental groups (using ttest2). 



%  INPUT:   -S.workingdir:         path where the data from MEG_SR_Beam_LBPD.m is stored
%           -S.plot_nifti_name:    character with name for nifti files (and saved mat files) (it may be useful if you run separate analysis)
%           -S.list1:              list of subjects (group 1) (list with paths to mat files)
%           -S.list2:              list of subjects (group 2) (list with paths to mat files)
%           -S.Aarhus_clust:       1 to use paralle computing (Aarhus University, contact me, Leonardo Bonetti, for more information)
%                                  leonardo.bonetti@clin.au.dk

%  OUTPUT:  -statistics and nifti images saved on disk (in S.workingdir)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Oxford, UK, 08/04/2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    









if S.Aarhus_clust == 1
    %LBPD_startup_D
    pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
    addpath(pathl);
    LBPD_startup_D(pathl);
end

%loading data (group 1)
load([S.list1(1).folder '/' S.list1(1).name]) %loading one subject to extract information 
data1 = zeros(size(OUT.sources_ERFs,1),size(OUT.sources_ERFs,2),size(OUT.sources_ERFs,3),length(S.list1)); %brain voxels, time-windows conditions (baby cry and other cry together), subjects
for ii = 1:length(S.list1) %over subjects
    load([S.list1(ii).folder '/' S.list1(ii).name]) %loading subject ii (sources)
    data1(:,:,:,ii) = OUT.sources_ERFs;
    disp(['loading subj ' num2str(ii) ' / ' num2str(length(S.list1)) ' - group 1'])
end
%loading data (group 2)
data2 = zeros(size(OUT.sources_ERFs,1),size(OUT.sources_ERFs,2),size(OUT.sources_ERFs,3),length(S.list2)); %brain voxels, time-windows conditions (baby cry and other cry together), subjects
for ii = 1:length(S.list2) %over subjects
    load([S.list2(ii).folder '/' S.list2(ii).name]) %loading subject ii (sources)
    data2(:,:,:,ii) = OUT.sources_ERFs;
    disp(['loading subj ' num2str(ii) ' / ' num2str(length(S.list2)) ' - group 2'])
end

%t-tests
P = zeros(size(OUT.sources_ERFs,1),size(OUT.sources_ERFs,2),size(OUT.sources_ERFs,3)); %brain voxels, time-windows
T = zeros(size(OUT.sources_ERFs,1),size(OUT.sources_ERFs,2),size(OUT.sources_ERFs,3)); %brain voxels, time-windows
for cc = 1:size(OUT.sources_ERFs,3)
    for ii = 1:size(OUT.sources_ERFs,1) %over brain voxels
        for jj = 1:size(OUT.sources_ERFs,2) %over time-points
            [~,p,~,stats] = ttest2(squeeze(data1(ii,jj,cc,:)),squeeze(data2(ii,jj,cc,:))); %contrasting baby cry versus the average of the other crys..
            P(ii,jj,cc) = 1-p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp(['brain voxel ' num2str(ii) ' / ' num2str(size(OUT.sources_ERFs,1)) ' - condition ' num2str(cc) ' / ' num2str(size(OUT.sources_ERFs,3))])
    end
end

%printing nifti images
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for cc = 1:size(OUT.sources_ERFs,3)
    %P-values
    fnamenii = [S.workingdir '/' S.plot_nifti_name '_' OUT.S.inversion.conditions{cc} 'Pval.nii.gz']; %path and name of the image to be saved
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),size(OUT.sources_ERFs,2));
    for ii = 1:size(P,1) %over brain sources
        dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = P(ii,:,cc); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]);
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image'])
    save_nii(nii,fnamenii); %printing image
    %T-values
    fnamenii = [S.workingdir '/' S.plot_nifti_name '_' OUT.S.inversion.conditions{cc} 'Tval.nii.gz']; %path and name of the image to be saved
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),size(OUT.sources_ERFs,2));
    for ii = 1:size(T,1) %over brain sources
        dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = T(ii,:,cc); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]);
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image'])
    save_nii(nii,fnamenii); %printing image
end




end


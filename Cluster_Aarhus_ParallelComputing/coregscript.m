%% on the cluster

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')


% config of the cluster server see cluster_split.m for more deatails
clusterconfig('slot', 1);
clusterconfig('long_running', 1);
% clusterconfig('scheduler', 'cluster');

MRIsubj = dir('/projects/MINDLAB2017_MEG-LearningBach/raw/0*');
countINV2 = 0;
countMoCo = 0;
countMNI152 = 0;


countblock = 0;
for ii = 1:71
    countblock = countblock + 3;
    S = [];
    S.spm_files_recog_basen = spm_files_recog_basen;
    S.countblock = countblock;
    S.workingdir = workingdir;
    S.D = [workingdir '/dff' spm_files_recog_basen{countblock}]; %check this if you want 'epoched' or 'continuous'
    %D = spm_eeg_load(spm_files{ii});
    %IT SHOULD BE FINE, BUT MAYBE CHECK!!
    if ~isempty(dir(['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/NIFTI_INV2/' MRIsubj(ii).name '/*INV2.nii']))
        MRI_Nifti = dir(['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/NIFTI_INV2/' MRIsubj(ii).name '/*INV2.nii']);
        countINV2 = countINV2 + 1;
        S.mri = [MRI_Nifti.folder '/' MRI_Nifti.name];
    elseif ~isempty(dir(['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/NIFTI_INV2/' MRIsubj(ii).name '/MoCo*.nii']))
        MRI_Nifti = dir(['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/NIFTI_INV2/' MRIsubj(ii).name '/MoCo*.nii']);
        countMoCo = countMoCo + 1;
        S.mri = [MRI_Nifti.folder '/' MRI_Nifti.name];
    else
        S.mri = fullfile('/scratch5/MINDLAB2017_MEG-LearningBach/MRI_MNI152/MNI152_T1_2mm.nii');
        countMNI152 = countMNI152 + 1;
    end
    %S.mri = ['/scratch3/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/NIFTI_INV2/' MRIsubj(ii).name '/t1_mp2rage_sag_p2_iso_INV2.nii'];
    %[datadir '/' 't1_mp2rage_sag_p2_iso_UNI_Images.nii'];
    %S.mri = fullfile('/scratch3/MINDLAB2017_MEG-LearningBach/MRI_MNI152/MNI152_T1_2mm.nii'); %averaged MRI of 152 standard people  %('/projects/MINDLAB2017_MEG-LearningBach/scratch/Leonardo/MUMUFE/maxfilter_preproc/MPRAGE_MGH_variant_.nii'); %need the file converted in NIFTI
    S.useheadshape = 1;
    S.use_rhino = 1; %set 1 for having rhino, 0 for not having it
    S.forward_meg = 'Single Shell';
    S.fid.label.nasion = 'Nasion';
    S.fid.label.lpa = 'LPA';
    S.fid.label.rpa = 'RPA';
    input = S;
    % clusterconfig('scheduler', 'none');
    % distribute
    jobid = job2cluster(@coregfunc, input); 
end









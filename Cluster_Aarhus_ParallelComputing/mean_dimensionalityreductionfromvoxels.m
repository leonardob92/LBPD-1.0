%% OSL and varius paths

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA'); %path to LEiDA_MEG_leonardo functions
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Islands_new'); %path to islands and violin functions
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach')

%% settings for cluster

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')

% config of the cluster server
clusterconfig('slot', 1); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue


%%

%old source reconstruction
% recon_freq = 'source_localsphere_Hz053_20TrialsRT.oat';
% recon_frequ = 'source_localsphere_Hz053_20TrialsRT_broadband.oat';
% recon_freq = 'source_localsphere_Hz01_Hz3_20TrialsRT_delta2.oat';
% recon_freq = '/source_localsphere_Hz4_Hz15_20TrialsRT_alpha.oat';
% recon_freq = '/source_localsphere_Hz16_Hz31_20TrialsRT_beta.oat';
% recon_freq = '/source_localsphere_Hz32_Hz74_20TrialsRT_gamma.oat';
% recon_freq = '/source_localsphere_Hz04_Hz15_AllTrials_Alpha.oat';
% recon_freq = '/source_localsphere_Hz04_Hz15_AllTrials_Alpha_PCA90.oat';

%new source reconstruction
% recon_freq = '4_15Hz_longbaseline_testlocal.oat';
% recon_freq = '8_12Hz_longbaseline_testlocal.oat';
% recon_freq = '11_12Hz_longbaseline_testlocal.oat';
% recon_freq = '8_15Hz_longbaseline_testlocal.oat';
% recon_freq = '8_20Hz_longbaseline_testlocal.oat';
% recon_freq = 'source_localsphere_01_1Hz_AllTrials.oat';
% recon_freq = 'source_localsphere_2_8Hz_AllTrials.oat';
% recon_freq = 'source_localsphere_2_8Hz_AllTrials_longbaseline.oat';
% recon_freq = 'source_localsphere_01_2Hz_AllTrials_longbaseline.oat';
% recon_freq = 'source_localsphere_12_32Hz_AllTrials_longbaseline.oat';
recon_freq = 'source_localsphere_01_40Hz_AllTrials_longbaseline.oat';

%averaging across conditions
for ii = 1:71
    if ii ~= 39
%         filepath = ['/scratch3/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/' recon_freq '/concatMefsession' num2str(ii) '_spm_meeg.mat'];
        filepath = ['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/' recon_freq '/concatMefsession' num2str(ii) '_spm_meeg.mat'];
        S = [];
        S.D = filepath;
        S.prefix = 'mfs';
        jobid = job2cluster(@mean_D, S);
    end
end


%% parcellation AAL


%this is for giving input such as parcellation, ROI names and template (otherwise it can also find by itself the proper template, in this case 8mm)
parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
ROIsfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_ROIs.txt';
templatefile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/MNI152_T1_8mm_brain_PROVA.nii.gz';

p = parcellation(parcelfile,ROIsfile,templatefile);
p = p.remove_parcels(91:116);



%% parcellation 38 ROIs (and related files for labels and coordinates)


% p = parcellation('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/parcellations/fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz'); %38-ROI parcellation (the one suggested also in the OSL website)

% data_dir = '/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA';
% [~,~,raw_parc] = xlsread([data_dir '/Parcellation_38ROI_label_coordinates.xlsx']); %loading coordinates, labels and some information about the 38-ROI parcellation
% 
% coord = cell2mat(raw_parc(2:end,5:7)); %extracting coordinates and converting to double
% save('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/parc_38_ROI_coord','coord'); %saving double file with coordinates
% 
% lab = char(raw_parc(2:end,8)); %extracting labels and converting to characters
% save('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/parc_38_ROI_labels','lab'); %saving character file with labels



% p2 = parcellation(parcelfile);

% p2 = parcellation('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_2mm.nii');


% m = nii.load('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try4.nii.gz');
% 
% m3 = nii.load('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz');
% 
% m2 = nii.load('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_1mm.nii');




%% REDUCING DIMENSIONALITY OF VOXELS USING PCA

%old source reconstruction
% recon_freq = 'source_localsphere_Hz053_20TrialsRT_broadband.oat'; %3-40 Hz
% recon_freq = 'source_localsphere_Hz053_20TrialsRT.oat';
% recon_freq = 'source_localsphere_Hz01_Hz3_20TrialsRT_delta2.oat';
% recon_freq = 'source_localsphere_Hz4_Hz15_20TrialsRT_alpha.oat';
% recon_freq = 'source_localsphere_Hz16_Hz31_20TrialsRT_beta.oat';
% recon_freq = 'source_localsphere_Hz32_Hz74_20TrialsRT_gamma.oat';
% recon_freq = '/source_localsphere_Hz04_Hz15_AllTrials_Alpha.oat';
% recon_freq = '/source_localsphere_Hz04_Hz15_AllTrials_Alpha_PCA90.oat';

%new source reconstruction
% recon_freq = '4_15Hz_longbaseline_testlocal.oat';
% recon_freq = '8_12Hz_longbaseline_testlocal.oat';
% recon_freq = '11_12Hz_longbaseline_testlocal.oat';
% recon_freq = '8_15Hz_longbaseline_testlocal.oat';
% recon_freq = '8_20Hz_longbaseline_testlocal.oat';
% recon_freq = 'source_localsphere_01_1Hz_AllTrials.oat';
% recon_freq = 'source_localsphere_2_8Hz_AllTrials.oat';
% recon_freq = 'source_localsphere_2_8Hz_AllTrials_longbaseline.oat';


% list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz32_74_rest10.oat/con*');
% list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz2_Hz8_minorfl.oat/con*');
% list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_12_32Hz_AllTrials_longbaseline.oat/mfs*');
% list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_01_2Hz_AllTrials_longbaseline.oat/mfs*');
% list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_01_40Hz_AllTrials_longbaseline.oat/con*');
% list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz01_40_rest10.oat/con*');
list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/8_12Hz_longbaseline_testlocal.oat/mcon*');

mont = 2; %specify the SPM object montage you are working on
    
%averaging across conditions
for ii = 2:2:length(list)
%     if ii ~= 2 && ii ~= 3 && ii ~= 8
%         filepath = ['/scratch3/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/' recon_freq '/concatMefsession' num2str(ii) '_spm_meeg.mat'];
        S = [];
        S.p = p;
        S.montage = mont;
%         S.filepath = ['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/' recon_freq '/mfsconcatMefsession' num2str(ii) '_spm_meeg.mat'];
        S.filepath = [list(ii).folder '/' list(ii).name];
        jobid = job2cluster(@red_dim, S);
%     end
end




%%

D = D.montage('switch',2);





                
                
                


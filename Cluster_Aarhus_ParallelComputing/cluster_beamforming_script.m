

%% 

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')


%%


datadir = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc'; %with movement compensation
workingdir = datadir;
basefifdatadir = datadir;
basefifdatadirraw = '/raw/sorted/MINDLAB2017_MEG-LearningBach';

tempfif = dir([basefifdatadir '/SUB*']); %try what happen when you have more subject!!
for ii = 1:length(tempfif)
    spm_files_basenames{ii,1} = ['spmeeg_' tempfif(ii).name(1:end-4) '.mat'];
    spm_files{ii,1} = [datadir '/' 'spmeeg_' tempfif(ii).name(1:end-4) '.mat'];
end

%xlsx directory on the server and creation of xlsx_basenmas (procedure mirroring the spm_files)
xlsx_dir = '/scratch5/MINDLAB2017_MEG-LearningBach/BachStimuli';
xlsx_dir_behav = '/scratch5/MINDLAB2017_MEG-LearningBach/BehavioralResponses';

%creates a list of only the recog files (spm objects)
k = 0;
for ii = 1:length(tempfif)
    if strcmp(tempfif(ii).name(9:12),'reco')
        k = k + 1;
        %spm_files_recog{k,1} = spm_files{ii,1};
        spm_files_recog_basen{k,1} = spm_files_basenames{ii,1};
        spm_files_recog{k,1}=[workingdir '/dff' spm_files_recog_basen{k}];
        xlsx_basenames{k,1} = [spm_files_recog_basen{k,1}(12:15) spm_files_recog_basen{k,1}(21:23) '.xlsx'];
    end
end


%%

% config of the cluster server
clusterconfig('slot', 1); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
clusterconfig('scheduler', 'none'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue



oat = [];
%SOURCE RECON
count3 = 0;
for jj = 1:71
    count3 = count3 + 3;
    spm_files{jj}=['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/e80ldff'  spm_files_recog_basen{count3}]; %check this in order not to mix up 'epoched' and 'continuous'
%     D_epoched = spm_eeg_load(spm_files{ii});
%     D_epoched = D_epoched.montage('switch',1); %switch the montage to 1 in order to be safer
%     D_epoched.save();
    processed_file_epoched = [spm_files{jj}];
    spm_files{jj}=['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/dff'  spm_files_recog_basen{count3}]; %check this in order not to mix up 'epoched' and 'continuous'
%     D_continuous = spm_eeg_load(spm_files{ii});
%     D_continuous = D_continuous.montage('switch',1); %switch the montage to 1 in order to be safe
%     D_continuous.save();
    processed_file_continuous = [spm_files{jj}];
%     Beamform
%     oat                                = [];
    oat.source_recon.D_continuous(jj)    = {processed_file_continuous};%the file just after AFRICA and before the epoching
    oat.source_recon.D_epoched(jj)      = {processed_file_epoched};
    oat.source_recon.session_names(jj)      = {['session' num2str(jj)]};
    oat.source_recon.results_fnames(jj)     = {['session' num2str(jj) '_recon']};
    disp(num2str(jj))
end

% pca_dim_1 = 50;
pca_dim_1 = 50;
oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_continuous),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects

% oat.source_recon.D_continuous(ii)    = {processed_file_continuous{ii}};%the file just after AFRICA and before the epoching
% oat.source_recon.D_epoched(ii)      = {processed_file_epoched{ii}};
oat.source_recon.modalities        = {'MEGMAG'; 'MEGPLANAR'}; %'MEGMAG'};
oat.source_recon.conditions        = {'Old_Correct';'New_Correct'};
oat.source_recon.gridstep          = 8; % in mm
% oat.source_recon.time_range        = [-0.1 3.4]; % time range in secs
oat.source_recon.time_range        = [-2 5]; % time range in secs
% oat.source_recon.freq_range        = [3 40]; % frequency range in Hz
% oat.source_recon.freq_range        = [0.1 3]; % frequency range in Hz (delta2)
% oat.source_recon.freq_range        = [0.5 3]; % frequency range in Hz (delta)
% oat.source_recon.freq_range        = [4 15]; % frequency range in Hz (alpha large range)
oat.source_recon.freq_range        = [8 20]; % frequency range in Hz (alpha large range)
% oat.source_recon.freq_range        = [16 31]; % frequency range in Hz (beta)
% oat.source_recon.freq_range        = [32 74]; % frequency range in Hz (gamma); it is not possible higher than 75 because the D.fsample is 150, and D.fsample/2 cannot be higher than oat.source_recon.freq_range(2)
%S.source_recon.pca_order         = 250;
oat.source_recon.type              = 'Scalar';
oat.source_recon.method            = 'beamform';
% oat.source_recon.method            = 'none'; %set this for having sensor level (only) setup
% oat.source_recon.normalise_method  = 'none'; %set this for having sensor level (only) setup
oat.source_recon.normalise_method  = 'mean_eig';
%oat.source_recon.forward_meg       = 'Single Shell'
oat.source_recon.forward_meg       = 'MEG Local Spheres';
%S.source_recon.prefix            = '';
oat.source_recon.report.do_source_variance_maps = 1;
% oat.source_recon.dirname           = [workingdir '/source_localsphere_RHINOCORRECT_Hz053']; % spm_files_recog_basen{ii}(12:25)];
% oat.source_recon.dirname           = [workingdir '/source_localsphere_RHINOCORRECT_Hz053']; % spm_files_recog_basen{ii}(12:25)];
% oat.source_recon.dirname           = [workingdir '/source_localsphere_Hz053_GMMask.oat']; % spm_files_recog_basen{ii}(12:25)];
% oat.source_recon.dirname           = [workingdir '/source_localsphere_Hz053_20TrialsRT']; % spm_files_recog_basen{ii}(12:25)];
% oat.source_recon.dirname           = [workingdir '/source_localsphere_Hz053_20TrialsRT_broadband']; % spm_files_recog_basen{ii}(12:25)];
% oat.source_recon.dirname           = [workingdir '/source_localsphere_Hz4_Hz15_20TrialsRT_alpha'];
% oat.source_recon.dirname           = [workingdir '/source_localsphere_Hz16_Hz31_20TrialsRT_beta'];
% oat.source_recon.dirname           = [workingdir '/source_localsphere_Hz32_Hz74_20TrialsRT_gamma'];
% oat.source_recon.dirname           = [workingdir '/source_localsphere_Hz04_Hz15_AllTrials_Alpha_PCA90'];
% oat.source_recon.dirname           = [workingdir '/source_localsphere_Hz04_Hz15_AllTrials_Alpha_PCA90'];
% oat.source_recon.dirname           = ['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/4_15Hz_longbaseline']; 
oat.source_recon.dirname           = ['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/8_20Hz_longbaseline_cluster']; 








%%


for ii = 1:71 %loop with the number of total cluster to made
    % distribute
    oat.source_recon.sessions_to_do = [];
    oat.source_recon.sessions_to_do = [ii]; %sessions to do among the file_list
    jobid = job2cluster(@cluster_beamforming, oat);% this is the command for send the job to the cluster, in the brackets you can find the name on the function to run (afeter the @) and the variable for the input (in this case input)
    % look the script for more details about the function work
end





%% WHOLE-BRAIN FUNCTIONAL CONNECTIVITY UNDERLYING SOUND ENCODING

%%

%BEFORE PROCEEDING, PLEASE NOTE:

%As follows, every analysis has been reported to ensure full disclosure.
%Please, note that in this Leading script, I use built-in functions,
%functions from well-known toolboxes (e.g. OSL, SPM, FieldTrip) and
%in-house-built functions, which are reported together with this script in the collection named LBPD.
%Data can be provided according to the Danish regulation, upon reasonable request.
%If interested, please contact Leonardo Bonetti, leonardo.bonetti@clin.au.dk
%More information is provided in the ReadMe.mat file that I strongly advise to read.
%This script was used for the following paper:
%https://www.biorxiv.org/content/10.1101/2021.10.14.464389v1


%% START UP FUNCTIONS AND SETTING UP PATHS.. (LBPD_startup_D)

%this function was made to add some paths related to required external functions
%if you want to use/test the package of functions that we provided, you
%need to donwload the functions and then store them in a path that in this
%example script is called "pathl"
%please specify the path in your own computer where you stored the provided data
pathdata = 'C:/Users/calum/Desktop/Codes_Data_NN1/PaperNN1';
%starting up some of the functions that I wrote for the following analysis
% pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
pathl = 'C:/Users/calum/Desktop/Codes_Data_NN1/Leonardo_FunctionsPhD';
addpath(pathl);
LBPD_startup_D(pathl);

%%

%% *** PRE-PROCESSING ***

%% Maxfilter

%OBS!!! before running maxfilter, you need to close matlab, open the terminal, write: 'use anaconda', then open matlab and run maxfilter script
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup
basefifdatadir = '/raw/sorted/MINDLAB2017_MEG-LearningBach'; %0001/20171220_000000/MEG';
basefifdatadirraw = '/raw/sorted/MINDLAB2017_MEG-LearningBach'; %0001/20171220_000000/MEG';
datadir = '/projects/MINDLAB2017_MEG-LearningBach/scratch/Leonardo/LearningBach/maxfilter_preproc/'; %0001';
%datadir =
%'/Users/au550322/Documents/MATLAB/MajorMinor_Mumufe/Leonardo_NielsData';
%previous datadir that now is not useful
workingdir = datadir;
subnum = '1'; %set the sub ID for starting from a specific chosen subject in the temp1fif list
temp1fif = dir([basefifdatadir '/0*']);
%run again because you need to close matlab before for the use anaconda thing
for ii = str2double(subnum):length(temp1fif) %this iterates over subjects (remember to set only the subjects that you really want)
    temp2fif = [basefifdatadirraw '/' temp1fif(ii).name];
    dummy3fif = dir([temp2fif '/20*']);
    dummy2_3fif = dir([temp2fif '/' dummy3fif(1).name '/MEG/0*']);
    if isempty(dummy2_3fif)
        dummy2_3fif = dir([temp2fif '/' dummy3fif(2).name '/MEG/0*']);
    end
    if isempty(dummy2_3fif)
        display('fucking problem in getting the raw MEG data!');
    end
    fif_files{ii,1} = dummy2_3fif;
    for k = 1:length(dummy2_3fif) %this iterates over files for each subject
        temp4fif = [dummy2_3fif(k).folder '/' dummy2_3fif(k).name '/files/' dummy2_3fif(k).name(5:end) '.fif'];
        a{ii,k} = temp4fif; %just a cell for storing the all paths
        S = [];
        S.dataset = temp4fif;
        S.outfile = ['spmeeg_SUBJ' temp1fif(ii).name dummy2_3fif(k).name(5:end)];
%       D = spm_eeg_convert(S);
    end
end
%new proper lines for maxfilter
maxfilter_path = '/neuro/bin/util/maxfilter';
maxDir = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach_Fsample1000Hz'; %with movement compensation
project = 'MINDLAB2017_MEG-LearningBach';
%OBS!!! YOU NEED TO HAVE 'SUBJ0009' AND THEN 'SUBJ0010'.. THINK HOW TO DO THAT!!
for ii = str2double(subnum):length(temp1fif) %iterates over subjects
    fif_files_dummy = fif_files{ii,1};
    spm_files_dummy = fif_files{ii,1};
    for j = 1:length(spm_files_dummy) %iterates over experimental blocks
        if ~isempty(strfind(spm_files_dummy(j).name,'rest10')) || ~isempty(strfind(spm_files_dummy(j).name,'learminor')) || ~isempty(strfind(spm_files_dummy(j).name,'recogminor'))  %here for now I want only learminor, recogminor and rest10
            rawName{j} = [fif_files_dummy(j).folder '/' fif_files_dummy(j).name '/files' '/' fif_files_dummy(j).name(5:end) '.fif'];
            [~,n,~] = fileparts(rawName{j});
            maxfName = ['SUBJ00' num2str(ii) rawName{j}(regexp(rawName{j},'files/')+6:regexp(rawName{j},'.fif')-1)];
            badchans = [];
            cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ',num2str(project),' "',maxfilter_path,' -f ',[fif_files_dummy(j).folder '/' fif_files_dummy(j).name '/files' '/' fif_files_dummy(j).name(5:end) '.fif'],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 3 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
            system(cmd);
        end
    end
end

%% setting directories

clear
datadir = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc'; %with movement compensation in maxfilter
workingdir = datadir;
basefifdatadir = datadir;

%% 

%please note that these preprocessing steps have been run for a large
%number of data files. Some of them were not used in this specific paper
%but in other projects, however these analysese have been run together.

%% creating some paths (not the best codes ever.. but they work fine)

%convert fif file into SPM object
temp1fif = dir([basefifdatadir '/*.fif']);
for ii = 1:length(temp1fif) %this iterates over blocks through all subjects (remember to set only the subjects that you really want)
    temp2fif = [basefifdatadir '/' temp1fif(ii).name];
    fif_files{ii,1} = temp2fif;
    S = [];
    S.dataset = fif_files{ii,1};
%     D = spm_eeg_convert(S); 
end
%create the spm object files path and name
for ii = 1:length(temp1fif)
    spm_files_basenames{ii,1} = ['spmeeg_' temp1fif(ii).name(1:end-4) '.mat'];
    spm_files{ii,1} = [datadir '/' 'spmeeg_' temp1fif(ii).name(1:end-4) '.mat'];
end
%creates a list of only the recog files (spm objects)
k = 0;
jj = 0;
clear spm_files_lear_basen
for ii = 1:length(temp1fif)
    if strcmp(temp1fif(ii).name(9:12),'lear')
        k = k + 1;
        %spm_files_recog{k,1} = spm_files{ii,1};
        spm_files_lear_basen{k,1} = spm_files_basenames{ii,1}; %this is not computationally elegant and the space for these variable should be preallocated but it does not really change anything for my purposes..
        spm_files_lear{k,1}=[workingdir '/dff' spm_files_lear_basen{k}];
        xlsx_basenames{k,1} = [spm_files_lear_basen{k,1}(12:15) spm_files_lear_basen{k,1}(21:23) '.xlsx'];
    end
    if strcmp(temp1fif(ii).name(9:12),'rest')
        jj = jj + 1;
        spm_files_rest10_basen{jj,1} = spm_files_basenames{ii,1};
        spm_files_rest10{jj,1}=[workingdir '/dff' spm_files_rest10_basen{jj}];
    end
end

%% highpass-filtering
%data will get a 'f' before the name
for ii = 1:length(spm_files), % iterates over blocks !OBS!
    spm_files{ii} = [workingdir '/' spm_files_basenames{ii}];
    S2 = [];
    S2.D = spm_files{ii};
    S2.band = 'high';
    S2.freq = 0.1;
    D = spm_eeg_filter(S2);
end
%notch filter
%data will get a 'ff' before the name
for ii = 1:length(spm_files)
    spm_files{ii} = [workingdir '/f' spm_files_basenames{ii}];
    S2 = [];
    S2.D = [spm_files{ii}];
    S2.band = 'stop';
    S2.freq = [48 52]; % This defines the notch filter frequency range, i.e. around 50Hz. Because this is the frequency of the electric current or something like that
    D = spm_eeg_filter(S2);
end
%donwsampling
for ii = 1:length(spm_files) % iterates over experimental blocks %OBS!!!!!
    spm_files{ii} = [workingdir '/ff' spm_files_basenames{ii}];
    S = [];
    S.D = spm_files{ii};
    S.fsample_new = 150; % in Hz
    D = spm_eeg_downsample(S);    
end
%removing bad trials (if any) from the data by performing a visual
%inspection
for ii = 1:length(spm_fils)
    spm_files{ii}=[workingdir '/dff' spm_files_basenames{ii}];
    D = spm_eeg_load(spm_files{ii});
    D = oslview(D);
    D.save(); %it is needed to save the marked bad events (and/or channels) in oslview
end
% removing of eyeblink and heartbeat artefacts by computing ICA
for ii = 1:length(spm_files)
    spm_files{ii}=[workingdir '/dff' spm_files_basenames{ii}];
    D = spm_eeg_load(spm_files{ii});
    D = osl_africa(D,'do_ident',false,'do_remove',false,'used_maxfilter',true); %qui direi che decompone in independent components
    D = osl_africa(D,'do_remove',false); %questo apre la GUI e permette di fare la visual inspection
    D = osl_africa(D,'do_ident',false,'do_remove',true); %questo direi che serve per rimuovere visivamente le componenti
    D.save();
    disp(['just done subject number ' num2str(ii)])
end

%% epoching - (baseline = (-)100ms)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Francesco/Functions');
prefix_tobeadded = '2e2420'; %prefix to be added
note = ('/scratch5/MINDLAB2017_MEG-LearningBach/BachStimuli/learminonset.xlsx'); %load xlsx file with the onset of the notes 
ons = xlsread(note,'c:c'); %take all the value in the C row with the  onset of the notes
%creates a list of only the recog files (subgroup of the whole set of data files)
k = 0;
clear spm_files_recog
for ii = 1:length(temp1fif)
    if strcmp(temp1fif(ii).name(9:12),'lear')
        k = k + 1;
        %spm_files_recog{k,1} = spm_files{ii,1};
        spm_files_lear_basen{k,1} = spm_files_basenames{ii,1};
        spm_files_lear{k,1}=[workingdir '/dff' spm_files_lear_basen{k}];
        xlsx_basenames{k,1} = [spm_files_lear_basen{k,1}(12:15) spm_files_lear_basen{k,1}(21:23) '.xlsx'];
    end
end
%iterates for the recognition files only
for ii = 3:length(spm_files_lear_basen) %indexing only the files wanted for this paper
    %load the spm object
    D = spm_eeg_load(spm_files_lear{ii,1});
    events = D.events;
    %takes the correct triggers sent during the recording
    count_evval = 0;
    for ieve = 1:length(events)
        if strcmp(events(ieve).type,'STI101_up')
            if events(ieve).value == 5
                count_evval = count_evval + 1;
                trigcor(count_evval,1) = events(ieve).time + 0.010 - 0.05; %this takes the correct triggers and add 10ms of delay of the sound travelling into the tubes and removes 50ms (since the onsets of the notes were a little bit imprecise)
            end
        end
    end   
    trl_sam = zeros(length(trigcor),3);
    trl_sec = zeros(length(ons)*length(trigcor),3);
    rr = 0;
    for pp = 1:length(trigcor)
        for kk = 1:length(ons)
            rr = rr + 1;
            trl_sec(rr,1) = trigcor(pp,1) + ons(kk,1) - 0.1; % add the notes onset timing to the trigcore timing and then remove 0.1sec for the baseline
            trl_sec(rr,2) = trigcor(pp,1) + ons(kk,1) + 0.250; % add to the trigcor timing the total duration of a single note
            trl_sec(rr,3) = trl_sec(rr,2) - trl_sec(rr,1); %range time-windows in s
            trl_sam(rr,1) = floor(trl_sec(rr,1) * 150); %beginning time-window epoch in samples
            trl_sam(rr,2) = trl_sam(rr,1) + 53; %end time-window epoch in samples
            trl_sam(rr,3) = -15; %sample before the onset of the stimulus (corresponds to 0.100ms)
        end
    end
    input = [];
    %creates the epochinfo structure that is required for the source reconstruction later
    epochinfo.trl = trl_sam;
    epochinfo.time_continuous = D.time;
    %switch the montage to 0 because for some reason OSL people prefer to do the epoching with the not denoised data
    D = D.montage('switch',0);
    %store the epochinfo structure inside the D object
    D.epochinfo = epochinfo;
    D.save();
    %build structure for spm_eeg_epochs
    input.D = D;
    input.trl = trl_sam;
    input.prefix = prefix_tobeadded;
    
    jobid = job2cluster(@epoching_cluster, input);
end

%% averaging and combining planar gradiometers (used for power spectra analysis)

%averaging over trials
input_prefix_to_be_set_avg = '2e2420'; %setting a part of input name
output_prefix_to_be_set = 'm'; %setting prefix for the output
for ii = 3:3:length(spm_files_lear) %over files
    S = [];
    S.D = [workingdir '/' input_prefix_to_be_set_avg spm_files_lear{ii}];
    S.prefix = output_prefix_to_be_set;
    D = spm_eeg_average(S);
end
%combining planar gradiometers
input_prefix_to_be_set_cmb = 'm2e2420';
for ii = 3:3:length(spm_files_lear)
    D = spm_eeg_load([workingdir '/' input_prefix_to_be_set_cmb spm_files_lear{ii}]);
    D = D.montage('switch',1);
    D.save();
    S = [];
    S.D = [workingdir '/' input_prefix_to_be_set_cmb spm_files_lear{ii}];
    D = spm_eeg_combineplanar(S);
end

%% *** ANALYSIS ***

%% N100 and frequency analysis

%% LBPD_startup_D

pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);

%% extracting MEG sensor data and plotting topoplot

%%% THIS SECTION WAS USED BOTH TO EXTRACT MEG SENSOR DATA FROM SPM OBJECTS AND TO PLOT TOPOPLOTS %%%

load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/BrainActivity/grad_clust.mat');
load_data = 1; %set 1 if you want to load the data instead of extracting it from SPM objects
subjnum = {'01';'02';'03';'04';'05';'06';'07';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'38';'39';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';'51';'52';'53';'54';'55';'56';'57';'58';'59';'60';'61';'62';'63';'64';'65';'66';'67';'68';'69';'70';'71'}; %Everybody!!!!!!
datadir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc'; %path to data
outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/BrainActivity'; %path to write output in
S = [];
%computing data
S.outdir = outdir;
S.data = [];
if load_data == 1 %if you already computed and saved on disk the t-tests you can load them here
    load([outdir '/sensor_data.mat']);
    S.data = data_mat;
    S.chanlabels = chanlabels;
    S.time_real = time_sel;
else %otherwise you can extract the data from SPM MEEG objects (one for each subject)
    S.spm_list = cell(1,length(subjnum));
    for ii = 1:length(subjnum)        
        S.spm_list(ii) = {[datadir '/Pm2e2420dffspmeeg_SUBJ00' subjnum{ii} 'learminor_tsssdsm.mat']};
    end
end
S.conditions = {'Undefined'};
S.timeextract = [5:54]; %time-points to be extracted (from -0.07 seconds before the onset of the stimulus since if we take more pre-stimulus time we would get the N100 of the previous note..)
S.centerdata0 = 1; %1 to make data starting at 0
S.save_data = 1; %only meaningfull if you read data from SPM objects saved on disk
S.save_name_data = 'sensor_data';
%individual waveform plotting
S.waveform_singlechannels_label = 0; %1 to plot single channel waveforms
S.wave_plot_conditions_together = 0; %1 for plotting the average of all
S.mag_lab = 2; %1 for magnetometers; 2 for gradiometers
S.x_lim_temp_wave = []; %limits for time (in secs) (E.g. [-0.1 3.4])
S.y_lim_ampl_wave = []; %limit for amplitude (E.g. [0 120] magnetometes, [0 6] gradiometers)
%averaged waveform plotting
S.waveform_average_label = 1;
% S.left_mag = [9 11 13 15 19 21 23 25 31 117 119 133]; %misleading name, but it actually refers to whatever group of channels you want to plot
if S.waveform_average_label == 1
    clustplot = GRAD_clust{1,3}; %slightly elaborated way to get the channel original IDs of channels forming the significant cluster that you are considering
    S.left_mag = [];
    for ii = 1:size(clustplot,1)
        S.left_mag(ii) = find(cellfun(@isempty,strfind(chanlabels,clustplot{ii,1})) == 0);
    end
end
% S.signtp(1) = {20:50}; S.signtp(2) = {70:90}; S.signtp(3) = {140:180};
S.signtp = {[0.053 0.167]}; %21 corresponds to minimum (1) + 20 that is the starting time-point for calculating MCS (I started at 20 (corresponding to 0.050sec) since it did not make sense to look for N100 earlier.. 45 reflects the same concept
S.legc = 0;
S.sr = 150; %sampling rate (Hz)
S.avewave_contrast = 1; %1 to plot the contrast between conditions (averaged waveform)
S.save_label_waveaverage = 0;
S.label_plot = 'block_minor';
%t-tests
S.t_test_for_permutations = 0;
S.cond_ttests_tobeplotted_topoplot = [1]; %this is for both topoplot and t-tests!! (here [1 2] means cond1 vs cond2!!!!!!!)
%topoplotting
S.topoplot_label = 1;
S.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External';
S.topocontr = 0;
S.topocondsing = [];
% S.xlim = [0.053 0.167]; %time topolot (reported in the paper)
S.xlim = [0.080 0.150]; %time topolot
S.zlimmag = []; %amplitude topoplot mag
S.zlimgrad = []; %amplitude topoplot grad
S.colormap_spec = 'jet';
S.topoplot_save_label = 0;

[out] = MEG_sensors_plotting_ttest_LBPD_D2(S);

%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))

% if load_data == 0 %this is done because one subject (number 71) was computed with a different epoch and right now was faster to do this work around.. note also that subject 39 has been excluded because of practical problems during the experiment
%     D = spm_eeg_load([datadir '/Pme80dffspmeeg_SUBJ0071recogminor_tsssdsm.mat']);
%     hg = D(14:217,1:526,[2,3]);
%     data_mat2 = zeros(204,526,70,2);
%     data_mat2(:,:,1:69,:) = data_mat;
%     data_mat2(:,:,70,:) = hg;
%     data_mat = data_mat2;
%     save([outdir '/' S.save_name_data '.mat'],'data_mat','time_sel','chanlabels');
% end
% 
% save([outdir '/sensor_data.mat'],'data_mat','time_sel','chanlabels');

%% frequency analysis

%here you need to load the provided "sensor_data.mat" that you stored in your own directory
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/BrainActivity/sensor_data.mat');
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/BrainstormDemo_CodesTimeFrequencyDecomposition/functions') %to Dimitrios Morlet transform function
%frequency decomposition by using complex Morlet wavelet
mag = 1:204;
f = 1:1:60;
PI = zeros(length(f),size(data_mat,2),size(data_mat,3),size(data_mat,4));
for ii = 1:size(data_mat,3) %over subjects
    for jj = 1:size(data_mat,4) %over condition
        data2 = data_mat(mag,:,ii,jj);
        P = morlet_transform(data2,time_sel,f); %frequency decomposition for each channel, subject and condition
        P = squeeze(mean(P,1));
        PI(:,:,ii,jj) = P;        
    end
    disp(ii)
end
%power spectra t-tests (old vs new)
%calculating t-tests for each time-point and each frequency band (in the range 1 - 60 Hz)
P = zeros(size(PI,1),size(PI,2));
T = zeros(size(PI,1),size(PI,2));
% x = squeeze(mean(mean(PItask,1),2));
df = squeeze(mean(mean(PI(:,1:find(time_sel==0),:),2),1)); %getting mean value for power spectra matrix, for each subject
for ii = 1:size(PI,1)
%     df = squeeze(mean(PI(ii,1:find(time_sel==0),:),2)); %getting mean value for power spectra matrix, for each subject
    for jj = 1:size(PI,2)
        [~,p,~,stats] = ttest(squeeze(PI(ii,jj,:)),df);
        P(ii,jj) = p;
        T(ii,jj) = stats.tstat;
    end
    disp(ii)
end
%plotting
ges = 0;
Pold = mean(mean(PI,4),3);
%plotting
figure
imagesc(time_sel(1:end),f,Pold)
set(gca,'YDir','normal') %plotting in descending order
xlabel('time (s)'); ylabel('f (Hz)');
colorbar
if ges ~= 0
    caxis([0 ges])
end
set(gcf,'color','w')
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))
%%plotting contrast (this figure has not been reported in the paper but if you like you can produce also this image)
% T2 = T;
% % T2(T<0) = 0;
% T2(P>1.0e-18) = 0;
% %plotting
% figure
% imagesc(time_sel(find(time_sel==0):end),f,T2)
% set(gca,'YDir','normal') %plotting in descending order
% xlabel('time (s)'); ylabel('f (Hz)');
% colorbar
%testing contrast results (correcting for multiple comparison) by performing Monte Carlo simulations
P2 = P;
P2(T<0) = 1; %only positive contrast (cond1 > cond2)
thresh = 1.0e-18;
permut = 1000;
threshMC = 0.001;
perm_max = 1;
t1 = f;
t2 = time_sel(20:42); %from 0.05 until 0.2 seconds (time-range when it makes sense to look for N100)

[ OUT ] = twoD_MCS_LBPD_D( P2(:,20:42), thresh, permut, threshMC, perm_max, t1 , t2 )

%% MCS ERF (mainly making sure that N100 looks fine)

%% loading t-test results and reshaping them

%time-points to be selected
%gradiometers
min_time_point = 20; %since this is the window where it makes sense to look for N100
max_time_point = 42;

clear DATAP2 TSTAT2
outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/BrainActivity'; %path where t-test results are stored in (HERE YOU NEED TO SPECIFY YOUR OWN DIRECTORY!!)
load([outdir '/sensor_data_OUT_Undefined_vs_baseline.mat']);
DATAP2(:,:,1) = OUT.TSTATP_mag(:,min_time_point:max_time_point);
DATAP2(:,:,2) = OUT.TSTATP_grad(:,min_time_point:max_time_point);
TSTAT2(:,:,1) = OUT.TSTAT_mag(:,min_time_point:max_time_point);
TSTAT2(:,:,2) = OUT.TSTAT_grad(:,min_time_point:max_time_point);
chanlab = chanlabels(1:2:204)';
label = zeros(length(chanlab),1);
for ii = 1:length(chanlab)
    label(ii,1) = str2double(chanlab{ii,1}(4:end));
end
%binarising and reshaping data
p_thresh = 1.0e-16; %for binarising p-values matrices..
%individuating positive vs negative t-values
P = DATAP2; %trick for looking only into the positive t-values
P(P < p_thresh) = 1; %binarizing p-values according to threshold p_thresh
P(P < 1) = 0;
%old (pos) > new (neg)
P(TSTAT2 < 0) = 0; %deleting p-values for negative contrasts
TSTAT_mag_pos = P(:,:,1);
TSTAT_grad = P(:,:,2);
%negative vs positive p-values
P = DATAP2;
P(P < p_thresh) = 1; %binarizing p-values according to threshold 0.05
P(P < 1) = 0;
P(TSTAT2 > 0) = 0; %deleting p-values for positive contrasts
TSTAT_mag_neg = P(:,:,1);
%load a decent 2D approximation in a matrix of the MEG channels location
[~,~,raw_channels] = xlsread('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MatrixMEGChannelLayout_With0_2.xlsx');
% reshaping data for computational purposes
S = [];
S.label = label;
S.TSTAT_mag_pos = TSTAT_mag_pos;
S.TSTAT_mag_neg = TSTAT_mag_neg;
S.TSTAT_grad = TSTAT_grad;
S.raw_channels = raw_channels;

[MAG_data_pos, MAG_data_neg, GRAD_data] = MEG_sensors_MCS_reshapingdata_LBPD_D(S);

%actual Monte Carlo simulations
S = [];
%actual gradiometers data
S.data(:,:,:,1) = GRAD_data;
%zeros.. if you do not want to run the function for magnetometers
S.data(:,:,:,2) = zeros(size(MAG_data_pos,1),size(MAG_data_pos,2),size(MAG_data_pos,3));
S.data(:,:,:,3) = zeros(size(MAG_data_pos,1),size(MAG_data_pos,2),size(MAG_data_pos,3));
S.sensortype = [];
S.MEGlayout = cell2mat(raw_channels);
S.permut = 1000;
S.clustmax = 1;
S.permthresh = 0.001;

[MAG_clust_pos, MAG_clust_neg, GRAD_clust] = MEG_sensors_MonteCarlosim_LBPD_D(S);

%% extracting first and last significant time-point for each channel, converting them into seconds (from time-samples) and printing it as xlsx file 

clustnum = 1; %set the cluster number

%actual computation
hh = GRAD_clust{2,3};
time_sel2 = time_sel(1,min_time_point:max_time_point);
clear ff2
ff2(:,1) = hh(:,1); %extracting channel names
for ii = 1:size(ff2,1)
    ff2(ii,2) = {time_sel2(hh{ii,2}(1))}; %first significant time-point
    ff2(ii,3) = {time_sel2(hh{ii,2}(end))}; %last significant time-point
end
PDn = cell2table(ff2); %converting cell to table
% writetable(PDn,['Gradpos_clust2_tonevsbase.xlsx'],'Sheet',1) %saving xlsx file

%% plotting significant clusters (not reported in the paper)
% 
% S = [];
% S.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External';
% % S.PPall{1} = {[]};
% S.PPall{2} = {[]};
% S.PPall{3} = {[]};
% S.PPall{1} = GRAD_clust;
% % S.PPall{2} = MAG_clust_pos;
% % S.PPall{3} = MAG_clust_neg;
% S.clustnum_grad = 1;
% S.clustnum_magp = 1;
% S.clustnum_magn = 1;
% S.zlim = [];
% S.time = time_sel(1,min_time_point:max_time_point);
% 
% MEG_sensors_MCS_plottingclusters_LBPD_D(S)

%% SOURCE RECONSTRUCTION FOR GLM

%% OSL path

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA'); %path to LEiDA_MEG_leonardo functions
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Islands_new'); %path to islands and violin functions
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Islands_new/schemaball-master'); %path to schemaball functions
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/HMM'); %path to HMM l)eonardo
addpath(genpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/HMM/HMM-MAR')); %path to general HMM functions (%%% WRITE git pull I PRESUME IN THE DIRECTORY OF HMM-MAR OR SOMETHING LIKE THAT FOR GETTING AN UPDATE VRSION OF IT.. GET MORE INFORMATION ABOUT THAT %%%)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/altmany-export_fig-412662f'); %trying to get higher quality images
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Francesco/Functions'); %francesco's functions

%% CALCULATING BEAMFORMING FOR GLM (ACTIVATION AND NOT CONNECTIVITY)

input_prefix_to_be_set_cmb = '2e2420';
%getting already calculated rhino data (head model)
for ii = 12:3:213
    D = spm_eeg_load([workingdir '/dff' spm_files_lear_basen{ii}]);
    invh = D.inv;
    D2 = spm_eeg_load([workingdir '/' input_prefix_to_be_set_cmb 'dff' spm_files_lear_basen{ii}]);
    D2.inv = invh;
    D2.save();
    rhino_display(D2)
    disp(ii)
end

%% setting path

% workingdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Francesco';
% workingdir2 = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc';
% lista = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/e4l*');
%%% VERY IMPORTANT!! THESE EPOCHED FILES ARE THE ONES WITHOUT THE CORRECTION
%%% OF SUBTRACTION OF 50MS.. THEREFORE THAT CORRECTION IS DONE LATER..%%%
%SOURCE RECON
input_prefix_to_be_set_cmb = '2e2420';
oat = [];
cnt = 0;
for ii = 3:3:210 %:length(spm_files_recog_basen)
    cnt = cnt + 1;
    spm_files=[workingdir '/' input_prefix_to_be_set_cmb 'dff' spm_files_lear_basen{ii}];%check this in order not to mix up 'epoched' and 'continuous'
    D_epoched = spm_eeg_load(spm_files);
    D_epoched = D_epoched.montage('switch',1); %switch the montage to 1 in order to be safer
    D_epoched.save();
    processed_file_epoched{cnt} = spm_files;
    spm_files=[workingdir '/dff' spm_files_lear_basen{ii}]; %check this in order not to mix up 'epoched' and 'continuous'
    D_continuous = spm_eeg_load(spm_files);
    D_continuous = D_continuous.montage('switch',1); %switch the montage to 1 in order to be safe
    D_continuous.save();
    processed_file_continuous{cnt} = spm_files;
    oat.source_recon.session_names(cnt)      = {['session' num2str(cnt) '_recon.mat']};
    oat.source_recon.results_fnames(cnt)     = {['session' num2str(cnt) '_recon.mat']};
    disp(num2str(ii))
end
oat.source_recon.D_continuous    = processed_file_continuous; %the file just after AFRICA and before the epoching
oat.source_recon.D_epoched      = processed_file_epoched;
pca_dim_1 = 50;
oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_continuous),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects
oat.source_recon.sessions_to_do = [];
oat.source_recon.sessions_to_do = [9]; %sessions to do among the file_list
oat.source_recon.modalities        = {'MEGMAG'; 'MEGPLANAR'}; %'MEGMAG'};
oat.source_recon.conditions        = {'Undefined'};
oat.source_recon.gridstep          = 8; % in mm
oat.source_recon.time_range        = [-0.1 0.25]; % time range in secs
oat.source_recon.freq_range        = [2 8]; % frequency range in Hz (delta)
oat.source_recon.type              = 'Scalar';
oat.source_recon.method            = 'beamform';
oat.source_recon.normalise_method  = 'mean_eig';
oat.source_recon.forward_meg       = 'MEG Local Spheres';
oat.source_recon.report.do_source_variance_maps = 1;
oat.source_recon.dirname           = ['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_Hz2_Hz8_minorfl_forGLM'];

%% settings for cluster (parallel computing)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
% config of the cluster server
clusterconfig('slot', 1); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%% first level (each experimental block for each subject, independently)

for cnt = 1:71
    oat.source_recon.results_fnames(cnt)     = {['session' num2str(cnt) '_recon.mat']};
end
%FIRST LEVEL
design_matrix_summary = {};
% design_matrix_summary{1} = [1 0];design_matrix_summary{2} = [0 1]; %design_matrix_summary{3}=[0 0 1 0];design_matrix_summary{4}=[0 0 0 1];
design_matrix_summary{1} = [1]; %design_matrix_summary{3}=[0 0 1 0];design_matrix_summary{4}=[0 0 0 1];
oat.first_level.design_matrix_summary = design_matrix_summary;
% contrasts to be calculated:
oat.first_level.contrast = {};
% contrast design matrix
oat.first_level.contrast{1} = [1]'; %tone 
oat.first_level.contrast_name = {};
oat.first_level.contrast_name{1} = 'tone';
%original contrasts, try1 and try2 seem to give the exact same results
oat.first_level.report.first_level_cons_to_do = [1]; %[1 2 3]; %better to do 3 2 1 in order to get the information for the peak value of the contrast old-new
oat.first_level.time_range = [-0.1 0.25];
oat.first_level.post_tf_downsample_factor = 1;
oat.first_level.cope_type = 'coape'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
oat.first_level.name = ['wholebrain_first_level']; %REMEMBER TO CHECK THIS NAME!!
oat.first_level.bc = [0];
%to add if the oat has not been automatically saved
for ii = 1:71
    oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level.mat']};
end
oat.first_level.sessions_to_do = [1:2 4:70]; %this is probably not necessary but it may be useful to report if you run only subject_level analysis in oat..

%running first level on parallel computing
for ii = 2:70 %here remember that from subject 8 subject IDs are wrong (basically original subject ID 8 was missing) here it is not a problem since I do not do any group analysis (with GLM), at least for now: when I do group analysis with connectivity this problem is solved
    if ii ~= 3 %if subject is not 3
        oat.first_level.sessions_to_do = [];
        oat.first_level.sessions_to_do = [ii];
        jobid = job2cluster(@cluster_beamfirstlevel,oat);
    end
end

%% subject level (it simply puts together different experimental blocks for the same subject)

%SUBJECT LEVEL
%in this case it does not do anything since I have only one experimental
%block for each participant (however, for computational reasons, I have to
%run it)
% this is if you have a perfect correspondance between sessions and subjects
oat.subject_level.session_index_list = cell(1,71); %sarebbe length(subjects)
oat.subject_level.name = 'minor';
oat.subject_level.subjects_to_do = [];
oat.subject_level.subjects_to_do = [1:2 4:70];
% oat.subject_level.subjects_to_do = [3:3:114, 120:3:210]; %It might be necessary also for the group level analysis (to have the same subjects here that you also have in the group level analysis
%to update the names of the results files of subject level
for ii = 1:70
    oat.subject_level.session_index_list{ii} = ii;
    oat.subject_level.results_fnames(ii) = {['subject' num2str(ii) '_wholebrain_first_level_minor.mat']};
end

oat = osl_check_oat(oat);
oat.to_do = [0 0 1 0];
oat = osl_run_oat(oat);

%% GROUP LEVEL
oat.group_level = [];
oat.group_level.name = 'group_level_everybody_AAL'; %OBS!! REMEMBER TO UPDATE THE NAME!
oat.group_level.subjects_to_do = [];
oat.group_level.subjects_to_do = [1:2 4:69]; %this seems to be necessary since the function gets the indices for the group level from each progressive number of oat.subject_level.subjects_to_do (not very clever idea by the way..). Therefore not to mix up things, you should have here a vector from 1 to the total number of subjects and then specify the correct subject IDs that you want.. or at least, to me it looks like that..
%results name
oat.group_level.results_fnames = ['wholebrain_first_level_' oat.subject_level.name '_' oat.group_level.name '.mat'];
% Spatial and temporal averaging options
oat.group_level.time_range = [-0.1 0.24];
oat.group_level.space_average = 0;
oat.group_level.time_average = 0;
oat.group_level.time_smooth_std = 0; % secs
oat.group_level.use_tstat = 0;
%path to AAL template
oat.group_level.mask_fname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
% Spatial and temporal smoothing options
oat.group_level.spatial_smooth_fwhm = 0; % mm
oat.group_level.group_varcope_time_smooth_std = 100;
oat.group_level.group_varcope_spatial_smooth_fwhm = 100; % smooths the variance of the group copes. It is recommended to do this.
%store copes (useful for doing the permutation test later, otherwise it needs to compute again the group level analysis)
oat.group_level.store_lower_level_copes = 1;
% Set up design matrix and contrasts
%this is if you have only the general mean across the all participants
oat.group_level.group_design_matrix = ones(1,length(oat.group_level.subjects_to_do)); %if you want all of the participants
oat.group_level.group_contrast = [];
oat.group_level.group_contrast{1} = [1];
oat.group_level.group_contrast_name = {};
oat.group_level.group_contrast_name{1} = 'mean';
% Define which contrasts to perform for the report
oat.group_level.first_level_contrasts_to_do = [1]; % list of first level contrasts to run the group analysis on
oat.group_level.report.first_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.show_lower_level_copes = 0;
oat.group_level.report.show_lower_level_cope_maps = 0;

clust = 0; %parallel computing or local computation

if clust == 1
    jobid = job2cluster(@cluster_beamgrouplevel,oat);
else
    oat = osl_check_oat(oat);
    oat.to_do = [0 0 0 1];
    oat = osl_run_oat(oat);
end

%% producing nifti images of statistics

%group level
S2 = [];
S2.oat = oat;
S2.stats_fname = oat.group_level.results_fnames;
S2.first_level_contrasts = [1]; %remember that in this way you define the order of the output (this numbers refers to the order (numbers) defined in the contrasts; for example here tstat3 refers to the contrast old-new)
S2.group_level_contrasts = [1];
S2.resamp_gridstep = oat.source_recon.gridstep;
[statsdir,times,count] = oat_save_nii_stats(S2);

%% cluster-based permutation tests to get significant clusters of MEG sources

S = [];
S.oat = oat;
% S.cluster_stats_thresh = 0.7;
S.cluster_stats_thresh = 1.7;
S.cluster_stats_nperms = 5000; % we normally recommend doing 5000 perms
S.first_level_copes_to_do = [1];
S.group_level_copes_to_do = [1];
S.group_varcope_spatial_smooth_fwhm = S.oat.group_level.group_varcope_spatial_smooth_fwhm;
S.write_cluster_script = 0;
S.time_range = [0.05 0.2];
S.time_average = 1;
% Run the permutations
[ gstats ] = oat_cluster_permutation_testing(S);

%The outputted results were then used as mask for plotting only the
%significant voxels and their corresponding statistics. This operation was
%done by using fslmath.

%We used the following example line in the terminal:
% fslmaths tstat3_gc1_8mm.nii.gz.nii.gz -mas clustere_corrp_tstat3_gc1_8mm.nii.gz contr.nii.gz

%After that, we used Workbench for producing the
%images that we presented in the paper.
%In the following section we report the images that we obtained and used to
%produce a detailed statistical output concerning the significant voxels.
%Additional clarification about these algorithms can be found in the paper.
%Furthermore, you are very welcome to contatct us if you have specific
%questions.

%The resulting nifti image is "WB_N100.nii.gz" that has been then used in
%the subsequent section of this script to extract the statistics of each
%significant voxel

%% creating a file with information on significant clusters/voxels outputted by the permutation test (THINK TO MAKE THIS A FUNCTION)

fnn = 'WB_N100';
%path to file nifti
pathnii = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/BrainActivity'; %THIS MUST BE YOUR OWN DIRECTORY
%actual name plus path
fname = [pathnii '/' fnn '.nii.gz'];
%getting MNI coordinates of significant voxels within the provided image
[ mni_coords, xform ] = osl_mnimask2mnicoords(fname);
%loading the image
V = nii.load(fname);
%extracting statistics
VV = V(V~=0);
%indices of non-zero values of nifti image
VI = find(V~=0);
%path to AAL template
parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
%loading AAL labels
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat');
%extracting AAL coordinates information
K = nii.load(parcelfile);
%finding AAL non-zero coordinates values
% KI = find(K~=0);
%sorting results in order to have strongest voxels at the top (positive
%t-values) or at the bottom (negative t-values)
[VV2, II] = sort(VV,'descend');
VI = VI(II);
mni_coords = mni_coords(II,:);
%final cell
PD = cell(length(VV2),4);
%getting AAL indices
ROI = zeros(length(VI),1);
cnt = 0;
for ii = 1:length(VI)
    ROI(ii) = K(VI(ii));
    if ROI(ii) > 0 && ROI(ii) < 91
        cnt = cnt + 1;
        PD(cnt,1) = {lab(ROI(ii),3:end)}; %storing ROI
        PD(cnt,4) = {mni_coords(ii,:)}; %storing MNI coordinates
        if mni_coords(ii,1) > 0 %storing hemisphere
            PD(cnt,2) = {'R'};
        else
            PD(cnt,2) = {'L'};
        end
        PD(cnt,3) = {round(VV2(ii),2)}; %storing t-statistics
    end
end
PDn = cell2table(PD(~any(cellfun('isempty',PD),2),:)); %remove the possible empty cell
% writetable(PDn,[fnn '.xlsx'],'Sheet',1)

%% BEAMFORMING MUSIC LISTENING (FOR CONNECTIVITY ANALYSIS)

%setting path
workingdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Francesco';
workingdir2 = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc';
lista = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/e4l*');
%%% VERY IMPORTANT!! THESE EPOCHED FILES ARE THE ONES WITHOUT THE CORRECION
%%% OF SUBTRACTION OF 50MS.. THEREFORE THAT CORRECTION IS DONE LATER..%%%
pari = 0;
%SOURCE RECON
oat = [];
for ii = 1:length(lista)/2 %:length(spm_files_recog_basen)
    pari = pari + 2;
    spm_files=[workingdir2 '/' lista(pari).name];%check this in order not to mix up 'epoched' and 'continuous'
    D_epoched = spm_eeg_load(spm_files);
    D_epoched = D_epoched.montage('switch',1); %switch the montage to 1 in order to be safer
    D_epoched.save();
    processed_file_epoched{ii} = spm_files;
    spm_files=[workingdir2 '/' lista(pari).name(4:end)]; %check this in order not to mix up 'epoched' and 'continuous'
    D_continuous = spm_eeg_load(spm_files);
    D_continuous = D_continuous.montage('switch',1); %switch the montage to 1 in order to be safe
    D_continuous.save();
    processed_file_continuous{ii} = spm_files;
    if strcmp(lista(pari).name(20),'0')
        oat.source_recon.session_names(ii)      = {['session' lista(pari).name(21) '_recon.mat']};
        oat.source_recon.results_fnames(ii)     = {['session' lista(pari).name(21) '_recon.mat']};
    else
        oat.source_recon.session_names(ii)      = {['session' lista(pari).name(20:21) '_recon.mat']};
        oat.source_recon.results_fnames(ii)     = {['session' lista(pari).name(20:21) '_recon.mat']};
    end
    disp(num2str(ii))
end
oat.source_recon.D_continuous    = processed_file_continuous; %the file just after AFRICA and before the epoching
oat.source_recon.D_epoched      = processed_file_epoched;
pca_dim_1 = 50;
oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_continuous),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects
oat.source_recon.sessions_to_do = [];
oat.source_recon.sessions_to_do = [1:70]; %sessions to do among the file_list
oat.source_recon.modalities        = {'MEGMAG'; 'MEGPLANAR'}; %'MEGMAG'};
oat.source_recon.conditions        = {'Undefined'};
oat.source_recon.gridstep          = 8; % in mm
oat.source_recon.time_range        = [-2 152]; % time range in secs
oat.source_recon.freq_range        = [2 8]; % frequency range in Hz (delta)
oat.source_recon.type              = 'Scalar';
oat.source_recon.method            = 'beamform';
oat.source_recon.normalise_method  = 'mean_eig';
oat.source_recon.forward_meg       = 'MEG Local Spheres';
%S.source_recon.prefix            = '';
oat.source_recon.report.do_source_variance_maps = 1;
oat.source_recon.dirname           = [workingdir '/source_localsphere_Hz2_Hz8_minorfl'];

oat = osl_check_oat(oat);
oat.to_do = [1 0 0 0];
oat = osl_run_oat(oat);

%% BEAMFORMING RESTING STATE

%setting directories
clear
%basefifdatadir = '/scratch3/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter';
%datadir = '/scratch3/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter'; %without movement compensation
datadir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc'; %with movement compensation
workingdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Francesco';
basefifdatadir = datadir;
basefifdatadirraw = '/raw/sorted/MINDLAB2017_MEG-LearningBach';
tempfif = dir([basefifdatadir '/SUB*']); %try what happen when you have more subject!!
for ii = 1:length(tempfif)
    spm_files_basenames{ii,1} = ['spmeeg_' tempfif(ii).name(1:end-4) '.mat'];
    spm_files{ii,1} = [datadir '/' 'spmeeg_' tempfif(ii).name(1:end-4) '.mat'];
end
%xlsx directory on the server and creation of xlsx_basenmas (procedure mirroring the spm_files)
xlsx_dir = '/scratch7/MINDLAB2017_MEG-LearningBach/BachStimuli';
xlsx_dir_behav = '/scratch7/MINDLAB2017_MEG-LearningBach/BehavioralResponses';
%creates a list of only the recog and rest files (spm objects)
k = 0;
jj = 0;
for ii = 1:length(tempfif)
    if strcmp(tempfif(ii).name(9:12),'reco')
        k = k + 1;
        %spm_files_recog{k,1} = spm_files{ii,1};
        spm_files_recog_basen{k,1} = spm_files_basenames{ii,1};
        spm_files_lear{k,1}=[workingdir '/dff' spm_files_recog_basen{k}];
        xlsx_basenames{k,1} = [spm_files_recog_basen{k,1}(12:15) spm_files_recog_basen{k,1}(21:23) '.xlsx'];
    end
    if strcmp(tempfif(ii).name(9:12),'rest')
        jj = jj + 1;
        spm_files_rest10_basen{jj,1} = spm_files_basenames{ii,1};
        spm_files_rest10{jj,1}=[workingdir '/dff' spm_files_rest10_basen{jj}];
    end
end
%SOURCE RECON
for ii = 1:68 %:length(spm_files_recog_basen)
    spm_files{ii}=[datadir '/e1dff' spm_files_rest10_basen{ii}]; %check this in order not to mix up 'epoched' and 'continuous'
    D_epoched = spm_eeg_load(spm_files{ii});
    D_epoched = D_epoched.montage('switch',1); %switch the montage to 1 in order to be safer
    D_epoched.save();
    processed_file_epoched{ii} = [spm_files{ii}];
    spm_files{ii}=[datadir '/dff' spm_files_rest10_basen{ii}]; %check this in order not to mix up 'epoched' and 'continuous'
    D_continuous = spm_eeg_load(spm_files{ii});
    D_continuous = D_continuous.montage('switch',1); %switch the montage to 1 in order to be safe
    D_continuous.save();
    processed_file_continuous{ii} = [spm_files{ii}];
    oat.source_recon.D_continuous(ii)    = {processed_file_continuous{ii}};%the file just after AFRICA and before the epoching
    oat.source_recon.D_epoched(ii)      = {processed_file_epoched{ii}};
    oat.source_recon.session_names(ii)      = {['session' spm_files_rest10_basen{ii}(14:15)]}; %this is important to get the proper match between participants and sessions (right now it should be that the session that is being saved to disk corresponds to the actual number of the participant (NOTE THAT '%%%%%%%%%%%RUNNING OAT SOURCE RECON ON SESS = X%%%%%%%%%%%%%% TAKES ONLY AN ITERATIVE NUMBER SO X IS NOT EQUAL TO THE NUMBER OF PARTICIPANT.. BUT THAT SHOULD NOT BE A PROBLEM)
    oat.source_recon.results_fnames(ii)     = {['session' num2str(ii) '_recon']};
    disp(num2str(ii))
end
pca_dim_1 = 50;
oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_continuous),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects
oat.source_recon.sessions_to_do = [];
oat.source_recon.sessions_to_do = [1:68];
oat.source_recon.modalities        = {'MEGMAG'; 'MEGPLANAR'}; %'MEGMAG'};
oat.source_recon.conditions        = {'Undefined'};
oat.source_recon.gridstep          = 8; % in mm
oat.source_recon.time_range        = [0 550]; % time range in secs
oat.source_recon.freq_range        = [0.1 2]; % frequency range in Hz
oat.source_recon.type              = 'Scalar';
oat.source_recon.method            = 'beamform';
oat.source_recon.normalise_method  = 'mean_eig';
oat.source_recon.forward_meg       = 'MEG Local Spheres';
oat.source_recon.report.do_source_variance_maps = 1;
oat.source_recon.dirname           = [workingdir '/source_localsphere_Hz01_2_rest10']; % spm_files_recog_basen{ii}(12:25)];

oat = osl_check_oat(oat);
oat.to_do = [1 0 0 0];
oat = osl_run_oat(oat);

%% parcellation AAL and reducing dimensionality from voxel level to ROI level

%this is for giving input such as parcellation, ROI names and template (otherwise it can also find by itself the proper template, in this case 8mm)
parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
ROIsfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_ROIs.txt';
templatefile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/MNI152_T1_8mm_brain_PROVA.nii.gz';
p = parcellation(parcelfile,ROIsfile,templatefile);
p = p.remove_parcels(91:116);
lista = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz415_rest10_2.oat/concat*');
pari = 0; %set a pair counter
for ii = 1:length(lista)/2 %for loop from 1 until length of "lista" divided by 2
    pari = pari +2; %add 2 to pair counter
    D = spm_eeg_load([lista(pari).folder '/' lista(pari).name]); %load every spm file into D and estract the folder and the name from the "lista" with pair condition
    %to be run only the first time
    D = D.montage('switch',2);
    D3 = ROInets.get_node_tcs(D,p.parcelflag(true),'PCA');
    D3.save();
end               
   
%% extracting resting state and sound encoding (in some case called in thsi script free listening) data and save to a single file for later operations

%resting state
list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz2_8_rest10.oat/concat*meeg.mat'); %load dir file "concat" in lista
Drest = extracting_data_LBPD_D(list);
save -v7.3 restingstate.mat Drest

%sound encoding
list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz2_Hz8_minorfl.oat/concat*.mat'); %load concat files from dir into lista
Dfl = extracting_data_LBPD_D(list);
save -v7.3 freelistening.mat Dfl

%% preparing data - from task

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Francesco/Functions'); %francesco's functions
%getting onsets of tones
iww = '/scratch7/MINDLAB2017_MEG-LearningBach/BachStimuli/learminonset.xlsx';
col = 'c'; %column where you have the data
pre = 0; %pre-stimulus time in seconds (this is useful for removing later boundary artifcats)
post = 0; %post-stimulus time in seconds (same..)
d2 = Onset_length_tones_F(iww,col,pre,post); %onset of notes
%%% THE CORRECTION FOR THE 50MS IS DONE HERE.. SUBTRACTING 8 TIME-POINTS
%%% THAT CORRESPONDS TO 50MS.. %%% 300 time-samples are instead added since
%%% the original epoch started 2 seconds (300 time-samples) before the
%%% actual start of the musical piece (this may sound a bit confusing and
%%% could have done more smoothly, but it makes sense and it is ok)
d = d2 + 301 - 8;

%%% here it seems that concatenating some trials is necessary for later
%%% estimation of connectivity.
%%% for static FC that is easily understandable since we need to have more
%%% task data-points to calculate more robust Pearson's correlations. In
%%% the case of instantaneous phase of the envelope, it seems that it is
%%% still necessary to have enough task data-points for allowing the
%%% envelope to capture the actual oscillation associated to the task
%%% (approximately 2-8Hz). Too little data-points not allow it..
%%% Interestingly, if simply enlarging the epoch (so having more
%%% data-points) seems not to work since I guess the envelope tends to
%%% capture other oscillations (or pure noise) because the larger epoch
%%% has data-points that are not connected to the task

n_notes = 605; %number of single notes that you want to epoch
n_trials = 7; %number of trials that you want to concatenate
%epoching, concatenating trials, averaging, removing source leakage correction, calculating envelope
%task data
load('/scratch7/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz2_Hz8_minorfl.oat/freelistening.mat');
Dfl2 = Dfl([1:32 34:end],:); %removing subject 3 that has no resting state (to be used as baseline and that was not fully preprocessed..)
por = [1:4]; %repetitions of the musical piece you want to consider
if por == 4
    Dfl3 = Dfl2([1:42 44:end],:); %that subject is missing last repetition
    Dfl2 = Dfl3;
    clear Dfl3
end
fl_ep_sn = Epoch_singletone_fl_F3(Dfl2,d(1,3) - d(1,2) + 1,n_notes,n_trials,d,1,por);
clear Dfl2 Dfl
% save -v7.3 /scratch7/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz2_Hz8_minorfl.oat/cfl.mat fl_ep_sn % saving the file in a specific path for avoid memory issue

%same for resting state (used as baseline)
load('/scratch7/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz2_8_rest10.oat/restingstate.mat');
%the 201 here looks like a typo, I recall was much bigger, about 520 or so..
n_notes = 201; %the number of notes here is different (slightly smaller) since we do not have enough resting state data
rs10 = Epoch_singletone_fl_F3(Drest,d(1,3) - d(1,2) + 1,n_notes,n_trials,d,0,por);

%combining task and baseline (with extra care not to mix up subject IDs..)
%writing subject IDs for task data
IDtask = cell(71,1);
IDtask2 = cell(71,1);
for ii = 1:71
    IDtask(ii,1) = {['concatMefsession' num2str(ii) '_spm_meeg.mat']};
    IDtask2(ii,1) = {['concatMefsession' num2str(ii) '_recon.mat']};
end
%reordering baseline data (if existent)
dDd = cell(71,2);
data2 = cell(71,2);
for ii = 1:71 %over subjects (original ID)
    %baseline (fro, resting state)
    a = 0;
    count = 0;
    while a == 0 && count < length(rs10(:,2)) %looking for correct subject ID.. until it is found.. if it is not found after looking in all subjects, you later assign dDd = []
        count = count + 1;
        a = double(strcmp(IDtask{ii,1},rs10{count,2}));
    end
    if a == 1
        dDd(ii,1) = rs10(count,1);
        dDd(ii,2) = rs10(count,2);
    else
        dDd(ii,1) = {[]};
        data2(ii,1) = {[]}; %doing that also in data2 since we need to completely discard the subject without baseline..
        data2(ii,2) = {[]};
    end
    %task
    a = 0;
    count = 0;
    while a == 0 && count < length(fl_ep_sn(:,2)) %looking for correct subject ID.. until it is found.. if it is not found after looking in all subjects, you later assign dDd = []
        count = count + 1;
        a = double(strcmp(IDtask2{ii,1},fl_ep_sn{count,2}));
    end
    if a == 1 && ~isempty(dDd{ii,1})
        data2(ii,1) = fl_ep_sn(count,1);
        data2(ii,2) = fl_ep_sn(count,2);
    else
        dDd(ii,1) = {[]};
        data2(ii,1) = {[]}; %doing that also in data2 since we need to completely discard the subject without baseline..
        data2(ii,2) = {[]};
    end
end
%actual data
data(:,1) = data2(:,1); %task
data(:,2) = dDd(:,1); %baseline
%storing labels
L(1) = {'dataset_1_task'};
L(2) = {'dataset_2_RSBaseline'}; 
%storing data for subsequent function
DATA = [];
DATA.prep_data = data;
DATA.label = L;

%% static FC - Pearson's correlations

load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/StaticFC/DATASFC.mat') %THIS MUST BE YOUR OWN DIRECTORY WHERE YOU STORED THE PROVIDED DATA

S = [];
S.contrast = [2 1];
S.resam_number = [];
S.fsample = 150;
S.plot_conditions = [1 1];
S.plot_contr = 1;
S.perm_degree = 100;
S.perm_max = 0;
S.p_perm = 0;
S.p_threshMC = 0.001;
S.clims_contr = [-8 8];
S.clims_cond = [-0.2 0.2];
S.parametric = 1;
S.schemaball_contrast_label = 1;
S.schemaball_cond_label = 0;
S.schemplth = 100;
S.colbacksch = [];
S.lab_parc_ph = '/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat';
S.outpath = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/TestingPhDFunctions_Scratch';
% S.name_sch = 'testschemaballstaticFC';
S.name_sch = [];

[H, P_binary, P_ROIs, matf ] = static_FC_MEG_LBPD_D(S, DATA);
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))

P_ROIs
% writetable(cell2table(P_ROIs),['staticFC_new_1_190_degree10000perm.xlsx'],'Sheet',1)

%% DYNAMIC FUNCTIONAL CONNECTIVITY (DFC)

%% phase synchrony 

S = [];
S.data = DATA;
S.save_path = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/DynamicFC';
S.ROIs = 90;
S.left_bound = 37; %discarding the first concatenated trial
S.right_bound = 36; %discarding the last concatenated trial

phasesynchrony_LBPD_D(S);

%concatenated trials, important for Hilbert estimation, are now averaged together
%task (task1)
load([S.save_path '/task_1.mat'])
iFC_time_4D = squeeze(mean(reshape(iFC_time_4D,[90,90,180/5,5,68]),4));
save([S.save_path '/task_1.mat'],'iFC_time_4D','-v7.3')
%baseline (task2).. this is actually not really necessary since the baseline will be then averaged over time (because it is a baseline..) but for clarity I do this..
load([S.save_path '/task_2.mat'])
iFC_time_4D = squeeze(mean(reshape(iFC_time_4D,[90,90,180/5,5,68]),4));
save([S.save_path '/task_2.mat'],'iFC_time_4D','-v7.3')

%contrasting Sound Encoding DFC vs resting state DFC (here the resting state is used as baseline and different contrasts have been run for different groups of participants)
S = [];
S.load_path = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz';
% S.save_path = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub';
S.save_path = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub_FAM1';
% S.subjs_list = {[1,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71]}; %cell matrix with the number corresponding to the subjects grouped into different subsample; example: {[1,2,3];[4,5,6];[7,8,9]}; have {[1,2,3,4,5,6,7,8,9]} for only one group;
% S.subjs_list = {[1,9,10,11,13,15,17,23,24,27,28,32,33,35,37,38,40,42,45,46,47,48,50,51,57,58,60,62,66,69,71]}; %cell matrix with the number corresponding to the subjects grouped into different subsample; example: {[1,2,3];[4,5,6];[7,8,9]}; have {[1,2,3,4,5,6,7,8,9]} for only one group;
S.subjs_list = {[4,5,6,7,12,14,16,18,19,20,22,25,26,29,30,31,34,36,39,41,43,44,49,52,53,54,55,56,63,64,65,67,68,70]};
S.group_comparison = [1]; %[1] for one group; [1 2] for two groups (referred to the subjs_list)
S.condition_contrast = [2 1];
S.cond_contr_meanbaseline = 1;
S.Num_window = []; %leave it empty [] for not doing any subaveraging

ps_statistics_LBPD_D(S);


%degree MCS over time (calculated for all subjects or for different subjects grouped according to their musical expertise or general working memory skills (WM) or auditory working memory skills (MET)
phase_dir1(1) = {'/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub'}; %general path everybody
phase_dir1(2) = {'/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub_WMH'}; %high working memory
phase_dir1(3) = {'/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub_WML'}; %low working memory
phase_dir1(4) = {'/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub_METH'}; %high MET test
phase_dir1(5) = {'/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub_METL'}; %low MET test
phase_dir1(6) = {'/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub_MUS1'}; %musicians
phase_dir1(7) = {'/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub_MUS0'}; %non-musicians
phase_dir1(8) = {'/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub_FAM1'}; %non-musicians
phase_dir1(9) = {'/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub_FAM0'}; %non-musicians
for opo = 8:length(phase_dir1)
    %settings
    matrn = 0; %matrix number; 0 for every matrix
    thr = 0; %0 for non-binarising the t-values matrix; otherwise threshold for binarising it
    contr = '1_2';
    posl = 1; %1 for cond1 > cond2; otherwise for cond2 > cond1
    thresh = 0;
    permut = 1000;
    threshMC = 0.001;
    ROI_label = '/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat';
    reshlab = 0;
    perm_max = 0;
    phase_dir = phase_dir1{opo};
    load([phase_dir '/Contr_' contr '_STAT_zval_fullcontr_group_comp_1_1.mat']);
    load([phase_dir '/Contr_' contr '_STAT_p_fullcontr_group_comp_1_1.mat']);
    if thr ~= 0 
        STAT_p_pos_bin = double((double(STAT_p_fullcontr < thr) + double(STAT_zval_fullcontr > 0)) == 2); %isolating conditions when you have both p_values significant at the threshold defined by the user and the zval that is positive
        STAT_p_neg_bin = double((double(STAT_p_fullcontr < thr) + double(STAT_zval_fullcontr < 0)) == 2);
    else
        STAT_p_pos_bin = STAT_zval_fullcontr;
        STAT_p_neg_bin = STAT_zval_fullcontr * (-1);
    end
    if matrn == 0
        x = 1:size(STAT_p_fullcontr,3);
        ROII = cell(x(end),1);
    else
        x = 1;
        ROII = cell(1);
    end
    PP = zeros(90,length(x));
    PBIN = zeros(90,length(x));
    HSEG = zeros(length(x),1);
    for ii = x
        if posl == 1
            if matrn == 0
                P = STAT_p_pos_bin(:,:,ii);
            else
                P = STAT_p_pos_bin(:,:,matrn);
            end
        else
            if matrn == 0
                P = STAT_p_neg_bin(:,:,ii);
            else
                P = STAT_p_neg_bin(:,:,matrn);
            end
        end
        [h_deg, P_binary_deg, P_ROIs_deg, P_val, h_seg, ROI_seg_label, ROI_seg, glob_seg] = degree_segregation_MCS_LBPD_D(P,thresh,permut,threshMC,ROI_label,reshlab,perm_max);
        ROII(ii) = {P_ROIs_deg};
        PP(:,ii) = P_val;
        HSEG(ii) = h_seg;
        PBIN(:,ii) = P_binary_deg;
        disp(ii)
        P_ROIs_deg
        writetable(cell2table(P_ROIs_deg),['Contr_' contr '_pos_' num2str(posl) '_note' num2str(ii) '.xlsx'],'Sheet',1)
    end
    mop = num2str(threshMC);
    save([phase_dir '/tonevsbaseline_' num2str(permut) 'perm_nosubave_maxdeg' num2str(perm_max) 'thresh' num2str(thr) '_threshMC' mop(3:end) '_' phase_dir(76:78)],'PP','ROII','PBIN','HSEG');
end

%% MCS for contrasting the two time-windows (1 - 110ms and 111 - 220ms)

schempercl = 100; %set the percentage of connections to be shown in the schemaball

NAM{1} = 'TimeWind1';
NAM{2} = 'TimeWind2';
hemisph = 0; %1 to test the combined hemispheres; 0 for separate hemispheres
permut = 10000;
thresh = 2.7e-04;
%first directory (to be set for contrast)
phase_dir1 = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub'; %YOUR OWN PATH WHERE YOU PUT THE PROVIDED DATA
%second directory (to be set for contrast)
phase_dir2 = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper2/IFC/2_8Hz/nosub'; %YOUR OWN PATH WHERE YOU PUT THE PROVIDED DATA
%elaborated way to get the data from first 18 samples and last 18 samples..
kjk = dir([phase_dir1 '/old*']);
load([kjk(1).folder '/' kjk(1).name]);
a = PBIN;
kjk = dir([phase_dir2 '/old*']);
load([kjk(1).folder '/' kjk(1).name]);
b = PBIN;
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat');
ROI_lab = cell(90,1);
for ii = 1:90
    ROI_lab(ii) = {lab(ii,:)};
end
a = a(:,5:18); b = b(:,19:32);
X = 0; %on the all time-points of a (first 18s) and b (last 18s)
if hemisph == 1
    %putting together same areas of different hemispheres
    [n,m] = size(a);
    a2 = zeros(n/2,m);
    b2 = zeros(n/2,m);
    ROI_label2 = cell(n/2,1);
    cnt = -1;
    for ii = 1:n/2
        cnt = cnt + 2;
        a2(ii,:) = a(cnt,:) + a(cnt + 1,:);
        b2(ii,:) = b(cnt,:) + b(cnt + 1,:);
        ROI_label2(ii,1) = {ROI_lab{cnt}(3:end)};
    end
    ROI_final_hem_tog = diff_conditions_matr_degree_MCS_LBPD_D(a2, b2, ROI_label2, X, permut, thresh );
    ROI_final_hem_tog.diff
else
    out = diff_conditions_matr_degree_MCS_LBPD_D(a, b, ROI_lab, X, permut, thresh );
    out.diff
end
p_thresh = 0.05;
%INT
%first time-window
S = size(a);
INT = a;
%lab_ab
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat');
Order = [1:2:S(1) S(1):-2:2];
INT = INT(Order,:);
lab_ab = lab(Order,:);
%COUP
%first time-window (to be set for contrast)
load([phase_dir '/Contr_1_2_STAT_p_fullcontr_group_comp_1_1.mat']);
load([phase_dir '/Contr_1_2_STAT_zval_fullcontr_group_comp_1_1.mat']);
STAT_old = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
STAT_old = STAT_old(Order,Order,5:18);
COUP = STAT_old;
%r
r = out.diff;
r{end,3} = []; r{end-1,3} = [];
R{1} = r;
COUP2{1} = COUP;
INT2{1} = INT;
%second time-window
INT = b;
INT = INT(Order,:);
STAT_old = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
STAT_old = STAT_old(Order,Order,19:32);
COUP = STAT_old;
%r
r = out.diff;
r{1,3} = []; r{2,3} = [];
r{end,2} = r{end,2}*(-1); r{end-1,2} = r{end-1,2}*(-1);
R{2} = r;
COUP2{2} = COUP;
INT2{2} = INT;
%actual plotting settings and function
for pol = 1:2
    [POS,NEG] = signROIs_degree_connotherROIs_LBPD_D(INT2{pol}, COUP2{pol}, lab_ab, R{pol}, X);
    label_path = '/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat';
    matf = signROIdegree_otherROI_prepareplotting_LBPD_D(POS, NEG, label_path);
    %parcellation AAL
    %this is for giving input such as parcellation, ROI names and template (otherwise it can also find by itself the proper template, in this case 8mm)
    parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
    ROIsfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_ROIs.txt';
    templatefile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/MNI152_T1_8mm_brain_PROVA.nii.gz';
    p = parcellation(parcelfile,ROIsfile,templatefile);
    p = p.remove_parcels(91:116); %removing the cerebellum
    %actual plotting function
    S = [];
    S.MATF = matf;
    S.symmetric_l = 0;
    S.outpath = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/DynamicFC/BrainTemplateFigures';
    S.p = p;
    S.thr_cortex = [0.993]; f = num2str(S.thr_cortex);
    S.name_gif = [NAM{pol}];
%     S.name_gif = [];
    S.frame_vect = [15,90,130];
    S.fr_spec = {'Posterior_L_R';'Frontal_R_L';'Hem_L'};
    S.schball_l = 1;
    S.extr = [];
    S.lab_parc_ph = '/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat';
    S.colbacksch = [];
    S.schemplth = schempercl;
%     S.name_sch = [NAM{pol}];
    S.name_sch = [];
    %actual function
    FC_Brain_Schemaball_plotting_LBPD_D(S)
end

%% DFC over the whole time-window

schempercl = 100; %set the percentage of connections to be shown in the schemaball

permut = 10000;
thresh = 2.7e-04;
%first directory (to be set for contrast)
phase_dir1 = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/DynamicFC'; %YOUR OWN PATH WHERE YOU PUT THE PROVIDED DATA
kjk = dir([phase_dir1 '/old*']);
load([kjk(1).folder '/' kjk(1).name]);
a = PBIN;
b = zeros(size(a,1),size(a,2));
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat');
ROI_lab = cell(90,1);
for ii = 1:90
    ROI_lab(ii) = {lab(ii,:)};
end
X = 0; %on the all time-points of a (first 18s) and b (last 18s)
out = diff_conditions_matr_degree_MCS_LBPD_D(a, b, ROI_lab, X, permut, thresh );
out.diff
p_thresh = 0.05;
%INT
S = size(a);
INT = a;
%lab_ab
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat');
Order = [1:2:S(1) S(1):-2:2];
INT = INT(Order,:);
lab_ab = lab(Order,:);
%COUP
load([phase_dir1 '/Contr_1_2_STAT_p_fullcontr_group_comp_1_1.mat']);
load([phase_dir1 '/Contr_1_2_STAT_zval_fullcontr_group_comp_1_1.mat']);
STAT_old = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
STAT_old = STAT_old(Order,Order,:);
COUP = STAT_old;
%r
r = out.diff;
r{end,3} = []; r{end-1,3} = [];
R{1} = r;
COUP2{1} = COUP;
INT2{1} = INT;
%actual plotting settings and function
[POS,NEG] = signROIs_degree_connotherROIs_LBPD_D(INT2{pol}, COUP2{pol}, lab_ab, R{pol}, X);
label_path = '/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat';
matf = signROIdegree_otherROI_prepareplotting_LBPD_D(POS, NEG, label_path);
%parcellation AAL
%this is for giving input such as parcellation, ROI names and template (otherwise it can also find by itself the proper template, in this case 8mm)
parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
ROIsfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_ROIs.txt';
templatefile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/MNI152_T1_8mm_brain_PROVA.nii.gz';
p = parcellation(parcelfile,ROIsfile,templatefile);
p = p.remove_parcels(91:116); %removing the cerebellum
%actual plotting function
S = [];
S.MATF = matf;
S.symmetric_l = 0;
S.outpath = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/DynamicFC/BrainTemplateFigures';
S.p = p;
S.thr_cortex = [0.993]; f = num2str(S.thr_cortex);
S.name_gif = ['whole_timewindow'];
%     S.name_gif = [];
S.frame_vect = [15,90,130];
S.fr_spec = {'Posterior_L_R';'Frontal_R_L';'Hem_L'};
S.schball_l = 1;
S.extr = [];
S.lab_parc_ph = '/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat';
S.colbacksch = [];
S.schemplth = schempercl;
%     S.name_sch = [NAM{pol}];
S.name_sch = [];
%actual function
FC_Brain_Schemaball_plotting_LBPD_D(S)


%% CREATING NIFTI FILES TO BE THEN EXPORTED FROM THE SERVER AND USED IN WORKBENCH

%these are the codes used to produce nifti files that have then been further elaborated in WorkBench to obtain some of the figures reported in the paper.
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/creatingAALnifti') %path to function
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat');
vector = zeros(1,90);
for ii = 1:90
    idx = 0;
    cnt = 0;
    while idx == 0 %looking for original index (non symmetric indices here since I provide non-symmetric AAL labels)
        cnt = cnt + 1;
        idx = strcmp(lab(ii,:),out.diff{cnt,1});
    end
    if ~isempty(out.diff{cnt,3}) %if it is significant the cnt ROI
        vector(1,ii) = out.diff{cnt,2}; %storing the corresponding value
    end
%     vector(1,ii) = out.cond1{cnt,2};
end

%mask with significant ROIs and then taking average statistics for those ROIs
namenii = ['whole_timewindow.nii']; %name for the figure (remember to make it ending in '.nii'
symmetric = 0; %if the AAL vector is submitted in symmetric order
create_AALnifti(vector,namenii,symmetric);

%% comparing brain areas degree centrality during sound encoding among different groups (high and low general working memory (WM); high and low auditory working memory (MET); high and low musical expertise

grouplabel = 4; %set 1 for WM; 2 for MET; 3 for musical expertise; 4 = familiarity with the Bach's prelude
schempercl = 100; %set the percentage of connections to be shown in the schemaball
nam = 'try'; %specify name to be assigned to brain template figures

hemisph = 0; %1 to test the combined hemispheres; 0 for separate hemispheres
X = 0;
permut = 10000;
thresh = 9.2e-05;
phase_dir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/DynamicFC'; %THIS MUST BE YOUR OWN DIRECTORY WHERE YOU STORED THE PROVIDED DATA
if grouplabel == 1
    %first group (to be set for contrast)
    kjk = dir([phase_dir '/*b_WMH.mat']);
    load([kjk(1).folder '/' kjk(1).name]);
    a = PBIN;
    %second group (to be set for contrast)
    kjk = dir([phase_dir '/*b_WML.mat']);
    load([kjk(1).folder '/' kjk(1).name]);
    b = PBIN;
elseif grouplabel == 2
    %first group (to be set for contrast)
    kjk = dir([phase_dir '/*b_METH.mat']);
    load([kjk(1).folder '/' kjk(1).name]);
    a = PBIN;
    %second group (to be set for contrast)
    kjk = dir([phase_dir '/*b_METL.mat']);
    load([kjk(1).folder '/' kjk(1).name]);
    b = PBIN;
elseif grouplabel == 3
    %first group (to be set for contrast)
    kjk = dir([phase_dir '/*b_MUS1.mat']);
    load([kjk(1).folder '/' kjk(1).name]);
    a = PBIN;
    %second group (to be set for contrast)
    kjk = dir([phase_dir '/*b_MUS0.mat']);
    load([kjk(1).folder '/' kjk(1).name]);
    b = PBIN;
elseif grouplabel == 4
    %first group (to be set for contrast)
    kjk = dir([phase_dir '/*b_FAM1.mat']);
    load([kjk(1).folder '/' kjk(1).name]);
    a = PBIN;
    %second group (to be set for contrast)
    kjk = dir([phase_dir '/*b_FAM0.mat']);
    load([kjk(1).folder '/' kjk(1).name]);
    b = PBIN;
else
    warning('you need to select either 1, 2, 3 or 4..')
end
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat');
ROI_lab = cell(90,1);
for ii = 1:90
    ROI_lab(ii) = {lab(ii,:)};
end
if hemisph == 1
    %putting together same areas of different hemispheres
    [n,m] = size(a);
    a2 = zeros(n/2,m);
    b2 = zeros(n/2,m);
    ROI_label2 = cell(n/2,1);
    cnt = -1;
    for ii = 1:n/2
        cnt = cnt + 2;
        a2(ii,:) = a(cnt,:) + a(cnt + 1,:);
        b2(ii,:) = b(cnt,:) + b(cnt + 1,:);
        ROI_label2(ii,1) = {ROI_lab{cnt}(3:end)};
    end
    ROI_final_hem_tog = diff_conditions_matr_degree_MCS_LBPD_D(a2, b2, ROI_label2, X, permut, thresh );
    ROI_final_hem_tog.diff
else
    out = diff_conditions_matr_degree_MCS_LBPD_D(a, b, ROI_lab, X, permut, thresh );
    out.diff
end
%writetable(cell2table(out.diff),['Diffconddegree_10000permutations_oldvsnew_sub5.xlsx'],'Sheet',1) %exporting results in excel file
p_thresh = 0.05;
%INT
S = size(a);
INT = zeros(S(1),S(2),2);
INT(:,:,1) = a;
INT(:,:,2) = b;
%lab_ab
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat');
Order = [1:2:S(1) S(1):-2:2];
INT = INT(Order,:,:);
lab_ab = lab(Order,:);
%COUP
%cond1
if grouplabel == 1
    %first group (to be set for contrast)
    load([phase_dir '/Contr_1_2_STAT_p_fullcontr_group_comp_1_1_WMH.mat']);
    load([phase_dir '/Contr_1_2_STAT_zval_fullcontr_group_comp_1_1_WMH.mat']);
    STAT_old = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
    STAT_old = STAT_old(Order,Order,:);
    %second group (to be set for contrast)
    load([phase_dir '/Contr_1_2_STAT_p_fullcontr_group_comp_1_1_WML.mat']);
    load([phase_dir '/Contr_1_2_STAT_zval_fullcontr_group_comp_1_1_WML.mat']);
    STAT_new = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
    STAT_new = STAT_new(Order,Order,:);
elseif grouplabel == 2
    %first group (to be set for contrast)
    load([phase_dir '/Contr_1_2_STAT_p_fullcontr_group_comp_1_1_METH.mat']);
    load([phase_dir '/Contr_1_2_STAT_zval_fullcontr_group_comp_1_1_METH.mat']);
    STAT_old = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
    STAT_old = STAT_old(Order,Order,:);
    %second group (to be set for contrast)
    load([phase_dir '/Contr_1_2_STAT_p_fullcontr_group_comp_1_1_METL.mat']);
    load([phase_dir '/Contr_1_2_STAT_zval_fullcontr_group_comp_1_1_METL.mat']);
    STAT_new = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
    STAT_new = STAT_new(Order,Order,:);
elseif grouplabel == 3
    %first group (to be set for contrast)
    load([phase_dir '/Contr_1_2_STAT_p_fullcontr_group_comp_1_1_MUS1.mat']);
    load([phase_dir '/Contr_1_2_STAT_zval_fullcontr_group_comp_1_1_MUS1.mat']);
    STAT_old = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
    STAT_old = STAT_old(Order,Order,:);
    %second group (to be set for contrast)
    load([phase_dir '/Contr_1_2_STAT_p_fullcontr_group_comp_1_1_MUS0.mat']);
    load([phase_dir '/Contr_1_2_STAT_zval_fullcontr_group_comp_1_1_MUS0.mat']);
    STAT_new = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
    STAT_new = STAT_new(Order,Order,:);
elseif grouplabel == 4
    %first group (to be set for contrast)
    load([phase_dir '/Contr_1_2_STAT_p_fullcontr_group_comp_1_1_FAM1.mat']);
    load([phase_dir '/Contr_1_2_STAT_zval_fullcontr_group_comp_1_1_FAM1.mat']);
    STAT_old = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
    STAT_old = STAT_old(Order,Order,:);
    %second group (to be set for contrast)
    load([phase_dir '/Contr_1_2_STAT_p_fullcontr_group_comp_1_1_FAM0.mat']);
    load([phase_dir '/Contr_1_2_STAT_zval_fullcontr_group_comp_1_1_FAM0.mat']);
    STAT_new = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
    STAT_new = STAT_new(Order,Order,:);    
else
    warning('you need to select either 1, 2 or 3..')
end
S2 = size(STAT_old);
COUP = zeros(S2(1),S2(2),S2(3),2);
COUP(:,:,:,1) = STAT_old;
COUP(:,:,:,2) = STAT_new;
%r
r = out.diff;
[POS,NEG] = signROIs_degree_connotherROIs_LBPD_D(INT, COUP, lab_ab, r, X);
label_path = '/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat';
matf = signROIdegree_otherROI_prepareplotting_LBPD_D(POS, NEG, label_path);
%parcellation AAL
%this is for giving input such as parcellation, ROI names and template (otherwise it can also find by itself the proper template, in this case 8mm)
parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
ROIsfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_ROIs.txt';
templatefile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/MNI152_T1_8mm_brain_PROVA.nii.gz';
p = parcellation(parcelfile,ROIsfile,templatefile);
p = p.remove_parcels(91:116); %removing the cerebellum
%actual plotting function
S = [];
S.MATF = matf;
S.symmetric_l = 0;
S.outpath = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/DynamicFC/BrainTemplateFigures';
S.p = p;
S.thr_cortex = [0.993]; f = num2str(S.thr_cortex);
S.name_gif = [nam];
% S.name_gif = [];
S.frame_vect = [15,90,130];
S.fr_spec = {'Posterior_L_R';'Frontal_R_L';'Hem_L'};
S.schball_l = 1;
S.extr = [-6 6];
S.lab_parc_ph = '/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat';
S.colbacksch = [];
S.schemplth = schempercl;
S.name_sch = [];
%actual function
FC_Brain_Schemaball_plotting_LBPD_D(S)

%% CREATING NIFTI FILES TO BE THEN EXPORTED FROM THE SERVER AND USED IN WORKBENCH

%these are the codes used to produce nifti files that have then been further elaborated in WorkBench to obtain some of the figures reported in the paper
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/creatingAALnifti') %path to function
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat');
vector = zeros(1,90);
for ii = 1:90
    idx = 0;
    cnt = 0;
    while idx == 0 %looking for original index (non symmetric indices here since I provide non-symmetric AAL labels)
        cnt = cnt + 1;
        idx = strcmp(lab(ii,:),out.diff{cnt,1});
    end
    if ~isempty(out.diff{cnt,3}) %if it is significant the cnt ROI
        vector(1,ii) = out.diff{cnt,2}; %storing the corresponding value
    end
%     vector(1,ii) = out.cond1{cnt,2};
end

%mask with significant ROIs and then taking average statistics for those ROIs
namenii = ['FAM_Degree_1_36_everybody.nii']; %name for the figure (remember to make it ending in '.nii'
symmetric = 0; %if the AAL vector is submitted in symmetric order
create_AALnifti(vector,namenii,symmetric);

%% plotting IFC

phase_dir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/SoundEncodingFC/DynamicFC'; %THIS MUST BE YOUR OWN DIRECTORY WHERE YOU STORED THE PROVIDED DATA
load([phase_dir '/Contr_1_2_STAT_zval_fullcontr_group_comp_1_1.mat']);
%first time-window
a = STAT_zval_fullcontr(:,:,:);
S = [];
S.STAT = a;
S.ordresh = 1;
S.cond_contr_meanbaseline = 1;
S.Num_window = 2;
S.imag_contr_lim = [-2 2];
S.indimagesnumber = [1:2];
IFC_plotting_LBPD_D(S)

%PLEASE REMEMBER THAT IF YOU WANT TO VISUALIZE THE MATRICES IN RED-WHITE-BLUE COLORS YOU NEED TO SELECT THE FIGURE AND THEN RUN THE FOLLOWING TWO LINES..
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))

%%

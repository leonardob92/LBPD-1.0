%% Temporal pattern recognition in the human brain: a dual simultaneous processing - LEONARDO BONETTI

%%


%BEFORE PROCEEDING, PLEASE NOTE:

%As follows, every analysis that we made has been reported to ensure full disclosure.
%Please, note that in this Leading script, I use built-in functions,
%functions from well-known toolboxes (e.g. OSL, SPM, FieldTrip) and
%in-house-built functions, which are reported together with this script in the collection named LBPD.
%Data can be provided according to the Danish regulation, upon reasonable request.
%If interested, please contact Leonardo Bonetti.
%More information is provided in the ReadMe.mat file that I strongly advise to read.
%Finally, in this script the memorized (M) and novel (N) patterns have been
%sometimes called old (corresponding to M) and new (corresponding to N).


%Leonardo Bonetti
%leonardo.bonetti@psych.ox.ac.uk
%leonardo.bonetti@clin.au.dk



%%

%% START UP FUNCTIONS.. (LBPD_startup_D)

%starting up some functions for LBPD toolbox.
%in this work, this should be required only for computation of univariate
%tests and Monte-Carlo simulations (MCS)

%starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl);

%%

%% *** PRE-PROCESSING ***

%% Maxfilter

%some of the following codes (mainly in the pre-processing) are not perfecty optimized, but they work and now it is not worthy to spend time writing them much better
%OBS!!! before running maxfilter, you need to close matlab, open the terminal, write: 'use anaconda', then open matlab and run maxfilter script
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup
basefifdatadir = '/raw/sorted/MINDLAB2017_MEG-LearningBach'; %0001/20171220_000000/MEG';
basefifdatadirraw = '/raw/sorted/MINDLAB2017_MEG-LearningBach'; %0001/20171220_000000/MEG';
datadir = '/projects/MINDLAB2017_MEG-LearningBach/scratch/Leonardo/LearningBach/maxfilter_preproc/'; %0001';
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
        display('problem in getting the raw MEG data!');
    end
    fif_files{ii,1} = dummy2_3fif;
    for k = 1:length(dummy2_3fif) %this iterates over files for each subject
        temp4fif = [dummy2_3fif(k).folder '/' dummy2_3fif(k).name '/files/' dummy2_3fif(k).name(5:end) '.fif'];
        a{ii,k} = temp4fif; %just a cell for storing the all paths
        S = [];
        S.dataset = temp4fif;
        S.outfile = ['spmeeg_SUBJ' temp1fif(ii).name dummy2_3fif(k).name(5:end)];
    end
end
%new proper lines for maxfilter
maxfilter_path = '/neuro/bin/util/maxfilter';
maxDir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach_Fsample1000Hz'; %with movement compensation
project = 'MINDLAB2017_MEG-LearningBach';
for ii = str2double(subnum):length(temp1fif) %iterates over subjects
    fif_files_dummy = fif_files{ii,1};
    spm_files_dummy = fif_files{ii,1};
    for j = 1:length(spm_files_dummy) %iterates over experimental blocks
        if ~isempty(strfind(spm_files_dummy(j).name,'rest10')) || ~isempty(strfind(spm_files_dummy(j).name,'learminor')) || ~isempty(strfind(spm_files_dummy(j).name,'recogminor'))  %here for now I want only learminor, recogminor and rest10
            rawName{j} = [fif_files_dummy(j).folder '/' fif_files_dummy(j).name '/files' '/' fif_files_dummy(j).name(5:end) '.fif'];
            [~,n,~] = fileparts(rawName{j});
            maxfName = ['SUBJ00' num2str(ii) rawName{j}(regexp(rawName{j},'files/')+6:regexp(rawName{j},'.fif')-1)];
            badchans = [];
            %commands for submitting the job to the cluster
            cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ',num2str(project),' "',maxfilter_path,' -f ',[fif_files_dummy(j).folder '/' fif_files_dummy(j).name '/files' '/' fif_files_dummy(j).name(5:end) '.fif'],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 3 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
            system(cmd);
        end
    end
end

%%

%please note that these preprocessing steps have been run for a large
%number of data files to optmize the procedures. Some of them were not used in this specific paper
%but in other projects, however these preprocessing steps have been run together.

%% setting directories

clear
datadir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc'; %with movement compensation in maxfilter
workingdir = datadir;
basefifdatadir = datadir;

%% converting fif files into SPM objects and creating some file names for later steps

%convert fif file into SPM object
temp1fif = dir([basefifdatadir '/*.fif']);
for ii = 1:length(temp1fif) %this iterates over blocks through all subjects (remember to set only the subjects that you really want)
    temp2fif = [basefifdatadir '/' temp1fif(ii).name];
    fif_files{ii,1} = temp2fif;
    S = [];
    S.dataset = fif_files{ii,1};
    D = spm_eeg_convert(S); 
end
%create the spm object files path and name
for ii = 1:length(temp1fif)
    spm_files_basenames{ii,1} = ['spmeeg_' temp1fif(ii).name(1:end-4) '.mat'];
    spm_files{ii,1} = [datadir '/' 'spmeeg_' temp1fif(ii).name(1:end-4) '.mat'];
end

%% creates a list of only the recog files (spm objects)

k = 0;
kk = 0;
jj = 0;
clear spm_files_lear_basen spm_files_recog_basen spm_files_lear spm_files_recog xlsx_basenames
for ii = 1:length(temp1fif)
    if strcmp(temp1fif(ii).name(9:12),'lear')
        k = k + 1;
        %spm_files_recog{k,1} = spm_files{ii,1};
        spm_files_lear_basen{k,1} = spm_files_basenames{ii,1}; %this is not computationally elegant and the space for these variable should be preallocated but it does not really change anything for my purposes..
        spm_files_lear{k,1}=[workingdir '/dff' spm_files_lear_basen{k}];
%         xlsx_basenames{k,1} = [spm_files_lear_basen{k,1}(12:15) spm_files_lear_basen{k,1}(21:23) '.xlsx'];
    end
    if strcmp(temp1fif(ii).name(9:12),'reco')
        kk = kk + 1;
        %spm_files_recog{k,1} = spm_files{ii,1};
        spm_files_recog_basen{kk,1} = spm_files_basenames{ii,1}; %this is not computationally elegant and the space for these variable should be preallocated but it does not really change anything for my purposes..
        spm_files_recog{kk,1}=[workingdir '/dff' spm_files_recog_basen{kk}];
        xlsx_basenames{kk,1} = [spm_files_recog_basen{kk,1}(12:15) spm_files_recog_basen{kk,1}(21:23) '.xlsx'];  %%%%%%%%%%CHANGE HERE QUESTE 3 LINEE CONTROLLA IL NOME!!!!!!!!! ****** %%%%%%%%%%%%%%%%%%%%%
    end
    if strcmp(temp1fif(ii).name(9:12),'rest')
        jj = jj + 1;
        spm_files_rest10_basen{jj,1} = spm_files_basenames{ii,1};
        spm_files_rest10{jj,1}=[workingdir '/dff' spm_files_rest10_basen{jj}];
    end
end

%% filtering and downsampling

%high-pass filter
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

%% epoching - one epoch per old/new excerpt (baseline = (-)100ms)

prefix_tobeadded = 'e80'; %prefix to be added
%creates a list of only the recog files (subgroup of the whole set of data files)
k = 0;
for ii = 1:length(temp1fif)
    if strcmp(temp1fif(ii).name(9:12),'reco')
        k = k + 1;
        %spm_files_recog{k,1} = spm_files{ii,1};
        spm_files_recog_basen{k,1} = spm_files_basenames{ii,1};
        spm_files_recog{k,1}=[workingdir '/dff' spm_files_recog_basen{k}];
        xlsx_basenames{k,1} = [spm_files_recog_basen{k,1}(12:15) spm_files_recog_basen{k,1}(21:23) '.xlsx'];
    end
end
%iterates for the recognition files only
for ii = 3:3:length(spm_files_recog_basen) %indexing only the files wanted for this paper
    %load the spm object
    D = spm_eeg_load(spm_files_recog{ii,1});
    events = D.events;
    %takes the correct triggers sent during the recording
    count_evval = 0;
    for ieve = 1:length(events)
        if strcmp(events(ieve).type,'STI101_up')
            if events(ieve).value == 5
                count_evval = count_evval + 1;
                trigcor(count_evval,1) = events(ieve).time + 0.010; %this takes the correct triggers and add 10ms of delay of the sound travelling into the tubes
            end
        end
    end
    if count_evval ~= 4 && count_evval ~= 80
        disp('warning.. there is something wrong with the triggers');
        disp(spm_files_reco_basen{ii}(8:end-12));
    end
    trl_sam = zeros(length(trigcor),3);
    trl_sec = zeros(length(trigcor),3);
    deftrig = zeros(length(trigcor),1);
    for k = 1:length(trigcor)
        deftrig(k,1) = 0.012 + trigcor(k,1); %adding a 0.012 seconds delay to the triggers sent during the experiment (this delay was due to technical reasons related to the stimuli)
        trl_sec(k,1) = deftrig(k,1) - 0.1000; %beginning time-window epoch in s (please note that we computed this operation two times, obtaining two slightly different pre-stimulus times.
        %this was done because for some computations was convenient to have a slightly longer pre-stimulus time
        trl_sec(k,2) = deftrig(k,1) + 6.000; %end time-window epoch in s
        trl_sec(k,3) = trl_sec(k,2) - trl_sec(k,1); %range time-windows in s
        trl_sam(k,1) = round(trl_sec(k,1) * 150); %beginning time-window epoch in samples
        trl_sam(k,2) = round(trl_sec(k,2) * 150); %end time-window epoch in samples
        trl_sam(k,3) = -15; %sample before the onset of the stimulus (corresponds to 0.100ms)
    end
    %creates the epochinfo structure that is required for the source reconstruction later
    epochinfo.trl = trl_sam;
    epochinfo.time_continuous = D.time;
    %switch the montage to 0 because for some reason OSL people prefer to do the epoching with the not denoised data
    D = D.montage('switch',0);
    %build structure for spm_eeg_epochs
    S = [];
    S.D = D;
    S.trl = trl_sam;
    S.prefix = prefix_tobeadded;
    D = spm_eeg_epochs(S);
    %store the epochinfo structure inside the D object
    D.epochinfo = epochinfo;
    D.save();
    %take bad segments registered in oslview and check if they overlap with the trials. if so, it gives the number of overlapped trials that will be removed later   
    count = 0;
    Bad_trials = zeros(length(trigcor),1);
    for kkk = 1:length(events) %over events
        if strcmp(events(kkk).type,'artefact_OSL')
            for k = 1:length(trl_sec) %over trials
                if events(kkk).time - trl_sec(k,2) < 0 %if end of trial is > than beginning of artifact
                    if trl_sec(k,1) < (events(kkk).time + events(kkk).duration) %if beginning of trial is < than end of artifact
                        Bad_trials(k,1) = 1; %it is a bad trial (stored here)
                        count = count + 1;
                    end
                end                  
            end
        end
    end
    %if bad trials were detected, their indices are stored within D.badtrials field
    if count == 0
        disp('there are no bad trials marked in oslview for');
        disp(spm_files_recog_basen(ii));
    else
        D = badtrials(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
        epochinfo = D.epochinfo;
        xcv = find(Bad_trials == 1);
        for jhk = 1:length(xcv)
            D = D.conditions(xcv(jhk),'Bad');
            epochinfo.conditionlabels(xcv(jhk)) = {'Bad'};
        end
        D.epochinfo = epochinfo;
        D.save(); %saving on disk
        disp('bad trials are ')
        D.badtrials
    end    
end
%define conditions - only 1 epoch for each old/new excerpt (baseline = (-)100ms)
basnam_con_oc = {'Old_Correct'}; %;'Old_Correct_2';'Old_Correct_3';'Old_Correct_4';'Old_Correct_5'};
basnam_con_nc = {'New_Correct'}; %'New_Correct_2';'New_Correct_3';'New_Correct_4';'New_Correct_5'};
basnam_con_ouc = {'Old_Uncorrect'}; %;'Old_Uncorrect_2';'Old_Uncorrect_3';'Old_Uncorrect_4';'Old_Uncorrect_5'};
basnam_con_nuc = {'New_Uncorrect'}; %;'New_Uncorrect_2';'New_Uncorrect_3';'New_Uncorrect_4';'New_Uncorrect_5'};
xlsx_dir_behav = '/scratch7/MINDLAB2017_MEG-LearningBach/BehavioralResponses'; %dir to MEG behavioral results
for ii = 3:3:length(spm_files_recog) %only over the data files I am interested in
    spm_files_recog{ii}=[workingdir '/e80dff' spm_files_recog_basen{ii}];
    D = spm_eeg_load(spm_files_recog{ii});
    [~,~,raw_recog] = xlsread([xlsx_dir_behav '/' xlsx_basenames{ii}]);
    for k = 1:length(D.trialonset)
        if strcmp(raw_recog{(k + 1),3},'No_response') %if there was no response
            D = D.conditions(k,'No_response');
        elseif raw_recog{(k + 1),2}(1) == raw_recog{(k + 1),3}(1) %if the response was correct
            if raw_recog{(k + 1),2}(1) == 'O' %if the response was old
                D = D.conditions(k,basnam_con_oc); %assign old correct
            else
                D = D.conditions(k,basnam_con_nc); %otherwise assign new correct
            end
        else %else the response was wrong
            if raw_recog{(k + 1),2}(1) == 'O' %if was old
                D = D.conditions(k,basnam_con_ouc); %assign old uncorrect
            else
                D = D.conditions(k,basnam_con_nuc); %assign new uncorrect
            end
        end
    end
    if ~isempty(D.badtrials) %overwriting badtrials (if any) on condition labels
        BadTrials = D.badtrials;
        for badcount = 1:length(BadTrials)
            D = D.conditions(BadTrials(badcount),'Bad_trial');
        end
    end
    D = D.montage('switch',1);
    D.epochinfo.conditionlabels = D.conditions; %to add for later use in the source reconstruction
    D.save(); %saving data on disk
    disp(num2str(ii))
end

%%% the preprocessing until here was done for all of the analysis steps %%%
%%% the beamforming sections started after the above steps %%%

%%% the indipendent univariate tests and successive Monte Carlo simulations (first analysis reported in the paper)
%%% required also the following computations:
%%% -averaging over trials
%%% -combining planar gradiometers

%% averaging and combining planar gradiometers (parallel computing on Aarhus University server)

%averaging
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Francesco/Functions'); %path to LEiDA_MEG_leonardo functions
output_prefix_to_be_set = 'm'; 
% config of the cluster server
clusterconfig('slot', 1);
clusterconfig('long_running', 1); 
list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/fe80ldff*');
for ii = 2:2:length(list) %over files
    %distribute 
    input = [];
    input.D = [list(ii).folder '/' list(ii).name];
    input.prefix = output_prefix_to_be_set;
    jobid = job2cluster(@sensor_average, input);
end
%combining planar gradiometers
list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/mfe80ldff*');
for ii = 2%:2:length(list) %over files
    input = [];
    input.D = [list(ii).folder '/' list(ii).name];
    D = spm_eeg_load(input.D);
    D = D.montage('switch',1);
    D.save();
    jobid = job2cluster(@combining_planar_cluster, input);
end

%% *** ANALYSIS ***

%%

%%

%The following section calls "MEG_sensors_plotting_ttest_LBPD_D", a
%function that is designed with 2 scopes:
%1)extracting data from SPM objects and calculating t-tests between experimental conditions
%2)plotting data (potentially taking into account also some statistical
%results to be overlayed in the plots).
%In this case, it was used to extract MEG sensor data and compute t-tests for eahc time-point and MEG sensor.

%% extracting MEG sensor data

%original settings to build the path to the data
subjnum = {'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'38';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';'51';'52';'53';'54';'55';'56';'57';'58';'59';'60';'61';'62';'63';'64';'65';'66';'67';'68';'69';'70';'71'}; %Everybody!!!!!!
datadir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc'; %path to data
outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach'; %path to write output in
S = [];
%computing data
S.outdir = outdir;
S.data = [];
S.spm_list = cell(1,length(subjnum));
for ii = 1:length(subjnum)
    S.spm_list(ii) = {[datadir '/Pme80dffspmeeg_SUBJ00' subjnum{ii} 'recogminor_tsssdsm.mat']};
end
S.conditions = {'Old_Correct','New_Correct'};
S.timeextract = []; %time-points to be extracted
S.save_data = 1; %only meaningfull if you read data from SPM objects saved on disk
S.save_name_data = 'sensor_data';
%plotting settings (WE SIMPLY PROVIDE THE INFORMATION THAT WE DO NOT WANT PLOTS AT THIS TIME)
S.waveform_singlechannels_label = 0; %1 to plot single channel waveforms
S.wave_plot_conditions_together = 0; %1 for plotting the average of all
S.mag_lab = 2; %1 for magnetometers; 2 for gradiometers
S.x_lim_temp_wave = []; %limits for time (in secs) (E.g. [-0.1 3.4])
S.y_lim_ampl_wave = []; %limit for amplitude (E.g. [0 120] magnetometes, [0 6] gradiometers)
%averaged waveform plotting
S.waveform_average_label = 0;
S.legc = 0; %set 1 for legend
S.left_mag = 2:2:204;
S.signtp = {[]};
S.sr = 150; %sampling rate (Hz)
S.avewave_contrast = 0; %1 to plot the contrast between conditions (averaged waveform)
S.save_label_waveaverage = 0;
S.label_plot = 'block_minor';
%t-tests
S.t_test_for_permutations = 1; %here we compute t-tests between conditions (Bach's original vs Bach variation) for each time-point and each MEG channel
S.cond_ttests_tobeplotted_topoplot = [1 2]; %this is for both topoplot and t-tests!! (here [1 2] means cond1 vs cond2!!!!!!!)
%topoplotting
S.topoplot_label = 0;
S.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External'; %this can be loaded from the folder with the provided functions
S.topocontr = 0;
S.topocondsing = [2];
S.xlim = [0.547 0.913]; %time topolot (cluster II)
S.zlimmag = [-4.2 4.2]; %magnetometers amplitude topoplot limits
S.zlimgrad = [-3.5 3.5]; %gradiometers amplitude topoplot limits
S.colormap_spec = 0;
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
S.colormap_spec_x = x;
S.topoplot_save_label = 0;

[out] = MEG_sensors_plotting_ttest_LBPD_D(S);

%% loading t-test results and reshaping them

%input information
grad = 1; %1 for gradiometers; 0 for magnetometers; time-points to be considered for MCS
contrn = 1; %set1 to load previously calculated t-tests Bach's original vs variation; set 0 for Bach variation vs original)
p_thresh = 0.01; %for binarising p-values matrices..

%actual computation
%time-points to be selected
if grad == 1
    %gradiometers
    min_time_point = 16; %16 = 0 seconds (first 15 points are pre-stimulus time)
    max_time_point = 399;
else
    %magnetometers on the basis of gradiometers results
    min_time_point = 98;
    max_time_point = 193;
end
clear DATAP2 TSTAT2
outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_sensordata_UnivariateTests_MonteCarloSimulations'; %path where t-test results are stored in
if contrn == 1
    load([outdir '/sensor_data_OUT_Old_Correct_vs_New_Correct.mat']);
else
    load([outdir '/sensor_data_OUT_New_Correct_vs_Old_Correct.mat']);
end
chanlabels = OUT.chanlabels;
time_sel = OUT.time_sel;
%here gradiometers and magnetometers are always extracted but in the
%following steps only the requested channels (either magnetometers or
%gradiometers) are used
DATAP2(:,:,1) = OUT.TSTATP_mag(:,min_time_point:max_time_point);
DATAP2(:,:,2) = OUT.TSTATP_grad(:,min_time_point:max_time_point);
TSTAT2(:,:,1) = OUT.TSTAT_mag(:,min_time_point:max_time_point);
TSTAT2(:,:,2) = OUT.TSTAT_grad(:,min_time_point:max_time_point);
chanlab = chanlabels(1:2:204)';
label = zeros(length(chanlab),1); %channels label
for ii = 1:length(chanlab)
    label(ii,1) = str2double(chanlab{ii,1}(4:end));
end

%individuating positive vs negative t-values
P = DATAP2; %trick for looking only into the positive t-values
P(P < p_thresh) = 1; %binarizing p-values according to threshold
P(P < 1) = 0;
%old (pos) > new (neg)
P(TSTAT2 < 0) = 0; %deleting p-values for negative contrasts
TSTAT_mag_pos = P(:,:,1);
TSTAT_grad = P(:,:,2);
%negative vs positive p-values
P = DATAP2;
P(P < p_thresh) = 1; %binarizing p-values according to threshold
P(P < 1) = 0;
P(TSTAT2 > 0) = 0; %deleting p-values for positive contrasts
TSTAT_mag_neg = P(:,:,1);
%load a 2D approximation in a matrix of the MEG channels location (IT IS CONTAINED IN THE FUNCTIONS FOLDER THAT WE PROVIDED)
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
if grad == 1
    %actual gradiometers data
    S.data(:,:,:,1) = GRAD_data;
    %zeros.. if you do not want to run the function for magnetometers
    S.data(:,:,:,2) = zeros(size(MAG_data_pos,1),size(MAG_data_pos,2),size(MAG_data_pos,3));
    S.data(:,:,:,3) = zeros(size(MAG_data_pos,1),size(MAG_data_pos,2),size(MAG_data_pos,3));
else
    %0s for gradiometers
    S.data(:,:,:,1) = zeros(size(MAG_data_pos,1),size(MAG_data_pos,2),size(MAG_data_pos,3));
    %actual magnetometers data
    S.data(:,:,:,2) = MAG_data_pos;
    S.data(:,:,:,3) = MAG_data_neg;
end
S.sensortype = [];
S.MEGlayout = cell2mat(raw_channels);
S.permut = 1000;
S.clustmax = 1;
S.permthresh = 0.001;

[MAG_clust_pos, MAG_clust_neg, GRAD_clust] = MEG_sensors_MonteCarlosim_LBPD_D(S);

%getting mean t-value over significant MEG channels and time-points
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/chansgrad.mat'); %loading combined planard gradiometers labels
t_vals = zeros(length(GRAD_clust{1,3}),1);
for ii = 1:length(GRAD_clust{1,3}) %over significant MEG channels
    id = find(strcmp(GRAD_clust{1,3}(ii,1),chansgrad)); %looking for index in TSTAT2 of the significant MEG channels, by using their labels
    timeid = GRAD_clust{1,3}{ii,2}; %extracting significant time-points
    t_vals(ii,1) = mean(TSTAT2(id,timeid,2)); %mean over time-points
end
disp('mean t-val of the significant cluster is ') %only one large significant cluster was detected, showing stronger brain activity for the recognition of 'memorized' vs 'novel' temporal patterns
mean(t_vals)

%% Computing power spectra

freqq = 0; %1 for 1-10Hz; 0 for 1-74Hz

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/BrainstormDemo_CodesTimeFrequencyDecomposition/functions') %path to Dimitrios function for Morlet transform
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_sensordata_UnivariateTests_MonteCarloSimulations/sensor_data.mat');
mag = 1:204;
if freqq == 1
    f = 1:0.1:10;
else
    f = 1:1:30; %frequencies..
end
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
%plotting mean conditions
ges = 0; %maximum limit.. 0 for automatic limits calculation
%mean over subjects (for plotting purposes) and removing pre-stimulus time
Pold = mean(mean(PI,4),3);
%plotting old and new together
figure
% imagesc(time_sel(16:end),f,Pold)
imagesc(time_sel,f,Pold)
set(gca,'YDir','normal') %plotting in descending order
xlabel('time (s)'); ylabel('f (Hz)');
colorbar
if ges ~= 0
    caxis([0 ges])
end
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))
set(gcf,'color','w')

%%

%%

%% *** SOURCE RECONSTRUCTION (BEAMFORMING AND STATISTICS) - BRAIN ACTIVITY ***

%Please remember that more information on the beamforming algorithm and subsequent statistics can be
%found in the OSL website/documentation: %https://ohba-analysis.github.io/osl-docs/

%Please, note that the source reconstruction has been run for both 0.1-1Hz and 2-8Hz.
%Please, also note that for 2-8Hz, the source reconstruction algorithms have been run twice:
% 1) with the absolute
%value to calculate more standard GLMs as suggested by OSL community (N100_nocoape = 0)
% 2) without absolute values so keeping the double polarity of
%the magnetic field (this was extremely important for the k-means
%clustering algorithm shown later in this script) (N100_nocoape = 1);

lab01 = 1; % 1 for 0.1-1Hz; 0 for 2-8Hz

if lab01 == 1
    % 0.1-1 Hz
    count1 = 0;
    oat = [];
    %SOURCE RECON
    for ii = 3:3:213  % = 1:length(spm_files_recog_basen)
        count1 = count1 + 1;
        spm_files = ['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/e80dff'  spm_files_recog_basen{ii}]; %check this in order not to mix up 'epoched' and 'continuous'
        %     D_epoched = spm_eeg_load(spm_files);
        %     D_epoched = D_epoched.montage('switch',1); %switch the montage to 1 in order to be safer
        %     D_epoched.save();
        processed_file_epoched = spm_files;
        spm_files = ['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/dff'  spm_files_recog_basen{ii}]; %check this in order not to mix up 'epoched' and 'continuous'
        %     D_continuous = spm_eeg_load(spm_files);
        %     D_continuous = D_continuous.montage('switch',1); %switch the montage to 1 in order to be safe
        %     D_continuous.save();
        processed_file_continuous = spm_files;
        oat.source_recon.D_continuous(count1) = {processed_file_continuous};%the file just after AFRICA and before the epoching
        oat.source_recon.D_epoched(count1) = {processed_file_epoched};
        oat.source_recon.session_names(count1) = {['session' num2str(count1)]};
        oat.source_recon.results_fnames(count1) = {['session' num2str(count1) '_recon']};
        disp(num2str(count1))
    end
    pca_dim_1 = 50;
    oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_continuous),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects
    oat.source_recon.sessions_to_do = [];
    oat.source_recon.sessions_to_do = [1:38 40:71]; %sessions to do among the file_list (subject 39 excluded because of problems during data collection..)
    oat.source_recon.modalities = {'MEGMAG'; 'MEGPLANAR'}; %'MEGMAG'};
    oat.source_recon.conditions = {'Old_Correct';'New_Correct'};
    oat.source_recon.gridstep = 8; % in mm
    oat.source_recon.time_range = [-0.1 3.4]; % time range in secs
    oat.source_recon.freq_range = [0.1 1]; % frequency range in Hz (alpha large range)
    oat.source_recon.type = 'Scalar';
    oat.source_recon.method = 'beamform';
    oat.source_recon.normalise_method = 'mean_eig';
    oat.source_recon.forward_meg = 'MEG Local Spheres';
    oat.source_recon.report.do_source_variance_maps = 1;
    oat.source_recon.dirname = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_01_1Hz_AllTrials_BC_2.oat.oat'; % spm_files_recog_basen{ii}(12:25)];
else
    %2-8 Hz
    count1 = 0;
    oat = [];
    %SOURCE RECON
    for ii = 3:3:213  % = 1:length(spm_files_recog_basen)
        count1 = count1 + 1;
        spm_files = ['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/e80dff'  spm_files_recog_basen{ii}]; %check this in order not to mix up 'epoched' and 'continuous'
        %     D_epoched = spm_eeg_load(spm_files);
        %     D_epoched = D_epoched.montage('switch',1); %switch the montage to 1 in order to be safer
        %     D_epoched.save();
        processed_file_epoched = spm_files;
        spm_files = ['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/dff'  spm_files_recog_basen{ii}]; %check this in order not to mix up 'epoched' and 'continuous'
        %     D_continuous = spm_eeg_load(spm_files);
        %     D_continuous = D_continuous.montage('switch',1); %switch the montage to 1 in order to be safe
        %     D_continuous.save();
        processed_file_continuous = spm_files;
        oat.source_recon.D_continuous(count1) = {processed_file_continuous};%the file just after AFRICA and before the epoching
        oat.source_recon.D_epoched(count1) = {processed_file_epoched};
        oat.source_recon.session_names(count1) = {['session' num2str(count1)]};
        oat.source_recon.results_fnames(count1) = {['session' num2str(count1) '_recon']};
        disp(num2str(count1))
    end
    pca_dim_1 = 50;
    oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_continuous),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects
    oat.source_recon.sessions_to_do = [];
    oat.source_recon.sessions_to_do = [1:38 40:71]; %sessions to do among the file_list (subject 39 excluded because of problems during data collection..)
    oat.source_recon.modalities = {'MEGMAG'; 'MEGPLANAR'}; %'MEGMAG'};
    oat.source_recon.conditions = {'Old_Correct';'New_Correct'};
    oat.source_recon.gridstep = 8; % in mm
    oat.source_recon.time_range = [-0.1 3.4]; % time range in secs
    oat.source_recon.freq_range = [2 8]; % frequency range in Hz (alpha large range)
    oat.source_recon.type = 'Scalar';
    oat.source_recon.method = 'beamform';
    oat.source_recon.normalise_method = 'mean_eig';
    oat.source_recon.forward_meg = 'MEG Local Spheres';
    oat.source_recon.report.do_source_variance_maps = 1;
    if N100_nocoape == 1
        %no acope nor coape, so keeping the double polarity (for k-means clustering)
        oat.source_recon.dirname = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_2_8Hz_AllTrials_BC_noacope_2.oat'; % spm_files_recog_basen{ii}(12:25)];
    else
        %absolute value after beamforming, so loosing double polarity (for GLMs)
        oat.source_recon.dirname = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_Hz2_Hz8_minorfl_forGLM.oat'; % spm_files_recog_basen{ii}(12:25)];
    end
end

%%
oat = osl_check_oat(oat);
oat.to_do = [1 0 0 0];
oat = osl_run_oat(oat);

%% FIRST LEVEL

%setting parallel computing
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
clusterconfig('slot', 2); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
% clusterconfig('scheduler', 'none'); % set automatically the long run queue
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%actual first level
lab01 = 0; % 1 for 0.1-1Hz; 0 for 2-8Hz
N100_nocoape = 0; % 1 for no coape; 0 for absolute value after beamforming

%loading a previously computed oat
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_2_8Hz_AllTrials_BC_noacope_2.oat/oat_wholebrain_first_level_BC_minor_group_level.mat');

%FIRST LEVEL
design_matrix_summary = {};
design_matrix_summary{1} = [1 0];design_matrix_summary{2} = [0 1]; %design_matrix_summary{3}=[0 0 1 0];design_matrix_summary{4}=[0 0 0 1];
oat.first_level.design_matrix_summary = design_matrix_summary;
% contrasts to be calculated:
oat.first_level.contrast = {};
% contrast design matrix
oat.first_level.contrast{1} = [1 0]'; % old 
oat.first_level.contrast{2} = [0 1]'; % new 
oat.first_level.contrast{3} = [1 -1]'; % old-new
oat.first_level.contrast_name = {};
oat.first_level.contrast_name{1} = 'old';
oat.first_level.contrast_name{2} = 'new';
oat.first_level.contrast_name{3} = 'old-new';
%original contrasts, try1 and try2 seem to give the exact same results
oat.first_level.report.first_level_cons_to_do = [3 2 1]; %[1 2 3]; %better to do 3 2 1 in order to get the information for the peak value of the contrast old-new
oat.first_level.time_range = [-0.1 3.4];
oat.first_level.post_tf_downsample_factor = 1;
if lab01 == 1
    oat.first_level.name = ['wholebrain_first_level_BC_minor']; %REMEMBER TO CHECK THIS NAME!!
    oat.first_level.bc = [1 1 0];
    oat.first_level.cope_type = 'coape'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
    %to add if the oat has not been automatically saved
    for ii = 1:71
        oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level_BC_minor.mat']};
    end
else
    if N100_nocoape == 1
        oat.first_level.name = ['wholebrain_first_level_BC']; %REMEMBER TO CHECK THIS NAME!!
        oat.first_level.cope_type = 'none'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
        oat.first_level.bc = [1 1 0];
        for ii = 1:71
            oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level_BC.mat']};
        end
    else
        oat.first_level.name = ['wholebrain_first_level_coape_2021']; %REMEMBER TO CHECK THIS NAME!!
        oat.first_level.cope_type = 'coape'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
        oat.first_level.bc = [1 1 0];
        for ii = 1:71
            oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level_coape_2021.mat']};
        end
    end
end
%running first level on parallel computing
% v = [71];
for ii = 1:71%length(v)
    if ii ~= 39
        oat.first_level.sessions_to_do = [];
        oat.first_level.sessions_to_do = [ii];
        jobid = job2cluster(@cluster_beamfirstlevel,oat);
    end
end

%% SUBJECT LEVEL

lab01 = 0; % 1 for 0.1-1Hz; 0 for 2-8Hz
N100_nocoape = 0; % 1 for no coape; 0 for absolute value after beamforming

%in this case it does not do anything since I have only one experimental
%block for each participant (however, for computational reasons, I have to
%run it)
% this is if you have a perfect correspondance between sessions and subjects
oat.subject_level.session_index_list = cell(1,71); %sarebbe length(subjects)
oat.subject_level.name = 'minor';
% oat.subject_level.subjects_to_do = [71];
oat.subject_level.subjects_to_do = [1:38 40:71];
%to update the names of the results files of subject level
for ii = 1:71
    if ii ~= 39
        oat.subject_level.session_index_list{ii} = ii;
        if lab01 == 1 || N100_nocoape == 1
            oat.subject_level.results_fnames(ii) = {['subject' num2str(ii) '_wholebrain_first_level_BC_minor.mat']};
        else
            oat.subject_level.results_fnames(ii) = {['subject' num2str(ii) '_wholebrain_first_level_coape_2021_minor.mat']};
        end
    end
end
%running on parallel computing
for ii = 1:71
    if ii ~= 39
        oat.subject_level.subjects_to_do = [];
        oat.subject_level.subjects_to_do = ii;
        jobid = job2cluster(@cluster_beamsubjlevel,oat);
    end
end

%% GROUP LEVEL

lab01 = 0; % 1 for 0.1-1Hz; 0 for 2-8Hz
N100_nocoape = 0; % 1 for no coape; 0 for absolute value after beamforming

%actual group level
oat.group_level = [];
oat.subject_level.subjects_to_do = [1:38 40:71];
oat.group_level.subjects_to_do = [1:70]; %this seems to be necessary since the function gets the indices for the group level from each progressive number of oat.subject_level.subjects_to_do (not very clever idea by the way..). Therefore not to mix up things, you should have here a vector from 1 to the total number of subjects and then specify the correct subject IDs that you want.. or at least, to me it looks like that..
if lab01 == 1 || N100_nocoape == 1
    oat.group_level.name = 'group_level_everybody_BC'; %OBS!! REMEMBER TO UPDATE THE NAME!
    %results name
    oat.group_level.results_fnames = 'wholebrain_first_level_BC_minor_group_level_everybody_BC';
    oat.group_level.time_range = [-0.1 3.4];
else
    oat.group_level.name = 'group_level_everybody';
    %results name
    oat.group_level.results_fnames = 'wholebrain_first_level_coape_2021_minor_group_level_everybody';
end
% Spatial and temporal averaging options
oat.group_level.space_average = 0;
oat.group_level.time_average = 0;
oat.group_level.time_smooth_std = 0; % secs
oat.group_level.use_tstat = 0;
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
% oat.group_level.first_level_contrasts_to_do = [1,2,3,4,5,6,7,8]; % list of first level contrasts to run the group analysis on
oat.group_level.first_level_contrasts_to_do = [1,2,3]; % list of first level contrasts to run the group analysis on
oat.group_level.report.first_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.show_lower_level_copes = 0;
oat.group_level.report.show_lower_level_cope_maps = 0;

jobid = job2cluster(@cluster_beamgrouplevel,oat);

%% PRODUCING NIFTI IMAGES FOR GROP LEVEL ANALYSIS

%The following lines allow you to create nifti files of the group level analysis that can be visualize in programs such as fslview.
S2 = [];
S2.oat = oat;
S2.stats_fname = oat.group_level.results_fnames; %group level
S2.first_level_contrasts = [3,1,2]; %remember that in this way you define the order of the output (this numbers refers to the order (numbers) defined in the contrasts; for example here tstat3 refers to the contrast old-new)
S2.group_level_contrasts = [1];
S2.resamp_gridstep = oat.source_recon.gridstep;
[statsdir,times,count] = oat_save_nii_stats(S2);

%% NIFTI IMAGES SUBJECT LEVEL

%% cluster settings (parallel computing)

%settings for cluster (parallel computing)
% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%%

lab01 = 2; % 1 = 0.1-1Hz; 2 = 2-8Hz

if lab01 == 1
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_01_1Hz_AllTrials_BC_2.oat/oat_wholebrain_first_level_BC_minor_group_level_everybody_BC.mat');
    oat.source_recon.dirname = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_01_1Hz_AllTrials_BC_2.oat';
else
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_2_8Hz_AllTrials_BC_noacope_2.oat/oat_wholebrain_first_level_BC_minor_group_level.mat');
    oat.source_recon.dirname = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_2_8Hz_AllTrials_BC_noacope_2.oat';
end
%actual computation
for ii = 2:71
    if ii ~= 39
        S2 = [];
        S2.oat = oat;
        S2.stats_fname = oat.subject_level.results_fnames{ii};
        S2.first_level_contrasts = [1,2,3]; %remember that in this way you define the order of the output (this numbers refers to the order (numbers) defined in the contrasts; for example here tstat3 refers to the contrast old-new)
        S2.resamp_gridstep = oat.source_recon.gridstep;
        jobid = job2cluster(@cluster_oat_save_nii_stats,S2);
%         [statsdir,times,count] = oat_save_nii_stats(S2);
        disp(ii)
    end
end

%% CLUSTER-BASED MCS

%% Defining clusters on 3D brain voxel statistics averaged over the 5 time-windows corresponding to the 5 tones forming the musical patterns

freq = 2; %1 = 0.1-1Hz; 2 = 2-8Hz
tvalbin = 2.3; %value to binarize the data with

for tt = 1:5 %over tones (objects of the temporal pattern)
    %loading p-values and t-values
    clear DATA
    tone = tt; %set to 1 2 3 4 or 5
    FREQ{1} = 'source_localsphere_01_1Hz_AllTrials_BC_2.oat'; FREQ{2} = 'source_localsphere_2_8Hz_AllTrials_BC_noacope_2.oat';
    if freq == 1
        T = load_nii(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/' FREQ{freq} '/wholebrain_first_level_BC_minor_group_level_everybody_BC_randomise_c3_dir_rev_1_tone_' num2str(tone) '/tstat3_gc1_8mm.nii.gz']);
    else
        T = load_nii(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/' FREQ{freq} '/wholebrain_first_level_coape_2021_minor_group_level_everybody_randomise_c3_dir_rev_1_tone_' num2str(tone) '/tstat3_gc1_8mm.nii.gz']);
    end
    %extracting matrix with statistics
    % P2 = P.img;
    T2 = T.img;
    
    %mask for non-0 voxels in brain imges (basically building a layout for actual brain voxels)
    mask = zeros(size(T2,1),size(T2,2),size(T2,3));
    mask(T2~=0) = 1; %assigning 1 when you have real brain voxels
    
    %removing non-significant results
    %preparing data
    data = T2;
    data(data==0) = NaN; %removing non-brain voxels
    data(T2>tvalbin) = 1; %
    data(T2<-tvalbin) = 1; %
    data(data~=1) = NaN; %assigning 1 to significant voxels
    % data(~isnan(data)) = 1; %assigning 1 to significant voxels
    DATA{1} = data; %storing data
    
    %getting MNI coordinates
    %OBS! you may get a warning since the skeletonized image is not exactly in MNI space, but close enough
    [ mni_coords, xform ] = osl_mnimask2mnicoords('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');
    
    %preparation of information and data for the actual function
    for ii = 1%:2 %over directions of the contrast (cond1>cond2 and cond2>cond1)
        S = [];
        if freq == 1
            S.T = ['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/' FREQ{freq} '/wholebrain_first_level_BC_minor_group_level_everybody_BC_randomise_c3_dir_rev_1_tone_' num2str(tone) '/tstat3_gc1_8mm.nii.gz'];
        else
            S.T = ['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/' FREQ{freq} '/wholebrain_first_level_coape_2021_minor_group_level_everybody_randomise_c3_dir_rev_1_tone_' num2str(tone) '/tstat3_gc1_8mm.nii.gz'];
        end
        S.outdir = ['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Cluster_Based_MCS/Freq_' num2str(freq)]; %output path
        S.parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz';
        S.labels = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat';
        S.MNIcoords = mni_coords; %MNI coordinates of 8mm MNI152T1 brain
        S.data = DATA{ii};
        S.mask = mask; %mask of the brain layout you have your results in
        S.permut = 1000; %number of permutations for Monte Carlo simulation
        S.clustmax = 1; %set 1 for only max cluster size of each permutation MCS (more strict); set 0 for every size of each cluster detected for each permutation MCS (less strict).
        S.permthresh = 0.001; %threhsold for MCS
        S.anal_name = ['Tval_tone_' num2str(tone) '_Cond1vsCond2']; %name for the analysis (used to identify and save image and results)
        
        %actual function
        PP = BrainSources_MonteCarlosim_3D_LBPD_D(S);
    end
end

%% HERE YOU NEED TO COMBINE IMAGE WITH MORE THAN ONE CLUSTERS!

%% Combining images (manual work in this case..)

% path1 = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Cluster_Based_MCS/Freq_1';
path2 = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Cluster_Based_MCS/Freq_2';
path = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Cluster_Based_MCS';

%list of images to be combined (1)
name1 = ['Tval_tone_1_Cond1vsCond2_SignClust_1_Tvals.nii.gz'];
% name1 = ['Tval_tone_2_Cond1vsCond2_SignClust_1_Tvals.nii.gz'];

%list of images to be combined (2)
name2 = ['Tval_tone_1_Cond1vsCond2_SignClust_2_Tvals.nii.gz'];
% name2 = ['Tval_tone_2_Cond1vsCond2_SignClust_2_Tvals.nii.gz'];

%list of images to be combined (3)
% name3 = ['Tval_tone_1_Cond1vsCond2_SignClust_3_Tvals.nii.gz'];

%list of output names
% output = ['Tone1_01_1Hz.nii.gz'];
output = ['Tone1_2_8Hz.nii.gz'];
% output = ['Tone2_2_8Hz.nii.gz'];

%actual command
% cmd = ['fslmaths ' path1 '/' name1 ' -add ' path1 '/' name2 ' ' path '/' output]; %OBS! be careful with the spacing
cmd = ['fslmaths ' path2 '/' name1 ' -add ' path2 '/' name2 ' ' path '/' output]; %OBS! be careful with the spacing
% cmd = ['fslmaths ' path2 '/' name1 ' -add ' path2 '/' name2 ' -add ' path2 '/' name3 ' ' path '/' output]; %OBS! be careful with the spacing
system(cmd)

%% creating a file with information on significant clusters/voxels outputted by the permutation test (THINK TO MAKE THIS A FUNCTION)

lab01 = 2; % 1 for 0.1-1Hz; 0 for 2-8Hz 

%current lines to get the proper images..
if lab01 == 1
    outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Cluster_Based_MCS/Freq1_Combined';
    path = dir([outdir '/Tone*']);
    freqq = '01_1';
else
    outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Cluster_Based_MCS/Freq2_Combined';
    path = dir([outdir '/Tone*']);
    freqq = '2_8';
end
for cc = 3%1:3 %over contrasts
    for tt = 1:5 %over musical tones
        %actual name plus path
        %current (only for contrast 3)
        fname = [path(tt).folder '/' path(tt).name];
        %previous (for all contrast; please, note that contrasts 1 and 2 were not changed since they were nearly equal and was pointless to do that here..
%         fname = [path(1).folder '/' path(1).name(1:15) num2str(cc) '_tone_' num2str(tt) '.nii.gz'];
        %getting MNI coordinates of significant voxels within the provided image
        [ mni_coords, xform ] = osl_mnimask2mnicoords(fname);
        %loading the image
        V = nii.load(fname);
        %extracting statistics
        VV = V(V~=0);
        %indices of non-zero values of nifti image
        VI = find(V~=0);
        %path to AAL template
        parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'; %load this from the provided codes folder
        %loading AAL labels
        load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat'); %load this from the provided codes folder
        %extracting AAL coordinates information
        K = nii.load(parcelfile);
        %sorting results in order to have strongest voxels at the top (positive t-values) or at the bottom (negative t-values)
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
        writetable(PDn,[outdir '/Freq_' freqq '_contr_' num2str(cc) '_tone_' num2str(tt) '.xlsx'],'Sheet',1) %printing excel file
    end
end

%%

%% *** FUNCTIONAL-SPATIAL K-MEANS CLUSTERING ***

%% Summary of the procedures related to the functional-spatial k-means clustering

% 1) STATISTICAL COMPARISON (SUBJECT-LEVEL) OF MEMORIZED VS NOVEL TEMPORAL PATTERNS USING K-MEANS CLUSTERING
%    i - loading (old school) group level statistics from "standard" GLM for memorized (old) and novel (new) temporal patterns
%    ii - fixing their sign if needed (this is only for 2-8Hz)
%    iii - averaging them together (to get brain areas that can reasonably be compared)
%    iv - running functional-spatial k-means clustering on such data
%    v - running again group-level statistics (from source-space subject-level), but only on the new parcels outputted by point (iv)

% 2) DEFINING BEST FUNCTIONAL BRAIN PARCELS AND TIMESERIES FOR MEMORIZED (OLD) TEMPORAL PATTERNS
%    This allows me to get the best approximation of the proper brain areas (and their activity) involved in the actual recognition of previously learned auditory patterns
%    i - loading (old school) group level statistics from "standard" GLM for memorized (old) and novel (new, only for 2-8Hz) temporal patterns
%    ii - fixing their sign if needed (this is only for 2-8Hz)
%    iii - running functional-spatial k-means clustering on such data
%    iv - fitting equations to the different brain parcels timeseris (this is done in python)

%% 1) MEMORIZED (OLD) (AND/OR NOVEL (NEW)) TEMPORAL PATTERNS USING K-MEANS CLUSTERING

%user inputs
slowneg = 0; %1 for slow negativity (0.1-1Hz); 0 for N100 (2-8Hz)
average_l = 0;  %1 for averaging together OLD and NEW; 0 for not averaging
old_l = 1; %1 for old (Bach's original); 2 for new (Bach variation). THIS IS MEANINGFUL ONLY IF average_l = 0, meaning that you actually choose between OLD and NEW and not average them together
brain_inv_each_parcel = 0; %set 1 to plot in brain templates each parcel independently
method_figure = 0; %1 to get the figure for the Method Figure of the paper (in this case, the previous three variables/labels become non-relevant)

%loading and defining inputs
if method_figure == 1
    slowneg = 1; %1 for slow negativity (0.1-1Hz); 0 for N100 (2-8Hz)
    average_l = 0;  %1 for averaging together OLD and NEW; 0 for not averaging
    old_l = 1;
end
VS = cell(1,2);
for cc = 1:2 %over conditions (old and new)
    %loading data
    %path to file nifti
    if slowneg == 1
        %0.1-1Hz (slow negativity)
        pathnii = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_01_1Hz_AllTrials_BC_2.oat/wholebrain_first_level_BC_minor_group_level_everybody_BC_dir';
    else
        %2 - 8Hz (N100)
        pathnii = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_2_8Hz_AllTrials_BC_noacope_2.oat/wholebrain_first_level_BC_minor_group_level_everybody_BC_dir'; %baseline corrected and nocope
    end
    fname = [pathnii '/tstat' num2str(cc) '_gc1_8mm.nii.gz']; %tstat name
    %loading the image
    V_stat4D = nii.load(fname);
    %removing zeros from the statistics and reshaping it.. here the idea is passing from the representation in the "brain layout" of the statistics into a 2D representation of the timeseries of the voxels..
    Vstat2 = reshape(V_stat4D(V_stat4D~=0),[3559,size(V_stat4D,4)]); %this is a correct representation of voxels (rows) and their values over time (columns)
    SS = size(Vstat2);
    % reversing the sign
    if slowneg == 0 %if N100 data, reversing the sign in oderd not to have ambiguity (i.e. to have negativity and positivity always consistent)
        Vdum = zeros(SS(1),SS(2)); %preallocating space
        vect = zeros(SS(1),1); %preallocating space
        for ii = 1:SS(1) %over brain voxels
            if mean(Vstat2(ii,30:40),2) < 0 %if the data in voxel ii is negative during N100 time
                Vdum(ii,:) = Vstat2(ii,:); %keeping the original data
                vect(ii,1) = 1; %storing a vector with 1 and -1 to be used for later statistics
            else
                Vdum(ii,:) = Vstat2(ii,:) .* (-1); %otherwise reversing..
                vect(ii,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
            end
        end
        Vstat2 = Vdum;
        if average_l == 1
            %saving vect with 1s and -1s only if we work on averaged 2-8Hz data (since this data will then be used for subject level statistics)
            save(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_' num2str(slowneg) '_Average_l_' num2str(average_l) '_Old_l_1/onesvector_negativereverse.mat'],'vect');
        end
    end
    VS{cc} = Vstat2;
end
%options for producing the nifti image in the k-means function
options = [];
options.mask_fname = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/wholebrain_first_level_minor_group_level_everybody_mask.nii.gz';
%previous settings for the OSL functions that I used to employ (no longer useful)
% options.interp = 'cubic';
% options.output_spat_res = 8;
% options.tres = 0.0067;
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/time.mat'); %loading time in seconds
%path to MNI template to get proper coordinates
caz = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz';
%getting MNI coordinates of significant voxels within the provided image
[ mni_coords_template, xform_template ] = osl_mnimask2mnicoords(caz);
if average_l == 1
    %averaging together OLD and NEW (just if needed..).. in this case you need first to load independently OLD and NEW and assign: OLD = Vstat2 and then NEW = Vstat2
    Vstat2 = zeros(SS(1),SS(2));
    for ii = 1:SS(1)
        for jj = 1:SS(2)
            Vstat2(ii,jj) = (VS{1}(ii,jj) + VS{2}(ii,jj)) / 2;
        end
    end
else %otherwise just assigning the OLD or NEW, depending on user's request
    Vstat2 = VS{old_l};
end
%additional settings for user (definitely less crucial than the ones reported above)
S = [];
if method_figure == 1
    plotl = 1; %set 1 for all plots but SSDs.. 0 otherwise
else
    plotl = 0; %set 1 for all plots but SSDs.. 0 otherwise
end
S.ROItimescol = 0; %0 for different colors; 1 for different shades (hues) of red
S.plot.brain_inv_each_parcel = brain_inv_each_parcel; %set 1 to plot in brain templates each parcel independently
%actual computation of the clustering
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Papers/GaussianBach/Functions');
if method_figure == 1
    %here and in another parts of this section I changes the parameters and did several runs of the clustering function to get different data and images that I needed for the Method Figure of the paper
    S.clustsol = [6]; %specify number of temporal clusters for different cluster solutions for kmeans-clustering (e.g. 3:5 means that you want to test cluster solutions with 3, 4 and 5 clusters)
    S.clspatvect = [1 1 4 1 1 1]; %ones(1,20)*2; %specify different spatial cluster solutions for different temporal clusters (NOTE that this makes sense if clustsol has only 1 cluster solution and if that solution is equal to the length of clspatvect!!)
    S.clustspatial = []; %spatial clusters.. set it to 0 if you do not want to run it
    S.percentmax = 100; %value (percentage) of maximum values to be considered for clustering purposes (e.g 100 uses all voxels; 5 only the top 5% of the voxels)
    S.max_value = 0; %set 1 for temporal clusterization on maximum value; set 0 for temporal clusterization on time index of maximum value
    S.timep = [1:391]; %time-samples to do the computation on (usually done for 1:200 or 1:295 or 1:391 (slow negativity))
else
    if slowneg == 1
        S.max_value = 0; %set 1 for temporal clusterization on maximum value; set 0 for temporal clusterization on time index of maximum value
        S.timep = [1:391]; %time-samples to do the computation on (usually done for 1:200 or 1:295 or 1:391 (slow negativity))
        if average_l == 1 %if OLD and NEW together
            S.clustsol = [6]; %specify number of temporal clusters for different cluster solutions for kmeans-clustering (e.g. 3:5 means that you want to test cluster solutions with 3, 4 and 5 clusters)
            S.clspatvect = [1 5 4 3 4 6]; %ones(1,20)*2; %specify different spatial cluster solutions for different temporal clusters (NOTE that this makes sense if clustsol has only 1 cluster solution and if that solution is equal to the length of clspatvect!!)
            S.clustspatial = [2:10]; %spatial clusters.. set it to 0 if you do not want to run it
            S.percentmax = 100; %value (percentage) of maximum values to be considered for clustering purposes (e.g 100 uses all voxels; 5 only the top 5% of the voxels)
        else %good parameters (previously defined) for OLD
            %%%%%%%%% PARAMETERS TO BE CHECKED
            S.clustsol = [6]; %specify number of temporal clusters for different cluster solutions for kmeans-clustering (e.g. 3:5 means that you want to test cluster solutions with 3, 4 and 5 clusters)
            S.clspatvect = [1 4 3 3 3 5]; %specify different spatial cluster solutions for different temporal clusters (NOTE that this makes sense if clustsol has only 1 cluster solution and if that solution is equal to the length of clspatvect!!)
            S.clustspatial = [0]; %spatial clusters.. set it to 0 if you do not want to run it
            S.percentmax = 100; %value (percentage) of maximum values to be considered for clustering purposes (e.g 100 uses all voxels; 5 only the top 5% of the voxels)
        end
    else
        S.max_value = 1; %set 1 for temporal clusterization on maximum value; set 0 for temporal clusterization on time index of maximum value
        S.timep = [1:200];
        if average_l == 1 %if OLD and NEW together
            S.clustsol = [4]; %specify number of temporal clusters for different cluster solutions for kmeans-clustering (e.g. 3:5 means that you want to test cluster solutions with 3, 4 and 5 clusters)
            S.clspatvect = [3 3 2 2]; %ones(1,20)*2; %specify different spatial cluster solutions for different temporal clusters (NOTE that this makes sense if clustsol has only 1 cluster solution and if that solution is equal to the length of clspatvect!!)
            S.clustspatial = [2:10]; %spatial clusters.. set it to 0 if you do not want to run it
            S.percentmax = 20; %value (percentage) of maximum values to be considered for clustering purposes (e.g 100 uses all voxels; 5 only the top 5% of the voxels)
        else %not averaged
            if old_l == 1 %ONLY OLD
                S.clustsol = [4]; %specify number of temporal clusters for different cluster solutions for kmeans-clustering (e.g. 3:5 means that you want to test cluster solutions with 3, 4 and 5 clusters)
                S.clspatvect = [2 3 2 2]; %ones(1,20)*2; %specify different spatial cluster solutions for different temporal clusters (NOTE that this makes sense if clustsol has only 1 cluster solution and if that solution is equal to the length of clspatvect!!)
                S.clustspatial = [2:10]; %spatial clusters.. set it to 0 if you do not want to run it
                S.percentmax = 20; %value (percentage) of maximum values to be considered for clustering purposes (e.g 100 uses all voxels; 5 only the top 5% of the voxels)
            else %ONLY NEW
                S.clustsol = [4]; %specify number of temporal clusters for different cluster solutions for kmeans-clustering (e.g. 3:5 means that you want to test cluster solutions with 3, 4 and 5 clusters)
                S.clspatvect = [2 3 2 2]; %ones(1,20)*2; %specify different spatial cluster solutions for different temporal clusters (NOTE that this makes sense if clustsol has only 1 cluster solution and if that solution is equal to the length of clspatvect!!)
                S.clustspatial = [2:10]; %spatial clusters.. set it to 0 if you do not want to run it
                S.percentmax = 20; %value (percentage) of maximum values to be considered for clustering purposes (e.g 100 uses all voxels; 5 only the top 5% of the voxels)
            end
        end
    end
end
S.data = Vstat2; %data
S.PCA = 0; %1 for MCS-based PCA; 0 for mean
S.PCA_method = 'occurrences'; %'occurrences' or 'max_abs' or 'average'
S.time = time; %time in seconds (its length can be bigger than S.timep, but not bigger!)
%output directory.. one different folder for each of the conditions examined
if method_figure == 1
    S.outdir = ['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_1_S3/SpatioTemporalParcels_TempClust3'];
else
    S.outdir = ['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_' num2str(slowneg) '_Average_l_' num2str(average_l) '_Old_l_' num2str(old_l)];
%     S.outdir = ['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/NII_TOOL']; %path used to test new nifti images creation
end
S.coords_template = mni_coords_template; %3D coordinates of the template used for the brain data (originally I used MNI152 8 mm brain..)
S.options = options; %parameters for plotting results in brain template (through the function: )
%NOT REALLY MEANINGFUL NOW, BUT THEY MAY BE RUN AS WELL
S.randsim = 0; %0 for actual data only; 1 for simulation on random numbers (rand function in Matlab; 2 for randomised order of actual data); FOR NOW WE USE ONLY 0
S.ROItimeseries_rand = 0; %1 to plot timeseries of randomly simulated data.. 0 not to do that
S.rperm = 1; %number of permutation for randomising the data rperm time..
%UNTIL HERE

if plotl == 1 %for now either having all plots or no plots.. just SSDs plots are exceptions..
    S.plot.scatterl = 1; %1 if you want the scatterplots.. 0 otherwise (only for actual data)
    S.plot.viol = 1; %if you set 1 for "scatterl', then you can choose to have violins and scatter plots together (viol = 1) or scatter plots only (viol = 0)
    S.plot.brainv = 0; %1 for plotting in the brain (nifti file).. 0 otherwise (only for actual data)
    S.plot.ROItimeseries = 1; %1 to get timeseries of clusters (for now by computing the mean over the voxels belonging to each different cluster.. THEN CONSIDER TO TRY PCA OR SOMETHING SIMILAR AND TEST HOW DIFFERE AND HOW MUCH DIFFER WITH THE MEAN)
    S.plot.ROItimeseries_rand = 0;%same but for randomised data
else
    S.plot.scatterl = 0; %1 if you want the scatterplots.. 0 otherwise (only for actual data)
    S.plot.viol = 0; %if you set 1 for "scatterl', then you can choose to have violins and scatter plots together (viol = 1) or scatter plots only (viol = 0)
    S.plot.brainv = 0; %1 for plotting in the brain (nifti file).. 0 otherwise (only for actual data)
    S.plot.ROItimeseries = 0; %1 to get timeseries of clusters (for now by computing the mean over the voxels belonging to each different cluster.. THEN CONSIDER TO TRY PCA OR SOMETHING SIMILAR AND TEST HOW DIFFERE AND HOW MUCH DIFFER WITH THE MEAN)
    S.plot.ROItimeseries_rand = 0;%same but for randomised data
end
%actual function
[idk3,kkk2,KJK] = FunctionalSpatialClustering_voxels2ROIs_LBPD_D(S);

%additional plotting for Method Figure (NOT really meaningful for understanding/running the clustering algorithm)
if method_figure ~= 1
    %saving outputted clusters
    save([S.outdir '/Clusters.mat'],'idk3','kkk2','KJK','S')
else
    for pp = 6:9
        data_m(1) = {Vstat2};
        %actual plotting
        S2 = [];
        S2.data = data_m; %subsets of dataset (or different datasets) provided in different cells
        %settings equal to the ones used for the function above..
        S2.clustsol = S.clustsol;
        S2.clspatvect = S.clspatvect;
        S2.time = S.time;
        S2.timep = 1:391;
        %getting only parcels with ID specified in vector v
        v = [pp];
        kkf = zeros(length(kkk2),1);
        for ii = 1:length(v)
            kkf(kkk2==v(ii)) = kkk2(kkk2==v(ii));
        end
        %until here
        %parcellations
        S2.idk3 = []; %temporal clustering parcellation
        S2.kkk2 = kkf; %functional-spatial clustering parcellation
        if slowneg == 1
            S2.ylim = [-1.5 4];
        else
            S2.ylim = [-25 20];
        end
        S2.dimc = [1.5]; %linewidth for temporal and functional-spatial clustering timeseries
        % S2.dimc = [1.5];
        % S2.col_ind = {[1 2],[3 4]}; %different datasets same colors.. leave empty [] for not doing that
        % S2.col_ind = {[1],[2]}; %different datasets same colors.. leave empty [] for not doing that
        S2.col_ind = [];
        S2.color = {'r'};
        S2.legs = 1;
        %actual plotting function
        STC_plottingtimeseries_LBPD(S2)
    end
    %path to MNI template to get proper coordinates
    caz = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz';
    %getting MNI coordinates of significant voxels within the provided image
    [ mni_coords_template, xform_template ] = osl_mnimask2mnicoords(caz);
    %loading amplitude values of parcels to create examples for Method Figure
    fname = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_1_S3/SpatioTemporalParcels_TempClust3/tempclust_k_6_spatialclust_clspatvect_time_1_391_actualdata.nii.gz';
    %loading the image
    V_parc = nii.load(fname);
    V_parc2 = reshape(V_parc(V_parc~=0),[3559,1]); %this is a correct representation of voxels (rows) and their values over time (columns)
    %getting coordinates of 3 parcels
    as = mni_coords_template((V_parc2==6),:);
    bs = mni_coords_template((V_parc2==7),:);
    cs = mni_coords_template((V_parc2==8),:);
    %actual plotting (3 parcels separate colours)
    LWW = 0.1;
    figure
    scatter3(mni_coords_template(:,1),mni_coords_template(:,2),mni_coords_template(:,3),'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',LWW)
    hold on
    scatter3(as(:,1),as(:,2),as(:,3),'MarkerEdgeColor','k','MarkerFaceColor','b','LineWidth',LWW)
    hold on
    scatter3(bs(:,1),bs(:,2),bs(:,3),'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',LWW)
    hold on
    scatter3(cs(:,1),cs(:,2),cs(:,3),'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',LWW)
    set(gcf,'Color','w')
    view(-126,23)
    axis off
    export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_1_S3/NewImagesFigure1/SpatialParcels_3colours.eps']);
    %actual plotting (3 parcels together, same colour)
    figure
    catc = cat(1,as,bs,cs);
    scatter3(mni_coords_template(:,1),mni_coords_template(:,2),mni_coords_template(:,3),'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',LWW)
    hold on
    scatter3(catc(:,1),catc(:,2),catc(:,3),'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',LWW)
    set(gcf,'Color','w')
    view(-126,23)
    axis off
    export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_1_S3/NewImagesFigure1/SpatialParcels_samecolour.eps']);
    %3 independent images (1 per parcel)
    figure
    scatter3(mni_coords_template(:,1),mni_coords_template(:,2),mni_coords_template(:,3),'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',LWW)
    hold on
    scatter3(as(:,1),as(:,2),as(:,3),'MarkerEdgeColor','k','MarkerFaceColor','b','LineWidth',LWW)
    set(gcf,'Color','w')
    view(-126,23)
    axis off
    export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_1_S3/NewImagesFigure1/SpatialParcel1.eps']);
    figure
    scatter3(mni_coords_template(:,1),mni_coords_template(:,2),mni_coords_template(:,3),'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',LWW)
    hold on
    scatter3(bs(:,1),bs(:,2),bs(:,3),'MarkerEdgeColor','k','MarkerFaceColor','r','LineWidth',LWW)
    set(gcf,'Color','w')
    view(-126,23)
    axis off
    export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_1_S3/NewImagesFigure1/SpatialParcel2.eps']);
    figure
    scatter3(mni_coords_template(:,1),mni_coords_template(:,2),mni_coords_template(:,3),'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',LWW)
    hold on
    scatter3(cs(:,1),cs(:,2),cs(:,3),'MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',LWW)
    set(gcf,'Color','w')
    view(-126,23)
    axis off
    export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_1_S3/NewImagesFigure1/SpatialParcel3.eps']);
    %temporal clustering - 6 parcels
    fname = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_1_S3/tempclust_k_6_spatialclust_clspatvect_time_1_391_actualdata.nii.gz'; %nifti image with 6 (temporal) parcels
    V_parc = nii.load(fname); %loading nifti image
    V_parc2 = reshape(V_parc(V_parc~=0),[3559,1]); %this is a correct representation of voxels (rows) and their values over time (columns)
    COL = [0 0.6 1; 1 0.4 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0]; %vector with colors
    for ooo = 1:5 %over parcels of temporal parcellation
        %getting coordinates of 3 parcels
        as = mni_coords_template((V_parc2==ooo),:);
        figure
        scatter3(mni_coords_template(:,1),mni_coords_template(:,2),mni_coords_template(:,3),'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',LWW)
        hold on
        scatter3(as(:,1),as(:,2),as(:,3),'MarkerEdgeColor','k','MarkerFaceColor',COL(ooo,:),'LineWidth',LWW)
        set(gcf,'Color','w')
        view(-126,23)
        axis off
        export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_1_S3/NewImagesFigure1/TemporalParcel_' num2str(ooo) '.eps']);
    end
    %loading and plotting scatter plot for temporal clusters and indices of maximum value for each brain voxel
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_1_S3/NewImagesFigure1/scatterstuff.mat');
    figure
    d2 = cell(1,clustnum);
    COL = [0 0.6 1; 1 0.4 0; 1 0 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0]; %vector with colors
    for jj = 1:clustnum %here there is probably a better way to do the scatter plot.. without doing the loop.. but for now this is fine..
        imclustI2 = imclustI{poiI(jj)}; %plotting clusters sorted by their size (descending order)
        scatter(imclustI2,ones(length(imclustI2),1)*jj,[],COL(jj,:)) %scatter plotting..
        hold on
        scatter(CO(jj),jj,'wx','LineWidth',1);
        hold on
    end
    set(gcf,'Color','w')
    grid minor
    %method figure.. functional part of the k-means, frequency 2-8Hz
    fname = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_0_Average_l_0_Old_l_1/tempclust_k_4_time_1_200_actualdata.nii.gz';
    V_parc = nii.load(fname); %loading nifti image
    V_parc2 = reshape(V_parc(V_parc~=0),[3559,1]); %this is a correct representation of voxels (rows) and their values over time (columns)
    COL = [0 0.6 1; 1 0.4 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0]; %vector with colors
    for ooo = 1:4 %over parcels of temporal parcellation
        %getting coordinates of 4 parcels
        as = mni_coords_template((V_parc2==ooo),:);
        figure
        scatter3(mni_coords_template(:,1),mni_coords_template(:,2),mni_coords_template(:,3),'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',LWW)
        hold on
        scatter3(as(:,1),as(:,2),as(:,3),'MarkerEdgeColor','k','MarkerFaceColor',COL(ooo,:),'LineWidth',LWW)
        set(gcf,'Color','w')
        view(-126,23)
        axis off
        export_fig(['Freq_2_8Hz_TemporalParcel_' num2str(ooo) '.eps']);
    end
    %loading and plotting scatter plot for temporal clusters and indices of maximum value for each brain voxel
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_0_Average_l_0_Old_l_1/Dudu_NuovaFigura1/maxval.mat');
    %loading temporal (functional) clustering indices
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_0_Average_l_0_Old_l_1/Clusters.mat');
    figure
    COL = [0 0.6 1; 1 0.4 0; 1 0 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0]; %vector with colors
    for jj = 1:4 %here there is probably a better way to do the scatter plot.. without doing the loop.. but for now this is fine..
        imclustI2 = im(idk3==jj); %plotting clusters sorted by their size (descending order)
        CO = mean(imclustI2);
        scatter(ones(1,length(imclustI2))*jj,imclustI2,[],COL(jj,:)) %scatter plotting..
        hold on
%         scatter(jj,CO,'x','LineWidth',1);
%         hold on
    end
    set(gcf,'Color','w')
    grid minor
    export_fig(['Freq_2_8Hz_Scatter.eps']);
end

%% EXTRACTING VOXELS FORMING THE PARCELS AND REPORTING THEM IN EXCEL FILES (one of the additional tables of the paper)

%indexing the previously calculated clusters solutions (in the different frequency bands and for averaged and non-averaged OLD and NEW)
list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Sl*');
outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach';
files1{1} = 'tempclust_k_4_spatialclust_clspatvect_time_1_200_parcel_';
files2{1} = '_of_9_parcels';
files1{2} = 'tempclust_k_4_spatialclust_clspatvect_time_1_200_parcel_';
files2{2} = '_of_9_parcels';
files1{3} = 'tempclust_k_4_spatialclust_clspatvect_time_1_200_parcel_';
files2{3} = '_of_10_parcels';
files1{4} = 'tempclust_k_6_spatialclust_clspatvect_time_1_391_parcel_';
files2{4} = '_of_19_parcels';
files1{5} = 'tempclust_k_6_spatialclust_clspatvect_time_1_391_parcel_';
files2{5} = '_of_23_parcels';

for iii = 1:length(list) %over parcellations
    list2 = dir([list(iii).folder '/' list(iii).name '/*parcels.nii.gz']); %getting the different parcels of the parcellation ii
    for jj = 1:length(list2) %over parcels of the parcellation ii
        clear PD
        %legend
        PD{1,1} = 'Brain region'; PD{1,2} = 'Hemisphere'; PD{1,3} = 'T-stat'; PD{1,4} = 'MNI coordinates'; %1st row
        PD{2,4} = 'X'; PD{2,5} = 'Y'; PD{2,6} = 'Z'; %2nd row
        %actual name plus path
        fname = [list2(iii).folder '/' files1{iii} num2str(jj) files2{iii} '.nii.gz'];
        %getting MNI coordinates of significant voxels within the provided image
        [ mni_coords, xform ] = osl_mnimask2mnicoords(fname);
        %loading the image
        V = nii.load(fname);
        %extracting statistics
        VV = V(V~=0);
        %indices of non-zero values of nifti image
        VI = find(V~=0);
        %path to AAL template
        parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'; %load this from the provided codes folder
        %loading AAL labels
        load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat'); %load this from the provided codes folder
        %extracting AAL coordinates information
        K = nii.load(parcelfile);
        %sorting results in order to have strongest voxels at the top (positive t-values) or at the bottom (negative t-values)
        [VV2, II] = sort(VV,'descend');
        VI = VI(II);
        mni_coords = mni_coords(II,:);
        %final cell
%         PD = cell(length(VV2),4);
        %getting AAL indices
        ROI = zeros(length(VI),1);
        cnt = 2;
        for ii = 1:length(VI)
            ROI(ii) = K(VI(ii));
            if ROI(ii) > 0 && ROI(ii) < 91
                cnt = cnt + 1;
                PD(cnt,1) = {lab(ROI(ii),3:end)}; %storing ROI
                PD(cnt,4) = {mni_coords(ii,1)}; %storing MNI coordinates
                PD(cnt,5) = {mni_coords(ii,2)}; %storing MNI coordinates
                PD(cnt,6) = {mni_coords(ii,3)}; %storing MNI coordinates
                if mni_coords(ii,1) > 0 %storing hemisphere
                    PD(cnt,2) = {'R'};
                else
                    PD(cnt,2) = {'L'};
                end
                PD(cnt,3) = {round(VV2(ii),2)}; %storing t-statistics
            end
        end
        PDn = cell2table(PD);
        %         PDn = cell2table(PD(~any(cellfun('isempty',PD),2),:)); %remove the possible empty cell
        if iii < 4
            writetable(PDn,[outdir '/Tables/Table_2/' list(iii).name '.xlsx'],'Sheet',jj) %printing excel file
        else
            writetable(PDn,[outdir '/Tables/Table_3/' list(iii).name '.xlsx'],'Sheet',jj) %printing excel file
        end
    end
end

%% PREPARING DATA (2-8 Hz) FOR PYTHON CURVE_FIT
%this is necessary since the clustering algorithms were used in a slightly shorter time-window, while curve_fit in python is used with a slightly larger time-window

VS = zeros(3559,526,2);
for cc = 1:2 %over conditions (old and new)
    %loading data
    %path to file nifti
    %2 - 8Hz (N100)
    pathnii = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_2_8Hz_AllTrials_BC_noacope_2.oat/wholebrain_first_level_BC_minor_group_level_everybody_BC_dir'; %baseline corrected and nocope
    fname = [pathnii '/tstat' num2str(cc) '_gc1_8mm.nii.gz']; %tstat name
    %loading the image
    V_stat4D = nii.load(fname);
    %removing zeros from the statistics and reshaping it.. here the idea is passing from the representation in the "brain layout" of the statistics into a 2D representation of the timeseries of the voxels..
    Vstat2 = reshape(V_stat4D(V_stat4D~=0),[3559,size(V_stat4D,4)]); %this is a correct representation of voxels (rows) and their values over time (columns)
    SS = size(Vstat2);
    % reversing the sign
    Vdum = zeros(SS(1),SS(2)); %preallocating space
    vect = zeros(SS(1),1); %preallocating space
    for ii = 1:SS(1) %over brain voxels
        if mean(Vstat2(ii,30:40),2) < 0 %if the data in voxel ii is negative during N100 time
            Vdum(ii,:) = Vstat2(ii,:); %keeping the original data
            vect(ii,1) = 1; %storing a vector with 1 and -1 to be used for later statistics
        else
            Vdum(ii,:) = Vstat2(ii,:) .* (-1); %otherwise reversing..
            vect(ii,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
        end
    end
    Vstat2 = Vdum;
    VS(:,:,cc) = Vstat2;
end
%loading parcels
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_0_Average_l_0_Old_l_1/Clusters.mat')
lko = zeros(max(kkk2),316,2);
for ii = 1:max(kkk2)
    lko(ii,:,1) = mean(VS(kkk2==ii,1:316,1),1); %getting time-series of the parcels (between -0.1 and 2 sec)
    lko(ii,:,2) = (std(VS(kkk2==ii,1:316,2),0,1))./sqrt(length(find(kkk2==ii))); %standard errors..
end
save('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/OLD.mat','lko') %saving new time-series
%loading parcels
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_0_Average_l_0_Old_l_2/Clusters.mat')
lkn = zeros(max(kkk2),316,2);
for ii = 1:max(kkk2)
    lkn(ii,:,1) = mean(VS(kkk2==ii,1:316,2),1); %getting time-series of the parcels (between -0.1 and 2 sec)
    lkn(ii,:,2) = (std(VS(kkk2==ii,1:316,2),0,1))./sqrt(length(find(kkk2==ii))); %standard errors..
end
save('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/NEW.mat','lkn') %saving new time-series

%%

%here the fitting is computed in Python..

%% PLOTTING CURVE_FIT FROM PYTHON RESULTS IN MATLAB

%after computing fitting using python (curve_fit), I take the results from
%python and I use them here in Matlab to produce the final plots (this is
%done to maintain coherence with the graphical display of all the figures)
LW = 1.5; %line-width (waveforms)

%loading time for plotting purposes
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/time.mat');
%0.1-1Hz
%loading standard errors over voxels
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_1_Average_l_0_Old_l_1/Clusters.mat');
%reading data created in python and saved as txts
dirg = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/'; %path to txt files
for ii = 1:19 %over timeseries
    KJK2 = KJK(ii,:,2); %extracting standard errors
    a = dir([dirg 'slowfreq_predicted_data_' num2str(ii) '.txt']); %checking whether file exists for parcel ii
    if ~isempty(a) %if exists
        fileID = fopen([dirg 'slowfreq_predicted_data_' num2str(ii) '.txt']);
        pr = fscanf(fileID,'%f'); %reading txt file for predicted data
        fileID = fopen([dirg 'slowfreq_real_data_' num2str(ii) '.txt']);
        rl = fscanf(fileID,'%f'); %reading txt file for real data
        figure
        plot(time(1:length(rl)),rl,'b','LineWidth',LW) %real data
        hold on
        plot(time(1:length(rl)),rl' + KJK2(1:length(rl)),'--','Color','b','LineWidth',0.5) %real data
        hold on
        plot(time(1:length(rl)),rl' - KJK2(1:length(rl)),'--','Color','b','LineWidth',0.5) %real data
        hold on
        plot(time(1:length(pr)),pr,'--','Color','r','LineWidth',LW) %predicted data (dash line)
        ylim([-2.9 3.8]) %amplitude limits
        xlim([time(1) time(length(rl))]) %time limits
        grid minor
        set(gcf,'color','w');
        %exporting images (png and eps)
%         export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/NewTimeseriesMatlabFromPython/SlowFreq/Parcel_' num2str(ii) '.png']);
%         export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/NewTimeseriesMatlabFromPython/SlowFreq/Parcel_' num2str(ii) '.eps']);
    end
end
%2-8Hz
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/OLD.mat');
KJKOLD = lko;
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/NEW.mat');
KJKNEW = lkn;
%reading data created in python and saved as txts
dirg = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/'; %path to txt files
for ii = 1:9 %over timeseries
    KJK2 = KJKOLD(ii,:,2); %extracting standard errors
    a = dir([dirg 'fastfreq_predicted_data_' num2str(ii) 'block_2.txt']); %checking whether file exists for parcel ii
    if ~isempty(a) %if exists
        %OLD
        fileID = fopen([dirg 'fastfreq_predicted_data_' num2str(ii) 'block_2.txt']);
        pr = fscanf(fileID,'%f'); %reading txt file for predicted data
        fileID = fopen([dirg 'fastfreq_real_data_' num2str(ii) 'block_2.txt']);
        rl = fscanf(fileID,'%f'); %reading txt file for real data
        figure
        plot(time(1:length(rl)),rl,'b','LineWidth',LW) %real data
        hold on
        plot(time(1:length(rl)),rl' + KJK2(1:length(rl)),'--','Color','b','LineWidth',0.5) %real data
        hold on
        plot(time(1:length(rl)),rl' - KJK2(1:length(rl)),'--','Color','b','LineWidth',0.5) %real data
        hold on
        plot(time(1:length(pr)),pr,'--','Color','r','LineWidth',LW) %predicted data (dash line)
        ylim([-24 19]) %amplitude limits
        xlim([time(1) time(length(rl))]) %time limits
        grid minor
        set(gcf,'color','w');
        %exporting images (png and eps)
        export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/NewTimeseriesMatlabFromPython/FastFreq/Parcel_' num2str(ii) '_block_2.png']);
        export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/NewTimeseriesMatlabFromPython/FastFreq/Parcel_' num2str(ii) '_block_2.eps']);
    end
    KJK2 = KJKNEW(ii,:,2); %extracting standard errors
    a = dir([dirg 'fastfreq_predicted_data_' num2str(ii) 'block_3.txt']); %checking whether file exists for parcel ii
    if ~isempty(a) %if exists
        %NEW
        fileID = fopen([dirg 'fastfreq_predicted_data_' num2str(ii) 'block_3.txt']);
        pr = fscanf(fileID,'%f'); %reading txt file for predicted data
        fileID = fopen([dirg 'fastfreq_real_data_' num2str(ii) 'block_3.txt']);
        rl = fscanf(fileID,'%f'); %reading txt file for real data
        figure
        plot(time(1:length(rl)),rl,'b','LineWidth',LW) %real data
        hold on
        plot(time(1:length(rl)),rl' + KJK2(1:length(rl)),'--','Color','b','LineWidth',0.5) %real data
        hold on
        plot(time(1:length(rl)),rl' - KJK2(1:length(rl)),'--','Color','b','LineWidth',0.5) %real data
        hold on
        plot(time(1:length(pr)),pr,'--','Color','r','LineWidth',LW) %predicted data (dash line)
        ylim([-24 19]) %amplitude limits
        xlim([time(1) time(length(rl))]) %time limits
        grid minor
        set(gcf,'color','w');
        %exporting images (png and eps)
        export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/NewTimeseriesMatlabFromPython/FastFreq/Parcel_' num2str(ii) '_block_3.png']);
        export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/NewTimeseriesMatlabFromPython/FastFreq/Parcel_' num2str(ii) '_block_3.eps']);
    end
end

%% CONTRAST ON THE NEW PARCELS (LOADING DATA FROM SUBJECT LEVEL)
% here I use the clusters computed by the functional-spatial k-means on the average of the OLD and NEW group-level statistics
%Please, note that this is done on copes, while in the paper I reported the t-values (e.g. taken from the images at subject-level)
%On the conceptual level, this is nearly the same

lab01 = 1; % 1 for 0.1-1Hz; 2 for 2-8Hz 

%loading outputted clusters (by functional-spatial kmeans) and other information
if lab01 == 1
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_1_Average_l_1_Old_l_1/Clusters.mat')
    subjdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_01_1Hz_AllTrials_BC_2.oat'; % spm_files_recog_basen{ii}(12:25)];
    list = dir([subjdir '/subj*']);
    outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_1_Average_l_1_Old_l_1';
else
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_0_Average_l_1_Old_l_1/Clusters.mat')
    subjdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_2_8Hz_AllTrials_BC_noacope_2.oat'; % spm_files_recog_basen{ii}(12:25)];
    list = dir([subjdir '/subj*BC_minor.mat']);
    %loading the previously calculated vector to fix the sign ambiguity issue
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_0_Average_l_1_Old_l_1/onesvector_negativereverse.mat');
    outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_0_Average_l_1_Old_l_1';
end
%extracting and preparing data
OLD = zeros(max(kkk2),526,length(list)); %hard-coding for 526 time-points.. anyway, this should not really change.. so it is ok..
NEW = zeros(max(kkk2),526,length(list)); %hard-coding for 526 time-points.. anyway, this should not really change.. so it is ok..
for ii = 1:length(list) %over subjets
    load([list(ii).folder '/' list(ii).name])
    if lab01 ~= 1 %if 2-8Hz data
        subjdum = oat_stage_results.cope(:,:,1:2) .* vect; %fixing the sign issue (only for first two contrast (i.e. OLD and NEW))
    else %otherwise you are working with absolute values (0.1-1Hz), thus you do not need to fix the sign issue (again)..
        subjdum = oat_stage_results.cope(:,:,1:2);
    end
    for pp = 1:max(kkk2) %over different parcels
        OLD(pp,:,ii) = mean(subjdum(kkk2==pp,:,1),1); %condition1 (OLD); storing each participant
        NEW(pp,:,ii) = mean(subjdum(kkk2==pp,:,2),1); %condition2 (NEW); storing each participant
    end
    disp(ii)
end
%saving data prepared for subsequent t-tests
save([outdir '/subjlevel.mat'],'OLD','NEW')

%% ACTUAL T-TESTS AND MCS

lab01 = 1; % 1 for 0.1-1Hz; 2 for 2-8Hz

%loading outputted clusters (by functional-spatial kmeans) and other information
if lab01 == 1
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_1_Average_l_1_Old_l_1/subjlevel.mat');
else
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_0_Average_l_1_Old_l_1/subjlevel.mat');
end
%actual computation of t-tests
P = ones(size(OLD,1),size(OLD,2)); %ones since later I binarize the vector with values < .05, so it is more convenient having here 1s instead of 0s (that otherwise would be binarized as significant as well..)
TVAL = zeros(size(OLD,1),size(OLD,2));
for ii = 1:size(OLD,1) %over parcels
    for jj = 16:size(OLD,2) %over time-points (starting from 16 which corresponds to 0ms from the onset of the stimulus, which is actuallu the onset of the stimulus..)
        a = squeeze(OLD(ii,jj,:)); %parcel ii, time-point jj and all subjects for condition1 (OLD)
        b = squeeze(NEW(ii,jj,:)); %parcel ii, time-point jj and all subjects for condition2 (NEW)
%         if lab01 == 2
%             a = abs(a); %just to make more easily readable the contrasts (so negative t-values always means NEW stronger than OLD and viceversa..)
%             b = abs(b);
%         end
        [~,p,~,stats] = ttest(a,b); %compute t-test between the two conditions (OLD vs NEW)
        P(ii,jj) = p; %store p-values in P
        TVAL(ii,jj) = stats.tstat; %store t-values in TVAL
    end
    disp(ii)
end
%Monte-Carlo simulation
%parameters for Monte-Carlo simulation on significant time-points emerged from t-tests
p = 0.05;
max_lab = 1;
permnum = 1000;
CL = cell(size(P,1),2); %1st column = 1.cluster size; 2.p-value; 3.significant time-points; 2nd column = mean t-value over significant time-points of the cluster(s)
MCS_thresh = 0.001;
for ii = 1:size(P,1) %over parcels
    %one parcel at the time
    P2 = P(ii,1:400); %here 400 is an arbitrary value indicating the time-points when the task was actually happening (the full epoch is a bit longer)
    %actual Monte-Carlo simulation
    sign_clust = oneD_MCS_LBPD_D(P2,p,max_lab,permnum,MCS_thresh);
    CL{ii,1} = sign_clust; %storing significant clusters for parcel ii
    dumcl = zeros(size(sign_clust,1),1);
    for cc = 1:size(sign_clust,1) %over significant clusters for parcel ii
        CL{ii,2}(cc) = mean(TVAL(ii,sign_clust{cc,3})); %storing mean t-value over significant time-points (for cluster cc)
    end
%     CL{ii,2}(cc) = mean(dumcl); %mean over t-values for significant clusters of parcel ii
    disp(ii)
end

%% SAVING SIGNIFICANT TIME-WINDOWS FOR CONTRASTS OLD-NEW OVER ALL PARCELS (BOTH 0.1-1Hz AND 2-8Hz) TO EXCEL FILES (reported in one of the supplementary tables of the paper)

%this requires the MCS run in the previous (above) section

if lab01 == 1
    outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/Table_7/ContrastAfterClustering_01_1Hz';
else
    outdir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/Table_5/ContrastAfterClustering_2_8Hz';
end
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/time.mat'); %loading time in seconds
%legend
clear PD
PD{1,1} = 'Parcel #'; PD{1,2} = 'k '; PD{1,3} = 'p '; PD{1,4} = 'mean t-value'; PD{1,5} = 'time 1'; PD{1,6} = 'time 2'; %1st row
for ii = 1:length(CL)
    if size(CL{ii},1) ~= 0 %if there is at least one actual significant time-window, highlight it/them
        cnt7 = -6;
        for jj = 1:size(CL{ii},1) %over significant clusters for parcel ii
            cnt7 = cnt7 + 7; %this is to get the significant time-windows of parcel ii on the same row (with one empty column separating the different significant time-windows)
            PD{ii+1,cnt7} = ii; %parcel progressive number
            PD{ii+1,cnt7+1} = CL{ii,1}{jj,1}; %cluster size (k)
            PD{ii+1,cnt7+2} = CL{ii,1}{jj,2}; %cluster p-value after MCS
            PD{ii+1,cnt7+3} = CL{ii,2}(jj); %mean t-value for cluster jj and parcel ii over the significant time-points
            PD{ii+1,cnt7+4} = time(CL{ii,1}{jj,3}(1)); %time1 of the significant time-window jj
            PD{ii+1,cnt7+5} = time(CL{ii,1}{jj,3}(end)); %time end (last time-point) of the significant time-window jj
        end
    else
        PD{ii+1,1} = ii; %parcel progressive number
    end
end

PDn = cell2table(PD); %remove the possible empty cell
writetable(PDn,[outdir '.xlsx'],'Sheet',1) %printing excel file
    
%% PLOTTING SUBJECT LEVEL DATA (SUBJECTS COMIBINED) WITH RESULTS OF T-TESTS (GREY AREAS) AND STANDARD DEVIATIONS

%this requires the MCS run two sections above

save_images = 0; %1 to save images; 0 otherwise

load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/time.mat'); %loading time in seconds
data = zeros(size(OLD,1),size(OLD,2),size(OLD,3),2); %preparing data
data(:,:,:,1) = OLD; data(:,:,:,2) = NEW;
S = [];
S.data = data; %data
S.condition_n = 1:2;
S.STE = 2; %1 for plotting standard errors (shadows); 2 (dot lines); 0 (nothing)
S.conds = {'old','new'};
S.colorline = {'r','b'}; %Select colorline for each group
if lab01 == 1
    S.y_lim = [-0.045 0.07]; %Set y limits
    S.x_lim = [time(1) time(391)]; % Set x limits
else
    S.y_lim = [-0.6 0.45]; %Set y limits
    S.x_lim = [time(1) time(316)]; % Set x limits
end
S.time = time(1:526);
S.legendl = 0;
S.signtp = {[]};
for ii = 1:length(CL) %over parcels
    S.chans_index = ii; %plot waveforms at given parcel
    S = rmfield(S,'signtp'); %clearing the S.signtp before assigning the significant time-windows (to avoid that previously significant time-windows (for parcel ii - 1) are used also for parcel ii)
    if size(CL{ii},1) ~= 0 %if there is at least one actual significant time-window, highlight it/them
        for jj = 1:size(CL{ii},1) %over significant clusters for parcel ii
            S.signtp{jj,1} = [time(CL{ii}{jj,3}(1)) time(CL{ii}{jj,3}(end))]; %quite elaborated way to get the significant time-points
        end
    else %otherwise just pass the empty vector
        S.signtp{1} = [2.495 2.5];
    end
    plot_wave_conditions(S) %actual function
    if save_images == 1
        %saving figures
        if lab01 == 1
            export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_S6/Timeseries/Parcel_' num2str(ii) '.png']);
            export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_S6/Timeseries/Parcel_' num2str(ii) '.eps']);
        else
            export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_4/Timeseries/Parcel_' num2str(ii) '.png']);
            export_fig(['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Figures/Figure_4/Timeseries/Parcel_' num2str(ii) '.eps']);
        end
    end
end

%% PLOTTING FREQ 2-8Hz PRIMARY AND SECONDARY AUDITORY CORTEX AND CINGULATE TOGETHER FOR OLD AND NEW..

%loading OLD and NEW (from k-means on averaged OLD and NEW data)
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_0_Average_l_1_Old_l_1/subjlevel.mat');
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/time.mat'); %loading time in seconds
OLDstde = std(OLD,0,3)./sqrt(size(OLD,3));
NEWstde = std(NEW,0,3)./sqrt(size(NEW,3));
col1 = 'b'; col2 = 'r'; col3 = 'k'; %colors of the different brain parcels
linwidmean = 2; %line width of mean
linwidstde = 0.5; %line width of standard errors
stde_l = 0; %1 for standard errors; 0 for not
contrast_l = 0; %1 for k-means on contrast; 0 for k-means on OLD

close all
if contrast_l == 1
    %OLD
    %left hemisphere and cingulate
    figure
    plot(time(1:526),mean(OLD(10,:,:),3),'Color',col1,'LineWidth',linwidmean) %A1
    hold on
    if stde_l == 1
        plot(time(1:526),mean(OLD(10,:,:),3) + OLDstde(10,:),'--','Color',col1,'LineWidth',linwidstde) %A1
        hold on
        plot(time(1:526),mean(OLD(10,:,:),3) - OLDstde(10,:),'--','Color',col1,'LineWidth',linwidstde) %A1
        hold on
    end
    plot(time(1:526),mean(OLD(8,:,:),3),'Color',col2,'LineWidth',linwidmean) %A2
    hold on
    if stde_l == 1
        plot(time(1:526),mean(OLD(8,:,:),3) + OLDstde(8,:),'--','Color',col2,'LineWidth',linwidstde) %A2
        hold on
        
        plot(time(1:526),mean(OLD(8,:,:),3) - OLDstde(8,:),'--','Color',col2,'LineWidth',linwidstde) %A2
        hold on
    end
    plot(time(1:526),mean(OLD(6,:,:),3),'Color',col3,'LineWidth',linwidmean) %cingulate
    hold on
    if stde_l == 1
        plot(time(1:526),mean(OLD(6,:,:),3) + OLDstde(6,:),'--','Color',col3,'LineWidth',linwidstde) %cingulate
        hold on
        plot(time(1:526),mean(OLD(6,:,:),3) - OLDstde(6,:),'--','Color',col3,'LineWidth',linwidstde) %cingulate
    end
    grid minor
    set(gcf,'color','w')
    xlim([time(1) time(300)])
    ylim([-0.55 0.4])
    title('OLD left hem and cing')
    %right hemisphere and cingulate
    figure
    plot(time(1:526),mean(OLD(9,:,:),3),'Color',col1,'LineWidth',linwidmean) %A1
    hold on
    if stde_l == 1
        plot(time(1:526),mean(OLD(9,:,:),3) + OLDstde(10,:),'--','Color',col1,'LineWidth',linwidstde) %A1
        hold on
        plot(time(1:526),mean(OLD(9,:,:),3) - OLDstde(10,:),'--','Color',col1,'LineWidth',linwidstde) %A1
        hold on
    end
    plot(time(1:526),mean(OLD(7,:,:),3),'Color',col2,'LineWidth',linwidmean) %A2
    hold on
    if stde_l == 1
        plot(time(1:526),mean(OLD(7,:,:),3) + OLDstde(8,:),'--','Color',col2,'LineWidth',linwidstde) %A2
        hold on
        plot(time(1:526),mean(OLD(7,:,:),3) - OLDstde(8,:),'--','Color',col2,'LineWidth',linwidstde) %A2
        hold on
    end
    plot(time(1:526),mean(OLD(6,:,:),3),'Color',col3,'LineWidth',linwidmean) %cingulate
    hold on
    if stde_l == 1
        plot(time(1:526),mean(OLD(6,:,:),3) + OLDstde(6,:),'--','Color',col3,'LineWidth',linwidstde) %cingulate
        hold on
        plot(time(1:526),mean(OLD(6,:,:),3) - OLDstde(6,:),'--','Color',col3,'LineWidth',linwidstde) %cingulate
    end
    grid minor
    set(gcf,'color','w')
    xlim([time(1) time(300)])
    ylim([-0.55 0.4])
    title('OLD right hem and cing')
    
    %NEW
    %left hemisphere and cingulate
    figure
    plot(time(1:526),mean(NEW(10,:,:),3),'Color',col1,'LineWidth',linwidmean) %A1
    hold on
    if stde_l == 1
        plot(time(1:526),mean(NEW(10,:,:),3) + NEWstde(10,:),'--','Color',col1,'LineWidth',linwidstde) %A1
        hold on
        plot(time(1:526),mean(NEW(10,:,:),3) - NEWstde(10,:),'--','Color',col1,'LineWidth',linwidstde) %A1
        hold on
    end
    plot(time(1:526),mean(NEW(8,:,:),3),'Color',col2,'LineWidth',linwidmean) %A2
    hold on
    if stde_l == 1
        plot(time(1:526),mean(NEW(8,:,:),3) + NEWstde(8,:),'--','Color',col2,'LineWidth',linwidstde) %A2
        hold on
        plot(time(1:526),mean(NEW(8,:,:),3) - NEWstde(8,:),'--','Color',col2,'LineWidth',linwidstde) %A2
        hold on
    end
    plot(time(1:526),mean(NEW(6,:,:),3),'Color',col3,'LineWidth',linwidmean) %cingulate
    hold on
    if stde_l == 1
        plot(time(1:526),mean(NEW(6,:,:),3) + NEWstde(6,:),'--','Color',col3,'LineWidth',linwidstde) %cingulate
        hold on
        plot(time(1:526),mean(NEW(6,:,:),3) - NEWstde(6,:),'--','Color',col3,'LineWidth',linwidstde) %cingulate
    end
    grid minor
    set(gcf,'color','w')
    xlim([time(1) time(300)])
    ylim([-0.55 0.4])
    title('NEW left hem and cing')
    %right hemisphere and cingulate
    figure
    plot(time(1:526),mean(NEW(9,:,:),3),'Color',col1,'LineWidth',linwidmean) %A1
    hold on
    if stde_l == 1
        plot(time(1:526),mean(NEW(9,:,:),3) + NEWstde(10,:),'--','Color',col1,'LineWidth',linwidstde) %A1
        hold on
        plot(time(1:526),mean(NEW(9,:,:),3) - NEWstde(10,:),'--','Color',col1,'LineWidth',linwidstde) %A1
        hold on
    end
    plot(time(1:526),mean(NEW(7,:,:),3),'Color',col2,'LineWidth',linwidmean) %A2
    hold on
    if stde_l == 1
        plot(time(1:526),mean(NEW(7,:,:),3) + NEWstde(8,:),'--','Color',col2,'LineWidth',linwidstde) %A2
        hold on
        plot(time(1:526),mean(NEW(7,:,:),3) - NEWstde(8,:),'--','Color',col2,'LineWidth',linwidstde) %A2
        hold on
    end
    plot(time(1:526),mean(NEW(6,:,:),3),'Color',col3,'LineWidth',linwidmean) %cingulate
    hold on
    if stde_l == 1
        plot(time(1:526),mean(NEW(6,:,:),3) + NEWstde(6,:),'--','Color',col3,'LineWidth',linwidstde) %cingulate
        hold on
        plot(time(1:526),mean(NEW(6,:,:),3) - NEWstde(6,:),'--','Color',col3,'LineWidth',linwidstde) %cingulate
    end
    grid minor
    set(gcf,'color','w')
    xlim([time(1) time(300)])
    ylim([-0.55 0.4])
    title('NEW right hem and cing')
else
    %loading time for plotting purposes
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/time.mat');
    COL(1) = col1; COL(2) = col2; COL(3) = col3; %barbaric way to get colours..
    %0.1-1Hz
    %loading standard errors over voxels
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_0_Average_l_0_Old_l_1/Clusters.mat');
    %reading data created in python and saved as txts
    dirg = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/'; %path to txt files
    %left hemisphere and cingulate
    v = [9 7 6]; %selected timeseries
    figure
    for ii = 1:length(v) %over timeseries
        KJK2 = KJK(v(ii),:,2); %extracting standard errors
        a = dir([dirg 'fastfreq_predicted_data_' num2str(v(ii)) 'block_2.txt']); %checking whether file exists for parcel ii
        if ~isempty(a) %if exists
            fileID = fopen([dirg 'fastfreq_predicted_data_' num2str(v(ii)) 'block_2.txt']);
            pr = fscanf(fileID,'%f'); %reading txt file for predicted data
            fileID = fopen([dirg 'fastfreq_real_data_' num2str(v(ii)) 'block_2.txt']);
            rl = fscanf(fileID,'%f'); %reading txt file for real data
            plot(time(1:length(rl)),rl,'Color',COL(ii),'LineWidth',2) %real data
            hold on
        end
    end
    grid minor
    set(gcf,'color','w')
    xlim([time(1) time(300)])
    ylim([-25 20])
    title('OLD left hem and cing')
    %right hemisphere and cingulate
    v = [8 5 6]; %selected timeseries
    figure
    for ii = 1:length(v) %over timeseries
        KJK2 = KJK(v(ii),:,2); %extracting standard errors
        a = dir([dirg 'fastfreq_predicted_data_' num2str(v(ii)) 'block_2.txt']); %checking whether file exists for parcel ii
        if ~isempty(a) %if exists
            fileID = fopen([dirg 'fastfreq_predicted_data_' num2str(v(ii)) 'block_2.txt']);
            pr = fscanf(fileID,'%f'); %reading txt file for predicted data
            fileID = fopen([dirg 'fastfreq_real_data_' num2str(v(ii)) 'block_2.txt']);
            rl = fscanf(fileID,'%f'); %reading txt file for real data
            plot(time(1:length(rl)),rl,'Color',COL(ii),'LineWidth',2) %real data
            hold on
        end
    end
    grid minor
    set(gcf,'color','w')
    xlim([time(1) time(300)])
    ylim([-25 20])
    title('OLD right hem and cing')
end

%% CREATING SOME TIME-SERIES WITH SOME RANDOM NOISE FOR METHODS FIGURES

%FIGURE 2 - Methods
LW = 1.5; %line-width (waveforms)
%loading time for plotting purposes
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/time.mat');
%0.1-1Hz
%loading standard errors over voxels
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_1_Average_l_0_Old_l_1/Clusters.mat');
%reading data created in python and saved as txts
cnt = 0;
dirg = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/'; %path to txt files
for ii = 2:4%19 %over timeseries
    KJK2 = KJK(ii,:,2); %extracting standard errors
    cnt = cnt + 1;
    a = dir([dirg 'slowfreq_predicted_data_' num2str(ii) '.txt']); %checking whether file exists for parcel ii
    if ~isempty(a) %if exists
        fileID = fopen([dirg 'slowfreq_predicted_data_' num2str(ii) '.txt']);
        pr = fscanf(fileID,'%f'); %reading txt file for predicted data
        fileID = fopen([dirg 'slowfreq_real_data_' num2str(ii) '.txt']);
        rl = fscanf(fileID,'%f'); %reading txt file for real data
        figure
        %plotting simulated contrast
        rlnoise = rl./2; %adding noise
        plot(time(1:length(rl)),rlnoise,'r','LineWidth',LW) %noisy data
        pr2 = (pr./5) - 0.01*rand(size(pr));
        hold on
        plot(time(1:length(pr)),pr2,'Color','b','LineWidth',LW) %predicted data (dash line)
        xlim([time(1) time(length(rl))]) %time limits
        grid minor
        set(gcf,'color','w');
        export_fig(['fake_contrast_01_1Hz_' num2str(cnt) '.eps']);
        %plotting simulated modelling
        figure
        rlnoise = rl./1.5; %adding noise
        plot(time(1:length(rl)),rlnoise,'r','LineWidth',LW) %noisy data
%         pr2 = (pr./2) + 0.6*rand(size(pr));
        hold on 
        plot(time(1:length(pr)),pr,'--','Color','k','LineWidth',LW) %predicted data (dash line)
%         ylim([-24 19]) %amplitude limits
        xlim([time(1) time(length(rl))]) %time limits
        grid minor
        set(gcf,'color','w');
        %exporting images (png and eps)
        export_fig(['fake_modelling_01_1Hz_' num2str(cnt) '.eps']);
    end
end
%2-8Hz
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/OLD.mat');
KJKOLD = lko;
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/NEW.mat');
KJKNEW = lkn;
%reading data created in python and saved as txts
dirg = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/'; %path to txt files
cnt = 0;
for ii = 7:9 %over timeseries
    KJK2 = KJKOLD(ii,:,2); %extracting standard errors
    cnt = cnt + 1;
    a = dir([dirg 'fastfreq_predicted_data_' num2str(ii) 'block_2.txt']); %checking whether file exists for parcel ii
    if ~isempty(a) %if exists
        %OLD
        fileID = fopen([dirg 'fastfreq_predicted_data_' num2str(ii) 'block_2.txt']);
        pr = fscanf(fileID,'%f'); %reading txt file for predicted data
        fileID = fopen([dirg 'fastfreq_real_data_' num2str(ii) 'block_2.txt']);
        rl = fscanf(fileID,'%f'); %reading txt file for real data
        figure
        %exponential decay to moderate amplitude of the wave
        tau1 = 1;
        ti   = 1;
        t = time(1:315);
        yinitial = exp(-tau1*t)*8;
        %plotting simulated contrast
        rlnoise = rl./yinitial'; %adding noise
        plot(time(1:length(rl)),rlnoise,'r','LineWidth',LW) %noisy data
        pr2 = (pr./yinitial') - 0.2*rand(size(pr2));
        hold on
        plot(time(1:length(pr)),pr2,'Color','b','LineWidth',LW) %predicted data (dash line)
        xlim([time(1) time(length(rl))]) %time limits
        grid minor
        set(gcf,'color','w');
        export_fig(['fake_contrast_2_8Hz_' num2str(cnt) '.eps']);
        %plotting simulated modelling
        figure
        plot(time(1:length(rl)),rlnoise,'r','LineWidth',LW) %noisy data
        pr2 = (pr./yinitial') + 0.6*rand(size(pr2));
        hold on 
        plot(time(1:length(pr)),pr2,'--','Color','k','LineWidth',LW) %predicted data (dash line)
%         ylim([-24 19]) %amplitude limits
        xlim([time(1) time(length(rl))]) %time limits
        grid minor
        set(gcf,'color','w');
        %exporting images (png and eps)
        export_fig(['fake_modelling_2_8Hz_' num2str(cnt) '.eps']);
    end
end
%FIGURE 1 - Methods
LW = 1.5; %line-width (waveforms)
%loading time for plotting purposes
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/time.mat');
%0.1-1Hz
%loading standard errors over voxels
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_1_Average_l_0_Old_l_1/Clusters.mat');
%reading data created in python and saved as txts
cnt = 0;
dirg = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/'; %path to txt files
for ii = 8:9%19 %over timeseries
    KJK2 = KJK(ii,:,2); %extracting standard errors
    cnt = cnt + 1;
    a = dir([dirg 'slowfreq_predicted_data_' num2str(ii) '.txt']); %checking whether file exists for parcel ii
    if ~isempty(a) %if exists
        fileID = fopen([dirg 'slowfreq_predicted_data_' num2str(ii) '.txt']);
        pr = fscanf(fileID,'%f'); %reading txt file for predicted data
        fileID = fopen([dirg 'slowfreq_real_data_' num2str(ii) '.txt']);
        rl = fscanf(fileID,'%f'); %reading txt file for real data
        figure
        %plotting simulated contrast
        rlnoise = rl./2; %adding noise
        plot(time(1:length(rl)),rlnoise,'r','LineWidth',LW) %noisy data
        ylim([-0.5 4])
        xlim([time(1) time(length(rl))]) %time limits
        grid minor
        set(gcf,'color','w');
        export_fig(['fake_slowfig1_01_1Hz_' num2str(cnt) '.eps']);
    end
end
%2-8Hz
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/OLD.mat');
KJKOLD = lko;
load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/NEW.mat');
KJKNEW = lkn;
%reading data created in python and saved as txts
dirg = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/'; %path to txt files
cnt = 0;
for ii = 3:4 %over timeseries
    KJK2 = KJKOLD(ii,:,2); %extracting standard errors
    cnt = cnt + 1;
    a = dir([dirg 'fastfreq_predicted_data_' num2str(ii) 'block_2.txt']); %checking whether file exists for parcel ii
    if ~isempty(a) %if exists
        %OLD
        fileID = fopen([dirg 'fastfreq_predicted_data_' num2str(ii) 'block_2.txt']);
        pr = fscanf(fileID,'%f'); %reading txt file for predicted data
        fileID = fopen([dirg 'fastfreq_real_data_' num2str(ii) 'block_2.txt']);
        rl = fscanf(fileID,'%f'); %reading txt file for real data
        figure
        %exponential decay to moderate amplitude of the wave
        tau1 = 1;
        ti   = 1;
        t = time(1:315);
        yinitial = exp(-tau1*t)*8;
        %plotting simulated contrast
        rlnoise = rl./yinitial'; %adding noise
        plot(time(1:length(rl)),rlnoise,'r','LineWidth',LW) %noisy data
        ylim([-2 2])
        xlim([time(1) time(length(rl))]) %time limits
        grid minor
        set(gcf,'color','w');
        export_fig(['fake_fastfig1_2_8Hz_' num2str(cnt) '.eps']);
    end
end

%%

%% SPLINES - NOT REPORTED IN THE PAPER(!).. BUT STILL POTENTIALLY INTERESTING..

%%

%%% *** FITTING *** %%%

%FITTING USING EQUATIONS HAS BEEN DONE IN PYTHON

%% FITTING USING SPLINES SOLUTION HAS BEEN IMPLEMENTED HERE USING THE SLMtools

%here, I used functions from:
%https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/26288/versions/6/previews/SLMTools/html/slm_tutorial.html?access_key=

%%% 1) DEFINING BEST NUMBER OF KNOTS (WINDOW) FOR EACH PARCEL
%%%    the goal here is to peak the minimum number of knots that nearly maximize the r-squared (find a "good-enough" fit)
%%%    (this concept is similar to the elbow rule used in the clustering world)

%%% 2) COMPUTE THE SPLINES FITTING BASED ON THE BEST NUMBER OF KNOTS (WINDOW) FOR EACH PARCEL

%% 1) DEFINING BEST NUMBER OF KNOTS (WINDOW) FOR EACH PARCEL

labblock = 3; % 1 = 0.1-1Hz OLD; 2 = 2-8Hz OLD; 3 = 2-8Hz NEW
% TO REPRODUCE IMAGEs IN THE FIGURE S8 USE:
% 1) labblock = 1; for ii = 2 (over parcels);
% 2) labblock = 2; for ii = 9 (over parcels)

if labblock == 1
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_1_Average_l_0_Old_l_1/Clusters.mat')
    knots_num = 20; %number of knots to be tested
elseif labblock == 2
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/OLD.mat')
    KJK = lko; %just not to change the script later..
    knots_num = 45; %number of knots to be tested
else
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/NEW.mat')
    KJK = lkn; %just not to change the script later..
    knots_num = 45; %number of knots to be tested
end
x = 1:size(KJK,2); %time-points
%%% 1) DEFINING BEST NUMBER OF KNOTS (WINDOW) FOR EACH PARCEL
RR = zeros(size(KJK,1),knots_num);
for ii = 1:size(KJK,1) %over parcels
    y = KJK(ii,:); %selecting timeseries for parcel ii
    for yy = 2:knots_num %over knots (must start from 2 since at least 2 points (knots) are required to define a window
        [slm,xp,yp] = slmengine(x,y,'plot','off','knots',yy,'deg','3','res','pp');
%         set(gcf,'color','w')
%         grid minor
        RR(ii,yy) = slm.stats.R2Adj; %storing adjusted R-squared (mean over the different windows of a specific knots-solution)
    end
    %plotting R-squared as a function of number of knots
    figure
    tim = 1:knots_num;
    plot(tim,RR(ii,:))
    set(gcf,'color','w')
    grid minor
    title(['Parcel ' num2str(ii)])
end

%% 2) COMPUTE THE SPLINES FITTING BASED ON THE BEST NUMBER OF KNOTS (WINDOW) FOR EACH PARCEL

%automatically creating an excel file with information
%option for images (if you want them)
splines_image = 'off'; %'on' for images after splines; 'off' for no images

for bb = 1:3
    labblock = bb; % 1 = 0.1-1Hz OLD; 2 = 2-8Hz OLD; 3 = 2-8Hz NEW
    %defining the block and some information/data
    if labblock == 1
        load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/K_means/Slowneg_1_Average_l_0_Old_l_1/Clusters.mat')
        %vector with "ideal" number of knots for each parcel isolated in the previous section
        vect = [2,8,8,5,4,4,4,8,6,4,6,6,4,8,6,7,6,7,7]; %0.1-1Hz OLD
    elseif labblock == 2
        load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/OLD.mat')
        KJK = lko; %just not to change the script later..
        %vector with "ideal" number of knots for each parcel isolated in the previous section
        vect = [27,28,27,28,28,26,27,28,28]; %2-8Hz OLD
    else
        load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/ToPython_2_8Hz/NEW.mat')
        KJK = lkn; %just not to change the script later..
        %vector with "ideal" number of knots for each parcel isolated in the previous section
        vect = [28,28,27,28,28,27,27,28,27]; %2-8Hz NEW
    end
    %actual computation of splines
    clear PDn
    x = 1:size(KJK,2); %time-points
    cnt7 = -6;
    for ii = 1:size(KJK,1) %over parcels
        y = KJK(ii,:); %selecting timeseries for parcel ii
        [slm,xp,yp] = slmengine(x,y,'plot',splines_image,'knots',vect(ii),'deg','3','res','pp');
        xlim([1 316])
        set(gcf,'color','w')
        grid minor
        title(['Parcel ' num2str(ii)])
        %storing values (r-squared and coefficients from splines) in a table to be then used to print an excel file
        cnt7 = cnt7 + 7; %1stparcel(time-windows,adjR,a,b,c,d,1spaceempty),THEN 2nd parcel
        PDn{1,cnt7} = ['Parcel ' num2str(ii)]; %PDn{1,2} = 'Hemisphere'; PDn{1,3} = 'T-stat'; PDn{1,4} = 'MNI coordinates'; %1st row
        PDn{2,cnt7} = 'Time-windows'; PDn{2,cnt7+1} = 'Adj R squared'; PDn{2,cnt7+2} = 'a'; PDn{2,cnt7+3} = 'b'; PDn{2,cnt7+4} = 'c'; PDn{2,cnt7+5} = 'd'; %2nd row
        for jj = 1:(vect(ii)-1) %over splines time-windows (-1 since vect contains the knots, while the time-windows correspond to number of knots -1)
            PDn{2+jj,cnt7} = jj; %time-window
            PDn{2+jj,cnt7+1} = slm.stats.R2Adj; %adjusted R-squared
            PDn{2+jj,cnt7+2} = slm.coefs(jj,1); %a
            PDn{2+jj,cnt7+3} = slm.coefs(jj,2); %b
            PDn{2+jj,cnt7+4} = slm.coefs(jj,3); %c
            PDn{2+jj,cnt7+5} = slm.coefs(jj,4); %d
        end
    end
    %saving files in different excel sheets (1 = 0.1-1Hz OLD; 2 = 2-8Hz OLD; 3 = 2-8Hz NEW)
    PDnn = cell2table(PDn); %remove the possible empty cell
    writetable(PDnn,['/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/Tables/Table_8/Splines_01_1HzOLD_1Sheet_2_8HzOLD_2Sheet_2_8HzNEW_3Sheet.xlsx'],'Sheet',bb) %printing excel file
end

%%
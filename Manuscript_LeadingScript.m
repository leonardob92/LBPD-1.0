%% PAPER "Spatiotemporal whole-brain dynamics of auditory patterns recognition" - LEONARDO BONETTI


%BEFORE PROCEEDING, PLEASE NOTE:

%#1
%the pre-proessing steps and some of the first analyses are reported to
%ensure full disclosure, however the raw (and in some cases the pre-processed) data is not provided because of
%their large sizes and of issues with potential disclosure of personal information.
%Data is provided when we were sure that it was completely anonymized and
%could not be linked to the original personal information of the
%participants, according to the regulations.

%#2
%Several steps of the pre-processing and data analysis were computed on an
%external server of Aarhus University (AU) to benefit from a strong
%computational power. Then, to make the scripts available, understandable
%and usable for everybody, I adapted the essential parts of them for running on individual computers.
%Since some of those operations are computationally very expensive, I have
%not rerun every step again, even if I have carefully checked every line of codes.
%However, if you notice a few typos, you are very welcome to let me know so that I can fix them.
%Of course, I am also available to provide more information about the codes
%and help to set them up, especially when they require external functions, if anybody is interested.
%Thanks.

%#3
%To use these functions and test/replicate the results/plottings reported
%in the paper, you need to donwload the LBPD functions and store them in a folder of your own computer (here called "pathl").
%Remember that some of the following scripts require the OSL toolbox. To donwload that
%follow this link: https://ohba-analysis.github.io/osl-docs/
%Remember also that you need to have iOS or Linux to use OSL.
%Once you download OSL, please copy it under 'pathl/LBPD/External/'.
%If you have Microsoft Windows you can only run the scripts that do not explicitely
%require OSL functions.
%Further, in that case, you should donwload FieldTrip and SPM.
%You find FielTrip (the current LBPD codes rely on the version: fieldtrip-20171231) here:
%http://www.fieldtriptoolbox.org/
%You find SPM here:
%https://www.fil.ion.ucl.ac.uk/spm/
%Then, please copy FieldTrip and SPM in the directory: 'pathl/LBPD/External/',
%Finally, to start LBPD up, you only need to:
%1)run: addpath('YOUR_OWN_PATH/LBPD')
%2)run: LBPD_startup_D
%After that, you will be able to run the parts of this script that do not
%require AU server and therefore reproduce the results and
%plots that we provided (read ReadMe.mat for additional information).
%Please, cite these codes and results referrring to the following paper:
% https://www.biorxiv.org/content/10.1101/2020.06.23.165191v2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Please note that LBPD functions rely on some external freely available functions/files
%that we provide but that we did not write. More information about this is provided in the ReadMe.mat file that I strongly advise to read.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% START UP FUNCTIONS AND SETTING UP PATHS.. (LBPD_startup_D)

%%% PLEASE REMEMBER TO RUN THIS SECTION IN ORDER TO BE ABLE TO RUN THE FOLLOWING PARTS OF THIS SCRIPT %%%

%this function was made to add some paths related to required external functions
%if you want to use/test the package of functions that I provided, you
%need to donwload the functions and then store them in a path that in this
%example script is called "pathl"

%please specify the path in your own computer where you stored the provided data
pathdata = '/Users/au550322/Documents/AarhusUniversitet/PhDThesis/PaperFinal1/Codes_Data_NN1/Data';
%starting up some of the functions that I wrote for the following analysis
% pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
pathl = '/Users/au550322/Documents/AarhusUniversitet/PhDThesis/PaperFinal1/PNAS/FirstRevision/Codes_27_05_2021/LBPD';
addpath(pathl);
LBPD_startup_D(pathl);

%%



%% *** PRE-PROCESSING ***

%the pre-processing steps are reported for full disclosure purposes but they cannot be run here since they were developed for the parallel computing of the external server of Aarhus University (AU)

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
maxDir = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach_Fsample1000Hz'; %with movement compensation
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
%number of data files to optmize the codes. Some of them were not used in this specific paper
%but in other projects, however these analysese have been run together.

%% setting directories

clear
datadir = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc'; %with movement compensation in maxfilter
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
% removing eyeblink and heartbeat artefacts by computing ICA
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
xlsx_dir_behav = '/scratch5/MINDLAB2017_MEG-LearningBach/BehavioralResponses'; %dir to MEG behavioral results
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
%additional low-pass filter
datadir = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc';
subjnum = {'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'38';'39';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';'51';'52';'53';'54';'55';'56';'57';'58';'59';'60';'61';'62';'63';'64';'65';'66';'67';'68';'69';'70';'71'}; %Everybody!!!!!!
for ii = 1:length(subjnum) % iterates over blocks !OBS!
    S2 = [];
    S2.D = [datadir '/e80ldffspmeeg_SUBJ00' subjnum{ii} 'recogminor_tsssdsm.mat'];
    S2.band = 'low';
    S2.freq = 40;
    jobid = job2cluster(@filtering,S2);
    disp(ii)
end

%%% the preprocessing until here was done for all of the analysis steps %%%
%%% the decoding and beamforming sections started after the above steps %%%

%%% the indipendent univariate tests and successive Monte Carlo simulations
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
list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/fe80ldff*');
for ii = 2:2:length(list) %over files
    %distribute 
    input = [];
    input.D = [list(ii).folder '/' list(ii).name];
    input.prefix = output_prefix_to_be_set;
    jobid = job2cluster(@sensor_average, input);
end
%combining planar gradiometers
list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/mfe80ldff*');
for ii = 2:2:length(list) %over files
    input = [];
    input.D = [list(ii).folder '/' list(ii).name];
    D = spm_eeg_load(input.D);
    D = D.montage('switch',1);
    D.save();
    jobid = job2cluster(@combining_planar_cluster, input);
end

%% *** ANALYSIS ***

%% decoding - support vector machine (SVM) - multivariate pattern analysis

%The decoding consisted of support vector machines (SVM) implemented in
%external functions that can be found at the following link:
%http://www.csie.ntu.edu.tw/~cjlin/libsvm/
%Those functions can be also used in Matlab with a few adjustments done to achieve a proper implementation
%as described in the paper and including crucial steps such as, for example, cross validation.
%In our case, we implemented them by using AU server, therefore the codes immediately below show how we provided the data to the AU server.
%After that, we provide the output of the SVM algorithm so that you can reproduce
%the statistics and the plotting solutions that we presented in the paper.
%If you wish to know more about this SVM implementation, you are
%very welcome to contact us.


%% settings for AU server and SVM algorithm

cd /projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab/external');
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab/external/libsvm-3.21')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab/external/libsvm-3.21/matlab')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/SpatialLocation')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/SpatialLocation/plotchannel')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/permutationlab')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/permutationlab/private1')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/permutationlab/developer')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
make

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
% config of the cluster server
clusterconfig('slot', 1); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
clusterconfig('scheduler','cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue

%% preparing data for decoding functions (for AU server)

time = [286 751]; %defining time of interest (this refers to 0.100sec before the onset of the stimuli)
numperm = 100; %number of permutations
kfold = 5; %number of categories for k-fold classification
subjnum = {'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'38';'39';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';'51';'52';'53';'54';'55';'56';'57';'58';'59';'60';'61';'62';'63';'64';'65';'66';'67';'68';'69';'70';'71'}; %Everybody!!!!!!
datadir = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/';
savedir = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper1/Decoding/Lowpass_40Hz/';
D1 = cell(length(subjnum),2);
subjs_real = cell(1,length(subjnum));
cty = 0;
for ii = 1:1:length(subjnum) %length(subjnum)
    if ii ~= 39 %number 39 is simply a non-existent subject, however for computational reasons we decided to not update the ID of all other subjects and we preferred to simply leave this subject empty and ignore it in the computations
        %loading data
        D = spm_eeg_load([datadir 'fe80ldffspmeeg_SUBJ00' subjnum{ii} 'recogminor_tsssdsm.mat']);
        old = find(strcmp('Old_Correct',D.conditions)); %getting condition old indices
        new = find(strcmp('New_Correct',D.conditions)); %getting condition new indices
        clear data condid ind
        condid(1:length(old)) = {'O'}; %creating condition labels (first old)
        condid(length(old):(length(old)+length(new))) = {'N'}; %(then new)
        ind = [old new]; %rearranging indices (first old and then new)
        idxMEGc1 = find(strcmp(D.chanlabels,'MEG0111')); %getting extreme channels indexes
        idxMEGc2 = find(strcmp(D.chanlabels,'MEG2643')); %same
        data = D(idxMEGc1:idxMEGc2,time(1):time(2),ind); %extracting data
        %structure with inputs for AU server
        S = [];
        S.data = data;
        S.condid = condid;
        S.numperm = numperm;
        S.kfold = kfold;
        S.savedir = savedir;
        S.ii = ii;
        disp(ii)
        %job to cluster
        jobid = job2cluster(@decoding, S);
    end
    if ii == 1 %saving D.time (time in seconds), only once
        time_sel = D.time(time(1):time(2));
        save([savedir 'time.mat'],'time_sel');
    end
end

%% permutation test on decoding accuracy time-series (pairwise-decoding and temporal generalization)

% *** YOU CAN RUN THIS SECTION AND REPRODUCE RESULTS/PLOTS REPORTED IN THE PAPER ***

%directory where decoding results are stored in
datadir = [pathdata '/MEG_sensordata_Decoding_Lowpass_40Hz/']; %path to decoding data
%loading and reshpaping data outputted by AU server (SVM algorithm) for statistics
load([datadir 'time.mat']);
%loading data
dd1 = zeros(length(time_sel),70);
ddpatt = zeros(306,length(time_sel),70);
ddTG = zeros(length(time_sel),length(time_sel),70);
cnt = 0;
for ii = 1:71
    if ii ~= 39 %number 39 is simply a non-existent subject, however for computational reasons we decided to not update the ID of all other subjects and we preferred to simply leave this subject empty and ignore it in the computations
        cnt = cnt + 1;
        load([datadir 'PD_SUBJ' num2str(ii) '.mat']);
        dd1(:,cnt) = d.d;
        ddpatt(:,:,cnt) = d.pattern;
        load([datadir 'TG_SUBJ' num2str(ii) '.mat'])
        ddTG(:,:,cnt) = d.d;
    end
    disp(ii)
end
%average over subjects of pairwise decoding accuracy time-series
dd12 = mean(dd1,2);
%reshaping pairwise decoding acuracy time-series for submitting them to statistical testing
dd2 = dd1 - 50; %subtracting 50 since the function tests significance against 0 (here represented by 50% chance level)
dat = cell(1,size(dd2,2));
for ii = 1:size(dd2,2)
   dat(1,ii) = {dd2(:,ii)'};
end
stat = pl_permtest(dat,'alpha',0.05); %actual function
%plotting p-values after FDR correction for multiple comparisons
xx = stat.FDR.criticalmap;
CCt = bwconncomp(xx); %getting significant 'islands'
%plotting average over subjects (pairwise decoding accuracy) and the corresponding significance
figure;hold on;grid on;grid minor;
plot(time_sel,dd12,'k','linewidth',2); xlabel('Time (s)');ylabel('Decoding accuracy (%)');
hold on
% Show significant time window
patch_color = [.85 .85 .85]; % Color of grey box marking time range
ylims = get(gca,'YLim');
sgf = CCt.PixelIdxList;
for ii = 1:length(sgf)
    sgf2 = sgf{ii}./150;
    patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.7)
    hold on
end
plot(time_sel,dd12,'k','linewidth',2); xlabel('Time (s)');ylabel('Decoding accuracy (%)'); %(slightly stupid..) trick to get the time-series in front..
set(gcf,'Color','w')
xlim([time_sel(1) time_sel(end)])
%plotting patterns (derived from weights) for pairwise decoding
signt = [0.5 2.0]; %significant time-point(s) you want to plot
%extracting magnetometers for plotting purposes
dp = mean(ddpatt,3); %average over subjects
avg = dp(1:3:end,:); %extracting magnetometers only
%plotting topoplot (fieldtrip function)
%creating the mask for the data
fieldtrip_example = load([pathl '/External/fieldmask.mat']);
label2 = fieldtrip_example.M_timb_lock_gradplanarComb.label;
label = label2(103:end);
%setting labels according to the ones in the layout
for ii = 1:102
    label{ii,1} = [label{ii,1}(1:3) label{ii,1}(5:end)];
end
cfgdummy = fieldtrip_example.M_timb_lock_gradplanarComb.cfg;
cfgdummy.previous = [];
data = [];
data.cfg = cfgdummy;
data.time(1,:) = time_sel(:);
data.label = label;
data.dimord = fieldtrip_example.M_timb_lock_gradplanarComb.dimord;
data.grad = fieldtrip_example.M_timb_lock_gradplanarComb.grad;
data.avg = avg;
%creating the cfg for the actual plotting (magnetometers)
cfg = [];
cfg.layout = [pathl '/External/neuromag306mag.lay'];
cfg.colorbar = 'yes';
cfg.xlim = signt; %set temporal limits (in seconds)
cfg.zlim = [];
cfg.colormap = 'jet';
figure
ft_topoplotER(cfg,data);
set(gcf,'Color','w')
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))
%plotting temoral generalization
ddTGm = mean(ddTG,3);
figure; imagesc(time_sel,time_sel,ddTGm); xlabel('Train time (s)'); ylabel('Test time (s)'); set(gca,'YDir','normal');
colorbar
set(gcf,'Color','w')
%reshaping temporal generalization decoding acuracy for submitting it to statistical testing
ddTG2 = ddTG - 50; %subtracting 50 since the function tests significance against 0 (here represented by 50% chance level)
dat = cell(1,size(ddTG2,3));
for ii = 1:size(ddTG2,3)
   dat(1,ii) = {ddTG2(:,:,ii)};
end
%plotting results
stat = pl_permtest(dat,'alpha',0.02); %actual function
figure; imagesc(time_sel,time_sel,stat.FDR.criticalmap); xlabel('Train time (s)'); ylabel('Test time (s)'); set(gca,'YDir','normal');
set(gcf,'Color','w')

%% MCS on sensor data

%The following section calls "MEG_sensors_plotting_ttest_LBPD_D2", a
%function that is designed with 2 scopes:
%1)extracting data from SPM objects and calculating t-tests between experimental conditions
%2)plotting data (potentially taking into account also some statistical
%results to be overlayed in the plots).
%Therefore, this function has been used TWICE.
%The first time to extract the data from SPM objects (this is not
%reproducible since SPM objects contained information that may be connected
%to the participants and are quite large in terms of size).
%The second time to plot the data and the correspondent statistics.
%Thus, you will find a first implementation of the function "MEG_sensors_plotting_ttest_LBPD_2" that you cannot
%use, then you will see the statistical functions (THAT YOU CAN RUN TO
%REPRODUCE OUR RESULTS). Finally, you will see again
%"MEG_sensors_plotting_ttest_LBPD_D2" and this time you will BE ABLE TO RUN
%IT AND PRODUCE THE PLOTS THAT WE REPORTED IN THE PAPER


%% extracting MEG sensor data

%original settings to build the path to the data
subjnum = {'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'38';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';'51';'52';'53';'54';'55';'56';'57';'58';'59';'60';'61';'62';'63';'64';'65';'66';'67';'68';'69';'70';'71'}; %Everybody!!!!!!
datadir = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc'; %path to data
outdir = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper1/MCS_Sensor/lowpass40Hz'; %path to write output in
S = [];
%computing data
S.outdir = outdir;
S.data = [];
S.spm_list = cell(1,length(subjnum));
for ii = 1:length(subjnum)
    S.spm_list(ii) = {[datadir '/fPme80dffspmeeg_SUBJ00' subjnum{ii} 'recogminor_tsssdsm.mat']};
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
S.waveform_average_label = 1;
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

%The first analysis has been conducted for gradiometers only. Then, we have
%computed a second analysis on magnetometers by using only the significant
%time-points emerged from gradiometers analysis.
%This was done to avoid sign ambiguity connected to magnetometers which are
%the MEG sensors that keep their polarity (differently from gradiometers
%that have been combined by sum-root square for this analysis).
%Further details are available in the paper.
%Therefore, in this script, you should reproduce AT FIRST the gradiometers
%and THEN the magnetometers. You can do that by simply selecting the proper
%label (then everything is already set up).
%Please, note that this will take some minutes to run.

%input information
grad = 1; %1 for gradiometers; 0 for magnetometers; time-points to be considered for MCS
contrn = 1; %set1 to load previously calculated t-tests memorized musical patterns vs novel musical patterns; set 0 for vice-versa)


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
outdir = [pathdata '/MEG_sensordata_UnivariateTests_MonteCarloSimulations']; %path where t-test results are stored in
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
[~,~,raw_channels] = xlsread([pathl '/External/MatrixMEGChannelLayout_With0_2.xlsx']);
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


%% extracting first and last significant time-point for each channel, converting them into seconds (from time-samples) and printing it as xlsx file 

%Here you can print the statistics in tables/excel files that can be
%convenient to read. You can simply select the significant cluster number
%that you want by specifying it in the line below.
%THIS RESULTS CAN BE FOUND IN THE PAPER IN SUPPLEMENTARY MATERIALS

%input parameters for user 
clustnum = 1;
typec = 1; %1 for gradiometers; 2 for magnetometers positive, 3 for magnetometers negative


%actual computation
load([pathdata '/MEG_sensordata_UnivariateTests_MonteCarloSimulations/Grad_clust.mat']);
load([pathdata '/MEG_sensordata_UnivariateTests_MonteCarloSimulations/MagPos_clust.mat']);
load([pathdata '/MEG_sensordata_UnivariateTests_MonteCarloSimulations/MagNeg_clust.mat']);
if typec == 1
    hh = GRAD_clust{clustnum,3};
elseif typec == 2
    hh = MAG_clust_pos{clustnum,3}; %defining cluster
else
    hh = MAG_clust_neg{clustnum,3}; %defining cluster
end

time_sel2 = time_sel(1,min_time_point:max_time_point);
clear ff2
ff2(:,1) = hh(:,1); %extracting channel names
for ii = 1:size(ff2,1)
    ff2(ii,2) = {time_sel2(hh{ii,2}(1))}; %first significant time-point
    ff2(ii,3) = {time_sel2(hh{ii,2}(end))}; %last significant time-point
end
PDn = cell2table(ff2); %converting cell to table
% writetable(PDn,['NegPos_clust1_OldvsNew.xlsx'],'Sheet',1) %saving xlsx file

%% plotting

%Here you can reproduce the plots repored in the paper.
%Please note that the waveforms reported in the figure and that you can
%obtain by using the following function and the already defined settings
%relate to gradiometers data.

contr_topoplot = 3; %1 to plot the topoplot showing the contrast memorized vs novel musical patterns; 2 to plot memorized musical patterns; 3 to plot novel musical patterns


%actual computation
outdir = [pathdata '/MEG_sensordata_UnivariateTests_MonteCarloSimulations']; %this must become your own directory to the data that we provided
load([outdir '/sensor_data.mat']);
load([outdir '/Grad_clust.mat']);

S = [];
%computing data
S.outdir = outdir;
S.data = [];
S.data = data_mat;
S.chanlabels = chanlabels;
S.time_real = time_sel;
S.conditions = {'Old_Correct','New_Correct'};
S.timeextract = []; %time-points to be extracted
S.centerdata0 = 0; %1 to make data starting at 0
S.save_data = 0; %only meaningfull if you read data from SPM objects saved on disk
S.save_name_data = 'sensor_data';
%individual waveform plotting (settings for individual channel plots are not used in the paper..)
S.waveform_singlechannels_label = 0; %1 to plot single channel waveforms
S.wave_plot_conditions_together = 0; %1 for plotting the average of all
S.mag_lab = 2; %1 for magnetometers; 2 for gradiometers
S.x_lim_temp_wave = []; %limits for time (in secs) (E.g. [-0.1 3.4])
S.y_lim_ampl_wave = []; %limit for amplitude (E.g. [0 120] magnetometes, [0 6] gradiometers)
%averaged waveform plotting
S.waveform_average_label = 1;
S.legc = 0; %set 1 for legend
clustplot = GRAD_clust{1,3}; %slightly elaborated way to get the channel original IDs of channels forming the significant cluster that you are considering
S.left_mag = [];
for ii = 1:size(clustplot,1)
    S.left_mag(ii) = find(cellfun(@isempty,strfind(chanlabels,clustplot{ii,1})) == 0);
end
S.signtp(1) = {[0.547 1.18]}; %in time samples
S.sr = 150; %sampling rate (Hz)
S.avewave_contrast = 1; %1 to plot the contrast between conditions (averaged waveform)
S.save_label_waveaverage = 0;
S.label_plot = 'block_minor';
%t-tests
S.t_test_for_permutations = 0;
S.cond_ttests_tobeplotted_topoplot = [1 2]; %this is for both topoplot and t-tests!! (here [1 2] means cond1 vs cond2!!!!!!!)
%topoplotting
S.topoplot_label = 1;
S.fieldtrip_mask = [pathl '/External'];
if contr_topoplot == 1
    S.topocontr = 1;
else
    S.topocontr = 0;
end
if contr_topoplot == 2
    S.topocondsing = [1];
elseif contr_topoplot == 3
    S.topocondsing = [2];
end
S.xlim = [0.547 1.18]; %time topolot
if contr_topoplot == 1
    S.zlimmag = [-4 4]; %magnetometers amplitude topoplot limits
    S.zlimgrad = [-4 4]; %gradiometers amplitude topoplot limits
else
    S.zlimmag = [-150 150]; %magnetometers amplitude topoplot limits
    S.zlimgrad = [0 4]; %gradiometers amplitude topoplot limits
end
S.colormap_spec = 0;
% x = []; x.bottom = [0 0 1]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 1 0.5]; x.top = [1 0.95 0]; %yellow - blue
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
S.colormap_spec_x = x;
S.topoplot_save_label = 0;
S.color_line = [1 0 0; 0 0 1];

[out] = MEG_sensors_plotting_ttest_LBPD_D2(S);

%% PLOTTING PARTICIPANTS T-TESTS BETWEEN OLD AND NEW ACCORDING TO MUSICIANSHIP, WM ABILITIES, LIKING AND FAMILIARITY WITH THE BACH'S MUSIC

%this requires the section above to be run first
%you can here reproduce the results reported in the supplementary materials showing the few significant associations between neural data
%underlying musical patterns recognition and behavioral measures indexing memory and musical skills
%please, note that the plot concering working memory (WM) is reported in the main text of the paper, while the other measures in supplementary materials
%select the plot that you want
plot_l = 4; %2)musicianship; 3)WM abilities; 4)Liking; 5)Familiarity; 6)Goldsmith - musical sophistication index;

%loading t-tests
load([pathdata '/MEG_sensordata_UnivariateTests_MonteCarloSimulations/sensor_data.mat']);
%reading excel file with: 1)IDs; 2)musicianship; 3)WM abilities; 4)Liking; 5)Familiarity; 6)Goldsmith - musical sophistication index; 7) and 8)Old and New RTs; 9)Old and New (together) correct responses (out of 80 total) ..only few participants (4-8) previously studied the piece and partially remember it..
[num,~,raw] = xlsread([pathdata '/MEG_sensordata_UnivariateTests_MonteCarloSimulations/FirstRevisionPNAS.xlsx']);
%MEG of SUBJ 39 is not available and therefore I erased it from the excel file containing the behavioural measures (important to remember)
% dataold = squeeze(nanmean(nanmean(data_mat(S.left_mag,:,:,1),3),1)); %old
% datanew = squeeze(nanmean(nanmean(data_mat(S.left_mag,:,:,2),3),1)); %new
dataold = squeeze(nanmean(data_mat(S.left_mag,:,:,1),1)); %old
datanew = squeeze(nanmean(data_mat(S.left_mag,:,:,2),1)); %new
datadiff = transpose(dataold-datanew);
if plot_l == 2
    %2)musicianship (PROBABLY NOT SIGNIFICANT)
    GROUP = cell(70,1);
    for ii = 1:size(num,1)
        if num(ii,2) == 1
            GROUP{ii} = 'pian';
        elseif num(ii,2) == 2
            GROUP{ii} = 'mus';
        else
            GROUP{ii} = 'nonmus';
        end
    end
    if ~exist('PANOVA','var')
        PANOVA = zeros(size(datadiff,2),1); C = zeros(3,6,size(datadiff,2));
        for tt = 1:size(datadiff,2)
            [p,~,stats] = anova1(datadiff(:,tt),GROUP,'off');
            [c,~,~,~] = multcompare(stats);
            PANOVA(tt) = p; %storing ANOVA significance (p-value) for each time point and each channel
            C(:,:,tt) = c; %storing post-hoc results - c is a matrix in which rows contain the result of one paired comparison test, col.1 and 2 contain the indices of the samples that are being compared, col 3 the lower CI, col.4 the estimate, col 5 the upper CI and col 6 the p-value
            disp(tt)
        end
    end
    P = PANOVA; P(P<0.05) = 2; P(P<2) = 0; P(P==2) = 1; %barbaric but effective solution..
    [ sign_clust ] = oneD_MCS_LBPD_D( P, 0, 1, 1000, 0.001 );
    COL(1) = 'b'; COL(2) = 'r'; COL(3) = 'k';
    figure
    for ii = 1:3 %over musicianship groups
        plot(time_sel,mean(datadiff(num(:,2)==ii,:),1),'Color',COL(ii),'LineWidth',2)
        hold on
        plot(time_sel,(mean(datadiff(num(:,2)==ii,:),1) + (std(datadiff(num(:,2)==ii,:),0,1)./sqrt(length(num(:,2)==ii)))),'-.','Color',COL(ii),'LineWidth',0.1)
        hold on
        plot(time_sel,(mean(datadiff(num(:,2)==ii,:),1) - (std(datadiff(num(:,2)==ii,:),0,1)./sqrt(length(num(:,2)==ii)))),'-.','Color',COL(ii),'LineWidth',0.1)
        hold on
    end
    %plotting the significant time-points with gray shades
    patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
    ylims = get(gca,'YLim');
    for ii = 1:size(sign_clust,1)
        sgf2 = time_sel(sign_clust{ii,3});
        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.5)
        hold on
    end
    for ii = 1:3 %over musicianship groups
        plot(time_sel,mean(datadiff(num(:,2)==ii,:),1),'Color',COL(ii),'LineWidth',2)
        hold on
        plot(time_sel,(mean(datadiff(num(:,2)==ii,:),1) + (std(datadiff(num(:,2)==ii,:),0,1)./sqrt(length(num(:,2)==ii)))),'-.','Color',COL(ii),'LineWidth',0.1)
        hold on
        plot(time_sel,(mean(datadiff(num(:,2)==ii,:),1) - (std(datadiff(num(:,2)==ii,:),0,1)./sqrt(length(num(:,2)==ii)))),'-.','Color',COL(ii),'LineWidth',0.1)
        hold on
    end
    set(gcf,'Color','w')
    grid minor
    xlim([time_sel(1) time_sel(end)])
%     export_fig(['musicianship_PNAS.eps'])
%     export_fig(['musicianship_PNAS.png'])
elseif plot_l == 3
    %3)WM
    dataold = dataold'; datanew = datanew';
    rwm = zeros(1,size(datadiff,2)); pwm = zeros(1,size(datadiff,2));
    for tt = 1:size(datadiff,2)
        [rwm(tt),pwm(tt)] = corr(num(:,3),datadiff(:,tt),'rows','complete');
    end
    P = pwm; P(P<0.05) = 2; P(P<2) = 0; P(P==2) = 1; %barbaric but effective solution..
    %Monte Carlo simulation function to correct for multiple comparisons
    [ sign_clust ] = oneD_MCS_LBPD_D( P, 0, 1, 1000, 0.001 );
    figure
    plot(time_sel,rwm,'Color','k','LineWidth',2)
    set(gcf,'Color','w'); xlim([time_sel(1) time_sel(end)])
    hold on
    %plotting the significant time-points with gray shades
    patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
    %     ylims = get(gca,'YLim');
    ylims = [-0.45 0.45];
    for ii = 1:size(sign_clust,1)
        sgf2 = time_sel(sign_clust{ii,3});
        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.5)
        hold on
    end
    plot(time_sel,rwm,'Color','k','LineWidth',2)
    set(gcf,'Color','w'); xlim([time_sel(1) time_sel(end)])
    hold on
    grid minor
    ylim([-0.45 0.45])
%     export_fig(['WM_PNAS.eps'])
%     export_fig(['WM_PNAS.png'])
elseif plot_l == 4
    %4)Liking
    rlike = zeros(1,size(datadiff,2)); plike = zeros(1,size(datadiff,2));
    for tt = 1:size(datadiff,2)
        [rlike(tt),plike(tt)] = corr(num(:,4),datadiff(:,tt),'rows','complete');
    end
    P = plike; P(P<0.05) = 2; P(P<2) = 0; P(P==2) = 1; %barbaric but effective solution..
    %Monte Carlo simulation function to correct for multiple comparisons
    [ sign_clust ] = oneD_MCS_LBPD_D( P, 0, 1, 1000, 0.001 );
    figure
    plot(time_sel,rlike,'Color','k','LineWidth',2)
    set(gcf,'Color','w'); xlim([time_sel(1) time_sel(end)])
    hold on
    %plotting the significant time-points with gray shades
    patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
    %     ylims = get(gca,'YLim');
    ylims = [-0.45 0.45];
    for ii = 1:size(sign_clust,1)
        sgf2 = time_sel(sign_clust{ii,3});
        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.5)
        hold on
    end
    plot(time_sel,rlike,'Color','k','LineWidth',2)
    set(gcf,'Color','w'); xlim([time_sel(1) time_sel(end)])
    hold on
    grid minor
    ylim([-0.45 0.45])
    %export_fig(['liking_PNAS.eps'])
    %export_fig(['liking_PNAS.png'])
elseif plot_l == 5
    %5)Familiarity
    rfam = zeros(1,size(datadiff,2)); pfam = zeros(1,size(datadiff,2));
    for tt = 1:size(datadiff,2)
        [rfam(tt),pfam(tt)] = corr(num(:,5),datadiff(:,tt),'rows','complete');
    end
    P = pfam; P(P<0.05) = 2; P(P<2) = 0; P(P==2) = 1; %barbaric but effective solution..
    %Monte Carlo simulation function to correct for multiple comparisons
    [ sign_clust ] = oneD_MCS_LBPD_D( P, 0, 1, 1000, 0.001 );
    figure
    plot(time_sel,rfam,'Color','k','LineWidth',2)
    set(gcf,'Color','w'); xlim([time_sel(1) time_sel(end)])
    hold on
    %plotting the significant time-points with gray shades
    patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
    ylims = [-0.45 0.45];
    for ii = 1:size(sign_clust,1)
        sgf2 = time_sel(sign_clust{ii,3});
        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.5)
        hold on
    end
    plot(time_sel,rfam,'Color','k','LineWidth',2)
    set(gcf,'Color','w'); xlim([time_sel(1) time_sel(end)])
    hold on
    grid minor
    ylim([-0.45 0.45])
    %export_fig(['famil_PNAS.eps'])
    %export_fig(['famil_PNAS.png'])
elseif plot_l == 6
    %6)Goldsmith - musical sophistication index
    rgol = zeros(1,size(datadiff,2)); pgol = zeros(1,size(datadiff,2));
    for tt = 1:size(datadiff,2)
        [rgol(tt),pgol(tt)] = corr(num(:,6),datadiff(:,tt),'rows','complete');
    end
    P = pgol; P(P<0.05) = 2; P(P<2) = 0; P(P==2) = 1; %barbaric but effective solution..
    %Monte Carlo simulation function to correct for multiple comparisons
    [ sign_clust ] = oneD_MCS_LBPD_D( P, 0, 1, 1000, 0.001 );
    figure
    plot(time_sel,rgol,'Color','k','LineWidth',2)
    set(gcf,'Color','w'); xlim([time_sel(1) time_sel(end)])
    hold on
    %plotting the significant time-points with gray shades
    patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
    ylims = [-0.45 0.45];
    for ii = 1:size(sign_clust,1)
        sgf2 = time_sel(sign_clust{ii,3});
        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.5)
        hold on
    end
    plot(time_sel,rgol,'Color','k','LineWidth',2)
    set(gcf,'Color','w'); xlim([time_sel(1) time_sel(end)])
    hold on
    grid minor
    ylim([-0.45 0.45])
    %export_fig(['goldsm_PNAS.eps'])
    %export_fig(['goldsm_PNAS.png'])
end

%%

%% *** SOURCE RECONSTRUCTION (BEAMFORMING AND STATISTICS) - BRAIN ACTIVITY ***

%Please note that the first section works on SPM objects containing some
%potentially personal information related to the participants and it is
%computationally quite expensive. Therefore, we provided the
%data/results/images to be reproduced after the first pre-processing and
%data analysis steps.
%Please remember that more information on the beamforming algorithm and subsequent statistics can be
%found in the OSL website/documentation: %https://ohba-analysis.github.io/osl-docs/

%starting up OSL and adding paths
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup
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

%beamforming
count1 = 0;
oat = [];
%SOURCE RECON
for ii = 3:3:213  % = 1:length(spm_files_recog_basen)
    count1 = count1 + 1;
    spm_files = ['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/e80dff'  spm_files_recog_basen{ii}]; %check this in order not to mix up 'epoched' and 'continuous'
    D_epoched = spm_eeg_load(spm_files);
    D_epoched = D_epoched.montage('switch',1); %switch the montage to 1 in order to be safer
    D_epoched.save();
    processed_file_epoched = spm_files;
    spm_files = ['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/dff'  spm_files_recog_basen{ii}]; %check this in order not to mix up 'epoched' and 'continuous'
    D_continuous = spm_eeg_load(spm_files);
    D_continuous = D_continuous.montage('switch',1); %switch the montage to 1 in order to be safe
    D_continuous.save();
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
oat.source_recon.sessions_to_do = [1:38 40:71]; %sessions to do among the file_list; number 39 is simply a non-existent subject, however for computational reasons we decided to not update the ID of all other subjects and we preferred to simply leave this subject empty and ignore it in the computations
oat.source_recon.modalities = {'MEGMAG'; 'MEGPLANAR'}; %'MEGMAG'};
oat.source_recon.conditions = {'Old_Correct';'New_Correct'};
oat.source_recon.gridstep = 8; % in mm
oat.source_recon.time_range = [-0.1 3.4]; % time range in secs
oat.source_recon.freq_range = [0.1 40]; % frequency range in Hz (alpha large range)
oat.source_recon.type = 'Scalar';
oat.source_recon.method = 'beamform';
oat.source_recon.normalise_method = 'mean_eig';
oat.source_recon.forward_meg = 'MEG Local Spheres';
oat.source_recon.report.do_source_variance_maps = 1;
oat.source_recon.dirname = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_01_40Hz_AllTrials'; % spm_files_recog_basen{ii}(12:25)];

%first level (each experimental block for each subject, independently)
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
oat.first_level.cope_type = 'coape'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
oat.first_level.name = ['wholebrain_first_level']; %REMEMBER TO CHECK THIS NAME!!
oat.first_level.bc = [0 0 0];
%to add if the oat has not been automatically saved
for ii = 1:71
    oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level.mat']};
end
oat.first_level.sessions_to_do = [1:38 40:71]; %this is probably not necessary but it may be useful to report if you run only subject_level analysis in oat..
%running first level on parallel computing
for ii = 1:71
    if ii ~= 39 %number 39 is simply a non-existent subject, however for computational reasons we decided to not update the ID of all other subjects and we preferred to simply leave this subject empty and ignore it in the computations
        oat.first_level.sessions_to_do = [];
        oat.first_level.sessions_to_do = [ii];
        jobid = job2cluster(@cluster_beamfirstlevel,oat);
    end
end

%subject level (it simply puts together different experimental blocks for the same subject)
%SUBJECT LEVEL
%in this case it does not do anything since I have only one experimental
%block for each participant (however, for computational reasons, I have to
%run it)
% this is if you have a perfect correspondance between sessions and subjects
oat.subject_level.session_index_list = cell(1,71); %sarebbe length(subjects)
oat.subject_level.name = 'minor';
oat.subject_level.subjects_to_do = [];
oat.subject_level.subjects_to_do = [1:38 40:71];
%to update the names of the results files of subject level
for ii = 1:71
    if ii ~= 39 %number 39 is simply a non-existent subject, however for computational reasons we decided to not update the ID of all other subjects and we preferred to simply leave this subject empty and ignore it in the computations
        oat.subject_level.session_index_list{ii} = ii;
        oat.subject_level.results_fnames(ii) = {['subject' num2str(ii) '_wholebrain_first_level_minor.mat']};
    end
end

%GROUP LEVEL
oat.group_level = [];
oat.group_level.name = 'group_level_everybody_AAL';
oat.group_level.subjects_to_do = [];
oat.group_level.subjects_to_do = 1:70;
oat.group_level.results_fnames = ['wholebrain_first_level_' oat.subject_level.name '_' oat.group_level.name '.mat']; %results name
% Spatial and temporal averaging options
oat.group_level.time_range = [-0.1 3.4];
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
oat.group_level.first_level_contrasts_to_do = [3,1,2]; % list of first level contrasts to run the group analysis on
oat.group_level.report.first_level_cons_to_do = [3]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.show_lower_level_copes = 0;
oat.group_level.report.show_lower_level_cope_maps = 0;

oat = osl_check_oat(oat);
oat.to_do = [1 1 1 1];
oat = osl_run_oat(oat);

%The following lines allow you to create nifti files of the group level
%analysis that can be visualized in programs such as fslview.
S2 = [];
S2.oat = oat;
S2.stats_fname = oat.group_level.results_fnames; %group level
S2.first_level_contrasts = [3,1,2]; %remember that in this way you define the order of the output (this numbers refers to the order (numbers) defined in the contrasts; for example here tstat3 refers to the contrast old-new)
S2.group_level_contrasts = [1];
S2.resamp_gridstep = oat.source_recon.gridstep;
[statsdir,times,count] = oat_save_nii_stats(S2);

%reversing sign of copes (contrasting the two conditions..)
%This is for reversing the sign (multiplying by -1) of
%lower_level_copes since the permutation function used later does not 
%work with negative values.. therefore here I work around this problem by
%calculating permutations for the reversed copes. Then I combine the two
%masks outputted for obtaining a final image with both positive and
%negative statistics (t-values) that formed the significant clusters..
%this is done only for the cases were cond2 had clusters with stronger
%activity than cond1
reverse = 1;
path = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_2_8Hz_AllTrials.oat'; %path to files
name = 'wholebrain_first_level_minor_group_level_everybody_AAL'; %name of specific file with statistics
pathful = [path '/' name '.mat']; %full path
if reverse == 1
    load(pathful); %loading data
    load([path '/oat_' name '.mat']) %loading oat mask with information on the analysis
    oat.group_level.results_fnames = [name '_reverse']; %creating new name for the new data that is being reversed
    ghi = oat_stage_results.lower_level_copes{3}; %extracting data of interest (here contrast number 3, the one with cond1 against cond2)
    ghi2 = ghi * (-1); %reversing the sign of contrast 3
    oat_stage_results.lower_level_copes{3} = ghi2; %storing data
    save([path '/' name '_reverse.mat'],'oat_stage_results','-v7.3'); %saving new oat mask with information on the analysis
    save([path '/oat_' name '_reverse.mat'],'oat'); %save new reversed data
end

%cluster-based permutation test
%The following lines implement a cluster-based permutation test that
%highlights the significant clusters of brain activity. This test has been
%run at first for the significant cluster emerged from MCS on MEG sensor
%data. Then, we run this test for each of the 5 tones forming the musical
%sequences (memorized sequecens and novel sequences). Details about the procedure (as
%well as correction for multiple comparisons) are provided in the paper.
S = [];
S.oat = oat;
S.cluster_stats_thresh = 1.7;
S.cluster_stats_nperms = 5000; % we normally recommend doing 5000 perms
S.first_level_copes_to_do = [3,2,1];
S.group_level_copes_to_do = [1];
S.group_varcope_spatial_smooth_fwhm = S.oat.group_level.group_varcope_spatial_smooth_fwhm;
S.write_cluster_script = 0;
%here the desired time-window has been selected
S.time_range = [0.547 1.180]; %time MEG sensor cluster I
% S.time_range = [0 0.250];
% S.time_range = [0.251 0.500];
% S.time_range = [0.501 0.750];
% S.time_range = [0.751 1.000];
% S.time_range = [1.001 1.250];
S.time_average = 1;
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
%Additional clarifications about these algorithms can be found in the paper.
%Furthermore, you are very welcome to contact us if you have specific
%questions.

%% creating a file with information on significant clusters/voxels outputted by the permutation test

%%% THIS SECTION REQUIRES OSL %%%

%YOU CAN USE THE FOLLOWING 3 LINES TO REPRODUCE SOURCE RESULTS THAT WE
%REPORTED IN THE PAPER (IN SUPPLEMENTARY MATERIALS)
tonebytone = 0; %number of tone to be considered; 0 if you want the time-windows related to the MEG sensor significant cluster
contrn = 1; %set the contrast that you want (1 = memorized musical patterns vs baseline; 2 = novel musical patterns vs baseline; 3 = memorized vs novel musical patterns)


%path to file nifti
% pathnii = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_sourcedata_BrainActivity';
pathnii = [pathdata '/MEG_sourcedata_BrainActivity'];
if tonebytone == 0
    hgf = 'SourceClusters_CorrespondentToMEGSensorClusters';
    klj = '0547_1180';
else
    hgf = 'MusicalToneByTone';
    klj = ['note' num2str(tonebytone)];
end
%actual name plus path
fname = [pathnii '/' hgf '/c' num2str(contrn) '_' klj '.nii.gz'];
%getting MNI coordinates of significant voxels within the provided image
[ mni_coords, xform ] = osl_mnimask2mnicoords(fname);
%loading the image
V = nii.load(fname);
%extracting statistics
VV = V(V~=0);
%indices of non-zero values of nifti image
VI = find(V~=0);
%path to AAL template
parcelfile = [pathl '/External/aal_8mm_try5.nii.gz']; %load this from the provided codes folder
%loading AAL labels
% load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat'); %load this from the provided codes folder
load([pathl '/External/AAL_labels.mat']); %load this from the provided codes folder
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
% writetable(PDn,['c' num2str(contrn) '.xlsx'],'Sheet',1) %printing excel file

% volumeInfo=spm_vol(fname)
% [intensityValues,xyzCoordinates ]=spm_read_vols(volumeInfo);

%%

%%% *** FUNCTIONAL CONNECTIVITY *** %%%

%%

%Data has been then analysed in terms of functional connectivity.
%To do that, we beamformed the data in 5 different frequency bands (delta, theta, alpha, beta, gamma), as
%described in the paper. This operation was conducted for both task data
%and resting state data that in this case has been used as baseline.
%Once again we provide the scripts that allowed us to estimate the neural
%sources by beamforming as well as the functions for extracting the data
%that has then been analysed in terms of functional connectivity. This
%extracted data is provided for you to allow you to explore it and
%replicate our results.


%% BEAMFORMING (same algorithm as before.. simply implemented for both task and resting state and 5 different frequency bands.. details are provided in the paper)

%starting up OSL and adding paths
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup
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
%actual beamforming
count1 = 0;
oat = [];
for ii = 3:3:213  % = 1:length(spm_files_recog_basen)
    count1 = count1 + 1;
    spm_files = ['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/e80ldff'  spm_files_recog_basen{ii}]; %check this in order not to mix up 'epoched' and 'continuous'
    D_epoched = spm_eeg_load(spm_files);
    D_epoched = D_epoched.montage('switch',1); %switch the montage to 1 in order to be safer
    D_epoched.save();
    processed_file_epoched = spm_files;
    spm_files = ['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/dff'  spm_files_recog_basen{ii}]; %check this in order not to mix up 'epoched' and 'continuous'
    D_continuous = spm_eeg_load(spm_files);
    D_continuous = D_continuous.montage('switch',1); %switch the montage to 1 in order to be safe
    D_continuous.save();
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
oat.source_recon.sessions_to_do = [1:38 40:71]; %sessions to do among the file_list; number 39 is simply a non-existent subject, however for computational reasons we decided to not update the ID of all other subjects and we preferred to simply leave this subject empty and ignore it in the computations
oat.source_recon.modalities = {'MEGMAG'; 'MEGPLANAR'}; %'MEGMAG'};
oat.source_recon.conditions = {'Old_Correct';'New_Correct'};
oat.source_recon.gridstep = 8; % in mm
oat.source_recon.time_range = [-2 6]; % time range in secs
%setting for different frequency bands
oat.source_recon.freq_range = [0.1 2]; % frequency range in Hz (delta)
% oat.source_recon.freq_range = [2 8]; % frequency range in Hz (theta)
% oat.source_recon.freq_range = [8 12]; % frequency range in Hz (alpha)
% oat.source_recon.freq_range = [12 32]; % frequency range in Hz (beta)
% oat.source_recon.freq_range = [32 74]; % frequency range in Hz (gamma)
oat.source_recon.type = 'Scalar';
oat.source_recon.method = 'beamform';
oat.source_recon.normalise_method = 'mean_eig';
oat.source_recon.forward_meg = 'MEG Local Spheres';
oat.source_recon.report.do_source_variance_maps = 1;
%example of output path for delta band
oat.source_recon.dirname = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_01_2Hz_AllTrials_longbaseline'; % spm_files_recog_basen{ii}(12:25)];

oat = osl_check_oat(oat);
oat.to_do = [1 0 0 0];
oat = osl_run_oat(oat);

%% averaging over conditions (only for task data)

recon_freq = 'source_localsphere_2_8Hz_AllTrials_longbaseline.oat'; %this is an examle for 1 frequency band only, but the operation has been conducted for the 5 frequency bands described in the paper
%averaging across conditions
for ii = 1:71
    if ii ~= 39 %number 39 is simply a non-existent subject, however for computational reasons we decided to not update the ID of all other subjects and we preferred to simply leave this subject empty and ignore it in the computations
        filepath = ['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/' recon_freq '/concatMefsession' num2str(ii) '_spm_meeg.mat'];
        S = [];
        S.D = filepath;
        S.prefix = 'mfs';
        D = spm_eeg_average(S);
    end
end

%% reducing dimensionality from voxel level to ROI level (by using AAL parcellation)

%this was done for both task and resting state data

%this is for giving input such as parcellation, ROI names and template (otherwise it can also find by itself the proper template, in this case 8mm)
parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'; %this and the following files can be found in the provided codes folder
ROIsfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_ROIs.txt';
templatefile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz';
p = parcellation(parcelfile,ROIsfile,templatefile); %given required inputs, creating a parcellation object for SPM MEEG objects
p = p.remove_parcels(91:116);
%folder
recon_freq = 'source_localsphere_2_8Hz_AllTrials_longbaseline.oat'; %example for theta band task data (even if this operation has been done for all frequency bands and both task and resting state)
mont = 2; %specify the SPM object montage you are working on
for ii = 1:71
    if ii ~= 39 %number 39 is simply a non-existent subject, however for computational reasons we decided to not update the ID of all other subjects and we preferred to simply leave this subject empty and ignore it in the computations
        D = spm_eeg_load(['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/' recon_freq '/mconcatMefsession' num2str(ii) '_spm_meeg.mat']);
        D = D.montage('switch',mont);
        D = ROInets.get_node_tcs(D,p.parcelflag(true),'PCA'); %using PCA
        D.save();
    end
end

%% resting state

%Here we simulated trials in the resting state data in order to make it
%comparable with the task data and use it as a proper baseline (details are
%provided in the paper)

%extracting resting state data and saving it to single files for later operations
%2-8Hz
%resting state
list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz2_8_rest10.oat/concat*meeg.mat'); %load dir file "concat" in lista
Drest = extracting_data_LBPD_D(list);
save -v7.3 restingstate.mat Drest
%0.1-2Hz
%resting state
list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz01_2_rest10.oat/concat*meeg.mat'); %load dir file "concat" in lista
Drest = extracting_data_LBPD_D(list);
save('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz01_2_rest10.oat/restingstate.mat','Drest','-v7.3')
%8-12Hz
%resting state
list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz8_12_rest10.oat/concat*meeg.mat'); %load dir file "concat" in lista
Drest = extracting_data_LBPD_D(list);
save('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz8_12_rest10.oat/restingstate.mat','Drest','-v7.3')
%12-32Hz
%resting state
list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz12_32_rest10.oat/concat*meeg.mat'); %load dir file "concat" in lista
Drest = extracting_data_LBPD_D(list);
save('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz12_32_rest10.oat/restingstate.mat','Drest','-v7.3')
%32-74Hz
%resting state
list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz32_74_rest10.oat/concat*meeg.mat'); %load dir file "concat" in lista
Drest = extracting_data_LBPD_D(list);
save('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz32_74_rest10.oat/restingstate.mat','Drest','-v7.3')
%0.1-40Hz
%resting state
list = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz01_40_rest10.oat/concat*meeg.mat'); %load dir file "concat" in lista
Drest = extracting_data_LBPD_D(list);
save('/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz01_40_rest10.oat/restingstate.mat','Drest','-v7.3')

%% preparing data - baseline from resting state

% x = 190; %trial length (for static FC)
x = 291; %trial length (for dynamic FC)
y = 80; %number of trials
% pathd = '/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz01_2_rest10.oat/restingstate.mat';
pathd = '/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz2_8_rest10.oat/restingstate.mat';
% pathd = '/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz8_12_rest10.oat/restingstate.mat';
% pathd = '/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz12_32_rest10.oat/restingstate.mat';
% pathd = '/scratch5/MINDLAB2017_MEG-LearningBach/Francesco/source_localsphere_Hz32_74_rest10.oat/restingstate.mat';

outt = preparing_baseline_from_restingstate_LBPD_D(x,y,pathd);

%% preparing data - task

%preparing data
S = [];
subjnum = {'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'38';'39';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';'51';'52';'53';'54';'55';'56';'57';'58';'59';'60';'61';'62';'63';'64';'65';'66';'67';'68';'69';'70';'71'}; %Everybody!!!!!!
datadir1 = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source';
datadir2 = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/after_maxfilter_mc/';
% recon_freq = 'source_localsphere_01_2Hz_AllTrials_longbaseline.oat';
recon_freq = 'source_localsphere_2_8Hz_RTsTrials_longbaseline.oat';
% recon_freq = '8_12Hz_longbaseline_testlocal.oat';
% recon_freq = 'source_localsphere_12_32Hz_AllTrials_longbaseline.oat';
% recon_freq = 'source_localsphere_Hz32_Hz74_20TrialsRT_gamma.oat';
S = [];
S.conditions = cell(1);
S.time_range = cell(1);
%this is just because the gamma data is stored in a different directory and has different names due to the 'history' of this analysis..
if strcmp(recon_freq,'source_localsphere_Hz32_Hz74_20TrialsRT_gamma.oat')
    datadir = datadir2;
    concas = 'mconcatMefsession';
    S.time_range(1) = {{[17 207]}};
    S.conditions(1) = {{'Old_Correct','New_Correct'}};
else
    datadir = datadir1;
    concas = 'mfsconcatMefsession';
%     S.time_range(1) = {{[301 490]}}; %this is for static FC
    S.time_range(1) = {{[276 567]}}; %this is for DFC (here we selected a
%     slightly larger time-window since the instantaneous phase estimation
%     gave rise to boundary artifacts. therefore we computed it on a larger
%     timeseries (lasting approximately 2000ms) and then we did not consider
%     the edges of it. Currently, there is no consensus about the ideal amount of
%     time-points to be removed, however testing different solutions
%     provided similar results).
    S.conditions(1) = {{'Old_Correct_Fast','New_Correct_Fast'}};
end
%preparing data
S.average_label = 1; %set to 1 for average across trials; set to 0 for single trial analysis; set to -1 if you want the single trial and the specific corresponding RT (RT for each trial and each subject); set -2 for averaging across small groups of trials (some sort of stabilization of similar trials)
for ii = 1:length(subjnum)
    if strcmp(recon_freq,'source_localsphere_Hz32_Hz74_20TrialsRT_gamma.oat')
        subdummy = 3*str2double(subjnum{ii});
    else
        subdummy = ii;
    end
    if str2double(subjnum{ii}) ~= 39 %number 39 is simply a non-existent subject, however for computational reasons we decided to not update the ID of all other subjects and we preferred to simply leave this subject empty and ignore it in the computations
        S.subjects_list{1,ii} = [datadir '/' recon_freq '/' concas num2str(subdummy) '_spm_meeg.mat'];
    end
end
S.ROIs = 90; %number of ROIs (AAL = 90, 38-ROI OSL = 38..)
S.subjnum = [1:38 40:71];
S.fsample = 150;
S.montage_numb(1,1) = 3;

data = ps_preparedata_spmobj_LBPD_D(S);

%% combining task and baseline (with extra care not to mix up subjects ID..)

if strcmp(recon_freq,'source_localsphere_2_8Hz_RTsTrials_longbaseline.oat')
    lab01 = 0; %1 if subject IDs are with 01 and not 1.. due to some kind of unexpected process in the naming system of a previous function in the resting state.. anyway it is not a big deal..
else
    lab01 = 1;
end    
if lab01 ~= 1
    %writing subject IDs for task data
    IDtask = cell(length(subjnum),1);
    for ii = 1:length(subjnum)
        IDtask(ii,1) = {['concatMefsession' num2str(ii) '_spm_meeg.mat']};
    end
else
    %writing subject IDs for task data
    IDtask = cell(length(subjnum),1);
    for ii = 1:length(subjnum)
        if ii < 10 %this is jsut because we want to have a structure with all subjects and the resting state data is only for 69 out of 71 subjects.. it is not a big deal either.. just a matter of reshaping matrices, etc..
            if strcmp(outt{1,2}(17),'0')
                IDtask(ii,1) = {['concatMefsession0' num2str(ii) '_spm_meeg.mat']};
            else
                IDtask(ii,1) = {['concatMefsession' num2str(ii) '_spm_meeg.mat']};
            end
        else
            IDtask(ii,1) = {['concatMefsession' num2str(ii) '_spm_meeg.mat']}; 
        end
    end
end
%storing real data task in data2
data2 = cell(length(subjnum),3);
data2(:,1) = data.prep_data(:,1);
data2(:,2) = data.prep_data(:,2);
%reordering baseline data (if existent)
dDd = cell(length(subjnum),1);
for ii = 1:length(subjnum)
    if ~isempty(data2{ii,1})  
        a = 0;
        count = 0;
        while a == 0 && count < length(outt(:,2)) %looking for correct subject ID.. until it is found.. if it is not found after looking in all subjects, you later assign dDd = []
            count = count + 1;
            a = double(strcmp(IDtask{ii,1},outt{count,2}));
        end
        if a == 1
            dDd(ii,1) = outt(count,1);
        else
            dDd(ii,1) = {[]};
            data2(ii,1) = {[]}; %doing that also in data2 since we need to completely discard the subject without baseline..
            data2(ii,2) = {[]};
        end
    end
end
%actual data
data2(:,3) = dDd;
%storing labels
L(1) = data.label(1);
L(2) = data.label(2);
L(3) = {'dataset_2_RSBaseline'};
%storing data for subsequent function
DATA = data;
DATA.prep_data = data2;
DATA.label = L;

%% combining memorized musical patterns (old) and novel musical patterns (new) conditions by averaging them together (this is used only in static Functional Connectivity (SFC))

DATAori = DATA;
prep_data2 = cell(size(DATA.prep_data,1),2);
for ii = 1:71
    hj1 = DATA.prep_data{ii,1}; %old
    hj2 = DATA.prep_data{ii,2}; %new
    pol = zeros(size(hj1,1),size(hj1,2));
    for jj = 1:size(hj1,1)
        for pp = 1:size(hj1,2)
            pol(jj,pp) = (hj1(jj,pp) + hj2(jj,pp))/2; %average for each ROI and each time-point
        end
    end
    prep_data2(ii,1) = {pol}; %storing average old and new
    prep_data2(ii,2) = DATA.prep_data(ii,3); %storing baseline
    disp(ii)
end

DATA.prep_data = prep_data2;

%% saving data

%The data computed so far has been saved on disk with the following lines and shared since it is
%completely anonymised and to our knowledge not linkable to any personal
%information of the participants. This data will be loaded and used in the
%following sections for computing SFC, DFC and the correspondent
%statistics.
save delta_SFC.mat DATA
save theta_SFC_Bachsoriginalandvariation.mat DATA
save theta_SFC.mat DATA
save alpha_SFC.mat DATA
save beta_SFC.mat DATA
save gamma_SFC.mat DATA

%% static FC - Pearson's correlations

%Here you can reproduce the results concerning static FC (memorized and novel musical patterns
%combined and contrasted vs baseline) for the 5 different frequency bands considered in the study.
%Results and figure are reported in the paper and in supplementary materials.

%specify the requested frequency band
fb = 2; %1 = delta; 2 = theta; 3 = alpha; 4 = beta; 5 = gamma

%loading requested data
bandf = {'delta';'theta';'alpha';'beta';'gamma'};
load([pathdata '/MEG_FC_Static/' bandf{fb} '_SFC.mat']);
%actual SFC computation
S = [];
S.contrast = [2 1]; %for detailed information about these input parameters have a look at the documentation of this function (the same concept applies to the other functions..)
S.resam_number = [];
S.fsample = 150;
S.plot_conditions = [];
S.plot_contr = 1;
S.perm_degree = 1000;
S.perm_max = 0;
S.p_perm = 0;
S.p_threshMC = 2.0e-04;
S.clims_contr = [-15 15];
S.clims_cond = [];
S.parametric = 1;
S.schemaball_contrast_label = 0;
S.schemaball_cond_label = 0;
S.schemplth = 100;
S.colbacksch = [];
S.lab_parc_ph = [pathl '/External/AAL_labels.mat']; %this can be found in the provided folder
% S.outpath = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/TestingPhDFunctions_Scratch';
S.outpath = [pathdata 'testing'];
S.name_sch = [];

[H, P_binary, P_ROIs, matf, cor_mat ] = static_FC_MEG_LBPD_D(S, DATA);

%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))

%% preparing data for calculating ANOVA between different frequency bands

j = 0.2944; %maximum value in averaged matrix calculated above for 2-8Hz
dop = zeros(size(cor_mat,4),1);
for ii = 1:size(cor_mat,4)
    bob = mean(mean(cor_mat(:,:,1,ii)))/j; %storing average old and new static FC
    obo = mean(mean(cor_mat(:,:,2,ii)))/j; %storing baseline static FC
    dop(ii,1) = bob - obo; %subtracting baseline from task
    disp(ii)
end
%loading that.. previously calculated matrix with data calculated as: meanSFCtask - meanSFCbaseline (mean over conditions and ROIs)
% load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper3/Final/StaticFCaveragedtaskminusaveragedbaseline_bandsfrequencyinprogression.mat');
% load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_FC_Static/ANOVA_SFC.mat');
load([pathdata '/MEG_FC_Static/ANOVA_SFC.mat']);
[p,t,stats] = anova1(fd)
[c,m,h,nms] = multcompare(stats)
%plotting data as violins..
figure
violin(fd,'xlabel',{'0.1-2Hz','2-8 Hz','8-12Hz','12-32Hz','32-74Hz'},'mc',[],'medc',[])
grid minor
%writing csv file to be exported for violin plot in R software
% csvwrite('violin_R.csv',fd);

%The violin plot reported in the paper has been further edited in R software but
%you can see that reprents the same results depicted in the figure produed in this section

%% plotting average activity over all MEG gradiometers channels (supplementary materials)

outdir = [pathdata '/MEG_sensordata_UnivariateTests_MonteCarloSimulations']; %this must become your own directory to the data that we provided
load([outdir '/sensor_data.mat']);
%computing data
S = [];
S.outdir = outdir;
S.data = [];
S.data = data_mat;
S.chanlabels = chanlabels;
S.time_real = time_sel;
S.conditions = {'Old_Correct','New_Correct'};
S.timeextract = []; %time-points to be extracted
S.centerdata0 = 0; %1 to make data starting at 0
S.save_data = 0; %only meaningfull if you read data from SPM objects saved on disk
S.save_name_data = 'sensor_data';
%individual waveform plotting (settings for individual channel plots are not used in the paper..)
S.waveform_singlechannels_label = 0; %1 to plot single channel waveforms
S.wave_plot_conditions_together = 1; %1 for plotting the average of all
S.mag_lab = 2; %1 for magnetometers; 2 for gradiometers
S.x_lim_temp_wave = []; %limits for time (in secs) (E.g. [-0.1 3.4])
S.y_lim_ampl_wave = []; %limit for amplitude (E.g. [0 120] magnetometes, [0 6] gradiometers)
%averaged waveform plotting
S.waveform_average_label = 1;
S.legc = 0; %set 1 for legend
S.left_mag = [2:2:204];
S.signtp(1) = {[]}; %in time samples
S.sr = 150; %sampling rate (Hz)
S.avewave_contrast = 0; %1 to plot the contrast between conditions (averaged waveform)
S.save_label_waveaverage = 0;
S.label_plot = 'block_minor';
%t-tests
S.t_test_for_permutations = 0;
S.cond_ttests_tobeplotted_topoplot = [1 2]; %this is for both topoplot and t-tests!! (here [1 2] means cond1 vs cond2!!!!!!!)
%topoplotting
S.topoplot_label = 0;
% S.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External';
S.fieldtrip_mask = [pathl '/External'];
S.topocontr = 0;
S.topocondsing = [];
S.xlim = []; %time topolot
S.zlimmag = []; %magnetometers amplitude topoplot limits
S.zlimgrad = []; %gradiometers amplitude topoplot limits
S.colormap_spec = 0;
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
S.colormap_spec_x = x;
S.topoplot_save_label = 0;

[out] = MEG_sensors_plotting_ttest_LBPD_D2(S);

%% time-frequency plotting (supplementary materials)

freqq = 2; %1 for 1-10Hz; 0 for 1-74Hz

% load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_sensordata_UnivariateTests_MonteCarloSimulations/sensor_data.mat');
load([pathdata '/MEG_sensordata_UnivariateTests_MonteCarloSimulations/sensor_data.mat']);
mag = 1:204;
if freqq == 1
    f = 0.1:0.1:10;
else
    f = 1:1:74; %frequencies..
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
%plotting old
figure
imagesc(time_sel(16:end),f,Pold)
set(gca,'YDir','normal') %plotting in descending order
xlabel('time (s)'); ylabel('f (Hz)');
colorbar
if ges ~= 0
    caxis([0 ges])
end
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))


%% Static FC - theta band - Bach's original and Bach variation

%%% TO REPRODUCE THE CONNECTIVITY WITHIN THE BRAIN YOU NEED OSL %%%

%Here we report codes for reproducing the figure depicting the static FC for
%theta band computed independently for memorized and  novel musical patterns.
%These figures are reported in the supplementary materials of the paper and include matrix,
%schemaball and brain layout representation of brain SFC.
cond = 1; %set 1 for memorized musical patterns; set 2 for novel musical patterns

close all
% load(['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_FC_Static/theta_SFC_Bachsoriginalandvariation.mat']);
load([pathdata '/MEG_FC_Static/theta_SFC_Bachsoriginalandvariation.mat']);
poi = {[2 0 1]; [0 2 1]};
pl = [1 3];
plo = {'Bachs original';'Bach variation'};
aosl = dir([pathl '/External/*osl']);
if ~isempty(aosl)
    %parcellation AAL
    %this is for giving input such as parcellation, ROI names and template (otherwise it can also find by itself the proper template, in this case 8mm)
    % parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'; %this and the following files can be found in the provided codes folder
    parcelfile = [pathl '/External/aal_8mm_try5.nii.gz']; %this and the following files can be found in the provided codes folder
    ROIsfile = [pathl '/External/aal_ROIs.txt'];
    templatefile = [pathl '/External/MNI152_T1_8mm_brain.nii.gz'];
    p = parcellation(parcelfile,ROIsfile,templatefile);
    p = p.remove_parcels(91:116); %removing the cerebellum
end
for ii = cond %this loop is of course not really necessary anymore but it does not affect the results
    %actual SFC computation
    S = [];
    S.contrast = poi{ii}; %for detailed information about these input parameters have a look at the documentation of this function (the same concept applies to the other functions..)
    S.resam_number = [];
    S.fsample = 150;
    S.plot_conditions = [0 0 0];
    S.plot_contr = 1;
    S.perm_degree = 100;
    S.perm_max = 0;
    S.p_perm = 0;
    S.p_threshMC = 2.0e-04;
    S.clims_contr = [-12.3 12.3];
    S.clims_cond = [];
    S.parametric = 1;
    S.schemaball_contrast_label = 1;
    S.schemaball_cond_label = 0;
    S.schemplth = 100;
    S.colbacksch = [];
    S.lab_parc_ph = [pathl '/External/AAL_labels.mat']; %this can be found in the provided folder
    S.outpath = pathdata; %specify your own output path
    S.name_sch = [];
    
    [H, P_binary, P_ROIs, matf, cor_mat ] = static_FC_MEG_LBPD_D(S, DATA);
    
    figure(ii) %simply adjusting colour of matrices..
    title(plo(ii))
    %colormap with white for 0 values
    x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
    colormap(bluewhitered_PD(0,x))
    
    %schemaball and connectivity in the brain (static FC)
    %this function saves on disk the plots of connectivity in the brain
    %layout.
    %The loop is built to save at first Bach's original plot and then
    %Bach variation one.
    %input parameters
    S = [];
    S.MATF = matf{4,1};
    S.symmetric_l = 1;
    S.outpath = pathdata; %specify your own output path
    if ~isempty(aosl) %if you have OSL
        S.p = p;
    else
        S.p = [];
    end
    S.thr_cortex = [0.98];
    if ~isempty(aosl) %if you have OSL
        S.name_gif = [plo{ii}];
    else
        S.name_gif = [];
    end
    S.frame_vect = [15,90,130];
    S.fr_spec = {'Posterior_L_R';'Frontal_R_L';'Hem_L'};
    %the following input parameters are not used but merely provided..
    S.schball_l = 0; %they are not used because of this label = 0
    S.extr = [];
    S.lab_parc_ph = [pathl '/External/AAL_labels.mat']; %this can be found in the provided folder
    S.colbacksch = [];
    S.schemplth = 100;
    S.name_sch = ['Schemaball_oldvsnew_note' num2str(ii)];
    %actual function
    FC_Brain_Schemaball_plotting_LBPD_D(S)
end

%% phase synchrony

%Here we computed dynamic FC. DATA has been extracted in the same way as
%done for the static FC (by using the sections of this script reported
%above). Then, through the function of this section
%("phasesynchrony_LBPD_D") we estimated the instantaneous phase and saved
%the outputted results that is provided to you. That will then be used for
%statistics and plotting purposes in the following sections.

S = [];
S.data = DATA;
S.save_path = [pathdata '/MEG_FC_Dynamic2'];
S.ROIs = 90;
S.left_bound = 26; %starting after 26 time-points
S.right_bound = 77; %removing last 77 time-points (so considering only the time while the notes were played , + a bit of more time (around 250ms) to detect the processing of the last tone)
%here we computed the analysis of a larger time-window than the duration of
%the musical sounds to be able to discard the edges and prevent potential
%boundary artifacts of instantaneous phase estimation

phasesynchrony_LBPD_D(S);

%The output of this section is saved in the provided data in the folder:
%"yourpath/MEG_FC_Dynamic"

%% contrasts

%Here we calculated statistics (that then will be corrected for multiple
%comparisons and tested for specific parameters such as degree centrality
%in the following sections) for every participant.
%In the current section we calculated statistics for IFC matrices in two
%cases: by subaveraging them in order to have 1 matrix for each musical
%tone and by keeping 1 matrix for each time-point. Details are provided in the paper.
%THe results of this procedure are reported in:
%[pathdata '/MEG_FC_Dynamic/Time_windows_5bis']

%contrbase = 1; %1 for contrasting memorized and novel musical patterns vs baseline; 2 for contrasting memorized vs novel musical patterns
%subavel = 1; %1 for subaveraging (1 IFC matrix for each tone); 0 for non-subaveraging (1 IFC matrix for each time-point)
S = [];
S.subjs_list = {[1,4:7,9:38,40:71]};
S.load_path = [pathdata '/MEG_FC_Dynamic'];

if subavel == 1
    S.save_path = ['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_FC_Dynamic/Time_windows_5']; %1 matrix for each tone forming the musical sequences (Bach's original and Bach variation)
    S.Num_window = [5]; %leave it empty [] for not doing any subaveraging
else
    S.save_path = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_FC_Dynamic/nosub'; %1 matrix for each time-point
    S.Num_window = []; %leave it empty [] for not doing any subaveraging
end
S.group_comparison = [1]; %[1] for one group; [1 2] for two groups (referred to the subjs_list)
if contrbase == 1
    S.condition_contrast = [2 0 1; 0 2 1];
    S.cond_contr_meanbaseline = 1;
else
    S.condition_contrast = [2 1 0];
    S.cond_contr_meanbaseline = 0;
end
%actual function
ps_statistics_LBPD_D(S);

%% plotting IFC matrices

%Here you can plot IFC matrices. In this provided section you can run the
%IFC matrices corresponding of each tone of memorized and novel musical patterns.
%These matrices are reported in the supplementary materials of the paper.

condit = 2; %set 1 for Bach's original; set 2 for Bach variation

phase_dir = [pathdata '/MEG_FC_Dynamic/Time_windows_5']; %general path
load([phase_dir '/Contr_' num2str(condit) '_3_STAT_zval_fullcontr_group_comp_1_1.mat']);
a = STAT_zval_fullcontr;
S = [];
S.STAT = a;
S.ordresh = 1;
S.cond_contr_meanbaseline = 1;
S.Num_window = 5;
S.imag_contr_lim = [-4 4];
S.indimagesnumber = [];

IFC_plotting_LBPD_D(S)
%changing colormap of figures (red-bue with white for 0 values)
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))

%% degree MCS over the different musical notes

%Here we computed the degree centrality of ROIs within the whole brain
%networks and test the statistical significance by MCS.
%Details about this methodology are described in the paper.
%In the supplementary materials we illustrated a graphical representation of this
%MCS algorithm.
%Please note that since MCS algorithms are based on the randomization of the data,
%they provide slightly different results each time they are computed. These differences
%are very small and they do NOT change the meaning of the results. However,
%it is important to specify it and this is also the reason why you will
%probably end up noticing small differences between the results reported in
%the paper and the ones that you can get when using this script.
%Additional information about the relationship between randomization,
%thresholds and small variability of the results are provided in the
%documentation of the function "degree_segregation_MCS_LBPD_D" that is used
%in this section of the script.
%If you want you can of course reproduce these results, otherwise you can simply load them from
%the folders that we provided and used them in following sections of this script.

%please change the following two lines to get the results reported in the paper
contr = '1_2'; %1_2 = memorized musical patterns (cond1) vs novel musical patterns (cond2); 1_3 = memorized musical patterns (cond1) vs baseline (cond2); 2_3 = novel musical patterns (cond1) vs baseline (cond2)
posl = 1; %1 for cond1 > cond2; otherwise for cond2 > cond1
subavel = 1; %1 for subaveraging (degree calculation for each tone of memorized and novel musical patterns); 0 for non-subaveraging (degree calculation for each time-point of memorized and novel musical patterns)


%actual computation
clear phase_dir1
if subavel == 1
    phase_dir1(1) = {[pathdata '/MEG_FC_Dynamic/Time_windows_5']}; %general path everybody
else
    phase_dir1(1) = {'/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_FC_Dynamic/nosub'};
end

%settings
matrn = 0; %matrix number; 0 for every matrix
thr = 0; %0 for non-binarising the t-values matrix; otherwise threshold for binarising it
if strcmp(contr,'1_2')
    contr2 = 'oldvsnew';
elseif strcmp(contr,'1_3')
    contr2 = 'oldvsbaseline';
else
    contr2 = 'newvsbaseline';
end
thresh = 0;
permut = 1000;
threshMC = 1.0e-04;
ROI_label = [pathl '/External/AAL_labels.mat'];
reshlab = 0;
perm_max = 0;
for opo = 1%:length((phase_dir1))
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
    if matrn == 0 %calculation on every time-point (or tone)
        x = 1:size(STAT_p_fullcontr,3);
        ROII = cell(x(end),1);
    else %or on specific time-point(s) (or tone(s))
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
%         writetable(cell2table(P_ROIs_deg),[phase_dir '/Contr_' contr '_pos_' num2str(posl) '_note' num2str(ii) '.xlsx'],'Sheet',1)
    end
    save([phase_dir '/' contr '_' num2str(posl) '_degree.mat'],'PP','ROII','PBIN','HSEG');
    mop = num2str(threshMC);
end

%% CREATING NIFTI FILES TO BE THEN EXPORTED FROM THE SERVER AND USED IN WORKBENCH

%memorized musical patterns and novel musical patterns

%This section created nifti files that were then used in Workbench for
%making the brain plots reported in the paper.
%Here you can reproduce those files, however, if you want to see them in
%the brain layout, you need to further elaborate them in Workbench.

contr = '1_3'; %1_3 for Bach's original vs baseline; 2_3 for Bach variation vs baseline

%set your own path in the following lines where you stored the data that you computed or that we provided
load(['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_FC_Dynamic/Time_windows_5/Contr_' contr '_STAT_zval_fullcontr_group_comp_1_1.mat'])
ST = squeeze(mean(STAT_zval_fullcontr,2)); %getting mean of connections bewteen one ROI and the other ones
ST2 = zeros(90,length(x));
load(['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_FC_Dynamic/Time_windows_5/' contr '_1_degree.mat'])
ST2(PBIN==1) = ST(PBIN==1); %getting mean zvalues for significant ROIs
%mask with significant ROIs and then taking average statistics for those ROIs
for ii = 1:length(x)
    % vector = cell2mat(IFC_new_alpha.P_pos_int_sum(:,2)); %specify vector in AAL space
    vector = ST2(:,ii)';
    namenii = ['Degree_WB_' contr '_note' num2str(ii) '.nii']; %name for the figure (remember to make it ending in '.nii'
    symmetric = 0; %if the AAL vector is submitted in symmetric order
    create_AALnifti(vector,namenii,symmetric);
end

%% CREATING NIFTI FILES TO BE THEN EXPORTED FROM THE SERVER AND USED IN WORKBENCH

%memorized vs novel musical patterns

%This section created nifti files that were then used in Workbench for
%making the brain plots reported in the paper.
%Here you can reproduce those files, however, if you want to see them in
%the brain layout, you need to further elaborate them in Workbench.

%set your own path in the following lines where you stored the data that you computed or that we provided
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat'); %loaded from the provided codes folder
%memorized musical patterns - getting indices of ROIs (non-symmetric, so LRLRLR) that were significantly central over time
%(over the 5 musical tones)
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_FC_Dynamic/Time_windows_5/1_2_1_degree.mat')
ROI2 = zeros(90,1,5);
for ii = 1:5
    RI = ROII{ii};
    for jj = 1:size(RI,1)
        a = RI(jj,1);
        cnt = 1;
        while strcmp(lab(cnt,:),a) ~= 1 
            cnt = cnt + 1;
        end
        ROI2(cnt,1,ii) = 1;
    end
end    
%same for novel musical patterns
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_FC_Dynamic/Time_windows_5/1_2_0_degree.mat')
ROI3 = zeros(90,1,5);
for ii = 1:5
    RI = ROII{ii};
    for jj = 1:size(RI,1)
        a = RI(jj,1);
        cnt = 1;
        while strcmp(lab(cnt,:),a) ~= 1
            cnt = cnt + 1;
        end
        ROI3(cnt,1,ii) = 1;
    end
end
%memorized vs novel musical patterns
%loading IFC matrix of required contrasts
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_FC_Dynamic/Time_windows_5/Contr_1_2_STAT_zval_fullcontr_group_comp_1_1.mat')
ST = squeeze(mean(STAT_zval_fullcontr,2)); %getting mean of connections bewteen one ROI and the other ones over the 5 tones
%significantly central ROIs for memorized and novel musical patterns used as mask for plotting their
%mean statistics within the brain
%here the mean statistics is for contrasting memorized vs novel musical patterns (no baseline)
for ii = 1:5
    v = zeros(90,1);
    v(ROI2(:,1,ii) == 1) = ST(ROI2(:,1,ii) == 1,ii); %here ROI2 correponds to old (and ST to old vs new contrast)
    v(ROI3(:,1,ii) == 1) = ST(ROI3(:,1,ii) == 1,ii); %here ROI3 correponds to new (and ST to old vs new contrast)
    V(:,ii) = v;
    namenii = ['Degree_WB_origvsvar_contraststat_note' num2str(ii) '.nii']; %name for the figure (remember to make it ending in '.nii'
    symmetric = 0; %if the AAL vector is submitted in symmetric order
    create_AALnifti(v,namenii,symmetric);
end
%memorized musical patterns vs baseline and novel musical patterns vs baseline
%loading IFC matrix of required contrasts (manually change which contrast you want..)
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_FC_Dynamic/Time_windows_5/Contr_1_3_STAT_zval_fullcontr_group_comp_1_1.mat')
STO = squeeze(mean(STAT_zval_fullcontr,2)); %getting mean of connections bewteen one ROI and the other ones over the 5 tones
load('/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/PaperNN1/MEG_FC_Dynamic/Time_windows_5/Contr_2_3_STAT_zval_fullcontr_group_comp_1_1.mat')
STN = squeeze(mean(STAT_zval_fullcontr,2)); %getting mean of connections bewteen one ROI and the other ones over the 5 tones
%significantly cenrtal ROIs for old and new used as mask for plotting their
%mean statistics within the brain
%here the mean statistics is for contrasting old vs baseline and new vs baseline (then put together in a unique figure with reversed sign for Bach variation (new))
for ii = 1:5
    v = zeros(90,1);
    v(ROI2(:,1,ii) == 1) = STO(ROI2(:,1,ii) == 1,ii); %here ROI2 correponds to old (and STO to old vs baseline contrast)
    v(ROI3(:,1,ii) == 1) = (STN(ROI3(:,1,ii) == 1,ii))*(-1); %here ROI3 correponds to new (and STO to new vs baseline contrast)
    V(:,ii) = v;
    namenii = ['Degree_WB_origvsvar_statvsbaseline_note' num2str(ii) '.nii']; %name for the figure (remember to make it ending in '.nii'
    symmetric = 0; %if the AAL vector is submitted in symmetric order
    create_AALnifti(v,namenii,symmetric);
end

%% calculating statistics for Kuramoto order parameter

%This sections reports the statistics and plot computed for Kuramoto order
%parameter. The plot is reported in the paper.
%Results from the statistics can be consulted in the structure
%"sign_clust" where the first column indicates the clusters number, the
%second the associated p-values and the third the significant time-points
%(here expressed in samples). If you want them expressed in seconds you
%should simply divide them by 150 (e.g. sign_clust{1,3}./150).
%Please note that also here, since we used MCS, the p-values may be
%slightly different than the one reported in the paper, because of the sthocastic nature of MCS.
%However those small differences do NOT affect the meaning of the results in any way.

%loading the previously calculated Kuramoto order parameter
load([pathdata '/MEG_FC_Dynamic/OP.mat']);
%calculating t-tests for each time-point
P = zeros(1,size(OP{1},2));
T = zeros(1,size(OP{1},2));
for ii = 1:size(OP{1},2)
    [~,p,~,stats] = ttest(OP{1}(:,ii),OP{2}(:,ii)); %contrasting cond1 vs cond2
    P(1,ii) = p;
    T(1,ii) = stats.tstat;
end
p_thresh = 0.05; %threshold for binarising p-values vector
Pbin = double((double(P < p_thresh) + double(T > 0)) == 2); %binarising p-values vector when the p-values are < p_thresh and t-values are positive (that means cond1 > cond2)
%Monte Carlo simulation function to correct for multiple comparisons
[ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 0, 1000, 0.05 )
%loading information
load([pathdata '/MEG_FC_Dynamic/S_original.mat']);
times = S.data.time_original(1,S.data.S_original_info.time_range{1}{1}(1) + S.left_bound - 1:S.data.S_original_info.time_range{1}{1}(2) - S.right_bound); %contort writing to get the proper time in seconds
stdeo = std(OP{1},0,1)/sqrt(size(OP{1},1)); %calculating standard errors
stden = std(OP{2},0,1)/sqrt(size(OP{1},1));
%plotting
figure
plot(times,squeeze(mean(OP{1},1)),'r','LineWidth',2)
hold on
plot(times,squeeze(mean(OP{1},1)) + stdeo,':r','LineWidth',0.5)
hold on
plot(times,squeeze(mean(OP{1},1)) - stdeo,':r','LineWidth',0.5)
hold on
plot(times,squeeze(mean(OP{2},1)),'b','LineWidth',2)
hold on
plot(times,squeeze(mean(OP{2},1)) + stden,':b','LineWidth',0.5)
hold on
xlim([times(1) times(end)])
%storing significant time-points
signtp = cell(1,3);
signtp(1,1) = sign_clust(1,3); signtp(1,2) = sign_clust(2,3); signtp(1,3) = sign_clust(3,3);
%plotting the significant time-points with gray shades
patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
ylims = get(gca,'YLim');
% ylims = 0.3;
for ii = 1:length(signtp)
    sgf2 = signtp{1,ii}./150;
    patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.5)
    hold on
end
plot(times,squeeze(mean(OP{1},1)),'r','LineWidth',2)
hold on
plot(times,squeeze(mean(OP{1},1)) + stdeo,':r','LineWidth',0.5)
hold on
plot(times,squeeze(mean(OP{1},1)) - stdeo,':r','LineWidth',0.5)
hold on
plot(times,squeeze(mean(OP{2},1)),'b','LineWidth',2)
hold on
plot(times,squeeze(mean(OP{2},1)) + stden,':b','LineWidth',0.5)
hold on
grid on
grid minor
xlim([times(1) times(end)])
ylim([0.12 0.28])
set(gcf,'color','w')

%% schemaball and connectivity in the brain (dynamic FC)

%%% FOR PLOTTING IN THE BRAIN YOU NEED OSL %%%

%Here we reshaped the results to be plotted into connectivity brain layout.
%You can choose the contrast that you want by indicating it in the line below.
%These plots are reported in the main results of the paper and in supplementary materials.

condit = '1_3'; %set 1_3 for memorized musical patterns vs baseline; set 2_3 for novel musical patterns vs baseline; set 1_2 for memorized vs novel musical patterns 

aosl = dir([pathl '/External/*osl']);
if ~isempty(aosl)
    %parcellation AAL
    %this is for giving input such as parcellation, ROI names and template (otherwise it can also find by itself the proper template, in this case 8mm)
    % parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'; %this and the following files can be found in the provided codes folder
    parcelfile = [pathl '/External/aal_8mm_try5.nii.gz']; %this and the following files can be found in the provided codes folder
    ROIsfile = [pathl '/External/aal_ROIs.txt'];
    templatefile = [pathl '/External/MNI152_T1_8mm_brain.nii.gz'];
    p = parcellation(parcelfile,ROIsfile,templatefile);
    p = p.remove_parcels(91:116); %removing the cerebellum
end
phase_dir = [pathdata '/MEG_FC_Dynamic/Time_windows_5']; %general path
load([phase_dir '/Contr_' condit '_STAT_zval_fullcontr_group_comp_1_1.mat']);
if strcmp(condit(3),'2') %scaling results for plotting purposes
    a = STAT_zval_fullcontr./4.5413; %this number represents the maximum (absolute) value of contrast 1_2
else
    a = STAT_zval_fullcontr./5.9905; %this number represents the maximum (absolute) value of contrasts 1_3 and 2_3
end
max(max(max(a)))
min(min(min(a)))
%actual plotting settings and function
for ii = 1:size(a,3)
    matf = a(:,:,ii);
    S = [];
    S.MATF = matf;
    S.symmetric_l = 0;
    S.outpath = pathdata;
    if ~isempty(aosl)
        S.p = p;
    else
        S.p = [];
    end
    S.p = p;
    S.thr_cortex = [0.993];
    S.spinlim = [-1 1];
    if ~isempty(aosl)
        S.name_gif = ['spinbrain_' condit '_note' num2str(ii)];
    else
        S.name_gif = [];
    end
%     S.name_gif = [];
    S.frame_vect = [15,90,130];
    S.fr_spec = {'Posterior_L_R';'Frontal_R_L';'Hem_L'};
    S.schball_l = 1;
    S.extr = [1.00001]; %matlab here considers 1 < 1.000.. I guess because of some rule that I am not aware of.. therefore I made this small trick to work around the problem
    S.lab_parc_ph = [pathl '/External/AAL_labels.mat']; %it can be found in the provided codes folder
    S.colbacksch = [];
    S.schemplth = 10;
    S.name_sch = [];
    %actual function
    FC_Brain_Schemaball_plotting_LBPD_D(S)
    title(num2str(ii))
end

%% coupling over time (reported in supplementary materials..)

%Here you can reproduce the results reported in the supplementary materials.
%Details about this analysis are also reported in supplementary materials.
plot_label = 1; %set 0 for computation only; set 1 for plotting only; set 2 for both computation and plotting

thr = 1; %set 1 for thresholding the data (with p_thresh); set 0 not to do that
Order = [1:2:90 90:-2:1];
p_thresh = 0.05;
load([pathdata '/MEG_FC_Dynamic/nosub/Contr_1_3_STAT_p_fullcontr_group_comp_1_1.mat']);
load([pathdata '/MEG_FC_Dynamic/nosub/Contr_1_3_STAT_zval_fullcontr_group_comp_1_1.mat']);
if thr == 1
    STAT_old = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
else
    STAT_old = STAT_zval_fullcontr;
end
load([pathdata '/MEG_FC_Dynamic/nosub/Contr_2_3_STAT_p_fullcontr_group_comp_1_1.mat']);
load([pathdata '/MEG_FC_Dynamic/nosub/Contr_2_3_STAT_zval_fullcontr_group_comp_1_1.mat']);
if thr == 1
    STAT_new = double((double(STAT_p_fullcontr < p_thresh) + double(STAT_zval_fullcontr > 0)) == 2);
else
    STAT_new = STAT_zval_fullcontr;
end
b = STAT_new(Order,Order,:);
a = STAT_old(Order,Order,:);
load([pathl '/External/AAL_labels.mat']); %it can be found in the provided codes folder
ROI_lab = cell(90,1);
for ii = 1:90
    ROI_lab(ii) = {lab(ii,:)};
end
X = 0;
permut = 10000;
threshMC = 6.2e-05;

[P_pos_fin, P_neg_fin, PM] = diff_conditions_matr_couplingROIs_MCS_LBPD_D(a, b, ROI_lab, X, permut, threshMC, plot_label);

%HERE YOU CAN GET THE EXACT COLORS REPORTED IN THE PAPER
%to get the red-blue-white color select the matrix that you want and then run the following lines
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))

%% plotting results of coupling over time in the brain and through schemaball

%%% FOR PLOTTING IN THE BRAIN YOU NEED OSL %%%

%%% you can reproduce here the schemaball and plots in the brain template reported in supplementary materials %%%

%loading results from previously calculated coupling analysis..
load([pathdata '/MEG_FC_Dynamic/nosub/COUP.mat']);
%reshaping matrices to prepare data for plotting
label_path = [pathl '/External/AAL_labels.mat'];
[ matf ] = generalCoupling_prepareplotting_LBPD_D( P_pos_fin, P_neg_fin, [], label_path, 0 );
aosl = dir([pathl '/External/*osl']);
if ~isempty(aosl)
    %parcellation AAL
    %this is for giving input such as parcellation, ROI names and template (otherwise it can also find by itself the proper template, in this case 8mm)
    % parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'; %this and the following files can be found in the provided codes folder
    parcelfile = [pathl '/External/aal_8mm_try5.nii.gz']; %this and the following files can be found in the provided codes folder
    ROIsfile = [pathl '/External/aal_ROIs.txt'];
    templatefile = [pathl '/External/MNI152_T1_8mm_brain.nii.gz'];
    p = parcellation(parcelfile,ROIsfile,templatefile);
    p = p.remove_parcels(91:116); %removing the cerebellum
end
%actual plotting function for schemaball and connectivity in the brain (dynamic FC)
%this provides the connectivity results of the coupling analysis displayed in the brain template and as schemaball
PM(1,3) = {matf}; %inserting in the coupling matrix only the significant values for the contrast
NAM{1} = 'Bachsoriginal'; NAM{2} = 'Bachsvariation'; NAM{3} = 'Bachsoriginalvsvariation';
%getting extreme values for plotting purposes
F(1) = (max(max(abs(PM{1})))); F(2) = (max(max(abs(PM{2}))));
EXTR{1} = [(max(F) * (-1)) max(F)]; %Bach's original and variation are scaled with the same numbers
EXTR{2} = [(max(F) * (-1)) max(F)];
EXTR{3} = [(max(max(abs(PM{3}))) * (-1)) max(max(abs(PM{3})))]; %Bach's original vs Bach variation is scaled independently
for ii = 1:3 %Bach's original, Bach variation, and their contrast
    nam = ['CouplingROIs_' NAM{ii}];
    S = [];
    S.MATF = PM{1,ii};
    if ii == 3 %the connectivity is stored differently in the matrices.. specifically Bach's original and variation are already symmetric, while their contrast is not..
        S.symmetric_l = 0;
    else
        S.symmetric_l = 1;
    end
    S.outpath = pathdata; %your OWN output directory
    if ~isempty(aosl)
        S.p = p;
    else
        S.p = [];
    end
    S.thr_cortex = [0.98]; f = num2str(S.thr_cortex);
    if ~isempty(aosl)
        S.name_gif = ['CouplingROIs_' NAM{ii}];
    else
        S.name_gif = [];
    end
    S.frame_vect = [15,90,130];
    S.fr_spec = {'Posterior_L_R';'Frontal_R_L';'Hem_L'};
    S.schball_l = 1;
    S.extr = EXTR{ii};
    S.lab_parc_ph = [pathl '/External/AAL_labels.mat']; %this is provided within the codes folder
    S.colbacksch = [];
    S.schemplth = 100;
    S.name_sch = [];
    %actual function
    FC_Brain_Schemaball_plotting_LBPD_D(S)
end

%%


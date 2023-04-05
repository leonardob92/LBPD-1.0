%%

%% SYSTEMATIC VARIATION OF AUDITORY SEQUENCES (AUDITORY PATTERN RECOGNITION 2020) - LEONARDO BONETTI

%%

%% *** START UP FUNCTIONS.. (LBPD_startup_D) ***

%starting up some functions for LBPD toolbox.

%starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl);

%%

%% *** PREPROCESSING ***

%%% OBS!! the preprocessing was computed for (nearly) all the data
%%% collected and not only for the data that was actually analyzed and
%%% reported in this paper.

%%

%% Maxfilter

%OBS! before running maxfilter you need to close matlab, open the terminal and write: 'use anaconda', then open matlab and run maxfilter script

maxfilter_path = '/neuro/bin/util/maxfilter';
project = 'MINDLAB2020_MEG-AuditoryPatternRecognition';
maxDir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2'; %output path

path = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition'; %path with all the subjects folders
jj = dir([path '/0*']); %list all the folders starting with '0' in order to avoid hidden files
for ii = 90:length(jj) %over subjects (STARTING FROM 13, SINCE THE PREVIOUS ONES WERE THE PILOTS (IF THEY DID ALSO THE PROPER EXPERIMENT, THE MAXFILTER COMPUTATION WILL BE DONE IN THE NEXT  SECTION))
    cart = [jj(ii).folder '/' jj(ii).name]; %create a path that combines the folder name and the file name
    pnana = dir([cart '/2*']); %search for folders starting with '2'
    for pp = 1:length(pnana) %loop to explore ad analyze all the folders inside the path above
        cart2 = [pnana(1).folder '/' pnana(pp).name];
        pr = dir([cart2 '/ME*']); %looks for meg folder
        if ~isempty(pr) %if pr is not empty, proceed with subfolders inside the meg path
            pnunu = dir([pr(1).folder '/' pr(1).name '/00*']);
            %if length(pnunu) > 1
            %warning(['subj ' num2str(ii) ' has files nuber = ' num2str(length(pnunu))]) %show a warning message if any subj has more thatn 1 meg sub-folder
            %end
            for dd = 1:length(pnunu)
                if strcmp(pnunu(dd).name(5:6),'re') || strcmp(pnunu(dd).name(5:6),'sa') || strcmp(pnunu(dd).name(5:6),'vi') || strcmp(pnunu(dd).name(5:6),'pd') %checks whether characters 5 to 6 are equal to 're', 'vi, 'sa' or 'pd'; the loop continues if this is true (1) and it stops if this is false (0)
                    %idx2 = strfind(pnunu(1).name,'Mus'); % search for musmelo folder in order to avoid other projects
                    %if ~isempty(idx2)
                    fpath = dir([pnunu(1).folder '/' pnunu(dd).name '/files/*.fif']); % looks for .fif file
                    rawName = ([fpath.folder '/' fpath.name]); %assigns the final path of the .fif file to the rawName path used in the maxfilter command
                    maxfName = ['SUBJ' jj(ii).name '_' fpath.name(1:end-4)]; %define the output name of the maxfilter processing
                    %movement compensation
                    cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
                    %no movement compensation (to be used if HPI coils did not work properly)
%                     cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
                    system(cmd);
                end
            end
        end
    end
end

%% Maxfilter - Subjects who did the experiment (version 2.0) after the first pilot (version 1.0)

%%% DONE SUBJECTS 2 10 11 7 %%%

%OBS! before running maxfilter you need to close matlab, open the terminal and write: 'use anaconda', then open matlab and run maxfilter script

%settings
maxfilter_path = '/neuro/bin/util/maxfilter';
project = 'MINDLAB2020_MEG-AuditoryPatternRecognition';
maxDir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2'; %output path

clear rawName maxfName
%path to raw files for these particular subjects
pathraw{1} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0002/20210413_000000/MEG'; %assigns the final path of the .fif file to the rawName path used in the maxfilter command
pathraw{2} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0010/20210416_000000/MEG';
pathraw{3} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0011/20210423_000000/MEG';
pathraw{4} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0007/20210511_000000/MEG';
pathraw{5} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0012/20210622_000000/MEG';

for ii = 5:length(pathraw) %over subjects
    list = dir([pathraw{ii} '/0*']);
    for jj = 1:length(list) %over experimental blocks
        if strcmp(list(jj).name(5:6),'re') || strcmp(list(jj).name(5:6),'sa') || strcmp(list(jj).name(5:6),'vi') || strcmp(list(jj).name(5:6),'pd') %checks whether characters 5 to 6 are equal to 're', 'vi, 'sa' or 'pd'; the loop continues if this is true (1) and it stops if this is false (0)
            
            fpath = dir([list(jj).folder '/' list(jj).name '/files/*.fif']); % looks for .fif file
            rawName = [fpath(1).folder '/' fpath(1).name];
            maxfName = ['SUBJ' pathraw{ii}(56:59) '_' fpath.name(1:end-4) '_bis']; %define the output name of the maxfilter processing
            %movement compensation
            cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
            %no movement compensation
%             cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
            system(cmd);
        end
    end
end

%% Converting the .fif files into SPM objects

%OBS! remember to run 'starting up OSL' first

%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

%%

fif_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/*.fif'); %creates a list with the .fif files

for ii = 340:341%:length(fif_list) %over the .fif files
    S = []; %structure 'S'                   
    S.dataset = [fif_list(ii).folder '/' fif_list(ii).name];
    D = spm_eeg_convert(S);
%     D = job2cluster(@cluster_spmobject, S); %actual function for conversion
end

%% Removing bad segments using OSLVIEW

%checks data for potential bad segments (periods)
%marking is done by right-clicking in the proximity of the event and click on 'mark event'
%a first click (green dashed label) marks the beginning of a bad period
%a second click indicates the end of a bad period (red)
%this will mean that we are not using about half of the data, but with such bad artefacts this is the best we can do
%we can still obtain good results with what remains
%NB: Push the disk button to save to disk (no prefix will be added, same name is kept)

%OBS! remember to check for bad segments of the signal both at 'megplanar' and 'megmag' channels (you can change the channels in the OSLVIEW interface)
%OBS! remember to mark the trial within the bad segments as 'badtrials' and use the label for removing them from the Averaging (after Epoching) 

spm_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/spmeeg*.mat'); %path to SPM objects

for ii = 1:3%:length(spm_list) %over experimental blocks %OBS!
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = oslview(D);
    D.save(); %save the selected bad segments and/or channels in OSLVIEW
    disp(ii)
end

%% AFRICA denoising (part I)

%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

%%

%ICA calculation
spm_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/spmeeg*.mat');

for ii = 339:length(spm_list) %OBS!
    S = [];
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    S.D = D;
    
    jobid = job2cluster(@cluster_africa,S);
%   D = osl_africa(D,'do_ica',true,'do_ident',false,'do_remove',false,'used_maxfilter',true); 
%   D.save();
end

%% AFRICA denoising (part II)

% v = [11 12 19 32];
%visual inspection and removal of artifacted components
%look for EOG and ECG channels (usually the most correlated ones, but check a few more just in case)
spm_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/spmeeg*.mat');

for ii = 339:length(spm_list) %OBS!%38:41
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = osl_africa(D,'do_ident','manual','do_remove',false,'artefact_channels',{'EOG','ECG'});
    %hacking the function to manage to get around the OUT OF MEMORY problem..
    S = [];
    S.D = D;
    jobid = job2cluster(@cluster_rembadcomp,S);
%   D.save();
    disp(ii)
end

%% Epoching: one epoch per old/new excerpt (baseline = (-)100ms)

prefix_tobeadded = 'e'; %adds this prefix to epoched files
spm_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/spmeeg*.mat');

for ii = 339:length(spm_list) %over .mat files
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]); %load spm_list .mat files
    dummy = D.fname; %OBS! D.fname does not work, so we need to use a 'dummy' variable instead
    %if strcmp(dummy(22:26), 'speed') %checks whether characters 22 to 26 are equal to 'speed'; the loop continues if this is true (1) and it stops if this is false (0)
    events = D.events; %look for triggers
    %takes the correct triggers sent during the recording
    clear trigcor
    count_evval = 0; %???
    for ieve = 1:length(events) %over triggers
        if strcmp(events(ieve).type,'STI101_up') %only triggers at the beginning of each stimuli
            if events(ieve).value ~= 103 && events(ieve).value ~= 104 && events(ieve).value ~= 128 && events(ieve).value ~= 8 && events(ieve).value ~= 132 && events(ieve).value ~= 48 && events(ieve).value ~= 32 && events(ieve).value ~= 64 %discard 104 and 128 for random triggers
                count_evval = count_evval + 1;
                trigcor(count_evval,1) = events(ieve).time; %+ 0.010; %this takes the correct triggers and add 10ms of delay of the sound travelling into the tubes
                %variable with all the triggers we need
            end
        end
    end
    trl_sam = zeros(length(trigcor),3); %prepare the samples matrix with 0's in all its cells
    trl_sec = zeros(length(trigcor),3); %prepare the seconds matrix with 0's in all its cells
    %deftrig = zeros(length(trigcor),1); %this is not useful
    for k = 1:length(trigcor) %over selected triggers
        %deftrig(k,1) = 0.012 + trigcor(k,1); %adding a 0.012 seconds delay to the triggers sent during the experiment (this delay was due to technical reasons related to the stimuli)
        trl_sec(k,1) = trigcor(k,1) - 0.1000; %beginning time-window epoch in s (please note that we computed this operation two times, obtaining two slightly different pre-stimulus times.
        %this was done because for some computations was convenient to have a slightly longer pre-stimulus time
        %remove 1000ms of baseline
        if strcmp(dummy(22:26), 'speed')
            trl_sec(k,2) = trigcor(k,1) + 5.5; %end time-window epoch in seconds
        else
            trl_sec(k,2) = trigcor(k,1) + 4.4; %end time-window epoch in seconds
        end
        trl_sec(k,3) = trl_sec(k,2) - trl_sec(k,1); %range time-windows in seconds
        trl_sam(k,1) = round(trl_sec(k,1) * 250) + 1; %beginning time-window epoch in samples %250Hz per second
        trl_sam(k,2) = round(trl_sec(k,2) * 250) + 1; %end time-window epoch in samples
        trl_sam(k,3) = -25; %sample before the onset of the stimulus (corresponds to 0.100ms)
    end
    dif = trl_sam(:,2) - trl_sam(:, 1); %difference between the end and the beginning of each sample (just to make sure that everything is fine)
    if ~all(dif == dif(1)) %checking if every element of the vector are the same (i.e. the length of the trials is the same; we may have 1 sample of difference sometimes because of different rounding operations..)
        trl_sam(:,2) = trl_sam(:,1) + dif(1);
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
    %take bad segments registered in OSLVIEW and check if they overlap with the trials. if so, it gives the number of overlapped trials that will be removed later   
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
    disp(spm_list(ii).name);
    if count == 0
        disp('there are no bad trials marked in oslview');
    else
        D = badtrials(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
%         D = conditions(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
        epochinfo = D.epochinfo;
        xcv = find(Bad_trials == 1);
        %this should be done only later.. in any case.. not a problem..
        for jhk = 1:length(xcv)
            D = D.conditions(xcv(jhk),'Bad');
            epochinfo.conditionlabels(xcv(jhk)) = {'Bad'};
            disp([num2str(ii) ' - ' num2str(jhk) ' / ' num2str(length(xcv))])
        end
        D.epochinfo = epochinfo;
        D.save(); %saving on disk
        disp('bad trials are ')
        length(D.badtrials)
    end
    D.save();
    disp(ii)
end


%% Defining the conditions - All blocks

%define conditions - 1 epoch for each old/new excerpt (baseline = (-)100ms)

xlsx_dir_behav = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/BehavioralTaskMEG/Version_2/Final_xlsx'; %dir to MEG behavioral results (.xlsx files)
epoch_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*.mat'); %dir to epoched files

for ii = 339:length(epoch_list) %over epoched data
    D = spm_eeg_load([epoch_list(ii).folder '/' epoch_list(ii).name]);
    dummy = D.fname;
    %barbaric solution.. to build the name to be read for the excel files with the MEG behavioral tasks performance
    if strcmp(dummy(18:23),'recogm')
        dumbloc = 'Block_3.xlsx';
        bl = 3;
    elseif strcmp(dummy(18:23),'recogs')
        dumbloc = 'Block_4.xlsx';
        bl = 4;
    elseif strcmp(dummy(18:23),'sameme')
        dumbloc = 'Block_5.xlsx';
        bl = 5;
    elseif strcmp(dummy(18:23),'visual')
        dumbloc = 'Block_6.xlsx';
        bl = 6;
    elseif strcmp(dummy(18:19),'pd')
        dumbloc = 'Project_PD.xlsx';
        bl = 7;
    end
    if strcmp(dummy((end-14):(end-14)+2),'bis')
        dumls = ['Subj_' dummy(13:16) 'bis_' dumbloc]; %getting subject ID directly from the SPM object to reduce the probability to make mistakes
    else
        dumls = ['Subj_' dummy(13:16) '_' dumbloc];
    end
    [~,~,raw_recog] = xlsread([xlsx_dir_behav '/' dumls]); %excel files
    %picking the current block
    if bl == 3 %block 3 (THIS IS THE ACTUAL BLOCK FOR THIS PAPER)
        for k = 1:length(D.trialonset)
            if raw_recog{(k + 1),3} == 0 %if there was no response
                D = D.conditions(k,'No_response');
            elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 1 %old correct
                D = D.conditions(k,'Old_Correct'); %assign old correct
            elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 2 %old incorrect
                D = D.conditions(k,'Old_Incorrect'); %otherwise assign new correct
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 2 %new t1 correct
                D = D.conditions(k,'New_T1_Correct');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 1 %new t1 incorrect
                D = D.conditions(k,'New_T1_Incorrect');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t2') && raw_recog{(k + 1),3} == 2 %new t2 correct
                D = D.conditions(k,'New_T2_Correct');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t2') && raw_recog{(k + 1),3} == 1 %new t2 incorrect
                D = D.conditions(k,'New_T2_Incorrect');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 2 %new t3 correct
                D = D.conditions(k,'New_T3_Correct');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 1 %new t3 incorrect
                D = D.conditions(k,'New_T3_Incorrect');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t4') && raw_recog{(k + 1),3} == 2 %new t4 correct
                D = D.conditions(k,'New_T4_Correct');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t4') && raw_recog{(k + 1),3} == 1 %new t4 incorrect
                D = D.conditions(k,'New_T4_Incorrect');
            end
        end
    elseif bl == 4 %block 4
        for k = 1:length(D.trialonset)
            if raw_recog{(k + 1),3} == 0 %if there was no response
                D = D.conditions(k,'No_response');
            elseif strcmp(raw_recog{(k + 1),2}(6:13),'fast_old') && raw_recog{(k + 1),3} == 1 %old fast correct
                D = D.conditions(k,'Old_Fast_Correct'); %assign old fast correct
            elseif strcmp(raw_recog{(k + 1),2}(6:13),'fast_old') && raw_recog{(k + 1),3} == 2 %old fast incorrect
                D = D.conditions(k,'Old_Fast_Incorrect'); %assign old fast incorrect
            elseif strcmp(raw_recog{(k + 1),2}(6:13),'slow_old') && raw_recog{(k + 1),3} == 1 %old slow correct
                D = D.conditions(k,'Old_Slow_Correct'); %assign old slow correct
            elseif strcmp(raw_recog{(k + 1),2}(6:13),'slow_old') && raw_recog{(k + 1),3} == 2 %old slow incorrect
                D = D.conditions(k,'Old_Slow_Incorrect'); %assign old slow incorrect
            elseif strcmp(raw_recog{(k + 1),2}(8:15),'fast_new') && raw_recog{(k + 1),3} == 2 %new fast correct
                D = D.conditions(k,'New_Fast_Correct'); %assign new fast correct
            elseif strcmp(raw_recog{(k + 1),2}(8:15),'fast_new') && raw_recog{(k + 1),3} == 1 %new fast incorrect
                D = D.conditions(k,'New_Fast_Incorrect'); %assign new fast incorrect
            elseif strcmp(raw_recog{(k + 1),2}(8:15),'slow_new') && raw_recog{(k + 1),3} == 2 %new slow correct
                D = D.conditions(k,'New_Slow_Correct'); %assign new slow correct
            else
                D = D.conditions(k,'New_Slow_Incorrect'); %assign new incorrect
            end
        end
    elseif bl == 5 %block 5
        for k = 1:20 %over the 20 encoding trials
            D = D.conditions(k,'Encoding');
        end
        for k = 1:(length(D.trialonset)-20) %over recognition trials (this would work even if the last trial was out of bonds and thus was not epoched..)
            if raw_recog{(k + 23),3} == 0 %if there was no response
                D = D.conditions(k + 20,'No_response');
            elseif strcmp(raw_recog{(k + 23),2}(9:11),'enc') && raw_recog{(k + 23),3} == 1 %encoding correct
                D = D.conditions(k + 20,'Old_Correct'); %assign encoding correct
            elseif strcmp(raw_recog{(k + 23),2}(9:11),'enc') && raw_recog{(k + 23),3} == 2 %encoding incorrect
                D = D.conditions(k + 20,'Old_Incorrect'); %otherwise assign encoding incorrect
            elseif strcmp(raw_recog{(k + 23),2}(9:11),'rec') && raw_recog{(k + 23),3} == 2 %recognition correct
                D = D.conditions(k + 20,'New_Correct'); %assign recognition correct
            else
                D = D.conditions(k + 20,'New_Incorrect'); %assign new incorrect
            end
        end
    elseif bl == 6 %block 6
        for k = 1:20 %over the 20 encoding trials
            D = D.conditions(k,'Encoding');
        end
        for k = 1:(length(D.trialonset)-20) %over recognition trials (this would work even if the last trial was out of bonds and thus was not epoched..)
            if raw_recog{(k + 23),3} == 0 %if there was no response
                D = D.conditions(k + 20,'No_response');
            elseif strcmp(raw_recog{(k + 23),2}(8:10),'old') && raw_recog{(k + 23),3} == 1 %encoding correct
                D = D.conditions(k + 20,'Old_Correct'); %assign encoding correct
            elseif strcmp(raw_recog{(k + 23),2}(8:10),'old') && raw_recog{(k + 23),3} == 2 %encoding incorrect
                D = D.conditions(k + 20,'Old_Incorrect'); %otherwise assign encoding incorrect
            elseif strcmp(raw_recog{(k + 23),2}(8:10),'new') && raw_recog{(k + 23),3} == 2 %recognition correct
                D = D.conditions(k + 20,'New_Correct'); %assign recognition correct
            else
                D = D.conditions(k + 20,'New_Incorrect'); %assign new incorrect
            end
        end
    elseif bl == 7
        for k = 1:length(D.trialonset) %over trials
            if strcmp(raw_recog{(k + 1),2}(3:7),'porco')
                D = D.conditions(k,'pd');
            else
                D = D.conditions(k,'pm');
            end
        end
    end
    %this is for every block
    if ~isempty(D.badtrials) %overwriting badtrials (if any) on condition labels
        BadTrials = D.badtrials;
        for badcount = 1:length(BadTrials) %over bad trials
            D = D.conditions(BadTrials(badcount),'Bad_trial');
        end
    end
    D = D.montage('switch',1);
    D.epochinfo.conditionlabels = D.conditions; %to add for later use in the source reconstruction
    D.save(); %saving data on disk
    disp(num2str(ii))
end

%% COMPUTING STATISTICS OF BEHAVIORAL TASKS IN MEG

xlsx_dir_behav = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/BehavioralTaskMEG/Version_2/Final_xlsx/'; %dir to MEG behavioral results (.xlsx files)
bb = 3;
list_beh = dir([xlsx_dir_behav 'Subj*' num2str(bb) '.xlsx']); %normal blocks
bl = bb;
if bl == 3
    Block_3 = cell(length(list_beh)+1,13);
elseif bl == 4
    Block_4 = cell(length(list_beh)+1,11);
elseif bl == 5
    Block_5 = cell(length(list_beh)+1,7);
elseif bl == 6
    Block_6 = cell(length(list_beh)+1,7);
end
for ii = 1:length(list_beh) %over subjects for block bb
    %barbaric solution.. to build the name to be read for the excel files with the MEG behavioral tasks performance
    [~,~,raw_recog] = xlsread([list_beh(ii).folder '/' list_beh(ii).name]); %excel files
    %picking the current block
    %legend
    Block_3{1,1} = 'Subject'; Block_3{1,2} = 'OLD_Cor'; Block_3{1,3} = 'OLD_Cor %'; Block_3{1,4} = 'New_T1_Cor'; Block_3{1,5} = 'New_T1_Cor %'; %1st row
    Block_3{1,6} = 'New_T2_Cor'; Block_3{1,7} = 'New_T2_Cor %'; Block_3{1,8} = 'New_T3_Cor'; Block_3{1,9} = 'New_T3_Cor %'; Block_3{1,10} = 'New_T4_Cor'; Block_3{1,11} = 'New_T4_Cor %'; Block_3{1,12} = 'No response'; Block_3{1,13} = 'No response %'; %1st row
    Block_3{1,14} = 'OLD_RT'; Block_3{1,15} = 'NEWT1_RT'; Block_3{1,16} = 'NEWT2_RT'; Block_3{1,17} = 'NEWT3_RT'; Block_3{1,18} = 'NEWT4_RT';
    nr = 0; old = 0; n1 = 0; n2 = 0; n3 = 0; n4 = 0;
    ort = []; n1rt = []; n2rt = []; n3rt = []; n4rt = [];
    for k = 1:135
        if raw_recog{(k + 1),3} == 0 %if there was no response
            nr = nr + 1;
        elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 1 %old correct
            old = old + 1;
            ort = cat(1,ort,raw_recog{(k + 1),4});
        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 2 %new t1 correct
            n1 = n1 + 1;
            n1rt = cat(1,n1rt,raw_recog{(k + 1),4});
        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t2') && raw_recog{(k + 1),3} == 2 %new t2 correct
            n2 = n2 + 1;
            n2rt = cat(1,n2rt,raw_recog{(k + 1),4});
        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 2 %new t3 correct
            n3 = n3 + 1;
            n3rt = cat(1,n3rt,raw_recog{(k + 1),4});
        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t4') && raw_recog{(k + 1),3} == 2 %new t4 correct
            n4 = n4 + 1;
            n4rt = cat(1,n4rt,raw_recog{(k + 1),4});
        end
    end
    disp(num2str(['Block ' num2str(bb) ' - Subject ' num2str(ii)]))
    Block_3{ii+1,1} = list_beh(ii).name(6:9); Block_3{ii+1,2} = old; Block_3{ii+1,3} = (old/27)*100; Block_3{ii+1,4} = n1; Block_3{ii+1,5} = (n1/27)*100;
    Block_3{ii+1,6} = n2; Block_3{ii+1,7} = (n2/27)*100; Block_3{ii+1,8} = n3; Block_3{ii+1,9} = (n3/27)*100; Block_3{ii+1,10} = n4; Block_3{ii+1,11} = (n4/27)*100; Block_3{ii+1,12} = nr; Block_3{ii+1,13} = (nr/135)*100;
    Block_3{ii+1,14} = mean(ort); Block_3{ii+1,15} = mean(n1rt); Block_3{ii+1,16} = mean(n2rt); Block_3{ii+1,17} = mean(n3rt); Block_3{ii+1,18} = mean(n4rt);
end
Block_3_t = cell2table(Block_3); %remove the possible empty cell

%statistics (ANOVAs)
data = cell2mat(Block_3(2:end,[2 4 6 8 10])); %correct responses
datart = cell2mat(Block_3(2:end,14:18)); %reaction times
%correct responses
% [p,t,stats] = anova1(data); %'off' for not showing the plot
[p,t,stats] = kruskalwallis(data); %'off' for not showing the plot
[c,m,h,nms] = multcompare(stats,'ctype','tukey-kramer'); %perform multiple comparison test based on anova1 output = post-hoc analysis (c needs to be saved for every time point and every channel)
%reaction times
% [p,t,stats] = anova1(datart); %'off' for not showing the plot
[p,t,stats] = kruskalwallis(datart); %'off' for not showing the plot
[c,m,h,nms] = multcompare(stats,'ctype','tukey-kramer'); %perform multiple comparison test based on anova1 output = post-hoc analysis (c needs to be saved for every time point and every channel)

%plotting correct responses
%some default color specifications for later plotting..
% cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
cl = [0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8;];
datadum = cell(1,size(data,2));
for ii = 1:size(data,2)
    datadum{1,ii} = data(:,ii);
end
figure
rm_raincloud2(datadum',cl)
grid minor
set(gcf,'color','w')
set(gcf,'Position',[200,200,400,550])
xlabel('Correct responses'); %set(gca,'YDir','normal');
xlim([-5 38])

%plotting reaction times
%some default color specifications for later plotting..
% cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
cl = [0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8;];
datadum = cell(1,size(datart,2));
for ii = 1:size(datart,2)
    datadum{1,ii} = datart(:,ii);
end
figure
rm_raincloud2(datadum',cl)
grid minor
set(gcf,'color','w')
set(gcf,'Position',[200,200,400,550])
xlabel('Time (ms)'); %set(gca,'YDir','normal');
xlim([1500 3500])

%% correlation between WM and MEG behavioral data
list_WM = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/gemma/AuditoryPatternRecognition/behavioral_data/WM_APR2020_Edited_SR.xlsx'); %dir to working memory excel file
[~,~,WM] = xlsread([list_WM(1).folder '/' list_WM(1).name]); %loads working memory data
%OLD
[rho,pval] = corr(cell2mat(WM(2:end,2)),cell2mat(Block_3(2:end,2)),'rows','complete')
%NEWT1
[rho,pval] = corr(cell2mat(WM(2:end,2)),cell2mat(Block_3(2:end,4)),'rows','complete')
%NEWT2
[rho,pval] = corr(cell2mat(WM(2:end,2)),cell2mat(Block_3(2:end,6)),'rows','complete')
%NEWT3
[rho,pval] = corr(cell2mat(WM(2:end,2)),cell2mat(Block_3(2:end,8)),'rows','complete')
%NEWT4
[rho,pval] = corr(cell2mat(WM(2:end,2)),cell2mat(Block_3(2:end,10)),'rows','complete')

%% Averaging and Combining planar gradiometers

%settings for cluster (parallel computing)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster'); %set automatically the long run queue
clusterconfig('long_running', 1); %set automatically the long run queue
clusterconfig('slot', 1); %set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)

%% averaging

output_prefix_to_be_set = 'm';

v = [1]; %only a selection of files

epoch_list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*mat'); %dir to epoched files (encoding)
for ii = 1:length(v)%1:length(epoch_list) %over epoched files
    %distribute 
    input = [];
    input.D = [epoch_list(v(ii)).folder '/' epoch_list(v(ii)).name];
    input.prefix = output_prefix_to_be_set;
    jobid = job2cluster(@sensor_average, input); % this is the command for send the job to the cluster, in the brackets you can find the name on the function to run (afeter the @) and the variable for the input (in this case input)
    % look the script for more details about the function work
end

%% combining planar gradiometers

average_list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/m*mat'); %dir to epoched files (encoding)
v = [339:342]; %only a selection of files

for ii = 1:length(v)%1:length(average_list) %over files
    input = [];
    input.D = [average_list(v(ii)).folder '/' average_list(v(ii)).name];
    D = spm_eeg_load(input.D);
    D = D.montage('switch',1);
    D.save();
    jobid = job2cluster(@combining_planar_cluster, input); % this is the command for send the job to the cluster, in the brackets you can find the name on the function to run (afeter the @) and the variable for the input (in this case input)
end

%% LBPD_startup_D

pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path where is the function for submit the jobs to the server

%% Extracting and plotting MEG sensor data

%%% This is currently configured for plotting MEG channel MEG0211 (= channels_plot = 9) 

block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd
channels_plot = [9]; %13;95;;9; %95, 13, 15, 11, 141, 101, 9 % empty for plotting single channels; otherwise number(s) of channels to be averaged and plotted (e.g. [13] or [13 18])
% channels_plot = []; % empty for plotting single channels; otherwise number(s) of channels to be averaged and plotted (e.g. [13] or [13 18])

%1321 1411
save_data = 0;
load_data = 1; %set 1 if you want to load the data instead of extracting it from SPM objects
% v = [1]; %subjects
%bad 8,9 and a bit 6 (recogminor)

S = [];
%computing data
if block == 3
    S.conditions = {'Old_Correct','New_T1_Correct','New_T2_Correct','New_T3_Correct','New_T4_Correct'};
%     S.conditions = {'Old_Correct','New_T1_Correct'};
elseif block == 4
    % S.conditions = {'Old_Fast_Correct','Old_Slow_Correct','New_Fast_Correct','New_Slow_Correct'};
%     S.conditions = {'Old_Fast_Correct','New_Fast_Correct'};
    S.conditions = {'Old_Slow_Correct','New_Slow_Correct'};
elseif block == 5
    S.conditions = {'Old_Correct','New_Correct'};
%     S.conditions = {'Old_Correct','Encoding'};
elseif block == 6
    S.conditions = {'Old_Correct','New_Correct'};
%     S.conditions = {'Old_Correct','Encoding'};
elseif block == 7
    S.conditions = {'pd','pm'};
end
if block == 3
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*_recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 4
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*recogspeed*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 5
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*samemel*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 6
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*visualpat*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 7
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/P*pd*_tsssdsm.mat'); %dir to epoched files (encoding)
end
v = 1:length(list); %subjects
% v = [2];
if ~exist('chanlabels','var')
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter/MEG_sensors/recogminor_all_conditions.mat', 'chanlabels')
end
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/MEG_sensors'; %path to write output in
S.outdir = outdir;
S.data = [];
if load_data == 1 %if you already computed and saved on disk the t-tests you can load them here
%     loadlead([outdir '/Block_' num2str(block) '.mat']);
    S.data = data_mat(:,:,v,:);
    S.chanlabels = chanlabels;
    S.time_real = time_sel;
else %otherwise you can extract the data from SPM MEEG objects (one for each subject)
%     S.spm_list = cell(1,length(list));
% v = 7;
    S.spm_list = cell(1,length(v));
    for ii = 1:length(v)
        S.spm_list(ii) = {[list(v(ii)).folder '/' list(v(ii)).name]};
    end
end

S.timeextract = []; %time-points to be extracted
S.centerdata0 = 0; %1 to make data starting at 0
S.save_data = save_data; %only meaningfull if you read data from SPM objects saved on disk
S.save_name_data = ['Block_' num2str(block)];

%individual waveform plotting
if isempty(channels_plot)
    S.waveform_singlechannels_label = 0; %1 to plot single channel waveforms
else
    S.waveform_singlechannels_label = 0; %1 to plot single channel waveforms
end
S.wave_plot_conditions_together = 0; %1 for plotting the average of all
S.mag_lab = 1; %1 for magnetometers; 2 for gradiometers
S.x_lim_temp_wave = [-0.1 3.4]; %limits for time (in secs) (E.g. [-0.1 3.4])
S.y_lim_ampl_wave = []; %limit for amplitude (E.g. [0 120] magnetometes, [0 6] gradiometers)

%averaged waveform plotting
if isempty(channels_plot)
    S.waveform_average_label = 0; %average of some channels
    S.left_mag = 95; %13 %37 (visual) %43 (visual) %199 (visual) %203 (visual) %channels for averaging
else
    S.waveform_average_label = 0; %average of some channels
    S.left_mag = channels_plot; %13 %37 (visual) %43 (visual) %199 (visual) %203 (visual) %channels for averaging
end
% S.left_mag = [2:2:204];
S.legc = 1; %set 1 for legend
% S.left_mag = 99;
S.signtp = {[]};
% S.sr = 150; %sampling rate (Hz)
S.avewave_contrast = 0; %1 to plot the contrast between conditions (averaged waveform)
S.save_label_waveaverage = 0;
S.label_plot = 'c';
%t-tests
S.t_test_for_permutations = 0;
S.cond_ttests_tobeplotted_topoplot = [1 2]; %this is for both topoplot and t-tests!! (here [1 2] means cond1 vs cond2!!!!!!!)

%topoplotting
S.topoplot_label = 1;
S.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External';
S.topocontr = 0;
S.topocondsing = [1]; %condition for topoplot
% S.xlim = [0.75 0.85]; %time topolot
% S.xlim = [1.1 1.2]; %time topolot
S.xlim = [0.25 0.25]; 
S.zlimmag = []; %magnetometers amplitude topoplot limits
S.zlimgrad = []; %gradiometers amplitude topoplot limits
S.colormap_spec = 0;
% x = []; x.bottom = [0 0 1]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 1 0.5]; x.top = [1 0.95 0]; %yellow - blue
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
S.colormap_spec_x = x;
S.topoplot_save_label = 0;
color_line = colormap(lines(5)); %extracting some colours from a colormap
S.color_line = color_line;
S.color_line(1,:) = color_line(2,:);
S.color_line(2,:) = color_line(1,:);
S.color_line(5,:) = [0.4 0.4 0.4];
 

[out] = MEG_sensors_plotting_ttest_LBPD_D2(S);

%%

%% *** DECODING (MULTIVARIATE PATTERN ANALYSIS ***

%%

%% decoding - support vector machine (SVM) - multivariate pattern analysis - FIGURE 3A

%The decoding consisted of support vector machines (SVM) implemented in
%external functions that can be found at the following link:
%http://www.csie.ntu.edu.tw/~cjlin/libsvm/
%Those functions can be also used in Matlab with a few adjustments done to achieve a proper imlplementation
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
make

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')
% config of the cluster server
clusterconfig('slot', 1); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
clusterconfig('scheduler','cluster'); % set automatically the long run queue
clusterconfig('long_running', 1); % set automatically the long run queue


%% preparing data for decoding functions (for AU server)

time = [1 1126]; %defining time of interest (this refers to 0.100sec before the onset of the stimuli)
numperm = 100; %number of permutations
kfold = 5; %number of categories for k-fold classification
pairl = 1; %1 = pairwise decoding; 2 = multiclass (confusion matrix)

list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*_recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
savedirdum = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Decoding';
D1 = cell(length(list),2);
subjs_real = cell(1,length(list));
for ii = 1:length(list) %length(subjnum)
        %loading data
        D = spm_eeg_load([list(ii).folder '/' list(ii).name]);
        old = find(strcmp('Old_Correct',D.conditions)); %getting condition old indices
        for nn = 1:4 %over NewTx conditions
            new = find(strcmp(['New_T' num2str(nn) '_Correct'],D.conditions)); %getting condition newTx indices
            savedir = [savedirdum '/Old_vs_NewT' num2str(nn)];
            mkdir(savedir);
            clear data condid ind
            condid(1:length(old)) = {'O'}; %creating condition labels (first old)
            condid(length(old):(length(old)+length(new))) = {'N'}; %(then new)
            ind = [old new]; %rearranging indices (first old and then new)
            idxMEGc1 = find(strcmp(D.chanlabels,'MEG0111')); %getting extreme channels indexes
            idxMEGc2 = find(strcmp(D.chanlabels,'MEG2643')); %same
            data = D(idxMEGc1:idxMEGc2,time(1):time(2),ind); %extracting data
            %structure with inputs for AU server
            S = [];
            S.pairl = pairl;
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

T_cond = 0; %1 = Old vs NewT1; 2 = Old vs NewT2; 3 = Old vs NewT3; 4 = Old vs NewT4; 0 = all together (only waveform plot);
tempgen_l = 0; %1 computing (and plotting) statistics for temporal generalization; 0 = not doing it (it takes time..)
topoplot_l = 0; %1 for topoplot; 0 for no topoplot (remember that you also have to specify the time for the topo-plot)


%directory where decoding results are stored in
datadirdum = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Decoding';
%loading and reshpaping data outputted by AU server (SVM algorithm) for statistics
load([datadirdum '/time.mat']);
if T_cond ~= 0
    datadir = [datadirdum '/Old_vs_NewT' num2str(T_cond)];
    listPD = dir([datadir '/PD*']);
    listGT = dir([datadir '/TG*']);
    %loading data
    dd1 = zeros(length(time_sel),length(listPD));
    ddpatt = zeros(306,length(time_sel),length(listPD));
    ddTG = zeros(length(time_sel),length(time_sel),length(listPD));
    cnt = 0;
    for ii = 1:length(listPD)
        cnt = cnt + 1;
        load([listPD(ii).folder '/' listPD(ii).name]);
        dd1(:,cnt) = d.d;
        ddpatt(:,:,cnt) = d.pattern;
        load([listGT(ii).folder '/' listGT(ii).name])
        ddTG(:,:,cnt) = d.d;
        disp(ii)
    end
    %average over subjects of pairwise decoding accuracy time-series
    dd12 = mean(dd1,2);
    dd12std = std(dd1,0,2)./sqrt(size(dd1,2)); %standard errors
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
    plot(time_sel,dd12,'k','linewidth',2,'DisplayName',['M_vs_NewT' num2str(T_cond)]);
    hold on
    plot(time_sel,dd12 + dd12std,':','color','k','linewidth',0.5);
    hold on
    plot(time_sel,dd12 - dd12std,':','color','k','linewidth',0.5);
    xlabel('Time (s)');ylabel('Decoding accuracy (%)');
    hold on
    % Show significant time window
    patch_color = [.85 .85 .85]; % Color of grey box marking time range
    ylims = get(gca,'YLim');
    sgf = CCt.PixelIdxList;
    for ii = 1:length(sgf)
        dumbum = sgf{ii};
        sgf2 = time_sel(dumbum);
%         patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.7)
        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[40 40 70 70],patch_color,'EdgeColor','none','FaceAlpha',.7)
        hold on
    end
    plot(time_sel,dd12,'k','linewidth',2,'DisplayName',['M_vs_NewT' num2str(T_cond)]);
    hold on
    plot(time_sel,dd12 + dd12std,':','color','k','linewidth',0.5);
    hold on
    plot(time_sel,dd12 - dd12std,':','color','k','linewidth',0.5);
    xlabel('Time (s)');ylabel('Decoding accuracy (%)'); %(slightly stupid..) trick to get the time-series in front..
    set(gcf,'Color','w')
    legend('show')
    xlim([time_sel(1) time_sel(end)])
    ylim([40 70]);
    if topoplot_l == 1
        %plotting patterns (derived from weights) for pairwise decoding
        warning('have you selected the proper time-window for topoplot..?')
        signt = [0.5 2.0]; %significant time-point(s) you want to plot
        %extracting magnetometers for plotting purposes
        dp = mean(ddpatt,3); %average over subjects
        avg = dp(1:3:end,:); %extracting magnetometers only
        %plotting topoplot (fieldtrip function)
        %creating the mask for the data
        fieldtrip_example = load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/fieldmask.mat');
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
        cfg.layout = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/neuromag306mag.lay';
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
    end
    %plotting temoral generalization
    ddTGm = mean(ddTG,3);
    figure; imagesc(time_sel,time_sel,ddTGm); xlabel('Train time (s)'); ylabel('Test time (s)'); set(gca,'YDir','normal');
    colorbar
    set(gcf,'Color','w')
    %reshaping temporal generalization decoding acuracy for submitting it to statistical testing
    ddTG2 = ddTG(1:776,1:776,:) - 50; %subtracting 50 since the function tests significance against 0 (here represented by 50% chance level)
    dat = cell(1,size(ddTG2,3));
    for ii = 1:size(ddTG2,3)
        dat(1,ii) = {ddTG2(:,:,ii)};
    end
    if tempgen_l == 1
        %computing statistics and plotting results
        stat = pl_permtest(dat,'alpha',0.05); %actual function
        STAT{T_cond} = stat; %storing statistics to compare images later and decide the limits for the colorbars
        time_selj = time_sel(1:776);
        figure; imagesc(time_selj,time_selj,stat.statmap); xlabel('Train time (s)'); ylabel('Test time (s)'); set(gca,'YDir','normal');
        colorbar
        set(gcf,'Color','w')
        caxis([-5 10])

        %testing contrast results (correcting for multiple comparison) by performing Monte Carlo simulations
        P2 = zeros(776,776);
        statt = stat.statmap;
        P2(statt>2.6) = 1; %threshold t-val = 2.6 corresponding to p-val < 0.01 (obtained by dividing 0.05 by the 4 comparisons employed here)
        thresh = 0;
        permut = 1000;
        threshMC = 0.001;
        perm_max = 1;
        t1 = time_sel; t2 = t1;
        
        [ OUT ] = twoD_MCS_LBPD_D( P2, thresh, permut, threshMC, perm_max, t1 , t2 )
        
        outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Statistics';
        PDn = cell2table(OUT); %table
        writetable(PDn,[outdir '/TemporalGeneralization_OldvsNewT' num2str(T_cond) 'Testing_Training.xlsx'],'Sheet',1); %printing excel file
    end
else %plotting together time series of decoding accuracy for all 4 decoding analyses (i.e. Old vs NewT1, Old vs NewT2, Old vs NewT3, Old vs NewT4)
    figure
    %     COL{1} = 'b'; COL{2} = 'r'; COL{3} = 'k'; COL{4} = 'g'; COL{5} = 'c'; COL{6} = 'm'; %colors (TO BE MADE THE SAME AS FOR WAVEFORMS OF DIFFERENT CONDITIONS!!)
    color_line = colormap(lines(5)); %extracting some colours from a colormap
    color_line2 = color_line;
    color_line2(1,:) = color_line(2,:);
    color_line2(2,:) = color_line(1,:);
    color_line2(5,:) = [0.4 0.4 0.4];
    for nn = 1:4 %over NewTX conditions
        datadir = [datadirdum '/Old_vs_NewT' num2str(nn)];
        listPD = dir([datadir '/PD*']);
        listGT = dir([datadir '/TG*']);
        %loading data
        dd1 = zeros(length(time_sel),length(listPD));
        ddpatt = zeros(306,length(time_sel),length(listPD));
        ddTG = zeros(length(time_sel),length(time_sel),length(listPD));
        cnt = 0;
        for ii = 1:length(listPD)
            cnt = cnt + 1;
            load([listPD(ii).folder '/' listPD(ii).name]);
            dd1(:,cnt) = d.d;
            ddpatt(:,:,cnt) = d.pattern;
            load([listGT(ii).folder '/' listGT(ii).name])
            ddTG(:,:,cnt) = d.d;
            disp([num2str(nn) ' - ' num2str(ii)])
        end
        %average over subjects of pairwise decoding accuracy time-series
        dd12 = mean(dd1,2);
        dd12std = std(dd1,0,2)./sqrt(size(dd1,2)); %standard errors
        %reshaping pairwise decoding acuracy time-series for submitting them to statistical testing
        dd2 = dd1 - 50; %subtracting 50 since the function tests significance against 0 (here represented by 50% chance level)
        dat = cell(1,size(dd2,2));
        for ii = 1:size(dd2,2)
            dat(1,ii) = {dd2(:,ii)'};
        end
        %plotting average over subjects (pairwise decoding accuracy)
        grid on;grid minor;
        plot(time_sel,dd12,'color',color_line2(nn+1,:),'linewidth',2,'DisplayName',['NewT' num2str(nn)]);
        hold on
        plot(time_sel,dd12 + dd12std,':','color',color_line2(nn+1,:),'linewidth',0.5,'DisplayName',['NewT' num2str(nn)]);
        hold on
        plot(time_sel,dd12 - dd12std,':','color',color_line2(nn+1,:),'linewidth',0.5,'DisplayName',['NewT' num2str(nn)]);
        xlabel('Time (s)');ylabel('Decoding accuracy (%)');
        hold on
    end
    legend('show')
    set(gcf,'Color','w')
    xlim([time_sel(1) time_sel(end)])
end
% %trick to export results in excel sheet.. I initialized dumbo in the command window
% dumbo(:,T_cond) = stat.FDR.criticalmap;
% dumbo = dumbo';
% dumbo2 = zeros(5,1126);
% dumbo2(1,1:1126) = time_sel;
% dumbo2(2:5,1:1126) = dumbo;
% outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Statistics';
% xlswrite([outdir '/PairwiseDecoding.xlsx'],dumbo2)

%%

%% *** SOURCE RECONSTRUCTION (LBPD) ***

%%

%% CREATING 8mm PARCELLATION FOR EASIER INSPECTION IN FSLEYES
%OBS!! This section is done only for better handling of some visualization purposes, but it does not affect any of the beamforming algorithm;
% it is just important not to mix up the MNI coordinates, thus I would recommend to use the following lines

%1) USE load_nii TO LOAD A PREVIOUS NIFTI IMAGE
imag_8mm = load_nii('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_T1_8mm_brain.nii.gz');
Minfo = size(imag_8mm.img); %get info about the size of the original image
M8 = zeros(Minfo(1), Minfo(2), Minfo(3)); %Initialize an empty matrix with the same dimensions as the original .nii image
cc = 0; %set a counter
M1 = imag_8mm.img;
for ii = 1:Minfo(1) %loop across each voxel of every dimension
    for jj = 1:Minfo(2)
        for zz = 1:Minfo(3)
            if M1(ii,jj,zz) ~= 0 %if we have an actual brain voxel
                cc = cc+1;
                M8(ii,jj,zz) = cc;
            end
        end
    end
end
%2) PUT YOUR MATRIX IN THE FIELD ".img"
imag_8mm.img = M8; %assign values to new matrix 
%3) SAVE NIFTI IMAGE USING save_nii
save_nii(imag_8mm,'/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_brain_diy.nii.gz');
%4) USE FSLEYES TO LOOK AT THE FIGURE
%Create parcellation on the 8mm template
for ii = 1:3559 %for each 8mm voxel
    cmd = ['fslmaths /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_brain_diy.nii.nii.gz -thr ' num2str(ii) ' -uthr ' num2str(ii) ' -bin /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/AAL_80mm_3559ROIs/' num2str(ii) '.nii.gz'];
    system(cmd)
    disp(ii)
end
%5) GET MNI COORDINATES OF THE NEW FIGURE AND SAVE THEM ON DISK
MNI8 = zeros(3559,3);
for mm = 1:3559 %over brain voxel
    path_8mm = ['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/' num2str(mm) '.nii.gz']; %path for each of the 3559 parcels
    [mni_coord,pkfo] = osl_mnimask2mnicoords(path_8mm);  %getting MNI coordinates
    MNI8(mm,:) = mni_coord; %storing MNI coordinates
end
%saving on disk
save('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_coord_dyi.mat', 'MNI8');

%% CONVERSION T1 - DICOM TO NIFTI

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/dicm2nii'); %adds path to the dcm2nii folder in osl
MRIsubj = dir('/projects/MINDLAB2020_MEG-AuditoryPatternRecognition/raw/0*');
MRIoutput = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/MRI_nifti';
MRIout_block{1} = 'Block_3'; MRIout_block{2} = 'Block_4'; MRIout_block{3} = 'Block_5'; MRIout_block{4} = 'Block_6'; MRIout_block{5} = 'Block_7';

for bb = 1:5 %over experimental blocks
    for ii = 1:length(MRIsubj) %over subjects
        asd = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/MRI_nifti/Block_' num2str(bb+2) '/' MRIsubj(ii).name];
        if ~exist(asd,'dir') %checking whether the directory exists
            mkdir(asd); %if not, creating it
        end
        if isempty(dir([asd '/*.nii'])) %if there are no nifti images.. I need to convert them
            flagg = 0;
            MRIMEGdate = dir([MRIsubj(ii).folder '/' MRIsubj(ii).name '/20*']);
            niiFolder = [MRIoutput '/' MRIout_block{bb} '/' MRIsubj(ii).name];
            for jj = 1:length(MRIMEGdate) %over dates of recording
                if ~isempty(dir([MRIMEGdate(jj).folder '/' MRIMEGdate(jj).name '/MR*'])) %if we get an MRI recording
                    MRI2 = dir([MRIMEGdate(jj).folder '/' MRIMEGdate(jj).name '/MR/*fatsat']); %looking for T1
                    if ~isempty(MRI2) %if we have it
                        flagg = 1; %determining that I could convert MRI T1
                        dcmSource = [MRI2(1).folder '/' MRI2(1).name '/files/'];
                        if ii ~= 68 || jj ~= 3 %this is because subject 0068 got two MRIs stored.. but the second one (indexed by jj = 3) is of another subject (0086); in this moment, subject 0086 is perfectly fine, but in subject 0068 there are still the two MRIs (for 0068 (jj = 2) and for 0086 (jj = 3))
                            dicm2nii(dcmSource, niiFolder, '.nii');
                        end
                    end
                end
            end
            if flagg == 0
                warning(['subject ' MRIsubj(ii).name ' has no MRI T1']);
            end
        end
        disp(ii)
    end
end

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% RHINO coregistration

%block to be run RHINO coregistrartion on
block = 6; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd

if block == 3
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 4
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogspeed*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 5
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*samemel*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 6
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*visualpat*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 7
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*pd*_tsssdsm.mat'); %dir to epoched files (encoding)
end

%running rhino
%OBS! check that all MEG data are in the same order and number as MRI nifti files!
a = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/MRI_nifti/Block_' num2str(block)]; %set path to MRI subjects' folders
for ii = 35%1:length(list) %OBS! change this depending on atonal vs. major
    S = [];
    S.ii = ii;
    S.D = [list(ii).folder '/' list(ii).name]; %path to major files
    D = spm_eeg_load(S.D);
    if ~isfield(D,'inv') %checking if the coregistration was already run
        dummyname = D.fname;
        if 7 == exist([a '/' dummyname(13:16)],'dir') %if you have the MRI folder
            dummymri = dir([a '/' dummyname(13:16) '/*.nii']); %path to nifti files (ending with .nii)
            if ~isempty(dummymri)
                S.mri = [dummymri(1).folder '/' dummymri(1).name];
                %standard parameters
                S.useheadshape = 1;
                S.use_rhino = 1; %set 1 for rhino, 0 for no rhino
                %         S.forward_meg = 'MEG Local Spheres';
                S.forward_meg = 'Single Shell'; %CHECK WHY IT SEEMS TO WORK ONLY WITH SINGLE SHELL!!
                S.fid.label.nasion = 'Nasion';
                S.fid.label.lpa = 'LPA';
                S.fid.label.rpa = 'RPA';
                jobid = job2cluster(@coregfunc,S); %running with parallel computing
            else
                warning(['subject ' dummyname(13:16) ' does not have the MRI'])
            end
        end
    else
        if isempty(D.inv{1}) %checking whether the coregistration was run but now it is empty..
            warning(['subject ' D.fname ' has an empty rhino..']);
        end
    end
    disp(ii)
end

%% checking (or copying) RHINO

copy_label = 0; % 1 = pasting inv RHINO from epoched data (where it was computed) to continuous data; 0 = simply showing RHINO coregistration
%block to be run RHINO coregistration on
block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd

if block == 3
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 4
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogspeed*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 5
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*samemel*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 6
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*visualpat*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 7
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*pd*_tsssdsm.mat'); %dir to epoched files (encoding)
end

for ii = 1:length(list)
    D = spm_eeg_load([list(ii).folder '/' list(ii).name]);
    if isfield(D,'inv')
        if copy_label == 0 %simply displaying RHINO coregistration
            if isfield(D,'inv') %checking if the coregistration was already run
                rhino_display(D)
            end
        else %pasting inv RHINO from epoched data (where it was computed) to continuous data
            inv_rhino = D.inv;
            D2 = spm_eeg_load([list(ii).folder '/' list(ii).name(2:end)]);
            D2.inv = inv_rhino;
            D2.save();
        end
    end
    disp(['Block ' num2str(block) ' - Subject ' num2str(ii)])
end

%%

%% BEAMFORMING

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% FUNCTION FOR SOURCE RECONSTRUCTION

%user settings
clust_l = 1; %1 = using cluster of computers (CFIN-MIB, Aarhus University); 0 = running locally
timek = 1:1026; %time-points
freqq = []; %frequency range (empty [] for broad band)
% freqq = [0.1 1]; %frequency range (empty [] for broad band)
% freqq = [2 8]; %frequency range (empty [] for broad band)
sensl = 1; %1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers (SUGGESTED 1!)
workingdir2 = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD'; %high-order working directory (a subfolder for each analysis with information about frequency, time and absolute value will be created)
block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd
invers = 1; %1-4 = different ways (e.g. mean, t-values, etc.) to aggregate trials and then source reconstruct only one trial; 5 for single trial independent source reconstruction

if isempty(freqq)
    absl = 0; % 1 = absolute value of sources; 0 = not
else
    absl = 0;
end
%actual computation
%list of subjects with coregistration (RHINO - OSL/FSL) - epoched
if block == 3
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Old_Correct','New_T1_Correct','New_T2_Correct','New_T3_Correct','New_T4_Correct'};
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
elseif block == 4
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogspeed*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Old_Fast_Correct','Old_Slow_Correct','New_Fast_Correct','New_Slow_Correct'};
elseif block == 5
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*samemel*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Encoding','Old_Correct','New_Correct'};
elseif block == 6
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*visualpat*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Encoding','Old_Correct','New_Correct'};
elseif block == 7
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*pd*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'pd','pm'};
end
if isempty(freqq)
    workingdir = [workingdir2 '/Block_' num2str(block) '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_broadband_invers_' num2str(invers)];
else
    workingdir = [workingdir2 '/Block_' num2str(block) '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_' num2str(freqq(1)) '_' num2str(freqq(2)) '_invers_' num2str(invers)];
end
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');
if ~exist(workingdir,'dir') %creating working folder if it does not exist
    mkdir(workingdir)
end
for ii = 1:length(list) %over subjects
    S = [];
    if ~isempty(freqq) %if you want to apply the bandpass filter, you need to provide continuous data
        %             disp(['copying continuous data for subj ' num2str(ii)])
        %thus pasting it here
        %             copyfile([list_c(ii).folder '/' list_c(ii).name],[workingdir '/' list_c(ii).name]); %.mat file
        %             copyfile([list_c(ii).folder '/' list_c(ii).name(1:end-3) 'dat'],[workingdir '/' list_c(ii).name(1:end-3) 'dat']); %.dat file
        %and assigning the path to the structure S
        S.norm_megsensors.MEGdata_c = [list(ii).folder '/' list(ii).name(2:end)];
    end
    %copy-pasting epoched files
    %         disp(['copying epoched data for subj ' num2str(ii)])
    %         copyfile([list(ii).folder '/' list(ii).name],[workingdir '/' list(ii).name]); %.mat file
    %         copyfile([list(ii).folder '/' list(ii).name(1:end-3) 'dat'],[workingdir '/' list(ii).name(1:end-3) 'dat']); %.dat file
    
    S.Aarhus_cluster = clust_l; %1 for parallel computing; 0 for local computation
    
    S.norm_megsensors.zscorel_cov = 1; % 1 for zscore normalization; 0 otherwise
    S.norm_megsensors.workdir = workingdir;
    S.norm_megsensors.MEGdata_e = [list(ii).folder '/' list(ii).name];
    S.norm_megsensors.freq = freqq; %frequency range
    S.norm_megsensors.forward = 'Single Shell'; %forward solution (for now better to stick to 'Single Shell')
    
    S.beamfilters.sensl = sensl; %1 = magnetometers; 2 = gradiometers; 3 = both MEG sensors (mag and grad) (SUGGESTED 3!)
    S.beamfilters.maskfname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz'; % path to brain mask: (e.g. 8mm MNI152-T1: '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz')
    
    S.inversion.znorml = 0; % 1 for inverting MEG data using the zscored normalized one; (SUGGESTED 0 IN BOTH CASES!)
    %                                 0 to normalize the original data with respect to maximum and minimum of the experimental conditions if you have both magnetometers and gradiometers.
    %                                 0 to use original data in the inversion if you have only mag or grad (while e.g. you may have used zscored-data for covariance matrix)
    %
    S.inversion.timef = timek; %data-points to be extracted (e.g. 1:300); leave it empty [] for working on the full length of the epoch
    S.inversion.conditions = condss; %cell with characters for the labels of the experimental conditions (e.g. {'Old_Correct','New_Correct'})
    S.inversion.bc = [1 26]; %extreme time-samples for baseline correction (leave empty [] if you do not want to apply it)
    S.inversion.abs = absl; %1 for absolute values of sources time-series (recommendnded 1!)
    S.inversion.effects = invers;
    
    S.smoothing.spatsmootl = 0; %1 for spatial smoothing; 0 otherwise
    S.smoothing.spat_fwhm = 100; %spatial smoothing fwhm (suggested = 100)
    S.smoothing.tempsmootl = 0; %1 for temporal smoothing; 0 otherwise
    S.smoothing.temp_param = 0.01; %temporal smoothing parameter (suggested = 0.01)
    S.smoothing.tempplot = [1 2030 3269]; %vector with sources indices to be plotted (original vs temporally smoothed timeseries; e.g. [1 2030 3269]). Leave empty [] for not having any plot.
    
    S.nifti = 1; %1 for plotting nifti images of the reconstructed sources of the experimental conditions
    S.out_name = ['SUBJ_' list(ii).name(13:16)]; %name (character) for output nifti images (conditions name is automatically detected and added)
    
    if clust_l ~= 1 %useful  mainly for begugging purposes
        MEG_SR_Beam_LBPD(S);
    else
        jobid = job2cluster(@MEG_SR_Beam_LBPD,S); %running with parallel computing
    end
end

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 0); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 3); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% STATISTICS OVER PARTICIPANTS (USING PARALLEL COMPUTING, PARTICULARLY USEFUL IF YOU HAVE SEVERAL CONTRASTS)

block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd
clust = 1; % 1 = using Aarhus cluster (parallel computing); 0 = run locally
analys_n = 2; %analysis number (in the list indexed below)

%building structure
asd = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_' num2str(block) '/Beam*']);
S = [];
% dumm = dir([asd(analys_n).folder '/' asd(analys_n).name '/SUBJ*mat']);
% S.list = dumm(1:10);
S.workingdir = [asd(analys_n).folder '/' asd(analys_n).name '/test']; %path where the data from MEG_SR_Beam_LBPD.m is stored
S.sensl = 1; % 1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers.
S.plot_nifti = 1; %1 to plot nifti images; 0 otherwise
S.plot_nifti_name = []; %character with name for nifti files (it may be useful if you run separate analysis); Leave empty [] to not  specify any name
% S.contrast = [1 0 0 0 0 0 -1; 0 1 0 0 0 0 -1; 0 0 1 0 0 0 -1; 0 0 0 1 0 0 -1; 0 0 0 0 1 0 -1; 0 0 0 0 0 1 -1; 1 1 1 1 1 1 -1]; %one row per contrast (e.g. having 3 conditions, [1 -1 0; 1 -1 -1; 0 1 -1]; two or more -1 or 1 are interpreted as the mean over them first and then the contrast. Leave empty [] for no contrasts. 
if block == 3
    S.contrast = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
elseif block == 4
        S.contrast = [1 0 -1 0; 0 1 0 -1];
elseif block == 5 || block == 6
    S.contrast = [0 1 -1; -1 1 0];
else
    S.contrast = [1 -1];
end
S.effects = 1; %mean over subjects for now
if clust == 1
    S.Aarhus_clust = 1; %1 to use paralle computing (Aarhus University, contact me, Leonardo Bonetti, for more information; leonardo.bonetti@clin.au.dk)
    %actual function
    jobid = job2cluster(@MEG_SR_Stats1_Fast_LBPD,S); %running with parallel computing
else
    S.Aarhus_clust = 0; %1 to use paralle computing (Aarhus University, contact me, Leonardo Bonetti, for more information; leonardo.bonetti@clin.au.dk)
    MEG_SR_Stats1_Fast_LBPD(S)
end

%% FIGURE - MEG0211 AND SOURCES OF ITS PEAK ACTIVITY (AND DIFFERENCES ACROSS EXPERIMENTAL CONDITIONS)

% 1) Based on decoding and inspection of the prototypical MEG0211 channel, we defined the following time-windows (and contrasts between conditions) and reconsturcted their brain sources:
%    -1) 0.5-0.7 Old vs NewT1 (0.5-0.6)
%    -2) 0.7-0.8 Old vs NewT1
%    -3) 0.95-1.05 Old vs NewT2 (0.98-1.02)
%    -4) 1.1-1.3 Old vs NewT1 (1.05-1.15)
%    -5) 1.1-1.3 Old vs NewT2 (1.05-1.15)
%    -6) 1.3-1.4 Old vs NewT3 (1.33-1.39)
%    -7) 1.45-1.55 Old vs NewT1
%    -8) 1.45-1.55 Old vs NewT2
%    -9) 1.45-1.55 Old vs NewT3
%    -10) 1.7-1.8 Old vs NewT4 (1.7-1.75)
%    -11) 1.8-1.9 Old vs NewT1 (1.75-1.85)
%    -12) 1.8-1.9 Old vs NewT2 (1.75-1.85)
%    -13) 1.8-1.9 Old vs NewT3 (1.75-1.85)
%    -14) 1.8-1.9 Old vs NewT4 (1.75-1.85)

clear names
names{1} = 'O_NT1_05_06'; names{2} = 'O_NT1_07_08'; names{3} = 'O_NT2_098_102'; names{4} = 'O_NT1_105_115'; names{5} = 'O_NT2_105_115';
names{6} = 'O_NT3_133_139'; names{7} = 'O_NT1_145_155'; names{8} = 'O_NT2_145_155'; names{9} = 'O_NT3_145_155';
names{10} = 'O_NT4_17_175'; names{11} = 'O_NT1_175_85'; names{12} = 'O_NT2_175_185'; names{13} = 'O_NT3_175_185'; names{14} = 'O_NT4_175-185';
% vecttime = [151 201; 202 226; 264 289; 301 351; 301 351; 352 376; 389 414; 389 414; 389 414; 451 476; 477 501; 477 501; 477 501; 477 501]; %time-windows, as reported above (here in time samples instead of seconds)
vecttime = [172 190; 202 226; 271 281; 289 314; 289 314; 359 375; 389 414; 389 414; 389 414; 451 464; 465 489; 465 489; 465 489; 465 489]; %time-windows, as reported above (here in time samples instead of seconds)

vectcond = [2,2,3,2,3,4,2,3,4,5,2,3,4,5]; %condition X to be contrasted against condition 1 (Old)
list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/SUBJ*.mat']);
data = zeros(3559,length(vecttime),5,length(list)); %sources x time-windows of intereste x conditions x subjects
for ii = 1:length(list) %over subjects
    load([list(ii).folder '/' list(ii).name])
    for jj = 1:length(vecttime) %over time-windows
        data(:,jj,:,ii) = mean(OUT.sources_ERFs(:,vecttime(jj,:),:),2); %averge over the selected jj time-window
    end
    disp(ii)
end
%computing t-tests
P = zeros(3559,length(vecttime));
T = zeros(3559,length(vecttime));
for ii = 1:size(P,1) %over sources
    for jj = 1:length(vecttime) %over time-windows
        [~,p,~,stats] = ttest(squeeze(data(ii,jj,1,:)),squeeze(data(ii,jj,vectcond(jj),:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
        P(ii,jj) = p;
        T(ii,jj) = stats.tstat;
    end
    disp(ii)
end
%creating nifti images of computed ROIs
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211';
mkdir(outdir)
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for iii = 1%:size(T,2)
    fnamenii = [outdir '/' names{iii} '.nii.gz']; %path and name of the image to be saved
    SO = T(:,iii);
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    for ii = 1:size(SO,1) %over brain sources
        dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = SO(ii,:); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - ROI ' num2str(iii)])
    save_nii(nii,fnamenii); %printing image
end

%cluster-based MCS to clean things up (either by removing random voxels automatically or by selecting only the biggest cluster)
%Cluster permutations test
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/*gz');
%OBS! you may get a warning since the skeletonized image is not exactly in MNI space, but close enough
[ mni_coords, xform ] = osl_mnimask2mnicoords('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');
MASK = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');
%preparing function for MCS
for ii = 1%:length(list) %over time-windows
    %getting MNI coordinates
    T = load_nii([list(ii).folder '/' list(ii).name]);
    %extracting matrix with statistics
    T2 = T.img(:,:,:); %extracting time-window ii
    %mask for non-0 voxels in brain imges (basically building a layout for actual brain voxels)
    mask = zeros(size(T2,1),size(T2,2),size(T2,3));
    data = mask;
    data(abs(T2)>3) = 1; %binarizing data (t-value threshold = absolute value of 3)
    mask(MASK.img~=0) = 1; %assigning 1 when you have real brain voxels
    %preparation of information and data for the actual function
    S = [];
    S.T = T2;
    S.outdir = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS']; %output path
    S.parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz';
    S.labels = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat';
    S.MNIcoords = mni_coords; %MNI coordinates of 8mm MNI152T1 brain
    S.data = data;
    S.mask = mask; %mask of the brain layout you have your results in
    %THE NUMBERS FOR THE PERMUTATIONS HAVE BEEN CHANGED IN DIFFERENT FOR THE LAST TWO CASES SINCE THERE WERE A LOT OF SIGNIFICANT VALUES, SO IT WAS NECESSARY EITHER TO USE A STRICTED THRESHOLD FOR THE T-VALUES
    %OR RELAZING THE THRESHOLD FOR THE MCS, THIS IS THE INTRINSIC LIMITATION OF THE MCS APPROACH.. BUT THIS DOES NOT CHANGE THE SIGNIFICANCE OF THE RESULTS IN ANY WAY..
    S.permut = 1000; %number of permutations for Monte Carlo simulation
    S.clustmax = 0; %set 1 for only max cluster size of each permutation MCS (more strict); set 0 for every size of each cluster detected for each permutation MCS (less strict).
    S.permthresh = 0.001; %threshold for MCS
    %UNTIL HERE
    S.anal_name = [list(ii).name(1:end-7) '_MCS']; %name for the analysis (used to identify and save image and results)
    
    %actual function
    PP = BrainSources_MonteCarlosim_3D_LBPD_D(S);
    
    disp(ii)
end

%Combining images with significant clusters
%here you need to combine images with more than one cluster (5 final images)
%general path
path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/ToBeCombined';
%figure 1
%cluster 1
name1 = ['O_NT2_098_102_MCS_SignClust_1_Tvals.nii.gz'];
%cluster 2
name2 = ['O_NT2_098_102_MCS_SignClust_2_Tvals.nii.gz'];
%output name
output = ['O_NT2_098_102.nii.gz'];
%command line and running
cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' ' path  '/' output]; %OBS! be careful with the spacing
system(cmd)
%figure 2
%cluster 1
name1 = ['O_NT2_175_185_MCS_SignClust_1_Tvals.nii.gz'];
%cluster 2
name2 = ['O_NT2_175_185_MCS_SignClust_2_Tvals.nii.gz'];
%output name
output = ['O_NT2_175_185.nii.gz'];
%command line and running
cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' ' path  '/' output]; %OBS! be careful with the spacing
system(cmd)
%figure 3
%cluster 1
name1 = ['O_NT3_145_155_MCS_SignClust_1_Tvals.nii.gz'];
%cluster 2
name2 = ['O_NT3_145_155_MCS_SignClust_3_Tvals.nii.gz'];
%output name
output = ['O_NT3_145_155.nii.gz'];
%command line and running
cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' ' path  '/' output]; %OBS! be careful with the spacing
system(cmd)
%figure 4
%cluster 1
name1 = ['O_NT3_175_185_MCS_SignClust_1_Tvals.nii.gz'];
%cluster 2
name2 = ['O_NT3_175_185_MCS_SignClust_2_Tvals.nii.gz'];
%output name
output = ['O_NT3_175_185.nii.gz'];
%command line and running
cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' ' path  '/' output]; %OBS! be careful with the spacing
system(cmd)
%figure 5
%cluster 1
name1 = ['O_NT4_175-185_MCS_SignClust_1_Tvals.nii.gz'];
%cluster 2
name2 = ['O_NT4_175-185_MCS_SignClust_2_Tvals.nii.gz'];
%output name
output = ['O_NT4_175_185.nii.gz'];
%command line and running
cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' ' path  '/' output]; %OBS! be careful with the spacing
system(cmd)

%after this, we created the final images using Workbench

%% Extracting information about the clusters at source level and reporting it in xlsx files

%here we obtain information about the brain regions forming the clusters
%the tables can be found in SUPPLEMENTARY MATERIALS

path = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench';
list = dir([path '/*gz']);

for ii = 1:length(list)
    fname = [path '/' list(ii).name]; %tone x cluster 1
    S = [];
    S.fname = fname;
    %actual function
    PDn = Extract_BrainCluster_Information_3D_LBPD_D(S);
    writetable(PDn,[path '/' list(ii).name(1:end-7) '.xlsx'],'Sheet',1); %printing excel file
end

%%

%% *** DEFINING FUNCTIONAL PARCELS ***

%%% Originally, I tried the functional k-means clustering, then I decided
%%% to proceed with a simpler solution, as described below

%%

%% *** SIMPLER SOLUTION TO GET THE FUNCTIONAL PARCELS (ROIs) ***

%%% 1) Getting key time-points of different activity between Old and New (and a short time-window around them).
%%%    This were based on the results obtained from the decdding and from the inspection of the prototypical MEG0211 channel.
%%% 2) Thresholding them (e.g. with a high t-value of 3 or even higher)
%%% 3) Running cluster-based MCS to remove spurious brain voxels who were randomly involved in steps 1) and 2)
%%% 4) In some cases, when the activity was left or right lateralized, getting the mirror ROI in the opposite hemisphere 

%%

%% ACTUAL COMPUTATION

% 1) Getting peaks
%%% The peaks of the different activity between Old and New corresponded to:
% - slow peak after the second (and third and fourth) tone of the sequence
% - slow peak after the last tone
% - later prediction error peak (NewT2) - hippocampus and VMPFC
% - first prediction error peak (NewT1) - auditory cortex

%0)output directory
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity';

%1)computing ROIs
%loading data
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/sources_contrast_1.mat');
Vstat2 = t_val_s;
names{1} = 'MC'; names{2} = 'HITR'; names{3} = 'HITL'; names{4} = 'VMPFC'; names{5} = 'ACL'; names{6} = 'ACR';
%1-MC
%defining peaks (using time-windows of + and - 5 time-points (+ and - 20ms)) and getting voxels above threshold in those peaks, independently for each time-point (taking all voxels with at least one significant time-point)
MC = zeros(size(Vstat2,1),11); %preparing dummy variable
MC(Vstat2(:,200:210) > 3) = 1; %getting ones when values are above the threshold (t-value = 3 here)
MC = sum(MC,2); %summing voxels to get how many times they were significant
MC(MC~=0) = 1; %if a voxel was significant at least one time, we take it
ROIs{1} = MC;
%2-HITR
%defining peaks (using time-windows of + and - 5 time-points (+ and - 20ms)) and getting voxels above threshold in those peaks, independently for each time-point (taking all voxels with at least one significant time-point)
HITR = zeros(size(Vstat2,1),11); %preparing dummy variable
HITR(Vstat2(:,483:493) > 3.2) = 1; %getting ones when values are above the threshold (t-value = 3 here)
HITR = sum(HITR,2); %summing voxels to get how many times they were significant
HITR(HITR~=0) = 1; %if a voxel was significant at least one time, we take it
ROIs{2} = HITR;
%5-ACL
%defining peaks (using time-windows of + and - 5 time-points (+ and - 20ms)) and getting voxels above threshold in those peaks, independently for each time-point (taking all voxels with at least one significant time-point)
ACL = zeros(size(Vstat2,1),11); %preparing dummy variable
ACL(Vstat2(:,159:169) < -3) = 1; %getting ones when values are above the threshold (t-value = 3 here)
ACL = sum(ACL,2); %summing voxels to get how many times they were significant
ACL(ACL~=0) = 1; %if a voxel was significant at least one time, we take it
ROIs{5} = ACL;
%4-VMPFC
%loading data
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/sources_contrast_2.mat');
Vstat2 = t_val_s;
%defining peaks (using time-windows of + and - 5 time-points (+ and - 20ms)) and getting voxels above threshold in those peaks, independently for each time-point (taking all voxels with at least one significant time-point)
VMPFC = zeros(size(Vstat2,1),11); %preparing dummy variable
VMPFC(Vstat2(:,269:279) < -3.7) = 1; %getting ones when values are above the threshold (t-value = 3 here)
VMPFC = sum(VMPFC,2); %summing voxels to get how many times they were significant
VMPFC(VMPFC~=0) = 1; %if a voxel was significant at least one time, we take it
ROIs{4} = VMPFC;

%2)creating nifti images of computed ROIs
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for iii = 1:length(ROIs)
    fnamenii = [outdir '/' names{iii} '.nii.gz']; %path and name of the image to be saved
    SO = ROIs{iii};
    if ~isempty(SO)
        %building nifti image
        SS = size(maskk.img);
        dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
        for ii = 1:size(SO,1) %over brain sources
            dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
            [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
            dumimg(i1,i2,i3,:) = SO(ii,:); %storing values for all time-points in the image matrix
        end
        nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
        nii.img = dumimg; %storing matrix within image structure
        nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
        disp(['saving nifti image - ROI ' num2str(iii)])
        save_nii(nii,fnamenii); %printing image
    end
end

%3)cluster-based MCS to clean things up (either by removing random voxels automatically or by selecting only the biggest cluster)
%Cluster permutations test
%loading main image 
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/*gz');
%OBS! you may get a warning since the skeletonized image is not exactly in MNI space, but close enough
[ mni_coords, xform ] = osl_mnimask2mnicoords('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');
MASK = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');
%preparing function for MCS
for ii = 1:length(list) %over time-windows
    %getting MNI coordinates
    T = load_nii([list(ii).folder '/' list(ii).name]);
    %extracting matrix with statistics
    T2 = T.img(:,:,:); %extracting time-window ii
    %mask for non-0 voxels in brain imges (basically building a layout for actual brain voxels)
    mask = zeros(size(T2,1),size(T2,2),size(T2,3));
    data = T2;
    mask(MASK.img~=0) = 1; %assigning 1 when you have real brain voxels
    %preparation of information and data for the actual function
    S = [];
    S.T = T2;
    S.outdir = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/MCS']; %output path
    S.parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz';
    S.labels = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat';
    S.MNIcoords = mni_coords; %MNI coordinates of 8mm MNI152T1 brain
    S.data = data;
    S.mask = mask; %mask of the brain layout you have your results in
    %THE NUMBERS FOR THE PERMUTATIONS HAVE BEEN CHANGED IN DIFFERENT BRAIN AREAS.. DOESN'T REALLY MATTER SINCE HERE I JUST WANT TO GET THE PROPER ROIs AND NOT DOING ANY STATISTICAL TEST
    S.permut = 1000; %number of permutations for Monte Carlo simulation
    S.clustmax = 1; %set 1 for only max cluster size of each permutation MCS (more strict); set 0 for every size of each cluster detected for each permutation MCS (less strict).
    S.permthresh = 0.001; %threshold for MCS
    %UNTIL HERE
    S.anal_name = [list(ii).name(1:end-7) '_MCS']; %name for the analysis (used to identify and save image and results)
    
    %actual function
    PP = BrainSources_MonteCarlosim_3D_LBPD_D(S);
    
    disp(ii)
end

%% GETTING FINAL ROIs (MAIN CLUSTER) AND ROIs MIRRORED IN THE OTHER HEMISPHERE WHEN THE ORIGINAL ROIs WERE LEFT- OR RIGHT-LATERALIZED

list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/MCS/*gz');
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_coord_dyi.mat'); %loading new cooridnate system..
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure

ROIs = zeros(3559,length(list)); %initializing ROIs
for mm = 1:length(list) %over ROIs
    roi = load_nii([list(mm).folder '/' list(mm).name]); %loading parcel (ROI) mm
    for ii = 1:size(maskk.img,1) %over dimension 1 of the 3D image
        for jj = 1:size(maskk.img,2) %over dimension 2 of the 3D image
            for zz = 1:size(maskk.img,3) %over dimension 3 of the 3D image
                if roi.img(ii,jj,zz) ~= 0 %if there is a voxel forming the ROI
                    ROIs(maskk.img(ii,jj,zz),mm) = 1; %assigning one to that voxel in ROIs (by taking number of voxel from the maskk.img)
                end
            end
        end
    end
end
%list of nifti images for the previously computed ROIs: 1)ACL; 2)HITR; 3)MC; 4)VMPFC
%3-HITL
ab = find(ROIs(:,2)==1);
HITL = zeros(3559,1);
for ii = 1:length(ab) %over voxels of HITR
    %[MIN,IDX] = min(abs(sum(MNI8-[MNI8(ab(ii),1)*(-1),MNI8(ab(ii),2),MNI8(ab(ii),3)],2)))
    dummy = zeros(3559,3);
    dummy((MNI8(ab(ii),1)*(-1)+4)==MNI8(:,1),1) = 1; %getting values in the other hemisphere (multiplying by -1 o reverse the sign and adding 4 since there is a very small (negligible) imperfection of the coordinate system (symmetric across hemispheres voxels are slightly misalligned.. nothing that really matters here in MEG)
    dummy((MNI8(ab(ii),2)==MNI8(:,2)),2) = 1;
    dummy((MNI8(ab(ii),3)==MNI8(:,3)),3) = 1;
    HITL((sum(dummy,2)==3),1) = 1;
end
%6-ACR
ab = find(ROIs(:,1)==1);
ACR = zeros(3559,1);
for ii = 1:length(ab) %over voxels of HITR
    %[MIN,IDX] = min(abs(sum(MNI8-[MNI8(ab(ii),1)*(-1),MNI8(ab(ii),2),MNI8(ab(ii),3)],2)))
    dummy = zeros(3559,3);
    dummy((MNI8(ab(ii),1)*(-1)+4)==MNI8(:,1),1) = 1; %getting values in the other hemisphere (multiplying by -1 o reverse the sign and adding 4 since there is a very small (negligible) imperfection of the coordinate system (symmetric across hemispheres voxels are slightly misalligned.. nothing that really matters here in MEG)
    dummy((MNI8(ab(ii),2)==MNI8(:,2)),2) = 1;
    dummy((MNI8(ab(ii),3)==MNI8(:,3)),3) = 1;
    ACR((sum(dummy,2)==3),1) = 1;
end

%reshaping everything in order
%order in ROIs: 1)ACL; 2)HITR; 3)MC; 4)VMPFC
%final order
names{1} = 'MC'; names{2} = 'HITR'; names{3} = 'HITL'; names{4} = 'VMPFC'; names{5} = 'ACL'; names{6} = 'ACR';
ROIs_DCM = zeros(3559,6);
ROIs_DCM(:,1) = ROIs(:,3); %MC
ROIs_DCM(:,2) = ROIs(:,2); %HITR
ROIs_DCM(:,3) = HITL; %HITL
ROIs_DCM(:,4) = ROIs(:,4); %VMPFC
ROIs_DCM(:,5) = ROIs(:,1); %ACL
ROIs_DCM(:,6) = ACR; %ACR
%saving on disk
save('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM','ROIs_DCM','names')

%% creating nifti images for the final ROIs

maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat');
for iii = 1:size(ROIs_DCM,2)
    fnamenii = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/MainROIs/' names{iii} '.nii.gz']; %path and name of the image to be saved
    SO = ROIs_DCM(:,iii);
    if ~isempty(SO)
        %building nifti image
        SS = size(maskk.img);
        dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
        for ii = 1:size(SO,1) %over brain sources
            dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
            [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
            dumimg(i1,i2,i3,:) = SO(ii,:); %storing values for all time-points in the image matrix
        end
        nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
        nii.img = dumimg; %storing matrix within image structure
        nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
        disp(['saving nifti image - condition ' num2str(iii)])
        save_nii(nii,fnamenii); %printing image
    end
end

%% FROM LBPD COORDINATES TO 3D NIFTI IMAGE
%same solution provided with a function

load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat');
S.data = ROIs_DCM; %data (voxels x ROIs (time-points))
S.fname = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/MainROIs/Test'; %path and name for the image to be saved
names{1} = 'MC'; names{2} = 'HITR'; names{3} = 'HITL'; names{4} = 'VMPFC'; names{5} = 'ACL'; names{6} = 'ACR';
S.names = names; %names for the different images
%actual function
FromCoordMatrix_2_3DNifti_8mm_LBPD_D(S);

%% FROM 3D NIFTI IMAGE (OR ALL ROIs OR MNI COORDINATES) TO LBPD COORDINATES
%not needed here but useful for future references

S = [];
S.input = 3; %1 = MNI coordinates; 2 = AAL ROIs; 3 = general image with non-zero values
S.coordd = [38 18 -16; -22 58 -16]; %coordinates in MNI space (x,y,z)
S.AAL_ROIs = [80]; %AAL ROIs numbers you want to use
% S.image = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/Contr_1_abs_0.nii.gz';
S.image = '/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/3000.nii.gz';
%actual function
idx_LBPD = From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D(S);

%%


%% GETTING FINAL ROIs IMAGES FILTERING AWAY THE WHITE MATTER

maskGM = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/osl/std_masks/lb_combinedGM8mm_mask_thr.nii.gz');
list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/MainROIs/*.gz');

for ii = 1:length(list)
    dum = load_nii([list(ii).folder '/' list(ii).name]);
    dum2 = dum.img;
    dum2(maskGM.img~=1) = 0;bdnf
    nii = make_nii(dum2,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dum2; %storing matrix within image structure
    nii.hdr.hist = maskGM.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - condition ' num2str(ii)])
    save_nii(nii,[list(ii).folder '/' list(ii).name(1:end-7) '_GM.nii.gz']); %printing image
end


%%


%% *** PLOTTING (AND DOING STATISTICS) ON THE TIME SERIES OF THE FUNCTIONALLY DERIVED ROIs ***

%%

%% PREPARING DATA FROM EACH SUBJECT USING THE ABOVE SELECTED ROIs (THIS SECTION WOULD ALSO PLOT DATA FROM OTHER DATASETS COLLECTED AT THE SAME TIME AS THE ONE USED IN THIS PAPER)

%task = 6 is the task reported in this paper
task = 6; %1 = elderly; 2 = fast/slow; 3 = learningbach; 4 = numbers (auditory); 5 = encoding block 5; 6 = Block 3 APR2020; 7 = MMN

%loading ROIs coordinatessou
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat')
if task == 1
    %list of source reconstructed data (for each subject)
    list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_1/SUBJ*.mat');
    %getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_1/sources_main_effects.mat');
elseif task == 2
    %list of source reconstructed data (for each subject)
    list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/SUBJ*.mat');
    %getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/sources_main_effects.mat');
elseif task == 3
    %list of source reconstructed data (for each subject)
    list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/source_LBPD/LearningBach/Beam_abs_0_sens_1_freq_broadband_time_1_466/SUBJ*.mat');
    %getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/source_LBPD/LearningBach/Beam_abs_0_sens_1_freq_broadband_time_1_466/sources_main_effects.mat');
elseif task == 4
    %list of source reconstructed data (for each subject)
    list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/SUBJ*.mat');
    %getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/sources_main_effects.mat');
elseif task == 5
    %list of source reconstructed data (for each subject)
    list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_5/Beam_abs_0_sens_1_freq_broadband_invers_1/SUBJ*.mat');
    %getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_5/Beam_abs_0_sens_1_freq_broadband_invers_1/sources_main_effects.mat');
elseif task == 6
    %list of source reconstructed data (for each subject)
    list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/SUBJ*.mat');
    %getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/sources_main_effects.mat');
elseif task == 7
    list = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_2/Beam_abs_0_sens_1_freq_broadband_invers_1/SUBJ*.mat');
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_2/Beam_abs_0_sens_1_freq_broadband_invers_1/MMNSubtracted_Average.mat');
end
if task == 3
    timex = 32:37;
elseif task == 5 || task == 4
    timex = 49:55;
elseif task == 7
    timex = 150:160;
else
    timex = 45:52;
end
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if task == 7
        if squeeze(mean(mean(data(jj,timex,1,:),4),2)) > 0 %if the data in voxel jj is positive during MMN time
            vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
        else
            vect(jj,1) = 1;
        end
    else
        if squeeze(mean(t_val_s(jj,timex,1),2)) > 0 %if the data in voxel jj is positive during N100 time
            vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
        else
            vect(jj,1) = 1;
        end
    end
end
ROIII = {1,2,3,4,5,6}; %selected ROIs; 1 = MC, 2 = HITR, 3 = HITL, 4 = VMPFC, 5 = ACL, 6 = ACR
load([list(1).folder '/' list(1).name])
if task == 7
    dum2 = zeros(length(ROIII),size(OUT.sources_ERFs,2),2,length(list));
else
    dum2 = zeros(length(ROIII),size(OUT.sources_ERFs,2),size(OUT.sources_ERFs,3),length(list));
end
for ii = 1:length(list) %over subjects
    disp(['loading source reconstructed data for subject ' num2str(ii) ' / ' num2str(length(list))])
    load([list(ii).folder '/' list(ii).name])
    if task == 7
        dum = zeros(size(OUT.sources_ERFs,1),size(OUT.sources_ERFs,2),2);
    else
        dum = zeros(size(OUT.sources_ERFs,1),size(OUT.sources_ERFs,2),size(OUT.sources_ERFs,3));
    end
    if task == 7
        dum3 = (OUT.sources_ERFs(:,:,3) - OUT.sources_ERFs(:,:,1));
        dum4 = (OUT.sources_ERFs(:,:,4) - OUT.sources_ERFs(:,:,2));
    end
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(OUT.sources_ERFs,1) %over brain voxels
            if task == 7
                if cc == 1 %here I subtract standard from deviant to obtain MMN (either global)
                    dum(jj,:,cc) = dum3(jj,:) .* vect(jj,1); %reversing (or not)..
                else %or local
                    dum(jj,:,cc) = dum4(jj,:) .* vect(jj,1); %reversing (or not)..
                end
            else
                dum(jj,:,cc) = OUT.sources_ERFs(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            end
            disp(['block ' num2str(task) ' - subject - ' num2str(ii) ' - condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for cc = 1:size(dum,3) %over conditions
        for pp = 1:length(ROIII) %over ROIs (with possibility of averaging together voxels of different ROIs (e.g. left and right auditory cortex (AC))
            dum2(pp,:,cc,ii) = mean(dum(sum(ROIs_DCM(:,ROIII{pp}),2)~=0,:,cc),1); %getting the indices of voxels of both ROIs (if you want two ROIs)
        end
    end
end
if task == 7
    condds = cell(1,2); condds{1} = 'GlobalMMN'; condds{2} = 'LocalMMN';
else
    condds = OUT.S.inversion.conditions;
end
save([list(1).folder '/ROIs_6.mat'],'dum2','condds')

%% PLOTTING MAIN ROIs AND CONDITIONS TOGETHER (ALL SUBJECTS AT THE SAME TIME)

conditions = [1 2]; %vector with conditions; tasks 1 and 4: 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4
%                      task 5: 1 = encoding; 2 = Old; 3 = NewT1); task 3: 1 = Old; 2 = New
%                      task 2: 1 = OldF; 2 = OldS; 3 = NewF; 4 = NewS
ROII = {1}; %selected ROIs; 1 = MC, 2 = HITR, 3 = HITL, 4 = VMPFC, 5 = ACL, 6 = ACR; ([2 3] averages together ROII 2 and 3; [2:4] averages together ROIs 2,3,4)
task = 3; %1 = block 3 TSA2021; 2 = fast/slow; 3 = learningbach; 4 = numbers (auditory); 5 = encoding block 5; 6 = block 3 APR2020
limmy = []; %limit for y-axis (amplitude of the signal)
% limmy = [];
leg_l = 1; %1 for legend; 0 otherwise
figexp = 0; %1 = export figures; 0 = not export figures


% close all
%defining colors
color_line = colormap(lines(5)); %extracting some colours from a colormap
color_line2 = color_line;
color_line2(1,:) = color_line(2,:);
color_line2(2,:) = color_line(1,:);
color_line2(5,:) = [0.4 0.4 0.4];
if task == 1
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_1/ROIs_6.mat');
    conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT3 ';
elseif task == 2
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/ROIs_6.mat');
    conds{1} = ' MemF '; conds{2} = ' MemS '; conds{3} = ' NewF '; conds{4} = ' NewS ';
elseif task == 3
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/source_LBPD/LearningBach/Beam_abs_0_sens_1_freq_broadband_time_1_466/ROIs_6.mat');
    conds{1} = ' Mem '; conds{2} = ' New ';
elseif task == 4
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/Block_4/Beam_abs_0_sens_1_freq_broadband_invers_1/ROIs_6.mat');
    conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT2 '; conds{4} = ' NewT3 '; conds{5} = ' NewT4 ';
elseif task == 5
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_5/Beam_abs_0_sens_1_freq_broadband_invers_1/ROIs_6.mat');
    conds{1} = ' Enc '; conds{2} = ' Mem '; conds{3} = ' NewT1 ';
elseif task == 6
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/ROIs_6.mat');
    conds{1} = ' Mem '; conds{2} = ' NewT1 '; conds{3} = ' NewT2 '; conds{4} = ' NewT3 '; conds{5} = ' NewT4 ';
end
if task == 3
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/source_LBPD/LearningBach/time.mat')
else
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/Source_LBPD/time_345.mat')
end
% COL{1} = 'b'; COL{2} = 'r'; COL{3} = 'k'; COL{4} = 'g'; COL{5} = 'c'; COL{6} = 'm'; %colors
if length(conditions) == 2 %line width
    LIN{1} = 2; LIN{2} = 1;
elseif length(conditions) == 3 %line width
    LIN{1} = 2.3; LIN{2} = 1.8; LIN{3} = 1;
elseif length(conditions) == 4 %line width
    LIN{1} = 2.8; LIN{2} = 2.3; LIN{3} = 1.8; LIN{4} = 1;
else
    LIN{1} = 3.2; LIN{2} = 2.8; LIN{3} = 2.3; LIN{4} = 1.8; LIN{5} = 1;
end
ROIN{1} = 'MC'; ROIN{2} = 'HITR'; ROIN{3} = 'HITL'; ROIN{4} = 'VMPFC'; ROIN{5} = 'ACL'; ROIN{6} = 'ACR';
% load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat')
figure
clear CC
for cc = 1:length(conditions) %over conditions
    for ii = 1:length(ROII) %over ROIs
        if length(ROII) == 1
            plot(time(1:size(dum2,2)),nanmean(dum2(ROII{ii},:,conditions(cc),:),4),'Color',color_line2(cc,:),'LineWidth',2,'DisplayName',[conds{conditions(cc)}])
            hold on
            plot(time(1:size(dum2,2)),nanmean(dum2(ROII{ii},:,conditions(cc),:),4) + (nanstd(dum2(ROII{ii},:,conditions(cc),:),0,4)./sqrt(size(dum2,4))),':','Color',color_line2(cc,:),'LineWidth',0.5,'DisplayName',[conds{conditions(cc)}])
            hold on
            plot(time(1:size(dum2,2)),nanmean(dum2(ROII{ii},:,conditions(cc),:),4) - (nanstd(dum2(ROII{ii},:,conditions(cc),:),0,4)./sqrt(size(dum2,4))),':','Color',color_line2(cc,:),'LineWidth',0.5,'DisplayName',[conds{conditions(cc)}])
        else
            if length(ROII{ii}) == 1
                plot(time(1:size(dum2,2)),nanmean(dum2(ROII{ii},:,conditions(cc),:),4),'Color',color_line(ii,:),'LineWidth',LIN{cc})
            else
                plot(time(1:size(dum2,2)),nanmean(nanmean(dum2(ROII{ii},:,conditions(cc),:),4),1),'Color',color_line(ii,:),'LineWidth',LIN{cc})
            end
        end
        hold on
    end
end
grid minor
if length(ROII) == 1
    if leg_l == 1
        legend('show')
    end
    title(ROIN{ROII{1}})
else %elaborated way (and probably barbaric) to create legends..
    dumR = cell(1,length(ROII)); %trick to get the legend..
    for ii = 1:length(ROII)
        if length(ROII{ii}) == 1
            dumR(ii) = ROIN(ROII{ii});
        else
            dumm = [];
            for pp = 1:length(ROII{ii})
                dumm = strcat(dumm,ROIN{ROII{ii}(pp)});
            end
            dumR(ii) = {dumm};
        end
    end
    if leg_l == 1
        legend(dumR)
    end
    s = [];
    for ii = 1:length(conditions) %trick to get the title..
        s = strcat(s,conds{conditions(ii)});
    end
    title(s)
end
% xlim([-0.1 time(size(dum2,2))])
xlim([-0.1 3.4])
if ~isempty(limmy)
    ylim(limmy)
end
set(gcf,'color','w')
if figexp == 1
    dumm = [];
    for pp = 1:length(ROII)
        dumm = strcat(dumm,ROIN{ROII{pp}});
    end
    if ~isempty(limmy)
        export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Figures_3_4_5/Cond_task' num2str(task) '_' conds{conditions} '_' dumm '_Scaled.png'])
        export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Figures_3_4_5/Cond_task' num2str(task) '_' conds{conditions} '_' dumm '_Scaled.eps'])
    else
        export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Figures_3_4_5/Cond_task' num2str(task) '_' conds{conditions} '_' dumm '.png'])
        export_fig(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Figures_3_4_5/Cond_task' num2str(task) '_' conds{conditions} '_' dumm '.eps'])
    end
end

%% COMPUTING STATISTICS CONTRASTING OLD VERSUS THE FOUR CATEGORIES OF NEW, ONE AT A TIME (THIS IS DONE FOR EACH OF THE SIX ROIs)

%this is the file that you have to load and use to run the statistics
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/ROIs_6.mat');
clear ROIN
ROIN{1} = 'MC'; ROIN{2} = 'HITR'; ROIN{3} = 'HITL'; ROIN{4} = 'VMPFC'; ROIN{5} = 'ACL'; ROIN{6} = 'ACR'; % ROIs in dum2

%t-tests
P = zeros(size(dum2,1),size(dum2,2),(size(dum2,3)-1)); %ROIs x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dum2,1),size(dum2,2),(size(dum2,3)-1)); %ROIs x time-points x contrasts (every NewTX versus Old)
for ii = 1:size(dum2,1) %over ROIs
    for jj = 1:size(dum2,2) %ove time-points
        for cc = 1:(size(dum2,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dum2(ii,jj,1,:)),squeeze(dum2(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end
end
%MCS
p_thresh = 0.05; %threshold for binarising p-values vector
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Statistics';
for ii = 1:size(dum2,1) %over ROIs
    for cc = 1:(size(dum2,3)-1) %over the 4 NewTX
        Pbin = zeros(1,size(dum2,2));
        Pbin(P(ii,:,cc)<p_thresh) = 1;
        tvals = T(ii,:,cc);
        [ sign_clust ] = oneD_MCS_LBPD_D( Pbin, 0, 1, 1000, 0.001, time(1:size(dum2,2)), tvals ) %Monte Carlo simulation function to correct for multiple comparisons
        PDn = cell2table(sign_clust); %table
        writetable(PDn,[outdir '/' ROIN{ii} '_OldvsNewT' num2str(cc) '.xlsx'],'Sheet',1); %printing excel file
    end
end

%% THE FOLLOWING SECTION WAS NOT REPORTED IN THE PAPER

%% FOCUS ON PREDICTION ERROR - EXTRACTING PEAKS FOR EACH SUBJECT OF PREDICITON ERROR IN ACL AND VMPFC AND COMPUTING ANOVAs

peakk = [1 4]; %peaks that you want to use

load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat'); %loading time
%this is the file that you have to load and use to run the statistics
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/ROIs_6.mat');
ROIN{1} = 'MC'; ROIN{2} = 'HITR'; ROIN{3} = 'HITL'; ROIN{4} = 'VMPFC'; ROIN{5} = 'ACL'; ROIN{6} = 'ACR'; % ROIs in dum2
mat_anal_ACL = {[159:169] [246:256] [344:354] [434:444]}; %time-points (already with plus/minus 20ms) of the peaks for ACL
mat_anal_VMPFC = {[174:184] [266:276] [359:369] [451:461]}; %time-points (already with plus/minus 20ms) of the peaks for VMPFC
% mat_anal_ACL = {[164] [251] [349] [439]}; %time-points (already with plus/minus 20ms) of the peaks for ACL
% mat_anal_VMPFC = {[179] [271] [364] [456]}; %time-points (already with plus/minus 20ms) of the peaks for VMPFC


mat_anal_ACL = mat_anal_ACL(peakk);
mat_anal_VMPFC = mat_anal_VMPFC(peakk);

DUM{1} = mat_anal_ACL; DUM{2} = mat_anal_VMPFC;
roisdum = [5 4]; %indicies of ACL and VMPFC
data = zeros(2,length(mat_anal_ACL),size(dum2,4)); %ROIs (ACL and VMPFC) x time-points x subjects
for ii = 1:2 %over ROIs (ACL and VMPFC)
    dumm = DUM{ii};
    for tt = 1:length(mat_anal_ACL) %over peaks
        data(ii,tt,:) = mean(dum2(roisdum(ii),dumm{tt},peakk(tt)+1,:),2);
    end
end
%statistics (ANOVAs)
dataACL = squeeze(data(1,:,:))'; %data ACL
dataVMPFC = squeeze(data(2,:,:))'; %data VMPFC
%ACL
[p,t,stats] = anova1(dataACL); %'off' for not showing the plot
[c,m,h,nms] = multcompare(stats,'ctype','tukey-kramer'); %perform multiple comparison test based on anova1 output = post-hoc analysis (c needs to be saved for every time point and every channel)
%VMPFC
[p,t,stats] = anova1(dataVMPFC); %'off' for not showing the plot
[c,m,h,nms] = multcompare(stats,'ctype','tukey-kramer'); %perform multiple comparison test based on anova1 output = post-hoc analysis (c needs to be saved for every time point and every channel)
%plotting correct responses
%some default color specifications for later plotting..
% cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];
cl = [0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8;];
datadum = cell(1,size(dataACL,2));
for ii = 1:size(dataACL,2)
    datadum{1,ii} = dataACL(:,ii);
end
figure
rm_raincloud2(datadum',cl)
grid minor
set(gcf,'color','w')
set(gcf,'Position',[200,200,400,550])
xlabel('Amplitude'); %set(gca,'YDir','normal');
% xlim([-5 38])

cl = [0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8; 0 0 0.8;];
datadum = cell(1,size(dataVMPFC,2));
for ii = 1:size(dataVMPFC,2)
    datadum{1,ii} = dataVMPFC(:,ii);
end
figure
rm_raincloud2(datadum',cl)
grid minor
set(gcf,'color','w')
set(gcf,'Position',[200,200,400,550])
xlabel('Amplitude'); %set(gca,'YDir','normal');

%%

%% *** TIME-FREQUENCY ANALYSIS ***

%% COMPUTING TIME-FREQUENCY ANALYSIS ON SIGNLE VOXELS AND THEN AVERAGING THE RESULTS ('PROPER' INDUCED RESPONSES)

% A) here you provide a matrix with LBPD coordinates with 0s and 1s to index specific ROIs 

% addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/BrainstormDemo_CodesTimeFrequencyDecomposition/functions') %path to Dimitrios function for Morlet transform

%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'none');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 3); %slot in the queu

%% ON THE ROIs ACTUALLY USED IN THE PAPER

ROII{1} = 'MC'; ROII{2} = 'HITR'; ROII{3} = 'HITL'; ROII{4} = 'VMPFC'; ROII{5} = 'ACL'; ROII{6} = 'ACR';
for ii = 1%:6 %over ROIs
    S = [];
    S.ROI_n = ii; %1 = MC; 2 = HITLR; 3 = VMPFC; 4 = ACLR; 5 = ACL; 6 = ACR
    S.f = 1:1:60; %frequencies used
    S.mask = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/DCM/ROIs_MainActivity/FinalROIs/ROIsDCM.mat';
    list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/SUBJ*.mat']);
        S.subjlist = list(82:end);
%     S.subjlist = list;
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
    time(1027:end) = [];
    S.time = time;
    S.Aarhus_clust = 1; %0 = working locally; integer number (i.e. 1) = sending one job for each subject to the Aarhus cluster (the number you insert here corresponds to the slots of memory that you allocate for each job.
    S.single_subj = 1; %1 for saving single subject data; 0 for not saving it
    S.outdir = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Time_Frequency_Induced/' ROII{ii} '_SingleSubjects_SingleVoxels.mat'];
    %actual function
    jobid = job2cluster(@InducedResponses_Morlet_ROIs_LBPD_D,S);
end

%% COMPUTING TIME-FREQUENCY ANALYSIS ON SINGLE VOXELS OF SPECIFIC COORDINATES OR AAL ROI(s) AND THEN AVERAGING THE RESULTS ('PROPER' INDUCED RESPONSES)

%%% NOT REPORTED IN THE PAPER BUT USEFUL TO HAVE IT HERE.. %%%

% B) here you provide either coordinates in MNI space or the AAL ROIs (following the AAL order)

namee = 'Occ_Sup_R'; %if you use the AAL ROI(s), here you can specify the name to be used when saving the data


%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'none'); % 'none' or 'cluster'
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 3); %slot in the queu

S = [];
S.Aarhus_clust = 1; %0 = working locally; integer number (i.e. 1) = sending one job for each subject to the Aarhus cluster (the number you insert here corresponds to the slots of memory that you allocate for each job.
S.average_trials = 0; %1 = frequency decomposition after averaging the trials (evoked responses); 0 = frequency decomposition for single trials and then average of the time-frequency results (induced responses)
S.f = 1:1:60; %frequencies used
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
S.time = time(1:1026);
S.conds = [1 2]; %experimental conditions
% S.coordd = [-38 18 -16; -22 58 -16; -30 34 -16]; %coordinates (empty for AAL ROIs)
S.coordd = []; %coordinates (empty for AAL ROIs)
% S.AAL_ROIs =[5,6,9,10,15,16,25,26]; %AAL ROIs indices (meaningful only if S.coordd is empty)
S.AAL_ROIs =[50]; %AAL ROIs indices (meaningful only if S.coordd is empty)
if S.average_trials == 1
    S.subjlist = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/SUBJ*.mat']);
else
    S.subjlist = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/SUBJ*.mat']); 
end
if isempty(S.coordd)
    S.outdir = [S.subjlist(1).folder '/Time_Freq_AllSubjs_AAL_' namee]; %remember to specify the proper AAL ROI name that you want!!
else
    S.outdir = [S.subjlist(1).folder '/Time_Freq_AllSubjs_singlevoxels'];
end

jobid = job2cluster(@InducedResponses_Morlet_Coords_AALROIs_LBPD_D,S);

%%

% alternative solution

%% COMPUTATION OF TIME-FREQUENCY ANALYSIS FOR ALL VOXELS OF THE BRAIN (DEPRECATED SINCE IT REQUIRES A LOT OF MEMORY..)

%this is an alternative way to do that.. so first computing analysis on whole brain and then index the ROIs
%NOT SUGGESTED FOR MEMORY ISSUES.. BUT STILL WORTHY TO BE REPORTED HERE

%currently, it seems to work only locally.. then, the single jobs for the single subjects are fine on the cluster.. surprising.. anyway..
%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'none'); % 'none' or 'cluster'
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 3); %slot in the queu

S = [];
S.Aarhus_clust = 1; %1 = Aarhus cluster for parallel computing (i.e. running one job per subject); 0 = local computation
S.f = 1:1:60; %frequencies used
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
S.time = time(1:1026);
S.conds = [1:5]; %experimental conditions
S.subjlist = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/SUBJ*.mat']);
S.outdir = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Time_Frequency_Induced/WholeBrain_SingleSubjects_SingleVoxels.mat'];

jobid = job2cluster(@InducedResponses_Morlet_WholeBrain_LBPD_D,S);

%% plotting mean conditions

ROI = 5; % 1 = MC; 2 = HITLR; 3 = VMPFC; 4 = ACLR; 5 = VMPFC with Morlet computed independently on each voxel and each trial and then the output was averaged (twice)
condition = 1; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4
fff = 20:60; %frequencies to be plotted
average_trials = 1; %1 = frequency decomposition after averaging the trials (evoked responses); 0 = frequency decomposition for single trials and then average of the time-frequency results (induced responses)


%ROIs name
ROII{1} = 'MC'; ROII{2} = 'HITLR'; ROII{3} = 'VMPFC'; ROII{4} = 'ACLR'; ROII{5} = 'VMPFC_vox';
%loading data
if ROI == 5
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/DCM/VMPFC_SingleSubjects_SingleVoxels.mat') %average of trials and then time-frequency
else
    if average_trials == 1
        load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/DCM/Time_Frequency_AllSubjects_AveragedTrials.mat') %average of trials and then time-frequency
    else
        load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/DCM/Time_Frequency_AllSubjects_SingleTrials.mat') %time-frequency on single trials and then average
    end
end
%mean over subjects (for plotting purposes) and removing pre-stimulus time
if ROI == 5
    Pold = mean(P3,3);
else
    Pold = mean(P2,5); %mean over subjects
end
%plotting
figure
if ROI == 5
    imagesc(time,f(fff),squeeze(Pold(fff,:)))
else
    imagesc(time,f(fff),squeeze(Pold(ROI,fff,:,condition)))
end
set(gca,'YDir','normal') %plotting frequencies in descending order
xlabel('time (s)'); ylabel('f (Hz)');
colorbar
% if ges ~= 0
%     caxis([0 300])
% end
set(gcf,'color','w')
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))
if ROI == 5
    title(['ROI ' ROII{ROI} ' - cond Old'])
else
    title(['ROI ' ROII{ROI} ' - cond ' num2str(condition)])
end

%% PLOTTING - LOADING INDEPENDENT SUBJECTS (COMPUTED INDEPENDENT VOXELS AND THEN AVERAGED IN THE ROIs)

condition = 0; % 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT5; 0 = Old - NewT1
contr = [1 2 3 4 5]; %conditions to be contrasted (NOW IT'S FROM 1 TO 5 BECAUSE I CREATED A LOOP FOR DOING ALL CONTRASTS), meaningul only if condition == 0
fff = 1:60; %frequencies to be plotted
ROI = 6; % 1 = MC; 2 = HITL; 3 = HITR; 4 = VMPFC; 5 = ACL; 6 = ACR with Morlet computed independently on each voxel and each trial and then the output was averaged (twice)
cacs = []; %limits for plots; leave empty for not having them
loadl = 1; %1 for loading; 0 for not loading, in case you already loaded the data and want to plot different conditions

outdir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband/Sources_Decoding_MEG0211/MCS/FinalImagesToWorkBench/Statistics/Time_Frequency_MCS';
ROII{1} = 'MC'; ROII{2} = 'HITL'; ROII{3} = 'HITR'; ROII{4} = 'VMPFC'; ROII{5} = 'ACL'; ROII{6} = 'ACR';
for ss = 1:length(ROII)
    ROI = ss;
    if loadl == 1
        if ROI == 1
            list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Time_Frequency_Induced/MC*.mat']);
        elseif ROI == 2
            list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Time_Frequency_Induced/HITL*.mat']);
        elseif ROI == 3
            list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Time_Frequency_Induced/HITR*.mat']);
        elseif ROI == 4
            list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Time_Frequency_Induced/VMPFC*.mat']);
        elseif ROI == 5
            list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Time_Frequency_Induced/ACL*.mat']);
        elseif ROI == 6
            list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Time_Frequency_Induced/ACR*.mat']);
        end
        load([list(1).folder '/' list(1).name]);
        P2 = zeros(size(Psubj,1),size(Psubj,2),size(Psubj,3),length(list));
        for ii = 1:length(list) %over subjects
            load([list(ii).folder '/' list(ii).name]);
            P2(:,:,:,ii) = Psubj;
            disp(ii)
        end
    end
    for cc = 1:4 %over contrasts (old versus NewTX)
        %plotting
        figure
        if condition ~= 0
            %mean over subjects (for plotting purposes) and removing pre-stimulus time
            Pold = squeeze(mean(P2,4)); %mean over subjects
            imagesc(time,f(fff),squeeze(Pold(fff,:,condition)))
        else
            Pold = P2;
            %     Pold = Pold - mean(Pold(:,1:26,:,:),2); %baseline correction
            %t-tests
            %computing t-tests
            P = zeros(size(Pold,1),size(Pold,2));
            T = zeros(size(Pold,1),size(Pold,2));
            for ii = 1:size(P,1) %over frequencies
                for jj = 1:size(P,2) %over time-points
                    [~,p,~,stats] = ttest(squeeze(Pold(ii,jj,contr(1),:)),squeeze(Pold(ii,jj,contr(cc+1),:))); %contrasting Old vs New
                    P(ii,jj) = p;
                    T(ii,jj) = stats.tstat;
                end
                disp(ii)
            end
            %     T(abs(T)<1.8) = 0;
            imagesc(time,f(fff),squeeze(T(fff,:)))
            %simple difference of mean
            %     Pold = mean(Pold,4); %average over subjects
            %     Pdum = Pold(:,:,1) - Pold(:,:,2);
            %     imagesc(time,f(fff),squeeze(Pdum(fff,:)))
        end
        set(gca,'YDir','normal') %plotting frequencies in descending order
        xlabel('time (s)'); ylabel('f (Hz)');
        colorbar
        % if ges ~= 0
%         if ~isempty(cacs)
%             caxis(cacs)   
        caxis([-6 6])
%         end
        % end
        set(gcf,'color','w')
        %colormap with white for 0 values
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        colormap(bluewhitered_PD(0,x))
        %export figures
        export_fig([outdir '/TF_Full_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ss} '.eps'])
        export_fig([outdir '/TF_Full_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ss} '.png'])
        %figure with significant values only
        figure
        T(abs(T)<2) = 0;
        imagesc(time,f(fff),squeeze(T(fff,:)))
        set(gca,'YDir','normal') %plotting frequencies in descending order
        xlabel('time (s)'); ylabel('f (Hz)');
        colorbar
        % if ges ~= 0
%         if ~isempty(cacs)
%             caxis(cacs)
%         end
        caxis([-6 6])
        % end
        set(gcf,'color','w')
        %colormap with white for 0 values
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        colormap(bluewhitered_PD(0,x))
        title(['Condition ' num2str(condition)])
        export_fig([outdir '/TF_Sign_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ss} '.eps'])
        export_fig([outdir '/TF_Sign_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ss} '.png'])
        
        %testing contrast results (correcting for multiple comparison) by performing Monte Carlo simulations
        P3 = zeros(size(T,1),size(T,2));
        P3(abs(T)>2) = 1; %threshold t-val = 2.6 corresponding to p-val < 0.01 (obtained by dividing 0.05 by the 4 comparisons employed here)
        thresh = 0;
        permut = 1000;
        threshMC = 0.001;
        perm_max = 1;
        t1 = f(fff); t2 = time;
        
        [ OUT ] = twoD_MCS_LBPD_D( P3, thresh, permut, threshMC, perm_max, t1 , t2 )
        
        PDn = cell2table(OUT); %table
        writetable(PDn,[outdir '/TF_cond' num2str(contr(1)) '_vs_cond' num2str(contr(cc+1)) '_' ROII{ss} '.xlsx'],'Sheet',1); %printing excel file
    end
end

%% not reported in the paper but potentially useful to have it as follows

%% plotting single ROIs from AAL (such as Heschl's gyrus)

condition = 0; % 1 = Old; 2 = NewT1; 0 = Old - NewT1
fff = 1:60; %frequencies to be plotted
single_subj = 1; %analysis with single subjects saved on disk 
listn = 1; %index of the analysis to be run (with reference to the variable "list)
loadl = 0; %1 for loading; 0 for not loading, in case you already loaded the data and want to plot different conditions


if loadl == 1
    %loading data
    if single_subj == 1
        list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/Time_*']);
        list = list([list(:).isdir]); %keeping only the folders (otherwise you would get .mat files as well..)
        list2 = dir([list(listn).folder '/' list(listn).name '/SUBJ*.mat']);
        load([list2(1).folder '/' list2(1).name]); %to get some information
%         load([list(listn).folder '/' list(listn).name]);
        P2 = zeros(size(Psubj,1),size(Psubj,2),size(Psubj,3),size(Psubj,4),length(list2));
        for ii = 1:length(list2)
            load([list2(ii).folder '/' list2(ii).name]);
            P2(:,:,:,:,ii) = Psubj;
            disp(['loading subject ' num2str(ii)])
        end
    else
        list = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/Time_*.mat']);
        load([list(listn).folder '/' list(listn).name]);
    end
end
%plotting
figure
if condition ~= 0
    %mean over subjects (for plotting purposes) and removing pre-stimulus time
    Pold = squeeze(mean(mean(P2,5),1)); %mean over subjects and over voxels
    imagesc(time,f(fff),squeeze(Pold(fff,:,condition)))
else
    Pold = squeeze(mean(P2,1)); %mean over voxels
    %t-tests
    %computing t-tests
    P = zeros(size(Pold,1),size(Pold,2));
    T = zeros(size(Pold,1),size(Pold,2));
    for ii = 1:size(P,1) %over frequencies
        for jj = 1:size(P,2) %over time-points
            [~,p,~,stats] = ttest(squeeze(Pold(ii,jj,1,:)),squeeze(Pold(ii,jj,2,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj) = p;
            T(ii,jj) = stats.tstat;
        end
        disp(ii)
    end
%     T(abs(T)<1.8) = 0;
    imagesc(time,f(fff),squeeze(T(fff,:)))
    %simple difference of mean
%     Pold = mean(Pold,4); %average over subjects
%     Pdum = Pold(:,:,1) - Pold(:,:,2);
%     imagesc(time,f(fff),squeeze(Pdum(fff,:)))
end
set(gca,'YDir','normal') %plotting frequencies in descending order
xlabel('time (s)'); ylabel('f (Hz)');
colorbar
% if ges ~= 0
% caxis([-3 3])
% end
set(gcf,'color','w')
%colormap with white for 0 values
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))
title(['Condition ' num2str(condition)])

%%

%% *** CROSS-CORRELATIONS ***

%% PLOTTING MAIN ROIs AND CONDITIONS TOGETHER AND CROSS-CORRELATIONS BETWEEN TIME SERIES OF AVERAGED DATA OVER SUBJECTS

% in the paper it is reported 5 (ACL) and 4 (VMPFC)

conditions = [2]; %vector with conditions; 1 = Old; 2 = NewT1; 3 = NewT2; 4 = NewT3; 5 = NewT4
ROII = {5 4}; %selected ROIs; 1 = MC, 2 = HITR, 3 = HITL, 4 = VMPFC, 5 = ACL, 6 = ACR; ([2 3] averages together ROII 2 and 3; [2:4] averages together ROIs 2,3,4)
xcorrl = 1; %1 for cross correlation (meaningful only if you have two timeseries..)
xtime = [101:526]; %time-points for cross correlation

% close all
if ~exist('DCM_data','var') %loading data, if not already loaded
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/DCMfinal_ROIsXtimeXtrXsubsXconds.mat');
end
conds{1} = ' Old '; conds{2} = ' NewT1 '; conds{3} = ' NewT2 '; conds{4} = ' NewT3 '; conds{5} = ' NewT4 ';
COL{1} = 'b'; COL{2} = 'r'; COL{3} = 'k'; COL{4} = 'm'; COL{5} = 'g'; COL{6} = 'c'; %colors
if length(conditions) == 2 %line width
    LIN{1} = 2; LIN{2} = 1;
elseif length(conditions) == 3 %line width
    LIN{1} = 2.3; LIN{2} = 1.8; LIN{3} = 1;
elseif length(conditions) == 4 %line width
    LIN{1} = 2.8; LIN{2} = 2.3; LIN{3} = 1.8; LIN{4} = 1;
else
    LIN{1} = 3.2; LIN{2} = 2.8; LIN{3} = 2.3; LIN{4} = 1.8; LIN{5} = 1;
end
ROIN{1} = 'MC'; ROIN{2} = 'HITR'; ROIN{3} = 'HITL'; ROIN{4} = 'VMPFC'; ROIN{5} = 'ACL'; ROIN{6} = 'ACR';
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat')
figure
clear CC
for cc = 1:length(conditions) %over conditions
    for ii = 1:length(ROII) %over ROIs
        if length(ROII{ii}) == 1
            plot(time(1:1026),nanmean(nanmean(DCM_data(ROII{ii},:,:,:,conditions(cc)),4),3),'Color',COL{ii},'LineWidth',LIN{cc})
            CC{ii} = nanmean(nanmean(DCM_data(ROII{ii},:,:,:,conditions(cc)),4),3);
        else
            plot(time(1:1026),nanmean(nanmean(nanmean(DCM_data(ROII{ii},:,:,:,conditions(cc)),4),3),1),'Color',COL{ii},'LineWidth',LIN{cc})
            CC{ii} = nanmean(nanmean(nanmean(DCM_data(ROII{ii},:,:,:,conditions(cc)),4),3),1);
        end
        
%         plot(time(1:1026),nanmean(nanmean(DCM_data(ROII(ii),:,:,:,conditions(cc)),4),3),'Color',COL{ii},'LineWidth',LIN{cc})
%         CC{ii} = nanmean(nanmean(DCM_data(ROII(ii),:,:,:,conditions(cc)),4),3);
        hold on
    end
end
grid minor
dumR = cell(1,length(ROII)); %trick to get the legend..
for ii = 1:length(ROII)
    if length(ROII{ii}) == 1
        dumR(ii) = ROIN(ROII{ii});
    else
        dumm = [];
        for pp = 1:length(ROII{ii})
            dumm = strcat(dumm,ROIN{ROII{ii}(pp)});
        end
        dumR(ii) = {dumm};
    end
end
legend(dumR)
% legend('show')
xlim([-0.1 3.4])
s = [];
for ii = 1:length(conditions) %trick to get the title..
    s = strcat(s,conds{conditions(ii)});
end
title(s)
set(gcf,'color','w')

%cross-correlation
if xcorrl == 1 && length(CC) == 2
%     [r,lags] = xcorr(CC{1}(25:800),CC{2}(25:800));
    [r,lags] = xcorr(CC{1}(xtime),CC{2}(xtime),'coeff');
    lags2 = lags.*0.004;
%     figure
%     stem(lags,r)
%     grid minor
    figure
    scatter(lags2,r)
    grid minor
    set(gcf,'color','w')
    xlim([lags2(1) lags2(end)])
end

%%

%% codes for additional brain figure

%loading MNI coordinates of AAL 2-mm centroids
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI_RC_AAL_2mm.mat')
%coordinates for centroids of ROIs
vcoord = [-42.3467 -20.0356 8.6756; 42.3467 -20.0356 8.6756; 0 -16.1010 40.2143; -25.2554 -21.9635 -11.3841; 25.2554 -21.9635 -11.3841; 0 52.5174 -8.8567];
vcol = [0 0.4470 0.7410; 0 0.4470 0.7410; 0.4000 0.4000 0.4000; 0.8500 0.0250 0.1980; 0.8500 0.0250 0.1980; 0.8500 0.0250 0.1980]; %color codes for ROIs
%setting connections
%GREY
SEG{1} = [-42.3467 -20.0356 8.6756; 0 -16.1010 40.2143]; %LAC-MC
SEG{2} = [42.3467 -20.0356 8.6756; 0 -16.1010 40.2143]; %RAC-MC
%BLUE
SEG{3} = [-42.3467 -20.0356 8.6756; -25.2554 -21.9635 -11.3841]; %LAC-LH
SEG{4} = [42.3467 -20.0356 8.6756; 25.2554 -21.9635 -11.3841]; %RAC-RH
%RED
SEG{5} = [-25.2554 -21.9635 -11.3841; 0 52.5174 -8.8567]; %LH-VMPFC
SEG{6} = [25.2554 -21.9635 -11.3841; 0 52.5174 -8.8567]; %RH-VMPFC
%GREY
SEG{7} = [-42.3467 -20.0356 8.6756; 0 -16.1010 40.2143]; %LAC-MC
SEG{8} = [42.3467 -20.0356 8.6756; 0 -16.1010 40.2143]; %RAC-MC
SEG{9} = [0 52.5174 -8.8567; 0 -16.1010 40.2143]; %VMPFC-MC
vcolcon = [0.4000 0.4000 0.4000; 0.4000 0.4000 0.4000; 0 0.4470 0.7410; 0 0.4470 0.7410; 0.8500 0.0250 0.1980; 0.8500 0.0250 0.1980; 0.4000 0.4000 0.4000; 0.4000 0.4000 0.4000; 0.4000 0.4000 0.4000]; %color codes for connections

%loading a brain template
openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
hold on
for ii = 1:length(vcoord) %over ROIs
    plot3(vcoord(ii, 1), vcoord(ii, 2), vcoord(ii, 3), ['.'], 'Color', vcol(ii,:), 'MarkerSize', 25); %centroid of ROI ii
    hold on
end
rotate3d on; axis off; axis vis3d; axis equal
for ii = 1:length(SEG)
    vdum = SEG{ii}; %vector with the two selected nodes (jj and zz)
    plot3(vdum(:,1),vdum(:,2),vdum(:,3),'Color',vcolcon(ii,:),'LineWidth',1.5)  %plotting line connecting the two points
    hold on
end
view([0 90])
export_fig('top.eps')
export_fig('top.png')
view([-90 0])
export_fig('left.eps')
export_fig('left.png')
view([90 0])
export_fig('right.eps')
export_fig('right.png')
view([0 -90])
export_fig('top2.eps')
export_fig('top2.png')
view([0 0])
export_fig('back.eps')
export_fig('back.png')
view([-180 0])
export_fig('front.eps')
export_fig('front.png')

%%
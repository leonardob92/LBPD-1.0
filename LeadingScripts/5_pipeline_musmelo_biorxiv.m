 
%% MUSMELO - PIPELINE

% Francesco Carlomagno under guidance and corrections of Leonardo Bonetti

%% maxfilter

%OBS!!! before running maxfilter, you need to close matlab, open the terminal, write: 'use anaconda', then open matlab and run maxfilter script

%new proper lines for maxfilter
maxfilter_path = '/neuro/bin/util/maxfilter';
project = 'MINDLAB2018_MEG-LearningBach-MemoryInformation';
maxDir = '/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F'; %output path

% cmd = ['submit_to_cluster -q maxfilter .q -n 4 -p ' ,num2str(project), ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 3 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
path = '/raw/sorted/MINDLAB2015_MEG-TunteetMM'; % path where are all the subjects folders
jj = dir([path '/0*']); % create a dir path ad search for all the folders starting with '0' in order to avoid hiden files
for ii = 6:length(jj) % loop to process all the subj folders
    cart = ([jj(2).folder '/' jj(ii).name]); % create a path combining the folder path name and the string with the name of the folder
    pnana = dir([cart '/2*']); %search for folder startint with a common detail
    for pp = 1:length(pnana) % loop to explore ad analyze all the folder inside the path above
        cart2 = ([pnana(1).folder '/' pnana(pp).name]); 
        pr = dir([cart2 '/M*']); % looks for meg folder
        idx = strfind(pr(1).name,'MEG');
        if ~isempty(idx) % if is not empty proceed with sub-folders inside the meg path
            pnunu = dir([pr(1).folder '/' pr(1).name '/00*']);
            if length(pnunu) > 1
                warning(['subj ' num2str(ii) ' has files nuber = ' num2str(length(pnunu))]) %show a warning message if any subj has more thatn 1 meg sub-folder
            end
            for dd = 1:length(pnunu)
                idx2 = strfind(pnunu(1).name,'Mus'); % search for musmelo folder in order to avoid other projects
                if ~isempty(idx2)
                    fpath = dir([pnunu(1).folder '/' pnunu(1).name '/' 'files' '/PRO*']); % combine the path with the common folder files and the PRO initial string of the .fif files present for every subj
                    rawName = ([fpath.folder '/' fpath.name]); %assigns the final path of the .fif file to the rawName path used in the maxfilter command
                    maxfName = [fpath.name(10:17) '_musmelo']; % define the output name of the maxfilter processing
                    %no movement compensation
                    cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,num2str(project), ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsss.fif'] ' -st 4 -corr 0.98 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsss.log"']];
                    system(cmd);
                end
            end
        end
    end
end

%% A FEW ADDITIONAL SUBJECTS RUN INDEPENDENTLY FOR MAXFILTER

maxfilter_path = '/neuro/bin/util/maxfilter';
project = 'MINDLAB2018_MEG-LearningBach-MemoryInformation';
maxDir = '/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F'; %output path


rawName = '/raw/sorted/MINDLAB2015_MEG-TunteetMM/0002/20150702_000000/MEG/001.MusMelo_2/files/PROJ0192_SUBJ0002_SER001_FILESNO001.fif';
maxfName = 'SUBJ0002_SER001_musmelo'; % define the output name of the maxfilter processing
%no movement compensation
cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,num2str(project), ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsss.fif'] ' -st 4 -corr 0.98 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsss.log"']];
system(cmd);
% 
% rawName = '/raw/sorted/MINDLAB2015_MEG-TunteetMM/0002/20150702_000000/MEG/002.MusMelo_1/files/PROJ0192_SUBJ0002_SER002_FILESNO001.fif';
% maxfName = 'SUBJ0002_SER002_musmelo'; % define the output name of the maxfilter processing
% %no movement compensation
% cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,num2str(project), ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsss.fif'] ' -st 4 -corr 0.98 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsss.log"']];
% system(cmd);

rawName = '/raw/sorted/MINDLAB2015_MEG-TunteetMM/0016/20150702_000000/MEG/001.MusMelo_2/files/PROJ0192_SUBJ0016_SER001_FILESNO001.fif';
maxfName = 'SUBJ0016_SER001_musmelo'; % define the output name of the maxfilter processing
%no movement compensation
cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,num2str(project), ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsss.fif'] ' -st 4 -corr 0.98 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsss.log"']];
system(cmd);

% rawName = '/raw/sorted/MINDLAB2015_MEG-TunteetMM/0016/20150702_000000/MEG/002.MusMelo_1/files/PROJ0192_SUBJ0016_SER002_FILESNO001.fif';
% maxfName = 'SUBJ0016_SER002_musmelo'; % define the output name of the maxfilter processing
% %no movement compensation
% cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,num2str(project), ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsss.fif'] ' -st 4 -corr 0.98 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsss.log"']];
% system(cmd);

rawName = '/raw/sorted/MINDLAB2015_MEG-TunteetMM/0073/20150702_000000/MEG/001.MusMelo_2/files/PROJ0192_SUBJ0073_SER001_FILESNO001.fif';
maxfName = 'SUBJ0073_SER001_musmelo'; % define the output name of the maxfilter processing
%no movement compensation
cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,num2str(project), ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsss.fif'] ' -st 4 -corr 0.98 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsss.log"']];
system(cmd);
% 
% rawName = '/raw/sorted/MINDLAB2015_MEG-TunteetMM/0073/20150702_000000/MEG/002.MusMelo_1/files/PROJ0192_SUBJ0073_SER002_FILESNO001.fif';
% maxfName = 'SUBJ0073_SER002_musmelo'; % define the output name of the maxfilter processing
% %no movement compensation
% cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,num2str(project), ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsss.fif'] ' -st 4 -corr 0.98 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsss.log"']];
% system(cmd);

rawName = '/raw/sorted/MINDLAB2015_MEG-TunteetMM/0094/20150702_000000/MEG/001.MusMelo_2/files/PROJ0192_SUBJ0094_SER001_FILESNO001.fif';
maxfName = 'SUBJ0094_SER001_musmelo'; % define the output name of the maxfilter processing
%no movement compensation
cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,num2str(project), ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsss.fif'] ' -st 4 -corr 0.98 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsss.log"']];
system(cmd);

% rawName = '/raw/sorted/MINDLAB2015_MEG-TunteetMM/0094/20150702_000000/MEG/002.MusMelo_1/files/PROJ0192_SUBJ0094_SER002_FILESNO001.fif';
% maxfName = 'SUBJ0094_SER002_musmelo'; % define the output name of the maxfilter processing
% %no movement compensation
% cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,num2str(project), ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsss.fif'] ' -st 4 -corr 0.98 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsss.log"']];
% system(cmd);
                  
%%
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %adding the path to OSL functions
osl_startup %starting the osl package

%%  CONVERTING THE .FIF FILES IN Spm OBJECTS

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster') %add the path where is the function for submit the jobs to the server
fif_list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/*fif'); %creating a list with the fif files outputted by Maxfilter
% fif_list ='/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/SUBJ0003_musmelo_tsss.fif'; %creating a list with the fif files outputted by Maxfilter
% fif_list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/SUBJ*SER*musmelo_tsss*'); %creating a list with the fif files outputted by Maxfilter

clusterconfig('long_running', 1);
clusterconfig('slot', 1);
for ii = 1:length(fif_list) %over the subjects
    S = [];                    
    S.dataset = [fif_list(ii).folder '/' fif_list(ii).name]
%     S.dataset = '/raw/sorted/MINDLAB2015_MEG-TunteetMM/0093/20150702_000000/MEG/001.MusMelo_1_2/files/PROJ0192_SUBJ0093_SER001_FILESNO001.fif'; %building the path for each subject
%     D = spm_eeg_convert(S);
    D = job2cluster (@cluster_spmobject, S); %actual function for conversion
end

%% REMOVING BAD TIME PERIODS USING OSLVIEW

spm_list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/spmeeg_SUBJ*SER*musmelo_tsss.mat*');

for ii = 1:length(spm_list) % iterates over experimental blocks %OBS!!!!!
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = oslview(D);
    D.save(); %it is needed to save the marked bad events (and/or channels) in oslview
    disp(ii)
end

%% AFRICA denoising

%clearvars -except datadir workingdir temp1fif spm_files spm_files_basenames

%%% OBS!! to be ran only the first time %%%
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 0)
clusterconfig('slot', 2);
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Francesco/Musmelo')
for jj = 1:length(spm_list)
    S = [];
    S.jj = jj;
    S.spm_list = spm_list;
    D = spm_eeg_load([spm_list(jj).folder '/' spm_list(jj).name]);
    S.D = D;
    jobid = job2cluster(@name_africa_cluster, S);  
end

% for ii = 1:length(v) 
for ii = 1:length(spm_list) %:length(spm_files) %OBS!!!!! 
    S = [];
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    S.D = D;
    jobid = job2cluster(@cluster_africa,S);
%     D = osl_africa(D,'do_ica',true,'do_ident',false,'do_remove',false,'used_maxfilter',true); 
%     D.save();
end
% v = [11 12 19 32];
%visual inspection and removal of artifacted components
for ii = 1:length(spm_list)%:2:length(spm_list) %:length(spm_files) %OBS!!!!!
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = osl_africa(D,'do_ident',true,'do_remove',false);
    %     hacking the function to manage to get around the OUT OF MEMORY problem..
    S = [];
    S.D = D;
    jobid = job2cluster(@cluster_rembadcomp,S);
%     D.save();
    disp(spm_list(ii))
end

%% epoching

spm_list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after max_filter/s*mat*');
subj = [26 38 42 77 110]; % list of subj with errors, the id correspond to the spm_list position and not to the subj num
prefix = 'e2';
% for ii = 2:2:length(spm_list)
% for ii = 1:length(spm_list)
for ii = 1:length(subj)
    S2 = [];
%     S2.D= [spm_list(ii).folder '/' spm_list(ii).name];
    S2.D= [spm_list(subj(ii)).folder '/' spm_list(subj(ii)).name]; % using the list of subj with errors
    D_continuous=spm_eeg_load(S2.D);
    D_continuous=D_continuous.montage('switch',0);
    
    pretrig = -100; % epoch start in ms
    posttrig = 600; % epoch end in ms   
    S2.timewin = [pretrig posttrig];
    S2.D_continuous = D_continuous;

    % event definitions
    S2.trialdef(1).conditionlabel = 'Standard_mr'; % standard for melody and rhythm mod. condition
    S2.trialdef(1).eventtype = 'STI 014_up';
    S2.trialdef(1).eventvalue = [2 3];
    S2.trialdef(2).conditionlabel = 'melody mod';
    S2.trialdef(2).eventtype = 'STI 014_up';
    S2.trialdef(2).eventvalue = [6 13];
    S2.trialdef(3).conditionlabel = 'rhythm mod';
    S2.trialdef(3).eventtype = 'STI 014_up';
    S2.trialdef(3).eventvalue = 14;
    S2.trialdef(4).conditionlabel = 'transposition';
    S2.trialdef(4).eventtype = 'STI 014_up';
    S2.trialdef(4).eventvalue = 5;
    S2.trialdef(5).conditionlabel = 'mistuning';
    S2.trialdef(5).eventtype = 'STI 014_up';
    S2.trialdef(5).eventvalue = 7;
    S2.trialdef(6).conditionlabel = 'timbre deviant';
    S2.trialdef(6).eventtype = 'STI 014_up';
    S2.trialdef(6).eventvalue = [8 9];
    S2.trialdef(7).conditionlabel = 'rhythm mistake';
    S2.trialdef(7).eventtype = 'STI 014_up';
    S2.trialdef(7).eventvalue = 12;
    S2.trialdef(8).conditionlabel = 'Standard_t'; % standard for transposition condition
    S2.trialdef(8).eventtype = 'STI 014_up';
    S2.trialdef(8).eventvalue = 1;
    S2.trialdef(9).conditionlabel = 'Standard_m'; % standard for mistuning condition
    S2.trialdef(9).eventtype = 'STI 014_up';
    S2.trialdef(9).eventvalue = 3;
    S2.trialdef(10).conditionlabel = 'Standard_tr'; % standard for timbre and rhythm mistake condition
    S2.trialdef(10).eventtype = 'STI 014_up';
    S2.trialdef(10).eventvalue = [3 4];
    
    %S2.trl = data.trialinfo;
    S2.prefix = prefix;
    S2.reviewtrials = 0;
    S2.save = 0;
    S2.epochinfo.padding = 0;
    S2.event = D_continuous.events;
    S2.fsample = D_continuous.fsample;
    S2.timeonset = D_continuous.timeonset;
    
    jobid = job2cluster(@cluster_epoch_osl,S2); %running with parallel
    disp(ii)
end

%% marking bad trials as 'bad trials'

spm_list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after max_filter/spmeeg_SUBJ*_musmelo_tsss.mat*');

for ii = 1:length(spm_list) %indexing only the files wanted for this paper
    D2 = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]); %loading continuous data with proper D2.event field
%     D2 = spm_eeg_load([spm_list(V(ii)).folder '/' spm_list(V(ii)).name]); %loading continuous data with proper D2.event field
    events = D2.events;
    %recreating proper event list..
    clear trl_sec
    ck = 0;
    for kkk = 1:length(events) %over events
        if strcmp(events(kkk).type,'STI 014_up') %with only STI 014_up
            ck = ck + 1;
            trl_sec(ck,1) = events(kkk).time;
            trl_sec(ck,2) = events(kkk).time + 0.3996;
        end
    end
    D = spm_eeg_load([spm_list(ii).folder '/e' spm_list(ii).name]); %loading epoched data
%     D = spm_eeg_load([spm_list2(V(ii)).folder '/' spm_list2(V(ii)).name]); %loading epoched data
    D = D.montage('switch',1);
    %take bad segments registered in oslview and check if they overlap with the trials. if so, it gives the number of overlapped trials that will be removed later   
    pd = 0;
    Bad_trials = zeros(size(trl_sec,1),1);
    for kkk = 1:length(events) %over events
        if strcmp(events(kkk).type,'artefact_OSL')
            for k = 1:length(trl_sec) %over trials
                if events(kkk).time - trl_sec(k,2) < 0 %if end of trial is > than beginning of artifact
                    if trl_sec(k,1) < (events(kkk).time + events(kkk).duration) %if beginning of trial is < than end of artifact
                        Bad_trials(k,1) = 1; %it is a bad trial (stored here)
                        pd = pd + 1;
                    end
                end                  
            end
        end
    end
    %if bad trials were detected, their indices are stored within D.badtrials field
    if pd == 0
        disp('there are no bad trials marked in oslview for');
        disp(ii);
    else
%         D = badtrials(D,1:size(D,3),0); %get the indices of the badtrials marked as '1' (that means bad)
        D = badtrials(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
%         D = conditions(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
        epochinfo = D.epochinfo;
        xcv = find(Bad_trials == 1);
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
end

%% MEG SENSOR LEVEL %%

%% averaging and combining planar gradiometers (parallel computing)

list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/espmeeg_SUBJ*_musmelo_tsss.mat*');
pd= 0;
for ii= 1:length(list)
    if list(ii).name(18) ~= 'S'
       pd = pd+1; 
       list_c(pd,1) = list(ii);
    end
end
list = list_c;
clear list_c;
%settings for cluster (parallel computing)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')
clusterconfig('slot', 2); % set manually the job cluster slots
% between 1 and 12 (n x 8gb of ram)
clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
% clusterconfig('scheduler', 'cluster'); % set automatically the long run queue
clusterconfig('long_running', 0); % set automatically the long run queue

%averaging
% addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Francesco/Functions'); %path to LEiDA_MEG_leonardo functions
% output_prefix_to_be_set = 'm'; 
output_prefix_to_be_set = 'm2'; 

for ii = 1:length(list) %over files
    %distribute 
    input = [];
    input.D = [list(ii).folder '/' list(ii).name];
    input.prefix = output_prefix_to_be_set;
    jobid = job2cluster(@sensor_average, input); % this is the command for send the job to the cluster, in the brackets you can find the name on the function to run (afeter the @) and the variable for the input (in this case input)
    % look the script for more details about the function work
end
%combining planar gradiometers
% list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after max_filter/me*.mat*');
list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after max_filter/m2e*.mat*');

for ii = 1:length(list) %over files
    input = [];
    input.D = [list(ii).folder '/' list(ii).name];
    D = spm_eeg_load(input.D);
    D = D.montage('switch',1);
    D.save();
    jobid = job2cluster(@combining_planar_cluster, input); % this is the command for send the job to the cluster, in the brackets you can find the name on the function to run (afeter the @) and the variable for the input (in this case input)
end

%% LBPD_startup_D

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster') %add the path where is the function for submit the jobs to the server
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);

%% extracting MEG sensor data

% label_lowfreq = 0; %set this if you want to take the filtered data below 3 Hz
% block = 0; %set the experimental block (1 = ato; 2 = maj; 3 = min; 0 = average!)
load_data = 1; %set 1 if you want to load the data instead of extracting it from SPM objects
% subjnum = {'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';'30';'31';'32';'33';'34';'35';'36';'37';'38';'40';'41';'42';'43';'44';'45';'46';'47';'48';'49';'50';'51';'52';'53';'54';'55';'56';'57';'58';'59';'60';'61';'62';'63';'64';'65';'66';'67';'68';'69';'70';'71'}; %Everybody!!!!!!
datadir = '/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after max_filter'; %path to data
outdir = '/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/megsensorresult'; %path to write output in
% outdir = '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/PhDThesis/Paper1/MCS_Sensor/lowpass40Hz/WM_RTs_Complexity'; %path to write output in
% blocklist = {'recogatonal';'recogmajor';'recogminor'}; %this is simply to get proper data among 3 different options
list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/Pm2e*.mat*');
S = [];
%computing data
S.outdir = outdir;
% S.spm_list = spm_list;
S.data = [];
if load_data == 1 %if you already computed and saved on disk the t-tests you can load them here
     load([outdir '/sensor_data_nostandard2.mat']);
     data_mat = PD;
%      load([outdir '/sensor_data_jan2021.mat']);
      
%     load([outdir '/sensor_data_also_uncorrect.mat']);
    S.data = data_mat;
    S.chanlabels = chanlabels;
    S.time_real = time_sel;
else %otherwise you can extract the data from SPM MEEG objects (one for each subject)
    S.spm_list = cell(1,length(list));
    for ii = 1:length(list)
        
        S.spm_list(ii) = {[list(ii).folder '/' list(ii).name]};
        
    end
end
% S.conditions = {'Old_Correct','New_Correct','Old_Uncorrect','New_Uncorrect'};
% S.conditions = {'Standard','melody mod','rhythm mod','transposition','mistuning','timbre deviant','rhythm mistake' };
S.conditions = {'melody mod','rhythm mod','transposition','mistuning','timbre deviant','rhythm mistake' };
% S.conditions = {'Standard_mr','melody mod','rhythm mod','transposition','mistuning','timbre deviant','rhythm mistake','Standard_t','Standard_m','Standard_tr'};
%assigning colors
color_line = colormap(lines(6)); %extracting some colours from a colormap
cl2 = color_line;
cl2(2,:) = color_line(1,:);
cl2(1,:) = color_line(2,:);
cl2(4,:) = color_line(3,:);
cl2(3,:) = color_line(4,:);
cl2(5,:) = [0 0.7 0];
S.color_line = cl2;

% S.timeextract = [286:751]; %time-points to be extracted
S.timeextract = []; %time-points to be extracted
S.centerdata0 = 0; %1 to make data starting at 0
S.save_data = 0; %only meaningfull if you read data from SPM objects saved on disk
% S.save_name_data = 'sensor_data';
%S.save_name_data = 'sensor_data_jan2021';
% S.save_name_data = 'sensor_data_lowpass3Hz_shortbaseline';
%individual waveform plotting
S.waveform_singlechannels_label = 0; %1 to plot single channel waveforms
S.wave_plot_conditions_together = 0; %1 for plotting the average of all
S.mag_lab = 1; %1 for magnetometers; 2 for gradiometers
S.x_lim_temp_wave = []; %limits for time (in secs) (E.g. [-0.1 3.4])
S.y_lim_ampl_wave = []; %limit for amplitude (E.g. [0 120] magnetometes, [0 6] gradiometers)
%averaged waveform plotting
S.waveform_average_label = 1;
% S.left_mag = [2:2:204];
S.legc = 1; %set 1 for legend

S.left_mag = 99;
S.signtp = {[0.05 0.22]};
% S.sr = 150; %sampling rate (Hz)
S.avewave_contrast = 0; %1 to plot the contrast between conditions (averaged waveform)
S.save_label_waveaverage = 0;
S.label_plot = 'block_minor';
%t-tests
S.t_test_for_permutations = 0;
S.cond_ttests_tobeplotted_topoplot = [2 1]; %this is for both topoplot and t-tests!! (here [1 2] means cond1 vs cond2!!!!!!!)
%topoplotting
S.xlim = [0.12 0.13]; %time topolot (cluster II)
S.topoplot_label = 0;
S.fieldtrip_mask = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External';
S.topocontr = 1;
S.topocondsing = [4];
% S.xlim = [0.913 1.160]; %time topolot (cluster I)
S.zlimmag = []; %magnetometers amplitude topoplot li  mits
S.zlimgrad = []; %gradiometers amplitude topoplot limits
S.colormap_spec = 0;
% x = []; x.bottom = [0 0 1]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 1 0.5]; x.top = [1 0.95 0]; %yellow - blue
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
S.colormap_spec_x = x;
S.topoplot_save_label = 0;

[out] = MEG_sensors_plotting_ttest_LBPD_D2(S);

%% BONFERRONI CORRECTION ON T-VALUES OBTAINED BY CONTRASTING DEVIANTS VS STANDARDS

min_time_point = 61; %16 = 0 seconds (first 15 points are pre-stimulus time)
max_time_point = 242;
time_point = max_time_point - min_time_point;
pvalue = 0.05;
ndeviant = 6; 

bonf_value = pvalue/(time_point* 6 * 102);

%% Bonferroni on p-value matrix

list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/megsensorresult/sensor_data_all*OUT*.mat');
load([list(1).folder '/' list(1).name]);
xx = size(OUT.TSTATP_grad);
Bonf = zeros(size(OUT.TSTATP_grad));
for ii = 1:length(list)
    load([list(1).folder '/' list(ii).name]);
    for jj = 1:xx(1)
        for aa = 1:xx(2)
            if OUT.TSTATP_grad(jj,aa) <= bonf_value
                Bonf(jj,aa) = 1;
            else
                Bonf(jj,aa) = 0;
            end
        end
        
    end
    %out = (['matxx' num2str(ii) '.mat']);
    save(['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/Paper/Bonf_dev_' list(ii).name(51:64)],'Bonf');
end

%% Removing standards from deviants to obtain MMN for later plotting purposes

load('sensor_data.mat');
D1 = data_mat;
PD = zeros(size(data_mat,1),size(data_mat,2),size(data_mat,3),size(data_mat,4)-4);
v = [1 1 8 9 10 10];
for ii = 1:6
    PD(:,:,:,ii) = D1(:,:,:,ii+1) - D1(:,:,:,v(ii));
end
save('sensor_data_nostandard2','PD')

%%

%% MEG SOURCE LEVEL %% 

%% getting co-registered data (MRI and MEG) from previously computed co-registrations for the same subjects

listsilvia = dir('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/es*mat');
listsbrency = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after max_filter/es*mat');
prc = 0;
for ii = 1:length(listsbrency)
    str = listsbrency(ii).name(1:16);
    for kk = 1:length(listsilvia)
        if strfind(listsilvia(kk).name, str)
            paths = [listsilvia(kk).folder '/' listsilvia(kk).name];
            Dsilvia = spm_eeg_load (paths);
            inv = Dsilvia.inv;
            pathf = [listsbrency(ii).folder '/' listsbrency(ii).name];
            Dsbrency = spm_eeg_load (pathf);
            Dsbrency.inv = inv;
            Dsbrency.save();
            prc = 1;
            break;
        else
            prc = 0;
        end
    end
    if prc == 0
        disp(listsbrency(ii).name)
    end
end

%% back to before AFRICA (with labels..)

%changing labels for compatibility reasons..
     
D2 = spm_eeg_load('/scratch5/MINDLAB2017_MEG-LearningBach/Portis/DIsp/spmeeg_SUBJ0002_mumufe_tsss.mat');
spm_list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/es*mat');

for ii = 20:length(spm_list)%60%4:2:length(spm_list)
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    S2 = [];
    S2.D = D; %file to be updated progressively
    S2.D2 = D2; %original labels
    jobid = job2cluster(@cluster_Dlabel,S2); %running with parallel
end

%% proper codes for source reconstruction using OAT (OHBA Oxford)
    
spm_list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/es*mat');

v = [2 10 11 12 14 33 61 62 80 81];
for ii = 1:length(spm_list) %6:2:length(spm_list)
    if isempty(find(ii == v))
        processed_file = [spm_list(ii).folder '/' spm_list(ii).name];    
        % Beamform
        oat = [];
        oat.source_recon.D_epoched(str2double(spm_list(ii).name(13:16)))         = {processed_file};
        %     oat.source_recon.D_epoched(str2double(spm_list(v(ii)).name(14:17)))         = {processed_file};
        pca_dim_1 = 50;
        oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_epoched),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects
        oat.source_recon.modalities = {'MEGMAG'; 'MEGPLANAR'};
        oat.source_recon.conditions        = {'Standard_mr','melody mod','rhythm mod','transposition','mistuning','timbre deviant','rhythm mistake','Standard_t','Standard_m','Standard_tr'};
        oat.source_recon.gridstep          = 8; % in mm
        oat.source_recon.time_range        = [-0.1 0.39]; % time range in secs
        oat.source_recon.freq_range        = [0.1 40]; % frequency range in Hz
        %     oat.source_recon.freq_range        = [0.1 10]; % frequency range in Hz
        %S.source_recon.pca_order         = 250;
        oat.source_recon.type              = 'Scalar';
        oat.source_recon.method            = 'beamform';
        oat.source_recon.normalise_method  = 'mean_eig';
        oat.source_recon.forward_meg       = 'Single Shell';
        %S.source_recon.prefix            = '';
        oat.source_recon.report.do_source_variance_maps = 1;
        oat.source_recon.sessions_to_do = [];
        oat.source_recon.sessions_to_do = [str2double(spm_list(ii).name(13:16))]; %sessions to do among the file_list (subject 39 excluded because of problems during data collection..)
        %     oat.source_recon.sessions_to_do = [str2double(spm_list(v(ii)).name(14:17))]; %sessions to do among the file_list (subject 39 excluded because of problems during data collection..)
        oat.source_recon.dirname           = [spm_list(1).folder '/source/' spm_list(ii).name(13:16)];
        %     oat.source_recon.dirname           = [spm_list(1).folder '/sourcetryindsubj600Hz/' spm_list(v(ii)).name(14:18)];
        
        jobid = job2cluster(@cluster_beamforming,oat); %running with parallel
    end
end

%% moving files computed independently for each subject into a common new folder

%this is a trick to be able to use the parallel computing of Aarhus cluster of computers even with OAT which would not be designed for such an environment

%%% TO BE RUN ONLY THE FIRST TIME %%%

list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/0*');
bs = '/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat';

for ii = 1:4%length(list)
    %concat file .dat
    as = ['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/' list(ii).name '/concatMfsession' num2str(str2double(list(ii).name(1:4))) '_spm_meeg.dat'];
    status = movefile(as,bs)
    %concat file .mat
    as = ['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/' list(ii).name '/concatMfsession' num2str(str2double(list(ii).name(1:4))) '_spm_meeg.mat'];
    status = movefile(as,bs)
    %session file
    as = ['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/' list(ii).name '/session' num2str(str2double(list(ii).name(1:4))) '_recon.mat'];
    status = movefile(as,bs)
end

%% codes to get a "normal" oat for source reconstruction (so all subjects in the same folder)..

spm_list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/es*mat');

v = [2 10 11 12 14 33 61 62 80 81];
oat = [];
cnt = 0;
for ii = 1:length(spm_list)
    if isempty(find(ii == v))
        cnt = cnt+1;
    processed_file = [spm_list(ii).folder '/' spm_list(ii).name];
    oat.source_recon.D_epoched(cnt)         = {processed_file};
    oat.source_recon.sessions_to_do(cnt) = {[str2double(spm_list(ii).name(13:16))]}; %sessions to do among the file_list (subject 39 excluded because of problems during data collection..)
    oat.source_recon.results_fnames(cnt) = {['session' num2str(str2double(spm_list(ii).name(13:16))) '_recon']};
    
    %     D_epoched = spm_eeg_load(processed_file);
    %     D_epoched = D_epoched.montage('switch',1); %switch the montage to 1 in order to be safer (so we have the AFRICA denoised data)
    %     D_epoched.save();
    end
end
% Beamform
% oat.source_recon.D_epoched(str2double(spm_list(ii).name(13:16)))         = {processed_file};

pca_dim_1 = 50;
oat.source_recon.pca_dim = pca_dim_1(ones(length(oat.source_recon.D_epoched),1)); %this seems to be necessary. It is for having the same number of values than the number of input D objects
oat.source_recon.modalities = {'MEGMAG'; 'MEGPLANAR'};
oat.source_recon.conditions        = {'Standard_mr','melody mod','rhythm mod','transposition','mistuning','timbre deviant','rhythm mistake','Standard_t','Standard_m','Standard_tr'};
oat.source_recon.gridstep          = 8; % in mm
oat.source_recon.time_range        = [-0.1 0.39]; % time range in secs
oat.source_recon.freq_range        = [0.1 40]; % frequency range in Hz
%     oat.source_recon.freq_range        = [0.1 10]; % frequency range in Hz
%S.source_recon.pca_order         = 250;
oat.source_recon.type              = 'Scalar';
oat.source_recon.method            = 'beamform';
oat.source_recon.normalise_method  = 'mean_eig';
oat.source_recon.forward_meg       = 'Single Shell';
%S.source_recon.prefix            = '';
oat.source_recon.report.do_source_variance_maps = 1;
%     oat.source_recon.sessions_to_do = [str2double(spm_list(v(ii)).name(14:17))]; %sessions to do among the file_list (subject 39 excluded because of problems during data collection..)
oat.source_recon.dirname           = [spm_list(1).folder '/source/Firstlevel'];

%% first level (each experimental block for each subject, independently)

%FIRST LEVEL
design_matrix_summary = {};
design_matrix_summary{1} = [1 0 0 0 0 0 0 0 0 0];design_matrix_summary{2} = [0 1 0 0 0 0 0 0 0 0]; design_matrix_summary{3}=[0 0 1 0 0 0 0 0 0 0];design_matrix_summary{4}=[0 0 0 1 0 0 0 0 0 0]; design_matrix_summary{5}=[0 0 0 0 1 0 0 0 0 0]; design_matrix_summary{6}=[0 0 0 0 0 1 0 0 0 0]; design_matrix_summary{7}=[0 0 0 0 0 0 1 0 0 0]; design_matrix_summary{8}=[0 0 0 0 0 0 0 1 0 0];design_matrix_summary{9}=[0 0 0 0 0 0 0 0 1 0];design_matrix_summary{10}=[0 0 0 0 0 0 0 0 0 1];
oat.first_level.design_matrix_summary = design_matrix_summary;
% contrasts to be calculated:
oat.first_level.contrast = {};
% contrast design matrix
oat.source_recon.conditions        = {'Standard_mr','melody mod','rhythm mod','transposition','mistuning','timbre deviant','rhythm mistake','Standard_t','Standard_m','Standard_tr'};
% oat.first_level.contrast{1} = [1 0 0 0 0 0 0]'; %Standard
% oat.first_level.contrast{2} = [0 1 0 0 0 0 0]'; %Pitch
% oat.first_level.contrast{3} = [0 0 1 0 0 0 0]'; %Timbre
% oat.first_level.contrast{4} = [0 0 0 1 0 0 0]'; %Localization
% oat.first_level.contrast{5} = [0 0 0 0 1 0 0]'; %Intensity
% oat.first_level.contrast{6} = [0 0 0 0 0 1 0]'; %Slide
% oat.first_level.contrast{7} = [0 0 0 0 0 0 1]'; %Rhythm
oat.first_level.contrast{1} = [-1 1 0 0 0 0 0 0 0 0]'; %melody mod - Standard_mr
oat.first_level.contrast{2} = [-1 0 1 0 0 0 0 0 0 0]'; %rhythm mod - Standard_mr
oat.first_level.contrast{3} = [0 0 0 1 0 0 0 -1 0 0]'; %transposition - Standard_t
oat.first_level.contrast{4} = [0 0 0 0 1 0 0 0 -1 0]'; %mistuning - Standard_m
oat.first_level.contrast{5} = [0 0 0 0 0 1 0 0 0 -1]'; %timbre deviant - Standard_tr
oat.first_level.contrast{6} = [0 0 0 0 0 0 1 0 0 -1]'; %rhythm mistake - Standard_tr

%contrast names
oat.first_level.contrast_name = {};
% oat.first_level.contrast_name{1} = 'Standard';
% oat.first_level.contrast_name{2} = 'Pitch';
% oat.first_level.contrast_name{3} = 'Timbre';
% oat.first_level.contrast_name{4} = 'Localization';
% oat.first_level.contrast_name{5} = 'Intensity';
% oat.first_level.contrast_name{6} = 'Slide';
% oat.first_level.contrast_name{7} = 'Rhythm';
oat.first_level.contrast_name{1} = 'melody mod - Standard_mr';
oat.first_level.contrast_name{2} = 'rhythm mod - Standard_mr';
oat.first_level.contrast_name{3} = 'transposition - Standard_t';
oat.first_level.contrast_name{4} = 'mistuning - Standard_m';
oat.first_level.contrast_name{5} = 'timbre deviant - Standard_tr';
oat.first_level.contrast_name{6} = 'rhythm mistake - Standard_tr';

%original contrasts, try1 and try2 seem to give the exact same results
oat.first_level.report.first_level_cons_to_do = [1 2 3 4 5 6]; %[1 2 3]; %better to do 3 2 1 in order to get the information for the peak value of the contrast old-new
oat.first_level.time_range = [-0.1 0.39];
oat.first_level.post_tf_downsample_factor = 1;
%slow negativity
oat.first_level.cope_type = 'coape'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
%N100
% oat.first_level.cope_type = 'none'; %this should calculate absolute values of glm copes before comparing different conditions (otherwise, by default the absolute value is done on the copes after comparing different conditions); the absolute value here seems to be necessary because of the sign ambiguity issue introduced by source reconstruction
oat.first_level.name = ['wholebrain_first_level_BC_coape']; %REMEMBER TO CHECK THIS NAME!!
oat.first_level.bc = ones(1,6);
%to add if the oat has not been automatically saved
list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/0*');
for ii = 1:length(list)
%     oat.first_level.results_fnames(ii) = {['session' num2str(str2double(list(ii).name(1:4))) '_wholebrain_first_level_BC.mat']};
    oat.first_level.results_fnames(ii) = {['session' num2str(str2double(list(ii).name(1:4))) '_wholebrain_first_level_BC_coape.mat']};
end
% oat.first_level.sessions_to_do = [1:38 40:71]; %this is probably not necessary but it may be useful to report if you run only subject_level analysis in oat..
% oat.first_level.sessions_to_do = [2 5 6 7 9 13 15 18:21 24 25 27:33 36:42 44:46 48 53 55 56 58 60 63 65 68 69 71]; %this is probably not necessary but it may be useful to report if you run only subject_level analysis in oat..
% a = [2 5 6 7 9 13 15 18:21 24 25 27:33 36:42 44:46 48 53 55 56 58 60 63 65 68 69 71]; %this is probably not necessary but it may be useful to report if you run only subject_level analysis in oat..

%% running first level on parallel computing

for ii = 1:length(list)
    oat.first_level.sessions_to_do = [];
    oat.first_level.sessions_to_do = [ii]; %here it seems that the session indexes oat.source_recon.results_fnames{ii} is directly related to the sequential 
%     oat.first_level.sessions_to_do = [str2double(list(ii).name(1:4))];
    jobid = job2cluster(@cluster_beamfirstlevel,oat);
end

%% SUBJECT LEVEL

%in this case it does not do anything since I have only one experimental
%block for each participant (however, for computational reasons, I have to
%run it)
%this is needed to read the proper subjects..
for ii = 1:140
%     oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level_BC.mat']};
    oat.first_level.results_fnames(ii) = {['session' num2str(ii) '_wholebrain_first_level_BC_coape.mat']};
end
% this is if you have a perfect correspondance between sessions and subjects
oat.subject_level.session_index_list = cell(1,140); %sarebbe length(subjects)
oat.subject_level.name = 'MMN';
oat.subject_level.subjects_to_do = [];
% oat.subject_level.subjects_to_do = [5:115];

%to update the names of the results files of subject level
for ii = 1:140
    oat.subject_level.session_index_list{ii} = ii;
    oat.subject_level.results_fnames(ii) = {['subject' num2str(ii) '_wholebrain_first_level_BC_coape_MMN.mat']};
end

for ii = 1:length(list)
    oat.subject_level.subjects_to_do = zeros(1);
    oat.subject_level.subjects_to_do(1) = str2double(list(ii).name(1:4));
    jobid = job2cluster(@cluster_subjlevel,oat);
end

%% GROUP LEVEL

for ii = 1:length(list)
    oat.subject_level.subjects_to_do(ii) = str2double(list(ii).name(1:4));
end
oat.group_level = [];
oat.group_level.name = 'group_level_everybody_BC'; %OBS!! REMEMBER TO UPDATE THE NAME!
oat.group_level.subjects_to_do = [];
% oat.group_level.subjects_to_do = 1:70; %this seems to be necessary since the function gets the indices for the group level from each progressive number of oat.subject_level.subjects_to_do (not very clever idea by the way..). Therefore not to mix up things, you should have here a vector from 1 to the total number of subjects and then specify the correct subject IDs that you want.. or at least, to me it looks like that..
oat.group_level.subjects_to_do = 1:104; %this seems to be necessary since the function gets the indices for the group level from each progressive number of oat.subject_level.subjects_to_do (not very clever idea by the way..). Therefore not to mix up things, you should have here a vector from 1 to the total number of subjects and then specify the correct subject IDs that you want.. or at least, to me it looks like that..
%results name
% oat.group_level.results_fnames = ['wholebrain_first_level_BC_MMN' '_' oat.group_level.name '.mat'];
oat.group_level.results_fnames = ['wholebrain_first_level_BC_coape_MMN' '_' oat.group_level.name '.mat'];
% Spatial and temporal averaging options
oat.group_level.time_range = [-0.1 0.39];
oat.group_level.space_average = 0;
oat.group_level.time_average = 0;
oat.group_level.time_smooth_std = 0; % secs
oat.group_level.use_tstat = 0;
%path to AAL template
% oat.group_level.mask_fname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/AAL_Nifti/aal_8mm_try5.nii.gz';
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
oat.group_level.glm_method='fixed_effects'; %ols or fixed-effects
% Define which contrasts to perform for the report
oat.group_level.first_level_contrasts_to_do = [1,2,3,4,5,6]; % list of first level contrasts to run the group analysis on
oat.group_level.report.first_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.group_level_cons_to_do = [1]; % purely used to determine which contrasts are shown in the report, the first one in list determines the contrast used to find the vox, time, freq with the max stat
oat.group_level.report.show_lower_level_copes = 0;
oat.group_level.report.show_lower_level_cope_maps = 0;

jobid = job2cluster(@cluster_beamgrouplevel,oat);

%% creatig nifti images with statistics..

load('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/oat_wholebrain_first_level_BC_coape_MMN_group_level_everybody_BC.mat');

S2 = [];
S2.oat = oat;
S2.stats_fname = oat.group_level.results_fnames;
S2.first_level_contrasts = 1:13; %remember that in this way you define the order of the output (this numbers refers to the order (numbers) defined in the contrasts; for example here tstat3 refers to the contrast old-new)
S2.group_level_contrasts = [1 2];
S2.resamp_gridstep = oat.source_recon.gridstep;

jobid = job2cluster(@cluster_oat_save_nii_stats,S2);

%% p (q) values through permutations

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster') %add the path where is the function for submit the jobs to the server
clusterconfig('slot', 3); % set manually the job cluster slots
clusterconfig('long_running', 1); % set automatically the long run queue

for ii = 1:6 %loop for all the condition
    %extracting the peaks of the MMNs and using it (plus/minus 20 ms) to decide the time-points to be used for our analyses
    [vl,ivl] = min(PD(99,61:242,ii));
    tm_s = time_sel(ivl+61);
    ts_x = tm_s - 0.025; % min
    ts_m = tm_s + 0.025; %max
    load('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/oat_wholebrain_first_level_BC_coape_MMN_group_level_everybody_BC.mat');
    
    S = [];
    S.oat = oat;
    S.cluster_stats_thresh = 3.3; % p-value = 0.01/6(deviants) = 0.0017 (corresponding t-value = 3.3 )
    S.cluster_stats_nperms = 5000; % we normally recommend doing 5000 perms
    S.first_level_copes_to_do = [ii];
    S.group_level_copes_to_do = [1];
    S.group_varcope_spatial_smooth_fwhm = S.oat.group_level.group_varcope_spatial_smooth_fwhm;
    S.write_cluster_script = 0;
    S.time_range = [ts_x ts_m];
    S.time_average=1;
    
    % Run the permutations
    % [ gstats ] = oat_cluster_permutation_testing(S);
    % Run the permutations (on cluster.. parallel computing)
    jobid = job2cluster(@clusterbasedpermutation_osl,S);
end

%%
























%% AVERAGING TOGETHER THE 3 HIGH DEVIANTS AND THE 3 LOW DEVIANTS

%HERE SOME EXAMPLES FOR USING NIFTI IMAGES - FROM LEONARDO TO FRANCESCO

lmt = [134,171,159,175,143,109];% signal peaks for each deviants
list2 = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/subject*dir'); %creating a list with the fif files outputted by Maxfilter

% creating empty structures
dp = zeros(length(asd(:,1,1)),length(asd(1,:,1)),length(asd(1,1,:)),length(list2),3);
dp2 = zeros(length(asd(:,1,1)),length(asd(1,:,1)),length(asd(1,1,:)),length(list2),3);

for ii = 1:length(list2) % loop for all the subjects 
    dirn = ([list2(ii).folder '/' list2(ii).name]); % building path 
    path = dir([dirn '/' 'tstat*.gz']); %taking only the tstat files
    pd = 1; 
    for bb = 1:2:length(path) % taking only the non-mip files
        pathf(pd) = path(bb); % creating a temporary path
        pd = pd + 1; % increase pd for changing the path
    end
    for kk = 1:length(pathf)
        fname = ([pathf(kk).folder '/' pathf(kk).name]); % build the path name for the gz file
        %loading the nifti image (fname is the path plus the name of the nifti file)
        aalnii = load_nii(fname); %aalnii becomes a structure with several fields
        img = aalnii.img; %.img is a key field, the one describing the image (it is a matrix usually in 3D or 4D space)
        asd = mean(img(:,:,:,lmt(kk)-25:lmt(kk)+25),4); % reduction of the dimensionality from a 4D to a 3D space
        %to save a nifti image
        if kk <=3 % splitting the High and low Deviants in two different matrix for all the subject and deviants
            dp(:,:,:,ii,kk) = asd(:,:,:);
        else
            dp2(:,:,:,ii,kk-3) = asd(:,:,:);
        end
    end
    disp(ii)
end
save('mean_allsubj_devhigh.mat', 'dp')
save('mean_allsubj_devlow.mat', 'dp2')

% reduction of the dimensionality from a 5D to a 4D (all the 3D coordinates plus all the subjects)
dpm1 = mean(dp,5);
dpm2 = mean(dp2,5);

%% T-TESTS CONTRASTING AVERAGED HIGH VS LOW DEVIANTS

% creating 2 matrices for the t and p value 
sz = size(dpm1); %taking the same dimension from the saved matrix without the subj dimensionality
T = zeros(sz(1),sz(2),sz(3)); % t value matrix
P = zeros(sz(1),sz(2),sz(3)); % p value matrix

for ii = 1:sz(1) %making loop for every dimension
    for yy = 1:sz(2)
        for zz = 1:sz(3)
            if dpm1(ii,yy,zz,1) ~= 0 % do the test for the coordinate only if is different to 0
                a = squeeze(dpm1(ii,yy,zz,:)); 
                b = squeeze(dpm2(ii,yy,zz,:)); 
                [h,p,ci,stats] = ttest(a,b); % simple t-test with: hypothesis test result, p value, confidence interval and test statistics
                T(ii,yy,zz) = stats.tstat; % store all the t test stats inside a new matrix
                P(ii,yy,zz) = p; % store all the p value inside a new matrix
            end
        end
    end
end

%%
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %adding the path to OSL functions
osl_startup %starting the osl package


%% creatig nifti images with statistics..

%1) USE load_nii TO LOAD A PREVIOUS NIFTI IMAGE
%FOR EXAMPLE:
imag_example = load_nii('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/wholebrain_first_level_BC_coape_MMN_group_level_everybody_BC_randomise_c1_dir/tstat1_gc1_8mm.nii.gz');

%2) PUT YOUR MATRIX (FOR INSTANCE "T" OR "P") IN THE FIELD ".img"
%FOR EXAMPLE:
imag_example.img = P;

%3) SAVE NIFTI IMAGE USING save_nii
%FOR EXAMPLE:
save_nii(imag_example,'/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Pstats_img_results.nii.gz');

%4) USE FSLEYES TO LOOK AT THE FIGURE
imag_example = load_nii('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Pstats_img_results.nii.gz');

%% LEONARDO 27-28 JANUARY 2021

%% MUSMELO FRANCESCO - DEFINING CLUSTERS ON 3D BRAIN VOXELS STATISTICS

%loading p-values and t-values
clear DATA
P = load_nii('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Pstats_img_results.nii.gz');
T = load_nii('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Tstats_img_results.nii.gz');
%extracting matrix with statistics
P2 = P.img;
T2 = T.img;

%mask for non-0 voxels in brain imges (basically building a layout for actual brain voxels)
mask = zeros(size(T2,1),size(T2,2),size(T2,3));
mask(T2~=0) = 1; %assigning 1 when you have real brain voxels
%preapring data
data = T2;
data(data==0) = NaN; %removing non-brain voxels
data(P2>0.05) = NaN; %removing non-significant voxels
%cond1>cond2 (so t-values are positive)
data1 = data;
data1(T2<0) = NaN; %removing negative t-values (so cond1>cond2)
data1(~isnan(data1)) = 1; %assigning 1 to significant voxels
DATA{1} = data1; %storing data
%cond2>cond1 (so t-values are negative)
data2 = data;
data2(P2>0.000001) = NaN; %removing non-significant voxels (MORE STRICT THRESHOLD HERE FOR COND2>COND1)
data2(T2>0) = NaN; %removing positive t-values (so cond2>cond1)
data2(~isnan(data2)) = 1; %assigning 1 to significant voxels
DATA{2} = data2; %storing data

%getting MNI coordinates (you may get a warning since the skeletonized image is not exactly in MNI space, but close enough; I would not worry too much about that at the moment)
[ mni_coords, xform ] = osl_mnimask2mnicoords('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');

%prepation of information and data for the actual function
for ii = 2 %over directions of the contrast (cond1>cond2 and cond2>cond1)
    S = [];
    S.T = '/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Tstats_img_results.nii.gz'; %path to image with t-values
    S.outdir = '/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat';
    S.parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz';
    S.labels = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat';
    S.MNIcoords = mni_coords; %MNI coordinates of 8mm MNI152T1 brain
    S.data = DATA{ii};
    S.mask = mask; %mask of the brain layout youhave your results in
    S.permut = 1000; %number of permutations for Monte Carlo simulation
    S.clustmax = 0; %set 1 for only max cluster size of each permutation MCS (more strict); set 0 for every size of each cluster detected for each permutation MCS (less strict).
    S.permthresh = 0.001; %threhsold for MCS
    if ii == 1
        S.anal_name = ['Cond1vsCond2']; %name for the analysis (used to identify and save image and results)
    else
        S.anal_name = ['Cond2vsCond1']; %name for the analysis (used to identify and save image and results)
    end
    
    %actual function
    PP = BrainSources_MonteCarlosim_3D_LBPD_D(S);
end

%% COMBINING IMAGES (COND1>COND2 AND COND2>COND1) WITH SIGNIFICANT CLUSTERS

path = '/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat';
name1 = 'Cond1vsCond2_SignClust_1_Tvals.nii.gz'; % loading the images with the different conditions
name2 = 'Cond2vsCond1_SignClust_1_Tvals.nii.gz';
outname = 'HighVsLowDev_image';

cmd = ['fslmaths ' path '/' name1 ' -add ' path '/' name2 ' ' path '/' outname]; % create the sys command with the path
system(cmd) % run the sys command

%%

%% ANALYSIS CONSIDERING MUSICAL EXPERTISE

load('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/groups_leonardo_F.mat')

%% MUSICAL TRAINING CORRELATION (SIX DEVIANTS) - NOT IN THE PAPER

% load the matrix of all the subjects mean for the 2 deviants groups (high and low)
load('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/mean_allsubj_devlow.mat');
load('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/mean_allsubj_devhigh.mat')

DP = cat(5,dp,dp2); % dp = high deviants; dp2 = low deviants ; concatenate che 5th dimention fo dp2 to the 5th dim of dp

% DEVIANTS ORDER
% High deviants = melody mod, rhythm mod, transposition
% Low deviants = mistuning, timbre deviant, rhythm mistake

dev{1} = 'melmod'; dev{2} = 'rhymod'; dev{3} = 'trans'; dev{4} = 'mistu'; dev{5} = 'timbd'; dev{6} = 'rhytmis'; % creating a vector with all the deviants
list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/subject*dir'); %creating a list with the fif files outputted by Maxfilter
[nume,~,raw] = xlsread('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/ordinemusmelonuovo.xlsx'); % reading the orders from xls file
b = nume(:,2); %taking only the second value column, years of training
% preparing empty matrix for linear correlation 
sz = size(DP); 
R = zeros(sz(1),sz(2),sz(3),sz(5));
P = zeros(sz(1),sz(2),sz(3),sz(5));
for kk = 1:sz(5) % loop for all the six deviants
    for ii = 1:sz(1) % loop for all the 3D coordinates
        for yy = 1:sz(2)
            for zz = 1:sz(3)
                if DP(ii,yy,zz,1,kk) ~= 0
                    a = squeeze(DP(ii,yy,zz,:,kk)); % Squeeze of the subj dimension
                    [rho,p] = corr(a,b); % correlation function 
                    R(ii,yy,zz,kk) = rho; % pairwise linear correlation coefficient matrix
                    P(ii,yy,zz,kk) = p; % p value matrix                
                end
            end
        end
    end
    
    % Images production section
    
    %1) USE load_nii TO LOAD A PREVIOUS NIFTI IMAGE
    imag_example = load_nii('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/wholebrain_first_level_BC_coape_MMN_group_level_everybody_BC_randomise_c1_dir/tstat1_gc1_8mm.nii.gz');
    %2) PUT YOUR MATRIX (FOR INSTANCE "T" OR "P") IN THE FIELD ".img"
    imag_example.img = P(:,:,:,kk);
    %3) SAVE NIFTI IMAGE USING save_nii
    save_nii(imag_example,['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Pval_mustrain_' dev{kk} '.nii.gz']);
    imag_example.img = R(:,:,:,kk);
    %3) SAVE NIFTI IMAGE USING save_nii
    save_nii(imag_example,['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Rho_mustrain_' dev{kk} '.nii.gz']);
    disp(kk)
end


%% MUSICAL TRAINING CORRELATION (AVERAGED HIGH AND LOW DEVIANTS) - IN THE PAPER

sel = 0; % set 0 if you want to process high deviants or 1 for the low deviants

load('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/mean_allsubj_devlow.mat')
load('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/mean_allsubj_devhigh.mat')

% mean of the fifth dimension, the one with the low and high deviants
dpm1 = mean(dp,5); %high deviants
dpm2 = mean(dp2,5); %low deviants
list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/subject*dir'); %creating a list with the fif files outputted by Maxfilter
% reading again the excell file with the musiciaship training data (years training and years playing)
[nume,~,raw] = xlsread('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/ordinemusmelonuovo.xlsx');
b = nume(:,2); %taking only the years playing vector for all the subj
% Processing separately high and low deviants and proceed with the correlations
if sel == 0
    sz = size(dpm1);
    DP = dpm1; %single variable used in both the condition for the ttest
    R = zeros(sz(1),sz(2),sz(3));%creating double with zeroes for the t test (rho and p value) to be filled up with the results
    P = zeros(sz(1),sz(2),sz(3));
else
    sz = size(dpm2);
    DP = dpm2;
    R = zeros(sz(1),sz(2),sz(3));
    P = zeros(sz(1),sz(2),sz(3));
end
%loop accross all the 3 dimensions
for ii = 1:sz(1)
    for yy = 1:sz(2)
        for zz = 1:sz(3)
            if DP(ii,yy,zz,1) ~= 0 %if the value is different to 0 proceed with the squeeze of values on the fourth dimension and produce the ttest
                a = squeeze(DP(ii,yy,zz,:)); %fourth dimension with high or low deviants values
                [rho,p] = corr(a,b); % correlation function
                R(ii,yy,zz) = rho; % pairwise linear correlation coefficient matrix
                P(ii,yy,zz) = p; % p value matrix
            end
        end
    end
end
% images production section
if sel == 0
    %1) USE load_nii TO LOAD A PREVIOUS NIFTI IMAGE
    imag_example = load_nii('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/wholebrain_first_level_BC_coape_MMN_group_level_everybody_BC_randomise_c1_dir/tstat1_gc1_8mm.nii.gz');
    %2) PUT YOUR MATRIX (FOR INSTANCE "T" OR "P") IN THE FIELD ".img"
    imag_example.img = 1-P;
    imag_example.img(DP(:,:,:,1)==0) = 0;
    %3) SAVE NIFTI IMAGE USING save_nii
    save_nii(imag_example,['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Pval_mustrain_LOW_Dev.nii.gz']);
    imag_example.img = R;
    %3) SAVE NIFTI IMAGE USING save_nii
    save_nii(imag_example,['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Rho_mustrain_LOW_Dev.nii.gz']);
else
    %1) USE load_nii TO LOAD A PREVIOUS NIFTI IMAGE
    imag_example = load_nii('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/wholebrain_first_level_BC_coape_MMN_group_level_everybody_BC_randomise_c1_dir/tstat1_gc1_8mm.nii.gz');
    %2) PUT YOUR MATRIX (FOR INSTANCE "T" OR "P") IN THE FIELD ".img"
    imag_example.img = 1-P;
    imag_example.img(DP(:,:,:,1)==0) = 0;
    %3) SAVE NIFTI IMAGE USING save_nii
    save_nii(imag_example,['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Pval_mustrain_HIGH_Dev.nii.gz']);
    imag_example.img = R(:,:,:);
    %3) SAVE NIFTI IMAGE USING save_nii
    save_nii(imag_example,['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Rho_mustrain_HIGH_Dev.nii.gz']);
end

%% Defining clusters on 3D brain voxel statistics

freq = 2; % select 1 if you want to proceed with high or 2 for the low

if freq == 1
    T = load_nii(['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Rho_mustrain_HIGH_Dev.nii.gz']);
    P = load_nii(['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Pval_mustrain_HIGH_Dev.nii.gz']);
else
    T = load_nii(['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Rho_mustrain_LOW_Dev.nii.gz']);
    P = load_nii(['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Pval_mustrain_LOW_Dev.nii.gz']);

end
% extracting matrix with statistics
% P2 = P.img;
T2 = T.img;
P2 = P.img;
%mask for non-0 voxels in brain imges (basically building a layout for actual brain voxels)
mask = zeros(size(T2,1),size(T2,2),size(T2,3));
mask(T2 ~= 0) = 1; %assigning 1 when you have real brain voxels
%data
data = zeros(size(T2,1),size(T2,2),size(T2,3));
data(P2 > 0.99) = 1; %assigning 1 when you have real brain voxels
% data(~isnan(data)) = 1; %assigning 1 to significant voxels
DATA{1} = data; %storing data
%getting MNI coordinates
%OBS! you may get a warning since the skeletonized image is not exactly in MNI space, but close enough
[ mni_coords, xform ] = osl_mnimask2mnicoords('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');
%preparation of information and data for the actual function
for ii = 1%:2 %over directions of the contrast (cond1>cond2 and cond2>cond1)
    S = [];
    if freq == 1
        S.T = ['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Rho_mustrain_HIGH_Dev.nii.gz'];
    else    
        S.T = ['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat/Rho_mustrain_LOW_Dev.nii.gz'];
    end    
    S.outdir = ['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/after_maxfilter/source/Firstlevel.oat']; %output path
    S.parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz';
    S.labels = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat';
    S.MNIcoords = mni_coords; %MNI coordinates of 8mm MNI152T1 brain
    S.data = DATA{ii};
    S.mask = mask; %mask of the brain layout you have your results in
    S.permut = 1000; %number of permutations for Monte Carlo simulation
    S.clustmax = 1; %set 1 for only max cluster size of each permutation MCS (more strict); set 0 for every size of each cluster detected for each permutation MCS (less strict).
    S.permthresh = 0.001; %threhsold for MCS
    %final names
    if freq == 1
        S.anal_name = ['High_dev_musiscian_cluster']; %name for the analysis (used to identify and save image and results)
    else
        S.anal_name = ['Low_dev_musiscian_cluster']; %name for the analysis (used to identify and save image and results)
    end
    
    %actual function
    PP = BrainSources_MonteCarlosim_3D_LBPD_D(S);
end

%% SUPPLEMENTARY TABLES FOR THE PAPER

%% Extracting information about the clusters at source level and reporting it in xlsx files

%SOURCE SPACE
list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/Paper/*gz');
for pp = 1:length(list)
    clear PD PDn
    %tone x cluster 1
    %getting MNI coordinates of significant voxels within the provided image
    [ mni_coords, xform ] = osl_mnimask2mnicoords([list(pp).folder '/' list(pp).name]);
    %loading the image
    V = nii.load([list(pp).folder '/' list(pp).name]);
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
    %setting legend
    PD{1,1} = 'Brain Region'; PD{1,2} = 'Hemisphere'; PD{1,3} = 't-value'; PD{1,4} = 'MNI Coordinates';
    PD{2,4} = 'x'; PD{2,5} = 'y'; PD{2,6} = 'z';
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
    PDn = cell2table(PD); %remove the possible empty cell
    writetable(PDn,['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/Paper/MusmeloSupplementaryTables_MEGSources.xlsx'],'Sheet',list(pp).name(1:10)) %printing excel file
end

%% MEG SENSOR SPACE

%TABLES AFTER BONFERRONI CORRECTION

%% BONFERRONI CORRECTION

% min_time_point = 61; %16 = 0 seconds (first 15 points are pre-stimulus time)
% max_time_point = 242;
% time_point = max_time_point - min_time_point;
% pvalue = 0.05;
% ndeviant = 6; 
% 
% bonf_value = pvalue/(time_point* 6 * 102);

load('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/Paper/chanlabels.mat');
load('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/Paper/time.mat');
list = dir('/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/Paper/Bonf*');

for ii = 1:length(list)
    clear PD PDn
    load([list(ii).folder '/' list(ii).name]);
    PD = cell(103,183);
    cnt = 0;
    for pp = 2:2:204 %over MEG channels (only combined planar gradiometers)
        cnt = cnt + 1;
        PD(cnt+1,1) = chanlabels(1,pp);
    end
    cn = 0;
    for tt = 61:242 %over time-points
        cn = cn + 1;
        PD{1,cn+1} = time_sel(tt);
    end
    for pp = 2:103 %over MEG channels (only combined planar gradiometers)
        cn = 0;
        for tt = 61:242 %over time-points
            cn = cn + 1;
            PD{pp,cn+1} = Bonf(pp-1,tt);
        end
    end
    PDn = cell2table(PD); %remove the possible empty cell
    writetable(PDn,['/scratch7/MINDLAB2018_MEG-LearningBach-MemoryInformation/musmelo_F/Paper/MusmeloSupplementaryTables_MEGSensors.xlsx'],'Sheet',list(ii).name(10:19)) %printing excel file
end

%%

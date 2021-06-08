function [ data ] = ps_preparedata_spmobj_LBPD_D( S )


% It extracts data from SPM objects (in one or more datasets) in the
% requested time-range(s).
% It computes the source leakage correction (by orthogonalisation, removing
% zero-lag correlations) and calculates the envelope.
% Currently, it is designed to extract averaged trials only, since the
% tests with single trials did not show any success.
% Future updates need to allow the user not to compute source leakage
% correction and envelope, as well as extracting data for multiple trials
% and for subaveraging of trials.



%   INPUT:  -S.subjects_list:    cell matrix with character paths to spm
%                                object files. If you have one dataset only
%                                the matrix is 1 x number of subjects. If
%                                you have two or more datasets (e.g. two blocks of
%                                the same subjects during the experiment)
%                                the matrix is number of datasets x
%                                number of subjects.
%                                Subjects must match across different
%                                datasets.
%           -S.subjnum:          subjects ID (mainly for later controlling reasons)
%           -S.ROIs:             number of ROIs
%           -S.conditions:       conditions to be looked for.
%                                Each row corresponds to the progressive
%                                dataset.
%                                Each row is a 1 x 1 cell with 1 x number
%                                of conditions cell.
%           -S.time_range:       Time range(s) of interest.
%                                Each row corresponds to the progressive
%                                dataset.
%                                Each row is a 1 x 1 cell with 1 x number
%                                of time-ranges cell.
%           -S.montage_numb:     montage numbers to be extracted (must be
%                                the same across all subjects).
%                                Each row corresponds to the progressive
%                                dataset.
%           -S.fsample:          sampling rate, not used in this function,
%                                but potentially used later.


%   OUTPUT: -data:               cell matrix with data
%                                rows: subjects
%                                columns: organized as:
%                                   -1st: cond1 time1 dataset1
%                                   -2nd: cond1 time2 dataset1
%                                   -3th: cond2 time1 dataset1
%                                   -4th: cond2 time2 dataset1
%                                   -1st: cond1 time1 dataset2
%                                   -2nd: cond1 time2 dataset2
%                                   -3th: cond2 time1 dataset2
%                                   -4th: cond2 time2 dataset2





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 10/11/2018
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 24/08/2019
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%check fields
if ~isfield(S,'subjects_list')
    error('Error.. you need to specify the subject list')
end

if ~isfield(S,'ROIs')
    error('Error.. you need to specify the number of ROIs of the parcellation that you used')
end

if ~isfield(S,'conditions')
    error('Error.. you need to specify the name of the conditions')
end

if ~isfield(S,'time_range')
    error('Error.. you need to specify the time range(s) that you want')
end

%getting some inputs
ROIs = S.ROIs; %just to check that user has control on the parcellation.. it sounds annoying but previous experience suggested me to force user to think about..
montage_numb = S.montage_numb;
%pre-allocating matrices
[ndataset,n_Subjects] = size(S.subjects_list);
FF = zeros(length(S.conditions));
for ii = 1:length(S.conditions)
    [~,n_cond] = size(S.conditions{ii,1});
    [~,n_tr] = size(S.time_range{ii,1});
    FF(ii) = n_cond.*n_tr;
end
prep_data = cell(n_Subjects,sum(sum(FF))); %subj x time-ranges and conditions and datasets
label_data = cell(1,sum(sum(FF))); %labels (organized in the same way as prep_data)

%extracting and preparing data
numdum = 0;
for sda = 1:ndataset
    ind = zeros(1,length(S.conditions{sda,:}));
    condfa = S.conditions{sda,1};
    trfa = S.time_range{sda,1};
    %iterates over subjects
    for s = 1:n_Subjects
        disp(['processing dataset ' num2str(sda) ' subject ' num2str(s)])
        if ~isempty(S.subjects_list{sda,s})
            D = spm_eeg_load(S.subjects_list{sda,s});
            D = D.montage('switch',montage_numb(sda));
            ROIs_D = size(D);
            if ROIs ~= ROIs_D(1)
                error(['for subj ' num2str(s) ' you are trying to run analysis in a different parcellation than the one of your data.. you might have indicated the wrong montage or forgot to update the S.ROIs field..'])
            end
            %find indexes.. necessary because of the randomization across subjects
            for jj = 1:length(condfa)
                ind(jj) = find(strcmp(D.conditions,condfa{1,jj}));
            end
            %iterates over conditions
            for task = 1:length(condfa)
                %iterates over time ranges
                for kk = 1:length(trfa(1,:))
                    time_range = trfa{1,kk}; 
                    %here there was some discussion on (1) extracting data and
                    %then correcting for source leakage and caculating the
                    %envelope or (2) first correcting for sl, calculating the
                    %envelope and then extracting data..
                    %originally I thought the second option was more
                    %correct, however empirical testing showed me that (1)
                    %was given me results that appeared less related to
                    %resting state FC and closer to areas and networks that
                    %we know were related to the task I was examining, when compared to
                    %(2).
                    %a possible guess (that is at this stage simply a
                    %guess) is that, at least in frequency range 4-15Hz,
                    %having a longer signal that is associated both to the
                    %task and some rest after that makes the envelope
                    %capturing the oscillation related to the task but also
                    %to the following rest and therefore the FC matrix
                    %coming out with later analysis tends to reflect both
                    %resting and task-related connectivity.
                    %this would explain why (2) was given results in the
                    %'middle' btween resting state and task-related
                    %networks.
                    %However, at the moment, this is simply a conclusion
                    %based on empirical observations and reasonings that are
                    %only partialy mathematically grounded and therefore
                    %I/we need to explore and think more about this issue.
                    %therefore, my advice is to extract the signal only in
                    %correspondance to the actual task plus some margines
                    %on both sides to prevent the data to be affected by
                    %later boundary artifacts, usually introduced by
                    %extracting the instantaneous phase.
                    %e.g. in one of my studies I was interested in a signal of around
                    %1.5 seconds (150Hz sampling rate) related to my task
                    %and I extracted the signal of about 230 time-points
                    %plus 25 tme-points before the onset of the task and 25
                    %time-points after then end of the task (these
                    %time-points have been later removed because affected
                    %by boundary artefacts).
                    ts = D(1:ROIs,time_range(1):time_range(2),ind(task));
                    tso_lk = remove_source_leakage_b(ts,'symmetric'); %leakage correction

%                     prep_data{s,kk+(task-1)*length(trfa)+numdum} = abs(ts); %running this and not the line above is for avoiding source leakage correction

                    prep_data{s,kk+(task-1)*length(trfa)+numdum} = hilbenv(tso_lk);
                    
                    %write the label only for the first subject available (since it is the same for everybody)
                    if isempty(label_data{1}) || isempty(label_data{2})
                        label_data(kk+(task-1)*length(trfa)+numdum) = {['dataset_' num2str(sda) '_' condfa{task} '_timerange_from_' num2str(time_range(1)) '_to_' num2str(time_range(2))]};
                    end                    
                    %write time only for the first available subject
                    if ~exist('times','var')
                        times = D.time;
                    end
                end
            end
        else
            disp(['dataset ' num2str(sda) ' subject ' num2str(s) ' is empty'])
        end
    end
    numdum = kk+(task-1)*length(trfa)+numdum; 
end

%preparing output
data.prep_data = prep_data;
data.label = label_data;
data.subjs_list = S.subjects_list;
data.subjnum = S.subjnum;
data.S_original_info = S;
data.time_original = times;

end


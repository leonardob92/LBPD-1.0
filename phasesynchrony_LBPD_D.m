function [] = phasesynchrony_LBPD_D(S)

% Calculates the instantaneous phase synchrony between each couple of brain
% areas by extracting the instantaneous phase of each brain area with Hilbert
% transform and then calculating the cosine similarity of the
% instantaneous phases of different brain areas.
% Calculates also the Kuramoto Order parameter for each functional
% connectivity matrix. This parameter can be considered a measure of global
% connectivity of the matrix.
% This procedure has been adapted to MEG from fMRI works by Enrico Glerean
% and Joana Cabral.
% Here data is saved automatically, therefore it is suggested to create a
% new folder each time this function is run.



%   INPUT:  -S.data:        data from ps_preparedata_spmobj_LBPD_D.m.
%                           n of subjects x n of tasks cell matrix.
%                           Each cell contains doubles (actual data).
%                           Here tasks refer to potentially different
%                           time-ranges/conditions/datasets.. user must take care of
%                           not mixing them up.
%                           Currently, the length of the tasks must be the
%                           same. In the future, it may be useful to make
%                           this more flexible.
%           -S.save_path:   path where you want to save data in.
%           -S.ROIs:        number of ROIs
%           -S.left_bound:  time-points to be discarded in the beginning of the
%                           time-course (due to Hilbert transform boundary artifacts)
%           -S.right_bound: time-point to be removed from the end of the time-course


%   OUTPUT: -it saves:      -instantaneous phase synchrony matrices (ROIs x
%                            ROIs x time x subjs).
%                            One .mat file for each task
%                           -Kuramoto order parameter (number of tasks cell
%                            array with double in every cell (subjs x time)
%                           -labels (labels of the tasks)
%                           -subjs_real: to be used for subsequent statistics
%                           -S structure with all information used by the
%                            algorithms (with exception of the actual data).





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 10/11/2018
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 04/10/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%initialising variables and checking some parameters
N_areas = S.ROIs; %getting number of parcels
%computing phase coherence and Kuramoto order parameter
data = S.data; %getting the actual data plus some information about it
%saving labels if any
if ~isfield(data,'label')
    label = [];
else
    label = data.label;
end
%creating the directory if it does not exist already
if ~exist(S.save_path,'dir')
    mkdir(S.save_path);
end
save([S.save_path '/label'],'label');
%checking that you have the actual data
if ~isfield(data,'prep_data')
    error('Error.. you need to have the prepared data from ps_preparedata_spmobj_leonardo_D.m function or data in the same configuration')
end
%getting the actual data
prep_data = data.prep_data;
disp('Processing Phase Synchrony of MEG data')
%getting number of subjects and tasks/conditions
[n_Subjects, n_Task] = size(prep_data);
%getting original number (ID) of actual subjects (this is slightly redundant.. but for now I prefer to keep this implementation)
countsubs = 0;
for ppp = 1:n_Subjects
    if isempty(find(cellfun('isempty',prep_data(ppp,:)) > 0)) %contort way to check if a subject is empty
        countsubs = countsubs + 1; %count real subjects
        subs_real(countsubs) = ppp; %this stores the original numbers (ID) of the subjects that are being processed (it is important for later statistics)
        if countsubs == 1 %doing that only for first subject (assuming each subject has the same kind of parcellation)
            jkl = size(prep_data{ppp});
            if jkl(1) ~= N_areas %checking if the actual parcellation is coherent with the ones indicated by the user (this is done since from my experience it can be a good idea to force the user to carefully check his/her parcellation indication)
                error(['you set number of ROIs = ' num2str(N_areas) ' while you actual ROIs number seems to be = ' num2str(jkl(1))])
            end
        end
    end
end
save([S.save_path '/subs_real'],'subs_real'); %save real subjects on disk (used in later statistics..)  
%pre-allocating space for Kuramoto order parameter cell
OP = cell(n_Task,1);

%actual computation
for task = 1:n_Task %iterates over tasks/conditions/time-ranges
    countsub_REAL = 0;
    %getting T_max
    Tmax = 0;
    cpd2 = 0;
    while Tmax == 0 %this while loop is useful if subject 1 has no data and therefore it looks for the first subject with data
        cpd2 = cpd2 + 1;
        [~,Tmax] = size(prep_data{cpd2,task}); %assuming that every subject has the same amount of data-points within the same task
    end
    %pre-allocating space for Kuramoto order parameter
    kop = zeros(countsubs,Tmax-(S.left_bound-1+S.right_bound));
    %pre-allocating space for IFC for each time-point and subject
    iFC_time_4D = zeros(S.ROIs,S.ROIs,Tmax-(S.left_bound-1+S.right_bound),countsubs);
    %iterates over subjects
    for s = 1:n_Subjects
        disp(['processing subject ' num2str(s) ' and task ' num2str(task)])
        if isempty(find(cellfun('isempty',prep_data(s,:)) > 0)) %if it cannot find empty cells in the data matrix.. (otherwise just skip the subject)
            countsub_REAL = countsub_REAL + 1;  
            MEG_ind = prep_data{s,task};
%             [~,Tmax] = size(MEG_ind); %getting the Tmax for each subject and each task (this is probably not necessary if you consider same-length trials, but for now I leave it like that)
            Phase_MEG = zeros(N_areas,Tmax);
            %getting the instantaneous phase by applying the Hilbert transform
            for seed = 1:N_areas %over ROIs
                Phase_MEG(seed,:) = angle(hilbert((MEG_ind(seed,:)) - mean((MEG_ind(seed,:))))); %subtracting mean to center the data on 0
            end
            t_new = 0; %initialise new counter for time
            for t = S.left_bound:Tmax-S.right_bound %removing some time-points in the beginning and end of the time-course since the Hilbert transform tends to produce boundary artifacts
                %from some empirical testing, at least in frequency range
                %4-15Hz, it seems that there is a linear relationship
                %linking signal length and amount of time-points (at both
                %beginning and end of the signal) that are affected by
                %boundary artifacts.
                %Future explorations are called for confirm and specifically define
                %this linear relationship (if possible) and to investigate
                %it in different frequency ranges.
                %further explorations seemed to show that lower frequency
                %bands tend to produce a boundary artifact that occurr for
                %longer time.
                %My advice here is to check the results and make some
                %testing for detecting the extent of the boundary artifact.
                t_new = t_new + 1; %increasing counter
                %Calculate the Instantaneous FC (cosine similarity of Phase Synchrony)
                iFC = zeros(N_areas);
                for n = 1:N_areas %this could be calculated only for the upper triangle and then mirrored since the matrix is square..
                    for p = 1:N_areas
                       iFC(n,p) = cos(Phase_MEG(n,t) - Phase_MEG(p,t)); %here 't' extract the proper time-points (e.g from 15 to whatever if you want to remove the first 15 time-points)
                    end
                end
                kop(countsub_REAL,t_new) = abs(mean(exp(1i*(Phase_MEG(:,t))))); %Kuramoto order parameter stored using 't_new' index
                iFC_time_4D(:,:,t_new,countsub_REAL) = iFC; %storing iFC matrix for each time-point and each subject
            end
        end
    end
    OP(task) = {kop};
    save([S.save_path '/task_' num2str(task)],'iFC_time_4D','-v7.3'); %otherwise saving them on disk
    clear iFC_time_4D %clearing it for memory limitations
end
save([S.save_path '/OP'],'OP'); %saving Kuramoto order parameter (even if this is not very large..)
S.data.prep_data = []; %removing actual data
save([S.save_path '/S_original'],'S'); %saving all information used by the algorithms

end





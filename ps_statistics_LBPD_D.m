function [] = ps_statistics_LBPD_D(S)

% Calculates t-tests (or sign-rank tests) comparing the phase synchrony of
% each couple of brain areas between two conditions (or one condition and two
% groups of subjects).
% To be used after: phasesynchrony_LBPD_D.m



%   INPUT:  -S.load_path:               path for loading data
%           -S.save_path:               path for saving new data computed by this function
%           -S.subjs_list:              list of subjects grouped into one or two groups (e.g. {[1,2,3,4]} or {[1,4],[2,3]})
%           -S.group_comparison:        [1] for one group; [1 2] for two groups (referred to the subjs_list, so the original IDs!!)
%           -S.condition_contrast:      contrasts to be done: remember that you always need to write numbers
%                                       consistently with data indicated in ps_preparedata_spmobj_LBPD_D.m.
%                                       Use [2 0 1 0] for contrasting task 1 against task 3.
%                                       Use [0 1 0 0] for same condition (task)
%                                       (and probably different groups..).
%                                       Each row corresponds to a different contrast between tasks.
%                                       If you want to compare different groups write:
%                                       [0 1 0 0; 0 0 1 0] (and S.group_comparison = [1 2]) for having different groups for
%                                       task2 and then task3.
%                                       Currently, if you want several contrasts between different groups of subjects
%                                       you need to run the function multiple times.
%                                       You also need to contrast tasks(conditions) with the same amount of data-points.
%           -S.cond_contr_meanbaseline: 1 if you want to average over time the task 2 of the contrast (usually done if task 2 is a baseline)
%           -S.Num_window:              number of windows for doing subaveraging of the data (this is implemented since the
%                                       instantaneous FC matrix is usually quite noisy at a single time-point, therefore doing
%                                       some sub-averaging can improve the reliability of the results, even if it reduces
%                                       the temporal resolution of the analysis.
%                                       Leave it empty [] for not doing any subaveraging.



%   OUTPUT: -it saves:                  -Contrasts for FC in .mat files (e.g. 'Contr_1_2_pval_etc.' and 'Contr_1_2_zval_etc.' means contrast
%                                        of task1 against task2).
%                                       -Contrasts for Kuramoto Order paramter in POP (p-val), ZOP (z- or t- val) (time-course x contrast).
%                                       -original S structure with input for the function.
%                                        Be careful not to overwrite it!! It is suggested to have one independent folder for each contrast you do!






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 10/11/2018
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 25/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%initialising variables and checking some parameters
%subjs list and groups comparison
if ~isfield(S,'save_path') %if save_path is not specified it is assumed that it is the same as load_path..
    S.save_path = S.load_path;
end
%creating the directory if it does not exist already
if ~exist(S.save_path,'dir')
    mkdir(S.save_path);
end
subjs_list = S.subjs_list;
groups_comp = S.group_comparison;
%conditions about groups of subjects
if length(groups_comp) ~= length(subjs_list)
   error('probably you have specified a comparison between two groups but you have only one group or viceversa..')
end

%checking for not mixing up baseline with groups contrast
if S.cond_contr_meanbaseline == 1 && length(S.group_comparison) > 1
    warning('you are probably comparing two groups and using one group as "baseline", so averaged over time and then compared to the other group.. it probably does not make sense.. are you sure you want to continue?')
end

%you need to have subs_real to properly define the subjects to be compared
load([S.load_path '/subs_real']);


%defining indexes for comparing (or not comparing) different groups of subjects)
if length(groups_comp) == 1
    [~,~,idx_subj1] = intersect(subjs_list{groups_comp(1)},subs_real); % this necessary otherwise you do not get the correct subjects..
    [~,~,idx_subj2] = intersect(subjs_list{groups_comp(1)},subs_real); %here it is the same since in this case you actually have only one group of subjects.. this is simply a computational trick for treating more easily some later analyses
else
   if length(groups_comp) > 2
       error('you cannot compare more than 2 groups')
   else
       %finding the indexes of subs_real of the original subject IDs provided by the user
       [~,~,idx_subj1] = intersect(subjs_list{groups_comp(1)},subs_real); %THAT MEANS THAT YOU NEED TO PROVIDE THE TWO GROUPS CONSIDERING THE ORIGINAL SUBJECT IDs!!!!!!!!
       [~,~,idx_subj2] = intersect(subjs_list{groups_comp(2)},subs_real); %THAT MEANS THAT YOU NEED TO PROVIDE THE TWO GROUPS CONSIDERING THE ORIGINAL SUBJECT IDs!!!!!!!!
   end
end


%calculating contrasts between conditions
if ~isempty(find(S.condition_contrast > 0)) %if contrasts have been specified
    SE = size(S.condition_contrast);
    for kk = 1:SE(1) %this is thought if there are more than just one contrast to be computed.. 
        disp(['processing contrast ' num2str(kk)])
        cond1 = find(S.condition_contrast(kk,:) == 2); %getting the condition 1 index
        cond2 = find(S.condition_contrast(kk,:) == 1); %getting the condition 2 index
        if isempty(cond1) %if you are using the same condition and comparing two groups here you need this trick for computing the proper contrasts later (REMEMBER TO PROPERLY SET THIS OPTION IN THE CONTRASTS.. e.g. [0 1 0 0]!!!)
            cond1 = cond2(1);
%             cond2 = cond2(1);
        end
        if length(cond2) ~= 1
            warning('you indicated the number "1" in the contrasts not only one time.. that is probably going to lead you towards terrible problems..consider to kill the process and check the inputs you gave..')
        end
        if S.cond_contr_meanbaseline == 0 %it is more meaningfull comparing the matrices at specific time-points with a baseline averaged across the whole time-window (for example, time-window = 100ms and you want to see matrices each 10ms.. it is better to average in the 10 10ms-timewindows the task and then comparing each of the 10 matrices with the averaged-over-100ms matrix of the baseline..)
            disp('loading Kuramoto order parameter..')
            %calculating Kuramoto order parameter (cond1)
            load([S.load_path '/OP']); %these loadings may be done only once, reshaping a bit the function, but for now it is ok..
            OP2 = OP{cond1};
            OP_cond1 = OP2(idx_subj1,:);
            %Kuramoto.. cond2
            OP2 = OP{cond2};
            OP_cond2 = OP2(idx_subj2,:);
            %IFC..
            disp('loading data.. cond 1')
            iFC_cond1 = load([S.load_path '/task_' num2str(cond1)]); %remember that this 'cond1' refers to the physical space occupied by the task/condition '1' outputted by LEiDA_MEG_leonardo_preparedata
            disp('loading data.. cond 2')
            iFC_cond2 = load([S.load_path '/task_' num2str(cond2)]);
            iFC_cond1_matr = iFC_cond1.iFC_time_4D;
            clear iFC_cond1
            iFC_cond2_base = iFC_cond2.iFC_time_4D;
            clear iFC_cond2
            SC = size(iFC_cond1_matr);
            SCdum = size(iFC_cond2_base);
            if SC(3) ~= SCdum(3)
                error('if you do not average the time-points the two conditions must have the same amount of data-points..')
            end
        else %calculating average over time for baseline
            disp('loading Kuramoto order parameter..')
            %calculating Kuramoto order parameter (cond1)
            load([S.load_path '/OP']);
            OP = OP{cond1};
            OP_cond1 = OP(idx_subj1,:);
            %Kuramoto.. cond2
            load([S.load_path '/OP']);
            OP = OP{cond2};
            OP_cond2 = mean(OP(idx_subj2,:),2); %calculating Kuramoto condisering condition 2 as baseline and therefore averaging it across the whole time-window
            %IFC..
            disp('loading data.. cond 1')
            iFC_cond1 = load([S.load_path '/task_' num2str(cond1)]);                        
            disp('loading data.. cond 2')
            iFC_cond2 = load([S.load_path '/task_' num2str(cond2)]);                        
            iFC_cond1_matr = iFC_cond1.iFC_time_4D;
            clear iFC_cond1
            iFC_cond2_base = squeeze(mean(iFC_cond2.iFC_time_4D,3));
            clear iFC_cond2
            SC = size(iFC_cond1_matr);
        end
        %subaveraging..
        if ~isempty(S.Num_window) %if user requires it
            disp('doing subaveraging..')
            gg = floor(SC(3)/S.Num_window); %getting number of time-points for each window (floored)
            ggm = mod(SC(3),S.Num_window); %getting discrepancies between actual data and maximum dividable data
            dum = iFC_cond1_matr(:,:,1:SC(3)-ggm,:); %extracting data cond1 (maximum time-points - ggm)
            clear iFC_cond1_matr %clearing previous data for space reasons
            d2 = reshape(dum,SC(1),SC(1),gg,S.Num_window,SC(4)); %reshaping data
            iFC_cond1_matr = squeeze(mean(d2,3)); %averaging data
            if S.cond_contr_meanbaseline == 0 %if the condition 2 is not a baseline, you need to do the same for it.. otherwise, of course, not..
                %here, I assume that the two conditions have the same
                %amount of data-points (as it should be..)
                dum = iFC_cond2_base(:,:,1:SC(3)-ggm,:); %extracting data cond2 (maximum time-points - ggm)
                clear iFC_cond2_base
                d2 = reshape(dum,SC(1),SC(1),gg,S.Num_window,SC(4)); %reshaping data
                iFC_cond2_base = squeeze(mean(d2,3)); %averaging data
            end
            SC = size(iFC_cond1_matr);
        end
        %preallocation of space for stats matrices
        STAT_p_fullcontr = zeros(SC(1),SC(2),SC(3));
        STAT_zval_fullcontr = zeros(SC(1),SC(2),SC(3));
        POP = zeros(SC(3),1);
        ZOP = zeros(SC(3),1);
        flag_ttest_general = 0;
        flag_ttest_op = 0;
        %PHASE SYNCHRONY
        disp('contrasting phase synchrony matrices and Kuramoto order parameter..')
        for zz = 1:SC(3)
            disp(['computing contrasts for time-point ' num2str(zz)])
            counttriangle = 1; %this is for calculating the statistics only in the upper triangle of the matrix
            for nn = 1:SC(1)
                counttriangle = counttriangle + 1; %this is for calculating the statistics only in the upper triangle of the matrix
                for mm = counttriangle:SC(2) %this is for calculating the statistics only in the upper triangle of the matrix
                    a = squeeze(iFC_cond1_matr(nn,mm,zz,idx_subj1));
                    if S.cond_contr_meanbaseline == 0
                        b = squeeze(iFC_cond2_base(nn,mm,zz,idx_subj2));
                    else
                        b = squeeze(iFC_cond2_base(nn,mm,idx_subj2));
                    end
                    [p,~,stats] = signrank(a,b); %computing contrasts                
                    STAT_p_fullcontr(nn,mm,zz) = p; %storing values
                    if isfield(stats,'zval') %this is not ideal since theoretically it could occur that some contrasts are done by using the sign rank test and other ones by using the t-test.. however, practically this is not going to happen and the tests to be used will be either one or the other one.. this is not computationally elegant, but it should be fine                              
                        STAT_zval_fullcontr(nn,mm,zz) = stats.zval;
                    else
                        flag_ttest_general = 1;
                        [~,p,~,stats] = ttest(a,b);
                        STAT_p_fullcontr(nn,mm,zz) = p;
                        STAT_zval_fullcontr(nn,mm,zz) = stats.tstat;
                    end                        

                end
            end
            STAT_p_fullcontr(:,:,zz) = (STAT_p_fullcontr(:,:,zz) + STAT_p_fullcontr(:,:,zz)') - eye(size(STAT_p_fullcontr(:,:,zz),1)).*diag(STAT_p_fullcontr(:,:,zz)); %computing the symmetric lower triangle in the matrix
            STAT_zval_fullcontr(:,:,zz) = (STAT_zval_fullcontr(:,:,zz) + STAT_zval_fullcontr(:,:,zz)') - eye(size(STAT_zval_fullcontr(:,:,zz),1)).*diag(STAT_zval_fullcontr(:,:,zz)); %computing the symmetric lower triangle in the matrix
            %KURAMOTO ORDER PARAMETER
            aop = OP_cond1(:,zz);
            if S.cond_contr_meanbaseline == 1
                bop = OP_cond2;
            else
                bop = OP_cond2(:,zz);
            end
            %actual statistics
            [pop,~,statsop] = signrank(aop,bop);
            %storing results
            POP(zz) = pop;
            if isfield(statsop,'zval')
                ZOP(zz) = statsop.zval;
            else
                flag_ttest_op = 1;
                [~,pop,~,statsop] = ttest(aop,bop);
                POP(zz) = pop;
                ZOP(zz) = statsop.tstat;
            end
        end
        if flag_ttest_general == 1
            warning('signrank test for general contrasts does not provide z-values.. it might not be working properly.. therefore t-test was calculated instead..')
        end
        if flag_ttest_op == 1
            warning('signrank test for OP does not provide z-values.. it might not be working properly.. therefore t-test was calculated instead..')
        end
        %characters showing with group/s have been compared.. this is used in the name of the data saved on the disk
        if length(S.group_comparison) == 1
            GP = [num2str(S.group_comparison(1)) '_' num2str(S.group_comparison(1))];
        else
            GP = [num2str(S.group_comparison(1)) '_' num2str(S.group_comparison(2))];
        end
        %saving statistics FC matrices
        save([S.save_path '/Contr_' num2str(cond1) '_' num2str(cond2) '_STAT_p_fullcontr_group_comp_' GP],'STAT_p_fullcontr');
        save([S.save_path '/Contr_' num2str(cond1) '_' num2str(cond2) '_STAT_zval_fullcontr_group_comp_' GP],'STAT_zval_fullcontr');
        %saving statistics Kuramoto Order parameter
        save([S.save_path '/Contr_' num2str(cond1) '_' num2str(cond2) '_POP_group_comp_' GP],'POP');
        save([S.save_path '/Contr_' num2str(cond1) '_' num2str(cond2) '_ZOP_group_comp_' GP],'ZOP');
    end
else
    error('you need to specify the contrast to be done!!')
end
%saving original S structure with inputs
save([S.save_path '/S_original'],'S');


end





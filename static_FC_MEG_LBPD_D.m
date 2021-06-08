function [ H, P_binary, P_ROIs, MATF, cor_mat ] = static_FC_MEG_LBPD_D( S, data )

% It calculates static FC by computing Pearson's correlations of the
% envelope of the data given in input. It does that independently for each
% condition.
% Additionally, it can perform contrast between two conditions, plot the FC
% matrix of the contrast and give back the result of Monte Carlo
% simulations for assessing whether some ROIs are significantly central
% within the brain network.
% Currently implemented for contrasting conditions only and not groups of
% subjects.




%   INPUT:  -data:                        cell matrix n x m with n = number of subjects and m =
%            (containing prep_data field) number of conditions (each cell has a double                                        
%                                         matrix ROIs x time-points with the actual
%                                         data).
%                                         Data is assumed to be sorted as:
%                                         LeftROI1, RightROI1, LeftROI2, etc. The
%                                         function will reshape it in the order:
%                                         LeftROI1, LeftROI2, LeftROI3.., RightROI1, RightROI2.., etc.
%                                         Data is outputted by ps_preparedata_spmobj_LBPD_D.m
%                                         (it can also be combined with the data outputted by
%                                         preparing_baseline_from_restingstate_LBPD_D.m).
%           -S.contrast:                  must be a double array of the same length of conditions
%                                         (column size of data.prep_data) expressing the conditions to be 
%                                         contrasted (e.g. [2 0 1 0] contrasts cond1 against cond3.
%                                         Currently, it is implemented for
%                                         contrasting conditions only and not
%                                         groups of subjects.
%           -S.resam_number:              resample parameter if you want to resample the
%                                         data (this can be a potential trick to increase
%                                         the strength of the correlation coefficients).
%                                         Leave empty [] for not resampling.
%                                         To me it does not make much sense because I am more interested
%                                         in spatial patterns that are not random more than absoluted
%                                         values, but since it does for others, it is worthy to mention it.
%           -S.fsample:                   sampling rate in Hz.
%                                         Leave empty [] for having the function looking for it
%                                         in 'data.S_original_info', since it is outputted by
%                                         'ps_preparedata_spmobj_LBPD_D.m' function.
%                                         Of course, this is possible only if you provide the specific data variable
%                                         outputted by 'ps_preparedata_spmobj_LBPD_D.m'.                                       
%           -S.plot_conditions:           array of lenght of conditions with 1 values for the conditions to
%                                         be plotted (e.g. [1 0 1 0] plots cond1 and cond3)
%                                         Remember to indicate with "0" the conditions that you do not want to plot.
%           -S.plot_contr:                set 1 to plot the contrast matrix
%           -S.perm_degree:               number of permutation for degree MCS
%           -S.perm_max:                  1 to calculate MCS considering only the maximum central ROI 
%                                          of each shuffled matrix (dimensionality: 1 x permutations).
%                                         0 for building the distribution for MCS with degree of all ROIs (dimensionality:
%                                          ROIs x permutations).
%           -S.p_perm:                    value for binarising P matrix after t-tests (or signrank tests).
%                                         Those matrix will then be used in MCS.
%                                         Here 5% has to be expressed as 0.05.
%                                         If set = 0 it takes the t-values (z-values previously calculated).
%                                         If empty [] the default value is 0.05.
%           -S.p_threshMC:                value for considering significant an ROI after estimating its probability
%                                         to not be random by using MCS.
%                                         Here 5% has to be expressed as 0.05.
%                                         If empty [] the default value is 0.05.
%           -S.clims_contr:               set limits for the figure (contrast, e.g. [-4 4])
%           -S.clims_cond:                set limits for the figure (conditions, e.g. [-0.3 0.3])
%           -S.parametric:                set 1 if you want parametric statistics and 0 for non-parametric statistics
%                                         when calculating contrasts between conditions.
%           -S.schemaball_contrast_label: 1 if you want schemaball for the contrast
%           -S.schemaball_cond_label:     1 for one schemaball per condition
%           -S.schemplth:                 set the value x for plotting the x percent strongest connections (default 100)
%                                         Here you need to provid 5% as 5,
%                                         not 0.05.
%           -S.colbacksch:                specify background color for schemaball.
%                                         Leavy empty [] for default.
%                                         Default is w (white).
%           -S.lab_parc_ph:               path to ROIs labels (assumed to be ordered as LRLRLR. The function will
%                                         reshape data and labels returning LLLRRR order).
%           -S.outpath:                   path where schemaball images have to be saved
%           -S.name_sch:                  name for schemaball image files (character).
%                                         Leave empty [] for not saving schemaball images.

%   OUTPUT: -H:                           1 if the Monte Carlo simulation was significant
%           -P_binary:                    1s = significantly central ROIs (reshaped LLLRRR) 
%           -P_ROIs_r:                    significantly central ROIs (labels) within the brain network according to Monte
%                                         Carlo simulations.
%           -MATF:                        cell containing connectivity matrices for the different conditions to be plotted,
%                                         plus two for the requested contrast (t-values and p-values) (if no contrast was
%                                         requested the last two cells are empty [].
%           -cor_mat:                     connectivity matrices (reshaped as LRLRLR) for each condition and subject
%                                         (ROI x ROI x condition x subject).
%           -plots..






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 16/01/2019
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 13/08/2019
% Leonardo Bonetti, Aarhus, DK, 27/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














%check some input parameters and set some variables for later usage
%here it would be better to increS.colbackschase the number of controls..
contrast = S.contrast;
resam_number = S.resam_number;
%checking contrast specification
% if length(find(contrast == 1)) ~= 1 || length(find(contrast == 2)) ~= 1
%     error('contrast array needs to have only one value = 2 and one value = 1 (corresponding to the conditions to be contrasted..)')
% end
%checking if there is the path for loading labels
if ~isfield(S,'lab_parc_ph')
    error('you need to specify the path of a file containing the ROI labels.. I do not have any by default for trying to avoid possible uncorrect uses of this function..')
else
    lab_parc_ph = S.lab_parc_ph;
end
%if it does not exist assigning 100% connections by default
if ~isfield(S,'schemplth')
    schemplth = 100;
else
    schemplth = S.schemplth;
end
%getting info about subjects and conditions
[n_subjs,n_conds] = size(data.prep_data);
%and about actual data (assuming every subject has the same amount of data)
ROIs = 0;
cpd2 = 0;
while ROIs == 0 %this while loop is useful if subject 1 has no data and therefore it looks for the first subject with data
    cpd2 = cpd2 + 1;
    [ROIs,~] = size(data.prep_data{cpd2,1});
end
%preallocating space for MATF cell, a cell containing the connectivity matrices for the different conditions and the requested contrast t-values and p-values
MATF = cell(length(S.plot_conditions) + 2,1);
%settings for resampling, if requested..
if ~isempty(S.resam_number)
    if isempty(S.fsample)
        if isfield(data,'S_original_info')
            fsample = S_original_info.sampling_rate; %this is assumed to be there since the current function is thought
        else
            error('if you want to resample you need to specify the current sampling rate..')
        end
    end
end
%if the threshold for binarising results after t-tests (or signrank tests)
%is not specified by the user, assigning 0.05 by default
if isempty(S.p_perm)
    S.p_perm = 0.05;
end
if isempty(S.p_threshMC) %same concept but for threshold after MCS
    S.p_threshMC = 0.05;
end
if isempty(S.colbacksch) %same concept for background color of schemaball
    S.colbacksch = 'w';
end
%getting information..
plot_conditions = S.plot_conditions;
plot_contr = S.plot_contr;
clims_contr = S.clims_contr;
clims_cond = S.clims_cond;
%to be used later for sorting data (first Left hemisphere than Right one)
order = [1:2:ROIs ROIs:-2:2];
%for properly initializing cor_mat
count = 0;
for ii = 1:n_subjs
    if isempty(find(cellfun('isempty',data.prep_data(ii,:)) > 0))
        count = count + 1;
    end
end
%if wanted setting some parameters for schemaball to be plotted later..
if S.schemaball_contrast_label == 1 || S.schemaball_cond_label == 1
    S2 = [];
    S2.name_gif = []; %for not plotting connectivity in the brain layout
    S2.schball_l = 1;
    S2.symmetric_l = 1;
    S2.lab_parc_ph = lab_parc_ph;
    S2.colbacksch = S.colbacksch;
    S2.schemplth = schemplth;
    if isempty(S.clims_contr)
        S2.extr = []; %leaving empty extreme values used to normalize values before schemaball
    else
        S2.extr = S.clims_contr; %using limits values for matrix as extreme values for standardizing schemaball
    end
    if ~isempty(S.name_sch)
        S2.outpath = S.outpath;
    end
end

%calculating static FC by using Pearson's correlation
cor_mat = zeros(ROIs,ROIs,n_conds,count);
cor_mat_ave = zeros(ROIs,ROIs,n_conds);
for cond = 1:n_conds
    disp(['processing condition ' num2str(cond)])
    actual_sub = 0;
    for s = 1:n_subjs
        if isempty(find(cellfun('isempty',data.prep_data(s,:)) > 0)) %if it cannot find empty cells in the data matrix.. (otherwise just skip the subject)
            actual_sub = actual_sub + 1;
            tso = data.prep_data{s,cond}; %selecting data
            %leonardo MIT stay
            Hen = tso; %this is because you already calculated the envelope in 'ps_preparedata_spmobj_LBPD_D.m' function
            %the above line is not meaningful now, but I leave it for chronological issues of this function
            %until here
            if isempty(resam_number) %if you want to resample or not
                cor_dummy = corr(Hen'); %Pearson's correlations
                cor_mat(:,:,cond,actual_sub) = cor_dummy(order,order); %here data is being reshaped                 
            else
                cor_dummy = corr(resample(Hen',resam_number,fsample));
                cor_mat(:,:,cond,actual_sub) = cor_dummy(order,order);
            end
        end
    end
    cor_mat_dummy = squeeze(cor_mat(:,:,cond,:)); %for not averaging the zeros of the initialization of cor_mat..
    cor_mat_ave(:,:,cond) = mean(cor_mat_dummy,3); %averaging across subjects
end

H = [];
P_ROIs = [];

%calculating contrast
if ~isempty(find(contrast > 0))    
    disp('processing contrast..')
    cond1 = find(S.contrast == 2); %getting the condition 1 index
    cond2 = find(S.contrast == 1); %getting the condition 2 index
    cond_1 = squeeze(cor_mat(:,:,cond1,:)); %getting data cond1
    cond_2 = squeeze(cor_mat(:,:,cond2,:)); %getting data cond2
    P = zeros(ROIs,ROIs);
    matrx_z = zeros(ROIs,ROIs);
    if S.parametric == 1 %here I give the possibility to choose parametric (t-test) or non-parametric (signrank test) statistics.. theoretically you would need at least to run a test to assess the level of 'normality' of the data.. however, this function does not provide any test like that
        for ii = 1:ROIs %here I could work only on the upper triangle, but the computation is anyway very fast..
            for jj = 1:ROIs
                if ii ~= jj %for not computing the statistics of the 1-value diagonal
                    a = squeeze(cond_1(ii,jj,:));
                    b = squeeze(cond_2(ii,jj,:));
                    [~,p,~,stats] = ttest(a,b);
                    P(ii,jj) = p;
                    matrx_z(ii,jj) = stats.tstat;
                end
            end
        end
    else
        for ii = 1:ROIs
            for jj = 1:ROIs
                if ii ~= jj %for not computing the statistics of the 1-value diagonal
                    a = squeeze(cond_1(ii,jj,:));
                    b = squeeze(cond_2(ii,jj,:));
                    [p,~,stats] = signrank(a,b);
                    if ~isfield('zval','stats')
                        warning('signrank test may have not worked properly.. consider to recompute contrasts by doing parametric t-tests..')
                    end
                    P(ii,jj) = p;
                    matrx_z(ii,jj) = stats.zval;
                end
            end
        end
    end
    stats = [];
    stats.cond_1 = cond_1;
    stats.cond_2 = cond_2;
    stats.matrx_z = matrx_z;
    P2 = P;
    if S.p_perm ~= 0 %if user wants to binarise the matrix
        %considering only positive contrasts (2 - 1 when you define the conditions)
        Pdum = double((double(P < S.p_perm) + double(matrx_z > 0)) == 2); %isolating conditions when you have both p_values significant at the threshold defined by the user and the zval that is positive
        P = Pdum;
    else %if not getting the t-values (z-values)
        P = matrx_z;
    end
    %calculating significance of degree of contrast matrix vertices
    [H,P_binary,P_ROIs,~,~,~] = degree_segregation_MCS_LBPD_D(P,0,S.perm_degree,S.p_threshMC,lab_parc_ph,1,S.perm_max);
    %plotting contrast
    disp('plotting contrast..')
    if plot_contr == 1
        figure
        colormap(jet)
        if isempty(clims_contr)
            imagesc(matrx_z+diag(zeros(ROIs,1)))
        else
            imagesc(matrx_z+diag(zeros(ROIs,1)),clims_contr)
        end
        axis square
        colorbar
        set(gcf,'Color',S.colbacksch)
        if ROIs == 90 %if AAL parcellation
            labelY = {'19','39','59','79','82','62','42','22','2'};
            labelX = labelY;
            ax = gca;
            ax.YTickLabel = labelY;
            ax.XTickLabel = labelX;
        elseif ROIs == 38 %if 38-ROI parcellation
            labelY = {'9','19','29','38','28','18','8'};
            labelX = labelY;
            ax = gca;
            ax.YTickLabel = labelY;
            ax.XTickLabel = labelX;
        end
        if isfield(data,'label')
            title(['Static FC - contrast between ' data.label{cond1} ' and ' data.label{cond2}])
        end    
    end
    
    %schemaball for contrast
    if S.schemaball_contrast_label == 1
        matrx_z = matrx_z.*~(eye(ROIs)); %assigning 0s to diagonal..
        S2.MATF = matrx_z;
        if ~isempty(S.name_sch)
            S2.name_sch = [S.name_sch '_contrast_' num2str(cond1) '_' num2str(cond2)];
        else
            S2.name_sch = [];
        end
        FC_Brain_Schemaball_plotting_LBPD_D(S2) %actual function
    end
end


%setting clims for conditions if they are not specified by user
if isempty(clims_cond)
    MAX = max(cor_mat_ave(cor_mat_ave < 1));
    MIN = min(cor_mat_ave(cor_mat_ave < 1));
    clims_cond = [MIN MAX];
end

%storing this data for outputting purposes
if ~isempty(plot_conditions)
    for ii = 1:length(plot_conditions)
        MATF(ii) = {cor_mat_ave(:,:,ii)+diag(nan(ROIs,1))};
    end
    idx23 = 3; %for later outputting procedures..
else
    idx23 = 0;
end

%plotting conditions (matrices)
if ~isempty(find(plot_conditions > 0)) > 0
    disp('plotting conditions..')
    idx2 = find(plot_conditions > 0);
    figure
    colormap(jet)
    for ii = 1:length(idx2)
        if mod(length(idx2),2) ~= 0 %this is for creating a proper grid for the subplot (if the number of matrices is odd then you have a grid floor(n/2) + 1,floor(n/2) + 1 otherwise you have floor(n/2),floor(n/2)
            grid = floor(length(idx2)/2) + 1;
        else
            grid = length(idx2)/2;
        end
        if length(idx2) == 2 %if the conditions are only 2 you need to have a grid 2,1
            subplot(grid+1,grid,ii) %grid 2 rows and 1 column
        else
            subplot(grid,grid,ii)
        end
        lkj = cor_mat_ave(:,:,ii).*~(eye(ROIs)); %assigning 0s to diagonal.. nicer when plotting..
        imagesc(lkj,clims_cond)
        axis square
        colorbar
        if ROIs == 90
            labelY = {'19','39','59','79','82','62','42','22','2'};
            labelX = labelY;
            ax = gca;
            ax.YTickLabel = labelY;
            ax.XTickLabel = labelX;
        elseif ROIs == 38
            labelY = {'9','19','29','38','28','18','8'};
            labelX = labelY;
            ax = gca;
            ax.YTickLabel = labelY;
            ax.XTickLabel = labelX;
        end    
        title(['Condition ' num2str(idx2(ii))])   
        set(gcf,'Color',S.colbacksch)
    end
    if S.schemaball_cond_label == 1
        %plotting one independent schemaball for each condition
        for ii = 1:length(idx2)
            cor_mat_ave2 = cor_mat_ave(:,:,ii);
            cor_mat_ave2 = cor_mat_ave2.*~(eye(ROIs)); %assigning 0s to diagonal..
            S2.MATF = cor_mat_ave2;
            if ~isempty(S.name_sch)
                S2.name_sch = [S.name_sch '_condition_' num2str(idx2(ii))];
            else
                S2.name_sch = [];
            end
            FC_Brain_Schemaball_plotting_LBPD_D(S2) %actual function
            title(['Condition ' num2str(ii)])
        end
    end    
end

if ~isempty(find(contrast > 0)) %if contrast was requested, storing statistics..
    MATF(idx23 + 1) = {matrx_z}; %storing t-values of contrast
    MATF(idx23 + 2) = {P2}; %storing p-values of contrast
end


end



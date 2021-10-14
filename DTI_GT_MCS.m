function [ p_final_pos, p_final_neg, node_pos_idx, node_neg_idx, diff1 ] = DTI_GT_MCS( S )

% It computes Monte-Carlo simulations (MCS) on distribution of differential
% data (contrasting two groups) derived from graph teorethical measures.
% If you do not want to use all your nodes (brain areas), but only a
% subset of them, simply provide only that subset for the two groups. This
% may be especially useful if you want to compare differences only between
% brain areas with a value higher (and lower) than a certain threshold.




%  INPUT:   -S.data:        data provided in a 1 x 2 cell array with one group in each cell.
%                           Specifically, each cell contain a double matrix = n. of subjects in the group x n. of nodes (brain areas).
%           -S.permtype:    choose the permutation type:
%                           1 for permuting the degrees of each couple of regions (high vs low PR) independently;
%                           2 for permuting all degrees together 
%           -S.permnum:     number of permutations for MCS (e.g. 1000)
%           -S.thr:         threshold for zeroing values that you do not want to take into account (i.e. keeping values
%                           smaller than -S.thr and bigger than S.thr).
%                           Leave empty [] for not using any threshold


%  OUTPUT:  -p_final_pos:   p-value for MCS (right side of the distribution)
%           -p_final_neg:   p-value for MCS (left side of the distribution)
%           -diff1:         difference of the ROIs between two groups






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 08/03/2021


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    







% 1) SOLUTION CONSIDERING ALSO STANDARD DEVIATIONS OF GROUPS AND NOT ONLY MEAN
%comparing the two groups (using t-values, but not really t-tests)
% diff1 = zeros(1,size(S.data{1},2));
% for ii = 1:size(S.data{1},2) %over nodes (brain areas)
%     a = S.data{1}(:,ii); %group 1
%     b = S.data{2}(:,ii); %group 2
%     [~,~,~,stats] = ttest2(a,b); %computing t-values
%     diff1(1,ii) = stats.tstat; %storing t-values
% end
% 2) SOLUTION CONSIDERING ONLY MEDIAN (MUCH FASTER EVEN THOUGH LESS ACCURATE; TO NOTE, FOR THIS PARTICULAR PURPOSE IS GOOD ENOUGH)
% SOLUTION 2 IS SUGGESTED
diff1 = median(S.data{1},1)-median(S.data{2},1); %Difference between mean of the two groups.. for this purpos  should be the same, even though slightly less accurate
% diff1 = mean(S.data{1},1)-mean(S.data{2},1); %Difference between mean of the two groups.. for this purpos  should be the same, even though slightly less accurate
if ~isempty(S.thr)
   diff1(diff1<S.thr & diff1>-S.thr) = 0; 
end
% figure
% scatter(1:size(S.data{1},2),diff1) %plotting the difference
% ylim([((-1) * (max(max(abs(diff1))) + max(max(abs(diff1)))*0.1)) (max(max(abs(diff1))) + max(max(abs(diff1)))*0.1)]) %minimum and maximum limits of the figure computed by adding 10% to ylim
distrdiff1 = length(find(diff1>0)); %number of times when diff1 is higher than 0
distrdiff2 = length(find(diff1<0)); %number of times when diff1 is smaller than 0

if ~isempty(S.thr) %if you asked for a thresholding
    permpos = zeros(1,S.permnum);
    permneg = zeros(1,S.permnum);
    for pp = 1:S.permnum
        vect = ones(1,length(diff1)); %vector of ones
        vect(1:round(length(diff1)/2)) = -1; %half vector of -1s
        idx_dummy = randperm(length(diff1)); %creating a permuted array with length of nodes (brain areas)
        perM3 = diff1 .* vect(idx_dummy); %shuffling sign of diff1 elements
        permpos(pp) = length(find(perM3>0)); %getting values for positive sign
        permneg(pp) = length(find(perM3<0)); %and negative sign
    end
%     figure
%     hist(permpos)
    p_final_pos = length(find(distrdiff1<permpos))/S.permnum;
    p_final_neg = length(find(distrdiff2<permneg))/S.permnum;
else
    %actual computations
    perM = zeros(size(S.data{1},2),2); %Initialize a new 2-column matrix for computing the permutations
    perM(:,1) = median(S.data{1},1); %degrees of subjects of group 1 in the 1st column
    perM(:,2) = median(S.data{2},1); %degrees of subjects of group 2 in the 2nd column
    %     perM(:,1) = mean(S.data{1},1); %degrees of subjects of group 1 in the 1st column
%     perM(:,2) = mean(S.data{2},1); %degrees of subjects of group 2 in the 2nd column
    diffpvec1 = zeros(1,S.permnum); %vector containing the differences between permutations of group 1 vs group 2
    diffpvec2 = zeros(1,S.permnum); %vector containing the differences between permutations of group 1 vs group 2
    for pp = 1:S.permnum %over permutations
        if S.permtype == 1 %permuting the degrees of each couple of regions independently
            idx_dummy = randperm(size(perM,1)*2); %creating a permuted array with length of nodes (brain areas) times number of groups (2)
            perM2 = perM(idx_dummy); %shuffling real data with randomized indices
            perM3 = reshape(perM2,[size(perM,1) size(perM,2)]); %reshaping the vector into a matrix
        else %permuting all degrees together
            perM3 = zeros(90,2);
            for qq = 1:size(S.data{1},2) %over nodes (brain areas)
                idx_lab = randperm(2); %permuting only 2 numbers (indices of the 2 groups)
                perM3(qq,:) = perM(qq,idx_lab); %randomly assigning real data using randomized indices (independently for each node (brain areas))
            end
        end
        diffperm = perM3(:,1) - perM3(:,2); %difference between permuted vectors
        if ~isempty(S.thr)
            diffperm(diffperm<S.thr & diffperm>-S.thr) = 0;
        end
        diffpvec1(pp) = length(find(diffperm>0)); %storing how many times the difference was positive (> 0); the distribution of this should tend to normality around 0, with increasing number of permutations
        diffpvec2(pp) = length(find(diffperm<0)); %storing how many times the difference was negaive (< 0); the distribution of this should tend to normality around 0, with increasing number of permutations
        disp(['permutation number ' num2str(pp) ' / ' num2str(S.permnum)])
    end
    %plotting histogram of computed differences
%     figure
%     hist(diffpvec1) %this distribution should tend to normality centerd on number of nodes (brain areas)/2; with increasing permutations number
    %computing final p-value as the ratio between how many times original numbers were higher than 0 compared to permuted ones and the total number of permutations
    p_final_pos = length(find(distrdiff1<diffpvec1))/length(diffpvec1);
    p_final_neg = length(find(distrdiff2<diffpvec2))/length(diffpvec2);
end
node_pos_idx = diff1>0;
node_neg_idx = diff1<0;




end


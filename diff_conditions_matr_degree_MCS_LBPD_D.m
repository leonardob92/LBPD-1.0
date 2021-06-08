function [ ROI_final_mat ] = diff_conditions_matr_degree_MCS_LBPD_D( a, b, ROI_lab, X, permut, thresh )

% It calculates the difference between the significant degree over time of brain
% areas (with the whole brain) in two conditions. Then, it calculates Monte
% Carlo simulations to assess which differences are significant.
% Please note that the assumption here is that if the ROIs have the same
% relevance in task x, then their significance occurrence over time
% (represented by 1s in a and b matrices) should be nearly equal. If that is the
% case, with simulated data you get no significance after this MCS. On the
% contrary, if your simulated data consist of random significant events
% (1s in a and b) that, when summed, return a normal distribution over the different
% ROIs, this MCS will very likely return some significant results.
% I strongly advice to carefully think about what you are doing before
% using this function. Moreover, I think that the p-values obtained from MCS
% should be divided according to how many differences you are looking at.
% This is because the single p-value is showing the probability after MCS of only one
% difference given the input data. Therefore, a good idea would be dviding
% your alpha level by the number of differences between ROIs that you are
% looking at. E.g. if you look at whole-brain level in 90-ROI AAL
% non-cerebellar parcellation, assuming an aplha level (thresh) = 0.05, you should
% probably adjust it as follows: 0.05/(2*90) = 2.7e-04. 2 because you look
% at both tales of the normal distribution of permuted data and 90 because
% of the number of ROIs in this specific parcellation. In this case it would
% be recommended to run something like at least 1000 permutations (in other words a number
% of permutations that is consistent with the adjusted alpha level). This looks very
% strict but from my tests it may be a more reliable approach. This
% function is conceived for reporting every p-value < thresh (variable passed
% by user), then the user should critically think if all of those
% "significant" values are to be kept or not, according to what I wrote few
% lines above.



%   INPUT:  -a:         matrix n x t, where n are ROIs and t are time-points of
%                       condition 1 (double)
%           -b:         matrix n x t, where n are ROIs and t are time-points of
%                       condition 2 (double)
%           -ROI_lab:   labels of the ROIs corresponding to a and b
%                       values (cell array)
%           -X:         set 0 for working on the whole matrices (double)
%                       set the time-window you want to work on (e.g. X = 1:180)
%           -permut:    number of permutations (double)
%           -thresh:    threshold for considering significant the
%                       difference between a and b after MCS.

%   OUTPUT: -ROI_final: cell with number of ROIs x 3 where:
%                       -first column: ROIs sorted according to the strongest
%                           differences between conditions
%                       -second column: number of times that the two conditions
%                           differed for the specific ROI
%                       -third column: p-value after Monte Carlo simulations (if
%                           the cell is empty it means that that specific ROI was
%                           not significant




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 10/11/2018
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 26/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%getting which time-points the user requested
if X == 0
    P1 = a;
    P2 = b;
else
    P1 = a(:,X);
    P2 = b(:,X);
end
SO = size(P1);

%calculating differences and sorting ROIs (cond1 - cond2)
PPP = sum((P1-P2),2);
[Val,II] = sort(PPP,'descend'); %cell2mat(IFC_sub(:,2)),'descend')
%     IFC_idx = ROI_lab(II,:);
IFC_idx = cell(SO(1),2);
for ii = 1:SO(1)
    IFC_idx(ii,1) = ROI_lab(II(ii));
    IFC_idx(ii,2) = {PPP(II(ii))};
end

%storing also the sorted vectors for the original conditions
P1s = sum(P1,2);
[~,IIp1] = sort(P1s,'descend'); %cell2mat(IFC_sub(:,2)),'descend')
IFC_idxp1 = cell(SO(1),2);
for ii = 1:SO(1)
    IFC_idxp1(ii,1) = ROI_lab(IIp1(ii));
    IFC_idxp1(ii,2) = {P1s(IIp1(ii))};
end
P2s = sum(P2,2);
[~,IIp2] = sort(P2s,'descend'); %cell2mat(IFC_sub(:,2)),'descend')
IFC_idxp2 = cell(SO(1),2);
for ii = 1:SO(1)
    IFC_idxp2(ii,1) = ROI_lab(IIp2(ii));
    IFC_idxp2(ii,2) = {P2s(IIp2(ii))};
end

%permuting matrices a and b, summing a and b over time (independently)
%and then calculating their differences; plotting hist(R3_diff) you
%obtain the (approximable to normal) distribution describing the differences
%of random dispositions of the elements of a and b
R3_diff = zeros(SO(1),permut);
for rr = 1:permut %over permutations
    %cond1
    f = SO(1)*SO(2);
    idx_2 = randperm(f); %create a permuted array from 1 to length(f)
    r_dummy = a(idx_2); %taking element in matrix a with index idx_2
    %rebuild 2D matrix with actual data
    r3 = reshape(r_dummy,[SO(1),SO(2)]);
%     [I_l,J_l] = ind2sub([SO(1),SO(2)],1:length(r_dummy));
%     %rebuild 2D matrix with actual data
%     r3 = zeros(SO(1),SO(2));
%     count = 0;
%     for ii = 1:length(r_dummy)
%         count = count + 1;
%         r3(I_l(ii),J_l(ii)) = r_dummy(count);
%     end
    R3_sum_cond1 = sum(r3,2); %cond1 sum
    %cond2
    idx_2 = randperm(f); %create a permuted array from 1 to length(f)
    r_dummy = b(idx_2); %taking element in matrix b with index idx_2
    %rebuild 2D matrix with actual data
    r3 = reshape(r_dummy,[SO(1),SO(2)]);
%     [I_l,J_l] = ind2sub([SO(1),SO(2)],1:length(r_dummy));
%     %rebuild 2D matrix with actual data
%     r3 = zeros(SO(1),SO(2));
%     for ii = 1:length(r_dummy)
%         r3(I_l(ii),J_l(ii)) = r_dummy(ii);
%     end
    R3_sum_cond2 = sum(r3,2); %cond2 sum
    R3_diff(:,rr) = R3_sum_cond1 - R3_sum_cond2; %difference between cond1 and cond2 that is stored for every permutation
end

%using distribution of R3_diff values for getting the significances for the original data
%positive differences
maxdum = cell2mat(IFC_idx(1,2)); %maximum difference
count = maxdum + 1; %counting variables..
p = 0;
count1 = 0;
while p < (thresh) %until the values are significant it proceeds.. 
    count1 = count1 + 1;
    count = count - 1;
    rrr = find(R3_diff > count);   
    p = length(rrr)/(SO(1)*permut);
    if p < (thresh) %this is not an elegant control, but for now I will keep it like that..
        P(count1,1) = count;
        P(count1,2) = p; %P becomes a theoretical table of p-values
    end
end
count1 = count1 - 1; %trick for not having a 0 value in P
%negative differences
mindum = cell2mat(IFC_idx(SO(1),2)); %maximum difference
count = mindum - 1; %counting variables..
p = 0;
while p < (thresh) %until the values are significant it proceeds..
    count1 = count1 + 1;
    count = count + 1;
    rrr = find(R3_diff < count);   
    p = length(rrr)/(SO(1)*permut);
    if p < (thresh) %this is not an elegant control, but for now I will keep it like that..    
        P(count1,1) = count;
        P(count1,2) = p;
    end
end
% P(count1,:) = []; %trick for not having a 0 value in P

%taking significant values stored in P and place them in the ROI_final
%cell matrix to be given as output of this function
ROI_final = cell(SO(1),3);
ROI_final(:,1:2) = IFC_idx;
if exist('P','var')
    for ii = 1:size(P,1) %this is not perfectly optimised but it should work well enough
        dummyI = find(P(ii,1) == Val);
        ROI_final(dummyI,3) = {P(ii,2)}; %here if it does not find anything (so if the theoretical p-value is not present in the real data) it simply does not index the matrix
    end
else
    ROI_final = [];
end
%preparing output
ROI_final_mat.diff = ROI_final;
ROI_final_mat.cond1 = IFC_idxp1;
ROI_final_mat.cond2 = IFC_idxp2;

end


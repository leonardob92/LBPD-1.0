function [ h_deg, P_binary_deg, P_ROIs_deg_r, P_val, h_seg, ROI_seg_label, ROI_seg, glob_seg ] = degree_segregation_MCS_LBPD_D( P, thresh, permut, threshMC, ROI_label, reshlab, perm_max )

% It calculates the degree of matrix vertices and modularity (by using the 
% 'modularity_und.m' function by Newman) in both original and
% randomised matrices.
% It assesses significance for vertices centrality and estimated modularity
% by performing Monte Carlo simulations.
% Simulations with uniformly distributed random numbers (rand built-in
% function in Matlab) showed that, if perm_max = 0, having a threshMC <= 0.001 prevents from
% having any falsely significant ROI when binarising the matrix (and nearly
% any when non-binarising it). A more common threshMC = 0.05 or 0.01
% usually results in 1-2 falsely significant ROIs when binarising the matrix (which may be not
% dramatic.. and maybe preferable than adopting more strict thresholds that
% could potentially result in loosing real information.. whatever people think, I believe
% it is important to know it).
% Otherwise, setting perm_max = 1 is another possibility, returning more
% strict results (working only on maximum degree permuted matrix) and therefore reducing the
% few false significant ROIs emerged with perm_max = 0 for more common threshMC = 0.05 or 0.01.


%   INPUT:  -P:             n x n of values (e.g. p-values obtained
%                           contrasting two matrices).
%                           The matrix cannot contain only 0s.
%           -thresh:        value for binarizing the matrix P (give 0 for working
%                           directly on the P matrix, e.g. in the case you already
%                           binarized it or you do not want to binarize it..)
%           -permut:        number of permutations to be run
%           -threshMC:      threshold for considering significant a value
%                           after Monte Carlo simulation (e.g. 5% must be
%                           given as 0.05).
%                           If left empty, default = 0.05.
%           -ROI_label:     path to the file with ROI labels for
%                           parcellation in use (must be a .mat file with a
%                           characted array named 'lab').
%                           Currently, it is not very flexible..
%                           Set 'null' if you do not want any label.
%           -reshlab:       1 to reshape the label from
%                           LeftROI1,RightROI1,LeftROI2,etc. to
%                           LeftROI1,LeftROI2,LeftROI3,etc.
%                           This is contort, but useful for how the
%                           different functions are concatenated together.
%                           0 not to reshape the label.
%                           Note that this reshapes ONLY the ROI labels and
%                           not the actual data!! It could be replaced by
%                           loading a file with already reshaped labels..
%           -perm_max:      -set 1 to calculate MCS considering only the maximum degree
%                            of each shuffled matrix (and build the reference distribution
%                            on those values only, dimensionality: 1 x permutations)
%                           -set 0 to build that distribution with degrees of all ROIs 
%                            for each shuffled matrix (dimensionality: ROIs x permutations)
%                           If left empty [], default = 0.

%   OUTPUT: -h:             1 if there is at least one ROI that is significantly central
%                           within the matrix; 0 if not.
%           -P_binary_deg:  binary vector n x 1 where 1 indicates the ROIs significantly
%                           central within the network.
%           -P_ROIs_deg_r:  name of the significantly central ROIs
%           -P_val:         p-values of each ROI degree
%           -h_seg:         1 if the matrix can be significantly divided in
%                           subgraphs
%           -ROI_seg_label: labels of ROIs and belonging to correspondent
%                           subgraph
%           -ROI_seg:       subdivision without labels
%           -glob_seg:      estimation of how much the matrix is "segregationable"..
%                           (values close to 1 means that the matrix can be reasonably
%                           divided in subgraphs).





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 10/11/2018
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 03/10/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%checking if matrix P contains only 0s.. (meaning that there are no
%significant connections between ROIs..)
if isempty(find(P~=0))
    error('matrix P cannot contain only 0s..')
end
%extracting size of P (for example it needs to be 90 if you have non-cerebellar AAL)
[~,m] = size(P);
%options for thresholding p-values within this function or outside
if thresh ~= 0
    P_dummy = double(P < thresh);
else
    P_dummy = P;
end
P_dummy = P_dummy.*~(eye(m)); %assigning 0s to diagonal..
% plotting binarised matrix
% figure
% imagesc(P_dummy);
% colorbar
%calculate degree of matrix P vertices
deg_all_P = sum(P_dummy); %calculate the sum of each ROI
% DEG_P = max(deg_all_P); %consider maximum number of the sums as the maximum degree of the matrix P
%calculate modularity (Newman)
[mod_ROI_ori,mod_glob_ori] = modularity_und(P_dummy);
%if not otherwise specified it works only on the maximum degree of the vertices for each shuffled matrix
if isempty(perm_max)
    perm_max = 0;
end
%if not otherwise specified threshMC is set as 0.05
if isempty(threshMC)
    threshMC = 0.05;
end

%actual computation
count_false_seg = 0;
if perm_max == 1
    DEG_R3_arr = zeros(permut,1);
else
    DEG_R3_arr = zeros(m,permut);
end
%permutations
for pp = 1:permut
    f = find(triu(ones(m),1) == 1); %find indexes of elements of the upper triangle of the square matrix
    idx_dummy = randperm(length(f)); %create a permuted array from 1 to length(f)
    idx_2 = f(idx_dummy); %shuffle the indexes in f according to the order in idx_dummy
    r_dummy = zeros(m*m,1); %initialise vector
    r_dummy(f) = P_dummy(idx_2); %taking element in matrix data with index idx_2 (shuffled indexes of the upper triangle)
%     [I_l,J_l] = ind2sub(m,1:length(r_dummy)); %convert indexes from 1D array to 2D matrix 
%     %rebuild 2D matrix with actual data
%     r3 = zeros(m);
%     for ii = 1:length(r_dummy)
%         r3(I_l(ii),J_l(ii)) = r_dummy(ii);
%     end
    %rebuild 2D matrix with actual data
    r3 = reshape(r_dummy,[m,m]);
    %ricreate the simmetric matrix for calculating more easily the sums
    r3 = (r3 + r3') - eye(size(r3,1)).*diag(r3);
    %mirroring the procedure above for P (degree)
    deg_all_r3 = sum(r3);
    if perm_max == 1
        DEG_R3_arr(pp) = max(deg_all_r3);
    else
        DEG_R3_arr(:,pp) = deg_all_r3;
    end
    %calculating segregation of the randomly rearranged matrix
    [~,mod_glob_r3] = modularity_und(r3); %only global value of segregation
    if mod_glob_r3 >= mod_glob_ori %count how many times the global value of the randomized matrix is >= than the original data one (if you store mod_glob_r3 and then plot hist you get a quasi-normal distribution)
        count_false_seg = count_false_seg + 1;
    end
end

%segregation MCS significance
ROI_seg = mod_ROI_ori;
glob_seg = mod_glob_ori;
ROI_seg_label = cell(m,2);
if count_false_seg/permut >= 0.05
    h_seg = 0;
else
    h_seg = 1;
end
%degree significance %this whole idea of doubl tails was wrong, but for now
%I decided to keep it here commented as a reminder
% if tail == 2 %if tail == 2 you want to see both tails of the normal distribution created by the MCS, therefore if your threshold was e.g. 0.05, I think you should look into 0.025 on both side (tails) of the normal distribution
%     threshMC = threshMC/2;
% end
%getting threshMC corresponding value in the randomized data
if perm_max == 1 %it works only with the maximum degree of the vertices of each shuffled matrix
    dummysort = sort(DEG_R3_arr,'descend');
    DEG_R3_array_number = dummysort(round((threshMC*100*permut/100) + 1)); %here I look for the highest threshMC (e.g. 5%) random number of vertices degrees calculated during the permutations
%     if tail == 2 %if you want 2 tails, looking at other side of sorted permuted data
%         DEG_R3_array_number_left = dummysort(round(length(dummysort) - ((threshMC*100*permut/100) + 1)));
%     end
else %it works with the degree of all vertices (ROIs) of every matrix
    INT_R3_all2 = reshape(DEG_R3_arr,[m*permut,1]);
    dummysort = sort(INT_R3_all2,'descend');
    DEG_R3_array_number = dummysort(round((m*threshMC*100*permut/100) + 1)); %here I look for the highest threshMC (e.g. 5%) random number of vertices degrees calculated during the permutations (defining the threshold!)
%     if tail == 2 %if you want 2 tails, looking at other side of sorted permuted data
%         DEG_R3_array_number_left = dummysort(round(length(dummysort) - ((m*threshMC*100*permut/100) + 1)));
%     end
end    
%comparing original with permuted data
P_binary_deg = double(deg_all_P > DEG_R3_array_number)'; %calculating which ROIs had a degree higher than the 95% random maximum vertices degree of the matrix calculated during permutations test
if ~isempty(find(P_binary_deg == 1))
    h_deg = 1;
else
    h_deg = 0;
end

%calculating p-values
P_val = zeros(m,1);
for ii = 1:length(deg_all_P) %over ROIs
    pk = find(deg_all_P(ii) == dummysort); %getting to what degree of the reference distribution the ROI degree ii is equal
    if ~isempty(pk)
        P_val(ii) = pk(end)/length(dummysort); %dividing that number (last one in case of multiple occurrences) by the number of elements of the reference distribution
    else
        if deg_all_P(ii) > dummysort(1) %this is because pk can be = [] both when the original ROI degree in deg_all_P is bigger or smaller than all of the permuted ROI degrees..
            P_val(ii) = 0; %storing 0 if the original ROI dergee was always bigger than all degrees in the reference distribution
        elseif deg_all_P(ii) < dummysort(end)
            P_val(ii) = 1; %otherwise if it is smaller than the smallest value in dummysort it stores 1
        else
            [~,olp] = min(abs(dummysort - deg_all_P(ii))); %otherwise you ned to find the closest original degree to the one in the reference distribution
            P_val(ii) = olp/length(dummysort); %and then calculate the p-value by dividing that index by the length of dummysort
        end
    end
end
        


% if tail == 2
%     bind = (-1)*double(deg_all_P < DEG_R3_array_number_left)';
%     P_binary_deg = P_binary_deg + bind;
% end
    
%getting labels (this will be updated in the future for increasing flexibility)
if ~strcmp(ROI_label,'null') 
    %reshaping AAL ROIS order to fit with the data    
    ROI_label_char = load(ROI_label);
    if isstruct(ROI_label_char)
        P_ROIs_dummy = ROI_label_char.lab;
    else
        P_ROIs_dummy = ROI_label_char;
    end
    if reshlab == 1 %if you want to reshape the labels from LRLRLR to LLLRRR (this is useful if you reshape the data before submitting it to this function (it could be otherwise solved by loading from disk a file with sorted labels..))
        order = [1:2:m m:-2:2];
        P_ROIs_dummy = P_ROIs_dummy(order,:);
    end
    fgf = find(P_binary_deg > 0); %getting significant ROI indices
    P_ROIs_deg_r = cell(length(fgf),2);
    for ioi = 1:length(fgf) %over significant ROI
        P_ROIs_deg_r(ioi,1) = {P_ROIs_dummy(fgf(ioi),:)}; %storing ROI label
        P_ROIs_deg_r(ioi,2) = {P_val(fgf(ioi))}; %storing p-value
    end
%     P_ROIs_deg_l = P_ROIs_dummy(P_binary_deg < 0,:);
    for pip = 1:m
        ROI_seg_label(pip,1) = {P_ROIs_dummy(pip,:)};
        ROI_seg_label(pip,2) = {ROI_seg(pip)};
    end
else
    P_ROIs_deg_r = [];
%     P_ROIs_deg_l = [];
    ROI_seg_label = [];    
end



end


function [ OUT ] = twoD_MCS_LBPD_D( P, thresh, permut, threshMC, perm_max, t1, t2 )

% It individuates clusters in a binary 2D matrix and tests their
% significance by Monte Carlo simulations.
% This function is designed to work with binary data.


%   INPUT:  -P:             n x n of values (e.g. p-values obtained
%                           contrasting two matrices).
%                           The matrix cannot contain only 0s.
%           -thresh:        value for binarizing the matrix P (give 0 for working
%                           directly on the P matrix, e.g. in the case you already
%                           binarized it).
%           -permut:        number of permutations to be run
%           -threshMC:      threshold for considering significant a value
%                           after Monte Carlo simulation (e.g. 5% must be
%                           given as 0.05).
%                           If left empty [], default = 0.05.
%           -perm_max:      -set 1 to calculate MCS considering only the
%                            maximum cluster of each shuffled matrix (and build the reference distribution
%                            on those values only, dimensionality: 1 x permutations).
%                           -set 0 to build that distribution with cluster sizes of all 
%                            clusters of each shuffled matrix.
%                            If left empty [], default = 1.
%           -t1:            vector of values describing first dimension (rows) of matrix P (e.g. frequencies)
%           -t2:            vector of values describing second dimension (columns) of matrix P (e.g. time in seconds)

%   OUTPUT: -OUT:           cell matrix with 7 columns corresponding to:
%                            -cluster #
%                            -cluster sizes
%                            -cluster p-values
%                            -minimum along first dimension (row)
%                            -maximum along first dimension (row)
%                            -minimum along second dimension (column)
%                            -maximum along second dimension (column)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 08/10/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











%checking if matrix P contains only 0s.. (meaning that there are no significant connections between ROIs..). If it is the case, returning error
if isempty(find(P~=0))
    error('matrix P cannot contain only 0s..')
end
%extracting sizes of P
[n,m] = size(P);
%options for thresholding p-values within this function or outside
if thresh ~= 0
    P_dummy = double(P < thresh);
else
    P_dummy = P;
end
%if they are not provided by user, assigning to t1 and t2 dimensions of matrix indices by default
if ~exist('t1','var')
    t1 = 1:n;
end
if ~exist('t2','var')
    t2 = 1:m;
end
% %plotting binarised matrix
% figure
% imagesc(P_dummy);
% colorbar
%identifying islands of 1s and 0s within matrix P
[AG,NSV,~] = islands(P_dummy);
NS = NSV(NSV(:,3)>0,:); %considering only the clusters of 1s
%sorting order according to cluster sizes
[~,kl] = sort(NS(:,2),'descend');
NS2 = NS(kl,:);
%if not otherwise specified it works only on the maximum cluster for each shuffled matrix
if isempty(perm_max)
    perm_max = 1;
end
%if not otherwise specified threshMC is set as 0.05
if isempty(threshMC)
    threshMC = 0.05;
end
%actual computation
NSVR = zeros(permut,1);
CH = [];
for pp = 1:permut %over total number of permutations
    f = find(ones(n,m)); %finding indexes of matrix elements
    idx_dummy = randperm(length(f)); %creating a permuted array from 1 to length(f)
    idx_2 = f(idx_dummy); %shuffling the indexes in f according to the order in idx_dummy
    r_dummy = zeros(n*m,1); %initialising vector
    r_dummy(f) = P_dummy(idx_2); %taking element in matrix data with index idx_2 (shuffled indexes of the original matrix)
    %rebuild 2D matrix with actual data
    r3 = reshape(r_dummy,[n,m]);
    %mirroring the procedure above for P (finding clusters)
    [~,NSVr,~] = islands(r3);
    hjk = NSVr(NSVr(:,3)>0,2); %considering only the clusters of 1s ad extracting only the cluster sizes
    if isempty(hjk)
        hjk = 0;
    end
    if perm_max == 1
        NSVR(pp) = max(hjk); %only maximum cluster size
    else
        CH = cat(1,CH,hjk); %or all cluster sizes
    end
end
if perm_max == 0
    NSVR = CH;
end
%sorting cluster sizes forming the reference distribution
dummysort = sort(NSVR,'descend');
gfd = dummysort(round(threshMC*length(NSVR) + 1)); %getting threshold value in the reference distribution
%comparing original data with threshold
NS3 = NS2(NS2(:,2)>gfd,1:2);
%calculating p-values
P_val = zeros(size(NS3,1),1);
for ii = 1:size(NS3,1) %over significant clusters
    pk = find(NS3(ii,2) == dummysort); %getting to what cluster size of the reference distribution the cluster size ii is equal
    if ~isempty(pk)
        P_val(ii) = pk(end)/length(dummysort); %dividing that number (last one in case of multiple occurrences) by the number of elements of the reference distribution
    else
        if NS3(ii,2) > dummysort(1) %this is because pk can be = [] both when the original cluster size in NS3 is bigger or smaller than all of the permuted data cluster sizes..
            P_val(ii) = 0; %storing 0 if the original cluster had a size always bigger than all sizes in the reference distribution
        elseif NS3(ii,2) < dummysort(end)
            P_val(ii) = 1; %otherwise if it is smaller than the smallest value in dummysort it stores 1
        else
            [~,olp] = min(abs(dummysort - NS3(ii,2))); %otherwise you ned to find the closest original cluster size to the one in the reference distribution
            P_val(ii) = olp/length(dummysort); %and then calculate the p-value by dividing that index by the length of dummysort
        end
    end
end
%considering only p-values < than threshMC
P_val2 = P_val;
P_val2(P_val > threshMC) = 1;
%preparing output
OUT = cell(size(NS3,1),7);
for ii = 1:size(NS3,1)
    OUT(ii,1) = {ii}; %cluster #
    OUT(ii,2) = {NS3(ii,2)}; %cluster size
    OUT(ii,3) = {P_val2(ii)}; %cluster p-value
    [dI,dJ] = find(AG == NS3(ii,1)); %finding location of the cluster
    OUT(ii,4) = {t1(min(dI))}; %min for t1 dimension
    OUT(ii,5) = {t1(max(dI))}; %max for t1 dimension
    OUT(ii,6) = {t2(min(dJ))}; %min for t2 dimension
    OUT(ii,7) = {t2(max(dJ))}; %max for t2 dimension
end

end


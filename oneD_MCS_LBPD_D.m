function [ sign_clust ] = oneD_MCS_LBPD_D( P, p, max_lab, permnum, MCS_thresh )

% One-dimensional (1D) Monte Carlo simulation.
% It calculates whether a cluster of significant values in 1 dimension (for
% example a time-course of p-values) is significant or not, assuming that
% random significant values do not occurr in clusters.
% It uses the following not built-in function: bwconncomp.
% Simulations using uniformly distributed random numbers (rand built-in
% function in Matlab) showed that increasing the length of the time-series
% the probability of having significant clusters increases. This suggests
% that this function should be used with time-series that cover the specific
% event(s) you are interested in and that if the time-series is built by
% many time-points (such as 1000) it could be a good idea to lower the
% threshold for considering a cluster significant (E.g. if you are testing
% the association between a specific event occurring between time t and time
% t + 10, I would test a time-series starting something like 10-30
% time-points before t and finishing 10-30 time-points after t + 10. This is
% just a very pragmatical suggestion to give an idea of what I am
% describing.. the key concept is that if in that example you take a
% time-series ending 1000 time-points after t, you dramatically increase the
% probability of getting false significative clusters (in that case you
% should consider to at least considerably lower the threshold for MCS if max_lab = 0).
% Otherwise, setting max_lab = 1, so using only the maximum cluster size
% of the permuted data for building the reference distribution, should also
% provide a good solution, without requiring to dramatically lower the MCS
% threshold).



% INPUT:    -P:          vector of p-values
%           -p:          alpha level used for binarizing the p-value vector.
%                        set 0 if you do not want to binarise the vector P (e.g. if you
%                        already binarized it).
%           -max_lab:    1 for having null distribution built on maximum
%                        cluster size of each permutation only; 0 for null distribution
%                        built on length of all clusters of each permutation
%           -permnum:    number of permutations
%           -MCS_thresh: threshold for MCS (e.g. .001)

% OUTPUT:   -sign_clust: clusters of P that survived the Monte Carlo
%                        simulation (1st column: cluster size, 2nd column: p-value, 3th column: temporal extent of original cluster)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 25/12/2018
% Leonardo Bonetti, Aarhus, DK, 25/06/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 04/10/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%checking if user wants to binarize the vector or if he/she provided an
%already binarized one
if p ~= 0
    dummyt = double(P < p); %binarizing p-values vector
else
    dummyt = P;
end
CCt = bwconncomp(dummyt); %getting the 'islands'
tclusters = CCt.PixelIdxList; %clusters of 'islands'
[so,io] = sort(cellfun('length',tclusters),'descend'); %sorting and getting indexes
sign_clust = cell(length(so'),3);
% sign_clust(:,1) = so';

rclusters = [];
for ii = 1:permnum
    pd = length(P); %randomizing the vector..
    idx_pd = randperm(pd); %randomization of indexes of pd (length(P))
    r = P(idx_pd); %randomising the element in P
    if p ~= 0 %if P vector was not already binarized..
        r2 = double(r < p); %binarizing r
    else
        r2 = r; %otherwise just passing r to r2..
    end
    CC = bwconncomp(r2);  %function for getting the 'islands'
    if max_lab == 1
        rclusters = cat(1,rclusters,max(cellfun('length',CC.PixelIdxList))); %getting the maximum value of the simulated (randomized) cluster lengths
    else
        rclusters = cat(2,rclusters,cellfun('length',CC.PixelIdxList));
    end
end

so2 = sort(rclusters,'descend');
% so2(permnum/100*5); % 5% of the number of permutations indicated by the user
for ii = 1:length(so)
    ab = find(so(ii)==so2); %finding where the original ii cluster size is located in the distribution of sizes of randomized data
    sign_clust(ii,1) = {so(ii)}; %storing original ii cluster size
    sign_clust(ii,3) = tclusters(io(ii)); %storing temporal extent of original ii cluster
    if ~isempty(ab)
        sign_clust(ii,2) = {ab(end)/length(rclusters)}; %storing p-value after Monte Carlo simulation
    else
        if so(ii) > so2(1) %this is because ab can be = [] both when the original cluster size in so is bigger or smaller than all of the permuted data cluster sizes.. so if it is bigger
            sign_clust(ii,2) = {0}; %storing 0 if the original cluster had a size always bigger than all sizes in the reference distribution
        elseif so(ii) < so2(end)
            sign_clust(ii,2) = {1}; %otherwise if it is smaller it stores 1
        else
            [~,olp] = min(abs(so2 - so(ii))); %otherwise you need to find the closest original cluster size to the one in the reference distribution
            sign_clust(ii,2) = {olp/length(so2)}; %and then calculate the p-value by dividing that index by the length of dummysort
        end
    end
    
    if sign_clust{ii,2} >= MCS_thresh %if the p-value after MCS is equal or higher than the MCS threshold
        sign_clust(ii,2) = {[]}; %deleting p-value
    end
end
%removing the entire original cluster when the p-value was erased (because was equal or higher than the MCS threshold)
sign_clust = sign_clust(~any(cellfun('isempty',sign_clust),2),:);


end


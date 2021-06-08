function [ MAG_clust_pos, MAG_clust_neg, GRAD_clust ] = MEG_sensors_MonteCarlosim_LBPD_D( S )

% It identifies spatio-temporal clusters (by using islands3D_LBPD_D.m) on a binarized
% 3D matrix and then computes Monte Carlo simulation for testing the significance of the
% clusters in the original data. The clusters always refer to only one direction of the
% contrast (e.g. cond1 > cond2 or viceversa, but not both at the same
% time).
% As conceivable, here we have the problem of the sign when dealing with
% magnetometers. Since mag have double polarity there is a
% non-extremely-likely-but-definitely-real possibility of detecting an
% (e.g.) positive cluster that is the sum of distinct significant patterns (ones where
% cond1 > cond2 and the other one where cond2 > cond1). Dealing with this specific
% issue on a computational level is quite tricky, therefore I would
% strongly recommend to run these functions and then carefully check the
% output and look for consistency between mag and grad significant clusters.
% If your effect is real and decently powerful, you are very likely to find
% coherent cluster configurations between mag and grad. If you do not find
% that you should probably question the reliability of your results (and in
% that case, if you really want to go anyway further, I would suggest to
% focus on gradiometers only).
% If submitted both gradiometers and magnetometers, it does one computation
% for gradiometers and two for magnetometers, since magnetometers have both
% positive and negative polarities. In this latter case be careful to
% provide data with significant elements only for (at first)
% positive and (then) negative magnetometers.
% This function should be used after MEG_sensors_MCS_reshapingdata_LBPD_D.m.



%  INPUT:   -S.data:        actual data (binarized p-values) (double).
%                           Sensorhor x sensorvert x time x channeltype.
%                           The 4th dimension must contain:
%                           -1st: gradiometers (GRAD_data in the output of MEG_sensors_MCS_reshapingdata_LBPD_D.m)
%                           -2nd: magnetometers positive (MAG_data_pos in the output of MEG_sensors_MCS_reshapingdata_LBPD_D.m)
%                           -3th: magnetometers negative (MAG_data_neg in the output of MEG_sensors_MCS_reshapingdata_LBPD_D.m)
%                           If you want to submit only gradiometers or
%                           magnetometers data, you need to fill the
%                           missing data with zeros and then using
%                           S.sensortype correctly to use only the channels
%                           that you want.
%           -S.sensortype:  3 for both
%                           1 for gradiometers
%                           2 for magnetometers
%                           If empty [], default is 3.
%           -S.MEGlayout:   approximation of MEG layout (magnetometers) (double)
%           -S.permut:      number of permutations for Monte Carlo simulation
%           -S.clustmax:    set 1 for only max cluster size of each permutation MCS (more strict).
%                           set 0 for every size of each cluster detected for each permutation MCS (less strict).
%           -S.permthresh:  threshold for considering significant the size of the original clusters
%                           (expressed between 0 and 1; e.g. 5% = 0.05)

%  OUTPUT:  -MAG_clust_pos: clusters for positive magnetometers data
%           -MAG_clust_neg: clusters for negative magnetometers data
%           -GRAD_clust:    clusters for gradiometers data





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 13/08/2019
% Leonardo Bonetti, Aarhus, DK, 28/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    








%input variables (think if you want to check them later..) and setting some
%parameters
% pthresh = S.pthresh;
datafull = S.data;
MEGlayout = S.MEGlayout;
permut = S.permut;
clustmax = S.clustmax;
permthresh = S.permthresh;
MAG_clust_pos = [];
MAG_clust_neg = [];
GRAD_clust = [];
%if no specification, assigning 0 by default (therefore calculations on
%both gradiometers and magnetometers)
if isempty(S.sensortype)
    S.sensortype = 3;
end
if S.sensortype == 1
    dummm = 1;
elseif S.sensortype == 2
    dummm = [2:3];
else
    dummm = [1:3];
end
dfg{1} = 'gradiometers';
dfg{2} = 'magnetometers positive';
dfg{3} = 'magnetometers negative';
    
%actual computation
for dd = dummm %over gradiometers and magnetometers
    %getting data
    data = datafull(:,:,:,dd);
    dummyt = data; %this is not useful but I leave it for historical reasons of this function..
    dummyt(dummyt==0) = NaN; %assigning 'NaN' to 0 values (that corresponds to the non-significant channels)
    [~,k_info_t,Idx1t,Idx2t,Idx3t] = islands3D_LBPD_D(dummyt); %extracting binary patterns (islands)
    if k_info_t ~= 0
        disp(['clusters found in tvalues are n = ' num2str(length(k_info_t(:,1)))])
        %permutations
        SS = size(data); %barbaric way to extract the proper indices that will then be randomised (we cannot randomise values mixing them with the zero padding previously done for MEGlayout)
        MEGdum = zeros(SS(1),SS(2),SS(3));
        for ii = 1:SS(3)
            MEGdum(:,:,ii) = MEGlayout;
        end
        ind = find(MEGdum > 0); %getting indices of non-0s channels in MEGlayout
        clustperm = [];
        maxclust = zeros(permut,1);
        for pp = 1:permut
            disp(['computing permutation ' num2str(pp) ' - ' dfg{dd}])
            idx_dummy = randperm(length(ind)); %create a permuted array from 1 to length(ind)
            idx_2 = ind(idx_dummy); %shuffle the indexes in ind according to the order in idx_dummy
            r_dummy = zeros(SS(1)*SS(2)*SS(3),1); %initialise vector
            r_dummy(ind) = dummyt(idx_2); %taking element in matrix data with index idx_2 
            r_resh = reshape(r_dummy,[SS(1),SS(2),SS(3)]); %reshaping the vector into a matrix shaped as the original data
            r_resh(r_resh==0) = NaN; %removing 0 values as done above for the actual data
            [~,k_info_r,~,~,~] = islands3D_LBPD_D(r_resh); %calculating clusters (islands) for permuted data
            if k_info_r == 0 %checking if some clusters have been found
                maxclust(pp,1) = 0;
                k_info_r = [0, 0, 0];
            else
                maxclust(pp,1) = max(k_info_r(:,2)); %storing maximum size of permuted clusters
            end
            clustperm = cat(1,clustperm,k_info_r); %storing all of the sizes of the detected permuted clusters
        end
        %sorting parameters of permuted clusters
        if clustmax == 1 %maximum sizes of permuted clusters (this usually tends to a decente normal distribution a bit shifted towards left..); to my understanding it is fine
            dummysort = sort(maxclust,'descend');
        else %or all of the sizes of the permuted clusters (this usually tends to only the right side of a normal distribution with a very peaky mean..); to my understanding it is fine
            dummysort = sort(clustperm(:,2),'descend');
        end
        threshfinal = dummysort(floor((permthresh*100*length(dummysort)/100) + 1)); %getting the final threshold for significance
        %storing the output
        d = find(k_info_t(:,2) > threshfinal); %looking for significant clusters (according to their sizes)
        cc2 = zeros(length(d),1);
        for kk = 1:length(d) %over significant clusters
            cc = find(k_info_t(d(kk),2) > dummysort); %getting a specific p-value by.. finding how many times original significant cluster was larger than permuted ones
            cc2(kk) = (length(dummysort) - length(cc))/length(dummysort); %then getting false positive (amount of simulated clusters - times when original significant cluster was larger than permuted ones) and dividing them by amount of simulated clusters
        end
        PP = cell(length(d),4);
        for hh = 1:length(d)
            PP(hh,1) = {hh}; %new ID of the clusters
            PP(hh,2) = {k_info_t(d(hh),2)}; %sizes
            PP(hh,4) = {cc2(hh,1)}; %p-values of each cluster
            dj = cell(length(Idx2t{d(hh)}),2); %this is for getting a better output.. so having the number of channels and then the significant time-points
            dj2 = cell(length(Idx2t{d(hh)}),1); %here proper names are stored
            countdj = 0;
            for kkk = 1:length(Idx2t{d(hh)})
                fxd = find(S.MEGlayout(Idx1t{d(hh)}(kkk),Idx2t{d(hh)}(kkk)) == cell2mat(dj(:,1))); %looking if the channel has already been stored
                if isempty(fxd) %if not store it
                    countdj = countdj + 1;
                    dj(countdj,1) = {S.MEGlayout(Idx1t{d(hh)}(kkk),Idx2t{d(hh)}(kkk))}; %store the channel ID
                    %leonardo Aarhus 23/08/2019
                    if dd == 1 %this is for storing proper GRAD and MAG labels
                        if length(num2str((S.MEGlayout(Idx1t{d(hh)}(kkk),Idx2t{d(hh)}(kkk))) + 1)) == 3
                            dj2(countdj,1) = {['MEG0' num2str((S.MEGlayout(Idx1t{d(hh)}(kkk),Idx2t{d(hh)}(kkk))) + 1) '+0' num2str((S.MEGlayout(Idx1t{d(hh)}(kkk),Idx2t{d(hh)}(kkk))) + 2)]}; %store the channel
                        else
                            dj2(countdj,1) = {['MEG' num2str((S.MEGlayout(Idx1t{d(hh)}(kkk),Idx2t{d(hh)}(kkk))) + 1) '+' num2str((S.MEGlayout(Idx1t{d(hh)}(kkk),Idx2t{d(hh)}(kkk))) + 2)]}; %store the channel
                        end
                    else
                        if length(num2str(S.MEGlayout(Idx1t{d(hh)}(kkk),Idx2t{d(hh)}(kkk)))) == 3
                            dj2(countdj,1) = {['MEG0' num2str(S.MEGlayout(Idx1t{d(hh)}(kkk),Idx2t{d(hh)}(kkk)))]};
                        else
                            dj2(countdj,1) = {['MEG' num2str(S.MEGlayout(Idx1t{d(hh)}(kkk),Idx2t{d(hh)}(kkk)))]};
                        end
                    end
                    %leonardo Aarhus 23/08/2019
                    dj(countdj,2) = {Idx3t{d(hh)}(kkk)}; %store the first significant time-point of the channel
                else %if it is already there, then use its index fxd
                    dj(fxd,2) = {cat(1,dj{fxd,2},Idx3t{d(hh)}(kkk))}; %and store the other significant time-points                    
                end
            end
%             lsd = dj(:,1); %storing original channels IDs.. (temporary variable.. just because I like to have proper MEG channel names in the 1st column and the original IDs in the 3th one..)
            dj(:,1) = dj2; %overwriting proper MEG channel names over IDs.. (1st column)
%             dj(:,3) = lsd; %storing in output file the original channel IDs (3th column)
            djnew = dj(~any(cellfun('isempty',dj),2),:); %remove the possible (and very likely to exist) empty cell
%             for lkl = 1:size(djnew,1) %providing a 4-digit format to the original IDs.. it is going to be useful later..
%                 if length(num2str(djnew{lkl,3})) == 3 %if it is a 3-figure number
%                     djnew(lkl,3) = {['0' num2str(djnew{lkl,3})]}; %apply an additional 0
%                 else
%                     djnew(lkl,3) = {num2str(djnew{lkl,3})}; %otherwise simply convert the number into character..
%                 end
%             end
            PP(hh,3) = {(djnew(:,:))}; %storing channels and corresponding time-points for output structure
        end
    else
        PP = {[]};
    end
    %storing clusters in the proper channel type (and polarity) specification
    if dd == 1
        GRAD_clust = PP;
    elseif dd == 2
        MAG_clust_pos = PP;
    else
        MAG_clust_neg = PP;
    end  
end


end


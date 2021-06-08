function [ Imag, Smag, Igrad, Sgrad ] = MEG_sensors_combining_magclustsign_LBPD_D( GRAD_clust, MAG_clust_pos, MAG_clust_neg )

% It calculates the temporal similarities between gradiometers and positive and negative
% magnetometer clusters outputted by MEG_sensors_MonteCarlosim_LBPD_D.m.
% The temporal similarities is estimated by the difference between weighted
% means of temporal occurrence of mag positive and negative clusters (1)
% and between gradiometer clusters and mean of most temporal similar
% magnetometers clusters, with reference to the positive ones (2).
% This is a useful approximation of which gradiometer and magnetometers clusters (positive
% and negative) have been probably elicited by the same event and can be
% used as plotting input for MEG_sensors_MCS_plottingclusters_LBPD_D.m.



%   INPUT:  -GRAD_clust:        gradiometer clusters outputted by   
%                               MEG_sensors_MonteCarlosim_LBPD_D.m
%           -MAG_clust_pos:     positive magnetometer clusters outputted by   
%                               MEG_sensors_MonteCarlosim_LBPD_D.m
%           -MAG_clust_pos:     negative magnetometer clusters outputted by   
%                               MEG_sensors_MonteCarlosim_LBPD_D.m

%   OUTPUT: -Imag:              double matrix with 2 columns describing best
%                               approximation of coupled magnetometer clusters (pos and neg):
%                                   -1st clm: positive clusters
%                                   -2nd clm: negative clusters
%           -Smag:              difference between temporal weighted mean of
%                               the two clusters indicated by I.
%           -Igrad:             double matrix with 2 columns describing best
%                               approximation of coupled clusters
%                               (gradiometer clusters with mean of most
%                               temporal similar magnetometer clusters,
%                               here with reference to positive mag
%                               clusters):
%                                   -1st clm: gradiometer clusters
%                                   -2nd clm: averaged magnetometer
%                                       clusters (with # reference to positives
%                                       ones..)
%           -Sgrad:             same concept but values describing the
%                               corresponding difference




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 13/09/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%looking for simultaneous magnetometers clusters
%gradiometer clusters
wmgp = zeros(1,size(GRAD_clust,1));
% mp = zeros(1,size(MAG_clust_pos,1)); %calculating normal mean gave almost
% the same results as weighted mean.. however I think that calculating the
% weighted mean is a more correct approach
for ii = 1:size(GRAD_clust,1)
    tem = [];
    for jj = 1:size(GRAD_clust{ii,3},1)
        tem = cat(1,tem,GRAD_clust{ii,3}{jj,2}); %extracting and concatenating time-points
    end
    [a,b] = hist(tem,unique(tem)); %finding occurrencies of each time-point for the iith cluster
    wmgp(ii) = a*b/(sum(a)); %weighted mean
%     mp(ii) = mean(b);
end
%magnetometers positive clusters
wmp = zeros(1,size(MAG_clust_pos,1));
% mp = zeros(1,size(MAG_clust_pos,1)); %calculating normal mean gave almost
% the same results as weighted mean.. however I think that calculating the
% weighted mean is a more correct approach
for ii = 1:size(MAG_clust_pos,1)
    tem = [];
    for jj = 1:size(MAG_clust_pos{ii,3},1)
        tem = cat(1,tem,MAG_clust_pos{ii,3}{jj,2}); %extracting and concatenating time-points
    end
    [a,b] = hist(tem,unique(tem)); %finding occurrencies of each time-point for the iith cluster
    wmp(ii) = a*b/(sum(a)); %weighted mean
%     mp(ii) = mean(b);
end
%magnetometers negative clusters
wmn = zeros(1,size(MAG_clust_neg,1));
% mn = zeros(1,size(MAG_clust_neg,1));
for ii = 1:size(MAG_clust_neg,1)
    tem = [];
    for jj = 1:size(MAG_clust_neg{ii,3},1)
        tem = cat(1,tem,MAG_clust_neg{ii,3}{jj,2});
    end
    [a,b] = hist(tem,unique(tem)); %finding occurrencies of each element of tem
    wmn(ii) = a*b/(sum(a));
%     mn(ii) = mean(b);
end
%calculating difference between temporal similarities of each couple of
%clusters
dwm = zeros(length(wmp),length(wmn));
Smag = zeros(length(wmp),1);
Imag = zeros(length(wmp),2);
mmg = zeros(length(wmp),1);
for ii = 1:length(wmp) %over positive cluster temporal weighted mean
    for jj = 1:length(wmn) %over negative cluster temporal weighted mean
        dwm(ii,jj) = abs(wmp(ii) - wmn(jj)); %abs of the difference
    end
    [s,i] = sort(dwm(ii,:)); %sorting the difference
    Smag(ii,1) = s(1); %storing the smallest difference
    Imag(ii,1) = ii; %storing iith positive cluster
    Imag(ii,2) = i(1); %storing negative cluster characterized by smallest temporal difference with iith positive cluster
    mmg(ii,1) = (wmp(Imag(ii,1)) + wmn(Imag(ii,2)))/2; %temporal mean of most similar magnetometers clusters
end
%looking for similarities between gradiometer clusters and average of most similar
%magnetomer clusters (it is possible also to look for all combinations
%between grad and mag pos and mag neg, but I decided to proceed in this
%2-step way to avoid providing too many non-very-redable results to the user
dwm = zeros(length(wmgp),length(mmg));
Sgrad = zeros(length(wmgp),1);
Igrad = zeros(length(wmgp),2);
for ii = 1:length(wmgp) %over positive cluster temporal weighted mean
    for jj = 1:length(mmg) %over negative cluster temporal weighted mean
        dwm(ii,jj) = abs(wmgp(ii) - mmg(jj)); %abs of the difference
    end
    [s,i] = sort(dwm(ii,:)); %sorting the difference
    Sgrad(ii,1) = s(1); %storing the smallest difference
    Igrad(ii,1) = ii; %storing iith positive cluster
    Igrad(ii,2) = i(1); %storing negative cluster characterized by smallest temporal difference with iith positive cluster
end



end


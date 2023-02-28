function [ PDn, PDnchan ] = Extract_MEGSensor_Information_LBPD_D( stats )

% It extracts information about the clusters of the provided MEG sensors clusters.
% To be used in connection with 'MEG_sensors_MonteCarlosim_LBPD_D' (provide
% its output as input to the current function).




%  INPUT:   -stats: :       file outputted by 'MEG_sensors_MonteCarlosim_LBPD_D'
%                           (e.g. MAG_clust_pos).


%  OUTPUT:  -PDn:           table with information on the clusters:
%                           -cluster #
%                           -size
%                           -number of channels involved
%                           -MCS p-value
%                           -max T-val
%                           -time (1)
%                           -time (end)
%           -PDnchan:       table with significant MEG channels for each cluster and correspondent significant time-points






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Aarhus, DK, 28/02/2023


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    







PDn = cell2table(stats); %overall table
PDn2 = cell(104,((size(stats,1)-1)*3)+(size(stats,1)-2));
cnt = 0;
for ii = 1:4:(size(stats,1)-1)*4 %over significant clusters
    cnt = cnt + 1;
    dum = stats{cnt + 1,3}; %extracting new cluster cnt
    PDn2{1,ii} = ['Cluster' num2str(cnt)]; %providing the cluster name
    PDn2{2,ii} = ['MEG channel']; %providing the cluster name
    PDn2{2,ii+1} = ['time 1']; %providing the cluster name
    PDn2{2,ii+2} = ['time end']; %providing the cluster name
    PDn2(3:size(dum,1)+2,ii) = dum(:,1); %MEG channel name
    for jj = 1:size(dum,1) %over significant MEG channels for cluster cnt      
        PDn2(jj+2,ii+1) = {dum{jj,2}(1)}; %getting time 1
        PDn2(jj+2,ii+2) = {dum{jj,2}(end)}; %getting time 2 (end)
    end
end
PDnchan = cell2table(PDn2); %table for MEG channels and correspondent significant time-points

end


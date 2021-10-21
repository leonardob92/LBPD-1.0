function [ ] = STC_plottingtimeseries_LBPD( S )

% It plots timeseries for a provided parcellation in a number of datasets.
% This is thought to be used as follows:
% 1)Estimate best parcellation using "FunctionalSpatialClustering_voxels2ROIs_LBPD_D" function on the full dataset.
% 2)Load different part of the above datasets (or potentially different datasets).
% 3)Use the current function for plotting the different (sub)datasets timeseries.






%   INPUT:

%           -S.data:                        cell array with one (sub)dataset in each cell.
%                                           The (sub)datasets are matrices with data (voxels x time-points).
%           -S.timep:                       time-samples to do the computation on (useful if the data matrix has too large time-series)
%           -S.time:                        time in seconds (its length can be bigger than S.timep, but not bigger!)
%          
%           -S.idk3:                        parcellation obtained after temporal clustering (from "SpatioTemporalClustering_voxels2ROIs_LBPD" function)
%                                           Leave empty [] if you do not want to plot it.
%           -S.kkk2:                        parcellation obtained after spatio-temporal clustering (from "SpatioTemporalClustering_voxels2ROIs_LBPD" function)
%                                           Leave empty [] if you do not want to plot it.
%           -S.clustsol:                    number of temporal clusters for different kmeans clustering solutions
%                                           This should correspond to the number of temporal clusters used in "SpatioTemporalClustering_voxels2ROIs_LBPD" function
%           -S.clspatvect:                  number of spatio-temporal clusters for different kmeans clustering solutions
%                                           This should correspond to the number of spatio-temporal clusters used in "SpatioTemporalClustering_voxels2ROIs_LBPD" function
%           -S.xlim:                        limits for time in seconds
%           -S.ylim:                        limits for amplitude. It can also be a custom number corresponding to the number of parcels that you want to plot.
%                                           Leave empty for automatic calculation of the limits.
%           -S.dimc:                        values of line-width for temporal clustering timeseries (one value for each dataset)
%           -S.col_ind:                     if you want two or more dataset to have the same color (and I guess different linewidth) for the plot, you can specify that here.
%                                           Indicate which datasets should have the same color (e.g. S.col_ind = {[1 2],[3 4]}).
%                                           Leave it empty [] not to do that.
%           -S.color:                       Specify the colors (correspondent to S.col_ind.. e.g. S.color = {'r','b'}
%           -S.legs:                        1 to show legends; 0 otherwise..


%   OUTPUT: -plotting..






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 25/06/2020


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%preparing colormaps
cmp = colormap(lines(S.clustsol(end))); %extracting some colours from a colormap
cmpk = colormap(lines(sum(S.clspatvect))); %extracting some colours from a colormap
% %preparing width of the lines
% azz = S.azz;
% dimc = zeros(length(S.data),1); %producing a colormap with proressive darker shades of red
% q = azz/(length(S.data)+1);
% for uu = 1:length(S.data)
%     dimc(uu,1) = azz-(q*uu);
% end
dimc = S.dimc;
timesel = S.time(S.timep);


%temporal clustering plotting
if ~isempty(S.idk3)
    figure
    for iii = 1:S.clustsol
        if ~isempty(find(S.idk3==iii)) %this since the user may set some voxels to 0 to exclude them.. (usually useful to exclude a whole parcel)
            for dd = 1:length(S.data) %over different datasets to be plotted together for visual comparison purposes (maybe RTs fast and RTs slow..)
                P2 = S.data{dd}(:,S.timep(1):S.timep(end));
                %getting voxels of cluster ii and combining them into a single ROI
                jg = mean(P2(S.idk3==iii,:),1); %for now trying with mean..
                if isempty(S.col_ind) %this is to have the same color for some categories of datasets (e.g. cond1 RTs fast and slow and cond2 RTs fast and slow)
                    plot(timesel,jg,'Color',cmp(iii,:),'DisplayName',[num2str(iii)],'LineWidth',dimc(dd))
                else
                    cop = 1;
                    while isempty(find(S.col_ind{cop}==dd))
                        cop = cop + 1;
                    end
                    plot(timesel,jg,'Color',S.color{cop},'DisplayName',[num2str(iii)],'LineWidth',dimc(dd))
                end
                hold on
            end
        end
    end
    grid minor
    if S.legs == 1
        legend('show')
    end
    set(gcf,'color','w');
    if ~isfield(S,'xlim')
        xlim([timesel(1) timesel(end)])
    else
        xlim(S.xlim)
    end;
    if ~isempty(S.ylim)
        ylim(S.ylim)
    end
end

%spatio-temporal clustering plotting
if ~isempty(S.kkk2)
    if ~isempty(S.clspatvect)
        figure
        for iiik = 1:sum(S.clspatvect)
            if ~isempty(find(S.kkk2==iiik)) %this since the user may set some voxels to 0 to exclude them.. (usually useful to exclude a whole parcel)
                %getting voxels of cluster ii and combining them into a single ROI
                for dd = 1:length(S.data) %over different datasets to be plotted together for visual comparison purposes (maybe RTs fast and RTs slow..)
                    P2 = S.data{dd}(:,S.timep(1):S.timep(end));
                    jgj = mean(P2(S.kkk2==iiik,:),1); %for now trying with mean..
                    if isempty(S.col_ind) %this is to have the same color for some categories of datasets (e.g. cond1 RTs fast and slow and cond2 RTs fast and slow)
                        plot(timesel,jgj,'Color',cmpk(iiik,:),'DisplayName',[num2str(iiik)],'LineWidth',dimc(dd))
                    else
                        cop = 1;
                        while isempty(find(S.col_ind{cop}==dd))
                            cop = cop + 1;
                        end
                        plot(timesel,jgj,'Color',S.color{cop},'DisplayName',[num2str(iiik)],'LineWidth',dimc(dd))
                    end
                    hold on
                end
            end
        end
        grid minor
        if S.legs == 1
            legend('show')
        end
        set(gcf,'color','w');
        if ~isfield(S,'xlim')
            xlim([timesel(1) timesel(end)])
        else
            xlim(S.xlim)
        end;
        if ~isempty(S.ylim)
            ylim(S.ylim)
        end
    end
end


end


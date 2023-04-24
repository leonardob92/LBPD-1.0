function O = waveplot_groups_server_LBPD(S)
O = []; 
% It plots waveform for multiple groups of participants. Typically used for
% plotting time series of MEG source ROIs according to different groups.



%  INPUTS:
%                    -S.data:           Data that you want to plot
%                                       ROIs (sensors) x time-points x subjects x conditions
%                    -S.condition_n:    Plot condition(s) (e.g. S.condition = 3 or = [1 3])
%                    -S.ROIs_labels:    Labels for ROIs (channels)
%                    -S.ROI_n:          indices of ROIs to be plotted
%                    -S.groups:         group labels
%                    -S.gsubj:          cell containing the indices of the subjects. Each vector contains the
%                                       indices of the subjects belonging to one of the groups.
%                                       !! The order of the vectors must correspond to the order of the group labels!!
%                    -S.signtp:         cell array (maximum 1 x 1 because it makes no sense to plot significant time-windows if you have more than one condition or ROI within one single plot..)
%                                       The significant time points are in seconds (clusters x time (begin and end)).
%                                       Leave empty for not plotting any significant time-window.
%                    -S.legendl:        set 1 for having the legend; 0 otherwise
%                    -S.x_lim:          Set temporal (x) limits in s
%                    -S.y_lim:          Set amplitudel (y) limits in fT
%                    -S.colorline:      Select colorline for each group (provide 3-element array RGB color code)
%                    -S.time_real:      time in seconds
%                    -S.conds:          cell array with condition names (for title)
%                    -S.STE:            1 = dot lines for standard errors
%                                       2 = shadows.
%                                       Do not provide field will result in dot lines.
%                    -S.transp:         value between 0 and 1 to indicate the level of transparecny of the shadow for standard errors (meaningfull only if S.STE = 2)

%  OUTPUT:           Figure






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Aarhus, DK, 01/03/2023
% Leonardo Bonetti, Bologna, Italy, 05/04/2023


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 








%extracting data
COL = S.colorline;
time = S.time_real;
data = S.data;
condition = S.condition_n;
ROIN = S.ROIs_labels;
roin = S.ROI_n;
gsubj = S.gsubj;
glab = S.groups;
signtp = S.signtp;
if ~isfield(S,'STE')
    STE = 1;
else
    STE = S.STE;
    transp = S.transp;
end

if length(COL) < length(condition) * length(glab)* length(ROIN) %if user did not provide enough colors..
    COL = colormap(lines(length(condition) * length(glab)* length(ROIN))); %extracting some colours from a colormap
end

%creating figure
figure
patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
cnt = 0;
for cc = 1:length(condition) %over conditions
    for ii = 1:length(roin) %over ROIs
        for gg = 1:length(glab)
            cnt = cnt + 1;
            %mean
            ananmean = squeeze(nanmean(data(roin(ii),:,gsubj{gg},condition(cc)),3));
            %standard error (positive)
            stdey = squeeze(nanstd(data(roin(ii),:,gsubj{gg},condition(cc)),0,3)./sqrt(length(gsubj{gg})));
            if length(condition) == 1
                plot(time,ananmean,'Color',COL(cnt,:),'LineWidth',1.5,'DisplayName',[glab{gg} ' ' ROIN{ii}])
            else
                plot(time,ananmean,'Color',COL(cnt,:),'LineWidth',1.5,'DisplayName',[glab{gg} ' ' ROIN{ii} ' ' S.conds{condition(cc)}])
            end
            hold on
            if STE == 1
                plot(time,ananmean + stdey,':','Color',COL(cnt,:),'LineWidth',1,'HandleVisibility','off')
                hold on
                %standard error (negative)
                plot(time,ananmean - stdey,':','Color',COL(cnt,:),'LineWidth',1,'HandleVisibility','off')
                hold on
            elseif STE == 2
                fill([time fliplr(time)],[ananmean + stdey fliplr(ananmean - stdey)],COL(cnt,:), 'FaceAlpha', transp,'linestyle','none');
                hold on
            end
        end
    end
end
%extracting limits of y-axis
if isfield(S,'y_lim')
    if ~isempty(S.y_lim)
        ylims = S.y_lim;
%         ylim(S.y_lim)
    else
        ylims = get(gca,'YLim');
    end
else
    ylims = get(gca,'YLim');
end

close all
%plotting greay areas to show significant time-points
if ~isempty(signtp)
    sgf1 = signtp{1,1}; %extracting clusters for ROIs ii and condition cc
end
if ~isempty(signtp)
    for iii = 1:size(sgf1,1) %over significant clusters
        sgf2 = sgf1(iii,:); %extracting cluster iii
        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.5,'HandleVisibility','off')
        hold on
    end
end
%plotting all the waveform again to have them not covered by the grey area.. not the most elegant solution but it is still quick and fine..
cnt2 = 0;
for cc = 1:length(condition) %over conditions
    for ii = 1:length(roin) %over ROIs
        for gg = 1:length(glab)
            cnt2 = cnt2 + 1;
            %mean
            ananmean = squeeze(nanmean(data(roin(ii),:,gsubj{gg},condition(cc)),3));
            %standard error (positive)
            stdey = squeeze(nanstd(data(roin(ii),:,gsubj{gg},condition(cc)),0,3)./sqrt(length(gsubj{gg})));
            if length(condition) == 1
                plot(time,ananmean,'Color',COL(cnt2,:),'LineWidth',1.5,'DisplayName',[glab{gg} ' ' ROIN{ii}])
            else
                plot(time,ananmean,'Color',COL(cnt2,:),'LineWidth',1.5,'DisplayName',[glab{gg} ' ' ROIN{ii} ' ' S.conds{condition(cc)}])
            end
            hold on
            if STE == 1
                plot(time,ananmean + stdey,':','Color',COL(cnt2,:),'LineWidth',1,'HandleVisibility','off')
                hold on
                %standard error (negative)
                plot(time,ananmean - stdey,':','Color',COL(cnt2,:),'LineWidth',1,'HandleVisibility','off')
                hold on
            elseif STE == 2
                fill([time fliplr(time)],[ananmean + stdey fliplr(ananmean - stdey)],COL(cnt2,:), 'FaceAlpha', transp,'linestyle','none','HandleVisibility','off');
                hold on
            end
        end
    end
end
grid minor
box on
if S.legendl == 1
    legend('show')
end
if isfield(S,'x_lim')
    xlim(S.x_lim)
end
ylim(ylims)
if length(condition) == 1
    s = [];
    for ii = 1:length(condition) %trick to get the title..
        s = strcat(s,S.conds{condition(ii)});
    end
    title(s)
end
set(gcf,'color','w')

end

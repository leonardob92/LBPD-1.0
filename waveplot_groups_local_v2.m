function O = waveplot_groups_local_v2(S)
O = []; 
% It plots waveform for multiple groups of participants. Typically used for
% plotting time series of MEG source ROIs according to different groups.

% OBS!! FUNCTION TO BE RUN ON LOCAL COMPUTERS (NOT ON AARHUS SERVER)
% THIS HAS BEEN UPDATED FOR BETTER PLOTTING OF SIGNIFICANT TIME-POINTS (OXFORD, 03/05/2024) 


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
%                    -S.signtp:         cell array with significant time points (in seconds).
%                                       This must correspond to the ID coding specified in S.signtp_col
%                                       Leave empty for not plotting any significant time-window.
%                    -S.lineplot:       first number (positive) indicating where the lines showing significant time-windows should be placed.
%                                       second number is either 1 or 2, depending on the orientation of the line (upper or lower part of the plot). 
%                                       For this option, you need to provide the color codes in S.signtp_col.
%                                       Leave empty [] for plotting shaded colors instead. 
%                    -S.signtp_col:     vector with indices corresponding to the color code/marker for each significant time-windows.
%                                       The color code are the ones provided in S.colorsign (!). 
%                                       Do not provide or leave empty [] for assigning automatically grey color.
%                                       Unssupported for subplots.
%                    -S.legendl:        set 1 for having the legend; 0 otherwise
%                    -S.x_lim:          Set temporal (x) limits in s
%                    -S.y_lim:          Set amplitudel (y) limits in fT
%                    -S.colorline:      Select colorline for each group (provide 3-element array RGB color code)
%                    -S.colorsign:      Colors or markers for the lines showing the significant clusters (not to be confounded with the colors of the waveforms, provided in S.colorline)
%                                       Markers should be provided in a cell array e.g. S.colorsign = {'*';'+'}; Colors in double, e.g. S.colorsign = [0.9 0.1 1; 1 0 0.1]
%                    -S.time_real:      time in seconds
%                    -S.conds:          cell array with condition names (for title)
%                    -S.STE:            1 = dot lines for standard errors
%                                       2 = shadows.
%                                       Do not provide field will result in dot lines.
%                    -S.transp:         value between 0 and 1 to indicate the level of transparecny of the shadow for standard errors (meaningfull only if S.STE = 2)
%                    -S.ii:             number for the figure (useful if you want to plot several figures)



%  OUTPUT:           Figure






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Aarhus, DK, 01/03/2023
% Leonardo Bonetti, Bologna, Italy, 05/04/2023
% Leonardo Bonetti, Oxford, UK, 03/05/2024


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 











%extracting data
COL = S.colorline;
% COL2 = S.colorsign;
time = S.time_real;
data = S.data;
if ~isempty(find(data==0))
    data(data==0) = nan; %assigning NaN to 0 values to later ignore them when average and standard deviations are computed; of course, you should be aware that if you have 0s that may be a problem.. this is why I inserted the warning
    warning('your data contains 0s..')
end
condition = S.condition_n;
ROIN = S.ROIs_labels;
roin = S.ROI_n;
gsubj = S.gsubj;
glab = S.groups;
sgf1 = S.signtp;
if ~isfield(S,'STE')
    STE = 1;
else
    STE = S.STE;
    transp = S.transp;
end

if length(COL) < length(condition) * length(ROIN) * length(glab) %if user did not provide enough colors..
    COL = colormap(lines(length(condition) * length(ROIN) * length(glab))); %extracting some colours from a colormap
end
if isfield(S,'lineplot') %if you want subplots you need to provide a field for it
    if ~isempty(S.lineplot) %and specify numbers
        linep = S.lineplot;
    else %otherwise no subplots
        linep = [];
    end
else %same
    linep = [];
end
if isfield(S,'ii') %if no figure number is provided
    figii = S.ii;
else
    figii = 1; %assigning 1 by default
end

%creating figure
figure(1000)
% patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
cnt = 0;
for cc = 1:length(condition) %over conditions
    for ii = 1:length(roin) %over ROIs
%         if ~isempty(signtp)
%             sgf1 = signtp{ii,cc}; %extracting clusters for ROIs ii and condition cc
%         end
        for gg = 1:length(glab)
            cnt = cnt + 1;
            %mean
             %mean
            ananmean = squeeze(mean(data(roin(ii),:,gsubj{gg},condition(cc)),3,'omitnan'));
            %standard error (positive)
            stdey = squeeze(std(data(roin(ii),:,gsubj{gg},condition(cc)),0,3,'omitnan')./sqrt(length(gsubj{gg})));
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
                fill([time fliplr(time)],[ananmean + stdey fliplr(ananmean - stdey)],COL(cnt,:), 'FaceAlpha', transp,'linestyle','none','HandleVisibility','off');
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
close(figure(1000))
figure(figii)
%plotting greay areas to show significant time-points
if ~isempty(sgf1)
    for iii = 1:length(sgf1) %over significant clusters
        sgf2 = sgf1{iii}; %extracting cluster iii
        if isfield(S,'signtp_col') %if you provide field for different colors for significant time-points
            if ~isempty(S.signtp_col) %if the field is not empty
                if isempty(linep) %if you do not want the significant time-windows plotted as lines
                    if ~iscell(S.colorsign) %if you provide colors
                        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],S.colorsign(S.signtp_col(iii),:),'EdgeColor','none','FaceAlpha',.1,'HandleVisibility','off') %you get them with shaded colors
                    else %otherwise patch color
                        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],[.85 .85 .85],'EdgeColor','none','FaceAlpha',.1,'HandleVisibility','off') %you get them with shaded colors
                    end
                else
                    for ff = round(find(sgf2(1)==time)):round(find(sgf2(2)==time))
                        if min(ananmean) < 0 %if the minimum is negative, I subtract which reduces the space occupied by the lines in the plot)
                            if ~iscell(S.colorsign) %if you provide colors
                                plot([sgf2(1),sgf2(end)],[ylims(linep(2)) - (ylims(linep(2)))/linep(1)*S.signtp_col(iii),ylims(linep(2)) - (ylims(linep(2)))/linep(1)*S.signtp_col(iii)],'color',S.colorsign(S.signtp_col(iii),:),'LineStyle','-','HandleVisibility','off','LineWidth',1.5)
                            else %you provide markers
                                plot([sgf2(1),sgf2(end)],[ylims(linep(2)) - (ylims(linep(2)))/linep(1)*S.signtp_col(iii),ylims(linep(2)) - (ylims(linep(2)))/linep(1)*S.signtp_col(iii)],'color','k','LineStyle','-','Marker',S.colorsign{S.signtp_col(iii)},'MarkerSize',5,'HandleVisibility','off','LineWidth',1.5)
                            end
                        else %otherwise, for obtaining the same results, I have to add
                            if ~iscell(S.colorsign) %if you provide colors
                                plot([sgf2(1),sgf2(end)],[ylims(linep(2)) + (ylims(linep(2)))/linep(1)*S.signtp_col(iii),ylims(linep(2)) + (ylims(linep(2)))/linep(1)*S.signtp_col(iii)],'color',S.colorsign(S.signtp_col(iii),:),'LineStyle','-','HandleVisibility','off','LineWidth',1.5)
                            else %you provide markers
                                plot([sgf2(1),sgf2(end)],[ylims(linep(2)) + (ylims(linep(2)))/linep(1)*S.signtp_col(iii),ylims(linep(2)) + (ylims(linep(2)))/linep(1)*S.signtp_col(iii)],'color','k','LineStyle','-','Marker',S.colorsign{S.signtp_col(iii)},'MarkerSize',5,'HandleVisibility','off','LineWidth',1.5)
                            end
                        end
                    end
                    hold on
                end
            else %otherwise assigning classicl grey
                patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],[.85 .85 .85],'EdgeColor','none','FaceAlpha',.5,'HandleVisibility','off')
            end
        else %same
            patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],[.85 .85 .85],'EdgeColor','none','FaceAlpha',.5,'HandleVisibility','off')
        end
        hold on
    end
end

%plotting all the waveform again to have them not covered by the grey area.. not the most elegant solution but it is still quick and fine..
cnt = 0;
for cc = 1:length(condition) %over conditions
    for ii = 1:length(roin) %over ROIs
%         if ~isempty(signtp)
%             sgf1 = signtp{ii,cc}; %extracting clusters for ROIs ii and condition cc
%         end
        for gg = 1:length(glab)
            cnt = cnt + 1;
            %mean
             %mean
            ananmean = squeeze(mean(data(roin(ii),:,gsubj{gg},condition(cc)),3,'omitnan'));
            %standard error (positive)
            stdey = squeeze(std(data(roin(ii),:,gsubj{gg},condition(cc)),0,3,'omitnan')./sqrt(length(gsubj{gg})));
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
                fill([time fliplr(time)],[ananmean + stdey fliplr(ananmean - stdey)],COL(cnt,:), 'FaceAlpha', transp,'linestyle','none','HandleVisibility','off');
                hold on
            end
        end
    end
end
grid minor
if S.legendl == 1
    legend('show')
end
if isfield(S,'x_lim')
    xlim(S.x_lim)
end
ylim(ylims)
% if length(condition) == 1
%     s = [];
%     for ii = 1:length(condition) %trick to get the title..
%         s = strcat(s,S.conds{condition(ii)});
%     end
%     title(s)
% end
set(gcf,'color','w')
box on

end

function [] = waveform_plotting_local_v2(S)

% It plots waveform for multiple conditions (one group of participants).
% Typically used for plotting time series of MEG channel/source, locally.



%  INPUTS:
%                    -S.data:           Data that you want to plot
%                                       ROIs (sensors) x time-points x subjects x conditions
%                    -S.condition_n:    Plot condition(s) (e.g. S.condition = 3 or = [1 3])
%                    -S.ROIs_labels:    Labels for ROIs (channels)
%                    -S.ROI_n:          indices of ROIs to be plotted together 
%                    -S.signtp:         cell array 1 x significant time-windows.
%                                       The significant time points are in seconds (e.g. two time-windows = {1,2}; S.signtp{1} = [1.5 2.3]; S.signtp{2} = [3.4 3.6]).
%                                       Leave empty for not plotting any significant time-window.
%                    -S.signtp_col:     vector with indices corresponding to the color code for each significant time-windows.
%                                       The color code are the ones provided in S.colorline (!). 
%                                       Do not provide or leave empty [] for assigning automatically grey color.
%                                       Unssupported for subplots.
%                    -S.ii:             number for the figure (useful if you want to plot several figures)
%                    -S.lineplot:       first number (positive) indicating where the lines showing significant time-windows should be placed.
%                                       second number is either 1 or 2, depending on the orientation of the line (upper or lower part of the plot). 
%                                       For this option, you need to provide the color codes in S.signtp_col.
%                                       Leave empty [] for plotting shaded colors instead. 

%%%%%%%%%%%% NOT YET SUPPORTED %%%%%%%%%%%%
%                    -S.subplot:        array with two entries corresponding to rows and column of subplots.
%                                       Do not provide field or leave empty [] if you do not want subplots.
%                                       Please, note that subplot here is supporting only for plotting single ROIs in different subplots of the figure. 
%%%%%%%%%%%% UNTIL HERE %%%%%%%%%%%%

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
% Leonardo Bonetti, Oxford, UK, 24/04/2023


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 








%extracting data
COL = S.colorline;
time = S.time_real;
data = S.data;
if ~isempty(find(data==0))
    data(data==0) = nan; %assigning NaN to 0 values to later ignore them when average and standard deviations are computed; of course, you should be aware that if you have 0s that may be a problem.. this is why I inserted the warning
    warning('your data contains 0s..')
end
condition = S.condition_n;
ROIN = S.ROIs_labels;
roin = S.ROI_n;
sgf1 = S.signtp;
if ~isfield(S,'STE')
    STE = 1;
else
    STE = S.STE;
    transp = S.transp;
end

if length(COL) < length(condition) * length(ROIN) %if user did not provide enough colors..
    COL = colormap(lines(length(condition) * length(ROIN))); %extracting some colours from a colormap
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
%     
% else
    %creating figure
    figure(1000)
    cnt = 0;
    for ii = 1:length(roin) %over ROIs
        for cc = 1:length(condition) %over conditions
            cnt = cnt + 1;
            %mean
            ananmean = squeeze(mean(data(roin(ii),:,:,condition(cc)),3,'omitnan'));
            %standard error (positive)
            stdey = squeeze(std(data(roin(ii),:,:,condition(cc)),0,3,'omitnan')./sqrt(size(data,3)));
            plot(time,ananmean,'Color',COL(cnt,:),'LineWidth',1.5,'DisplayName',[S.conds{condition(cc)} ' ' ROIN{ii}])
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
                        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],COL(S.signtp_col(iii),:),'EdgeColor','none','FaceAlpha',.1,'HandleVisibility','off') %you get them with shaded colors
                    else 
                        for ff = round(find(sgf2(1)==time)):round(find(sgf2(2)==time))
                            if min(ananmean) < 0 %if the minimum is negative, I subtract which reduces the space occupied by the lines in the plot)
                                plot([sgf2(1),sgf2(end)],[ylims(linep(2)) - (ylims(linep(2)))/linep(1)*S.signtp_col(iii),ylims(linep(2)) - (ylims(linep(2)))/linep(1)*S.signtp_col(iii)],'color',COL(S.signtp_col(iii),:),'LineStyle','-','HandleVisibility','off','LineWidth',1.5)
                            else %otherwise, for obtaining the same results, I have to add
                                plot([sgf2(1),sgf2(end)],[ylims(linep(2)) + (ylims(linep(2)))/linep(1)*S.signtp_col(iii),ylims(linep(2)) + (ylims(linep(2)))/linep(1)*S.signtp_col(iii)],'color',COL(S.signtp_col(iii),:),'LineStyle','-','HandleVisibility','off','LineWidth',1.5)
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
    cnt2 = 0;
    for ii = 1:length(roin) %over ROIs
        for cc = 1:length(condition) %over conditions
            cnt2 = cnt2 + 1;
            %mean
            ananmean = squeeze(mean(data(roin(ii),:,:,condition(cc)),3,'omitnan'));
            %standard error (positive)
            stdey = squeeze(std(data(roin(ii),:,:,condition(cc)),0,3,'omitnan')./sqrt(size(data,3)));
            plot(time,ananmean,'Color',COL(cnt2,:),'LineWidth',1.5,'DisplayName',[S.conds{condition(cc)} ' ' ROIN{ii}])
            hold on
            if STE == 1
                plot(time,ananmean + stdey,':','Color',COL(cnt2,:),'LineWidth',1,'HandleVisibility','off')
                hold on
                %standard error (negative)
                plot(time,ananmean - stdey,':','Color',COL(cnt2,:),'LineWidth',1,'HandleVisibility','off')
                hold on
            elseif STE == 2
                fill([time fliplr(time)],[ananmean + stdey fliplr(ananmean - stdey)],COL(cnt2,:), 'FaceAlpha', transp,'linestyle','none');
                hold on
            end
        end
    end
    box on
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
    grid minor
    % fontname("Helvetica Neue")
% end
end

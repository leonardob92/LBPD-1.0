function [] = plot_wave_conditions(S)

% It plots waveforms for a set of different conditions, with standard deviations.



%  INPUT:       -S.data:           Data that you want to plot (4D matrix)
%                                  It has the following configuration :
%                                  sensors/parcels x time points x subjects x conditions
%               -S.chans_index:    index(es) of the timeseries (e.g. of the brain parcels).
%                                  A series of indices returns the average over those parcels/channels
%               -S.condition_n:    conditions to be plotted (double)
%               -S.conds:          conditions labels
%               -S.signtp:         Significant time points (in seconds).
%                                  e.g. S.signtp = {[0.1 0.3];[0.5 0.55];[1.1 1.35]}
%                                  S.signtp = {[]} not to highlight any time-window
%               -S.legendl:        set 1 for having the legend; 0 otherwise
%               -S.STE:            set 1 if you want to plot the standard error (shadows, heavier, suggested in case of very clear differences).
%                                  set 2 if you want to plot the standard error (dot lines, lighter, suggested for more subtle differences)
%                                  set 0 not to plot any standard error
%               -S.time:           vector with time (in seconds)
%               -S.x_lim:          time limits for plotting (in seconds) e.g. [-0.1 2.3]
%               -S.y_lim:          limit for amplitude (e.g. [0 2])
%               -S.colorline:      colors for each condition (e.g. S.colorline = {'r','b'}


%  OUTPUT:      -plots..






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 19/01/2021


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









% Make sure that all needed fields are filled
if isempty(S.chans_index)
    warning('Define the channels/parcels to be plotted!!');
end

if isempty(S.condition_n)
    warning('Define the conditions to be plotted!!');
end

if isempty(S.colorline) || ~isfield(S,'colorline')
    S.colorline = {'k', 'r', 'b', 'g', 'm', 'c'}; % If no colors are defined, set default colors
elseif length(S.colorline)<length(S.conds)
    warning('The n. of colors does not coincide with the n. of conditions: the same color might be used for two or more conditions.');
end

if isempty(S.signtp{1})
    warning('Significant time points were not defined!');
end




figure
%preparing grey areas for significant time-windows
%Add gray patch to highlight significant time points (if specified)
if ~isempty(S.signtp{1})
    %plotting the significant time-points with gray shades
    patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
    for ll = 1:length(S.signtp) %over the significant time-windows
        sgf2 = S.signtp{ll}; %Find the mean value across the significant time points and create the patch around it (before and after)
        if ~isempty(S.y_lim)
            if S.STE == 1 %highlighting more the gray shadow if the standard errors are also plotted (just for visualization purposes)
                patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[S.y_lim(1) S.y_lim(1) S.y_lim(2) S.y_lim(2)], ...
                    patch_color,'EdgeColor','none','FaceAlpha',.8, 'HandleVisibility', 'off')
            else
                patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[S.y_lim(1) S.y_lim(1) S.y_lim(2) S.y_lim(2)], ...
                    patch_color,'EdgeColor','none','FaceAlpha',.5, 'HandleVisibility', 'off')
            end
            hold on
        else
            %getting maximum and minimum values of the timeseries that will then be plotted
            caz = zeros(2,length(S.condition_n));
            for ab = 1:length(S.condition_n)
                caz(1,ab) = max(squeeze(mean(mean(mean(S.data(S.chans_index,:,:,S.condition_n(ab)),1),3),4)));
                caz(2,ab) = min(squeeze(mean(mean(mean(S.data(S.chans_index,:,:,S.condition_n(ab)),1),3),4)));
            end
            ylimsmax = max(max(caz));
            ylimsmin = min(min(caz));
            %adding and subtracting a percentage (suggested 10-15%) of the range amplitude (maximum and minimum values of the timeseries) from the maximum and minimum values of the same timeseries
            ylims= [ylimsmin-(((ylimsmax - ylimsmin)/100)*15) ylimsmax+(((ylimsmax - ylimsmin)/100)*15)];
            if S.STE == 1 %highlighting more the gray shadow if the standard errors are also plotted (just for visualization purposes)
                patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)], ...
                    patch_color,'EdgeColor','none','FaceAlpha',.8, 'HandleVisibility', 'off')
            else
                patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)], ...
                    patch_color,'EdgeColor','none','FaceAlpha',.5, 'HandleVisibility', 'off')
            end
            hold on
        end
    end
end
%plotting the data
for aa = 1:length(S.condition_n)
    %Calculating the standard error mean (SEM) for each experimental group
    data_SEM = nanstd(squeeze(mean(mean(S.data(S.chans_index,:,:,S.condition_n(aa)),4),1)),0,2)/sqrt(size(S.data,3));
    amean = squeeze(mean(mean(mean(S.data(S.chans_index,:,:,S.condition_n(aa)),1),3),4));
    if S.STE == 1
        fill([S.time fliplr(S.time)],[amean+data_SEM fliplr(amean-data_SEM)],S.colorline{aa}, 'FaceAlpha', 0.005,'linestyle','none');
        %plot actual mean
        hold on;
    elseif S.STE == 2
        plot(S.time,amean + data_SEM',':','Color',S.colorline{aa},'LineWidth',1); %upper std error
        hold on
        plot(S.time,amean - data_SEM',':','Color',S.colorline{aa},'LineWidth',1); %lower std error
        hold on
    end
    p1(aa) = plot(S.time,amean,'Color',S.colorline{aa},'linewidth',2.5); %% change color or linewidth to adjust mean line
    hold on
end
if S.legendl == 1
    %sanguinoso trick to get the proper legend..
    legend([p1],S.conds)
end
set(gcf,'Color','w')
if ~isempty(S.y_lim)
    ylim(S.y_lim);
elseif exist('ylims','var') == 1
    ylim(ylims);
end
xlim(S.x_lim);
grid minor
%Set axis labels
xlabel('ms')
%         ylabel('fT')


end






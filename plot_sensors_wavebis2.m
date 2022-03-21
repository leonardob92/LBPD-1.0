function [] = plot_sensors_wavebis2(S)


%%% TO BE CHECKED AND SISTEMATA!!



%% Plot waveforms at sensors

% Plot waveform at MEG sensors (magnetometers or gradiometers) of interest. The
% signal can be plotted for single conditions (to be specified) or as the average across the selected conditions.
% Waveforms can be plotted for single channels or for all channels (if S.chans_index=0).
% Topographies can also be plotted.


% INPUT SETTINGS
% The argument of the function is the structure S. Before launching the
% function, the following fields of S need to be specified:


% MANDATORY FIELDS:
%                    -S.data:           Data that you want to plot.
%                                       data_mat is a 4D structure containing the data extracted from the
%                                       spm object.
%                                       It has the following configuration :
%
%                                                sensors x time points x subjects x conditions
%
%                                       (info can be visualized with size(data_mat))
%                                       Please provide data as MAG - GRAD - MAG - GRAD - MAG - GRAD.. 
%
%                    -S.condition_n:    Plot single condition (e.g. S.condition=3)
%                                       or average waveform across conditions (e.g. S.condition=2:7)
%                    -S.chans_index:    Plot waveform at single channel or the average across multiple
%                                       channels (specify channel index - you can find the channel labels
%                                       in 'chanlabels'); set 0 to plot all channels
%                    -S.groups:         group labels
%                    -S.sensors:        -1 for magnetometer, 0 for gradiometers (gradiometers have already
%                                       been combined).
%                                       This is relevant only if you want to plot all channels.. otherwise, just submit the specific
%                                       channel indices that you want to plot..
%                    -S.gsubj:          cell containing the indices of the subjects. Each vector contains the
%                                       indices of the subjects belonging to one of the groups.
%                                       !! The order of the vectors must correspond to the order of the group
%                                       labels!!
%                    -S.signtp:          Significant time points (in seconds). If the field is filled, the significant time
%                                        points will be highlighted with a gray patch. Significant time points
%                                        are stored in a 1 x time-points array. If no significant time points are specified, define it as [].
%                                        (only in the case of single channel or average across selected channels)
%                    -S.chanlabels:     chanlabels provided as a 1 x channels cell array with the label of each channel in each cell
%                    -S.legendl:        set 1 for having the legend; 0 otherwise
%                    -S.STE:            set 1 if you want to plot the standard error (shadows, heavier, suggested in case of very clear differences).
%                                       set 2 if you want to plot the standard error (dot lines, lighter, suggested for more subtle differences)
%                                       set 0 not to plot any standard error

%

% OPTIONAL FIELDS:
%
%                    -S.x_lim:           Set temporal (x) limits in s
%                    -S.y_lim:           Set amplitudel (y) limits in fT
%                    -S.colorline:       Select colorline for each group
%
%                    -S.time_real:       time points - needed to plot the
%                                        standard error (only in the case of single channel or average across selected channels)


% !! When single channels (or average across specified channels) are
% plotted, fancy features such as standard error and a gray patch
% highlighting are added to the waveform


% -------------------------------------------------------------------------
% LEONARDO BONETTI, SILVIA E.P. BRUZZONE
% July 2020
% Center for Music in the Brain (MIB), Aarhus, Denmark
% leonardo.bonetti@clin.au.dk
% silviaelisabettaportis.bruzzone@studenti.units.it
% -------------------------------------------------------------------------


% INPUT SETTINGS - example

% S.data = data_mat; %data to be plotted
% S.chans_index = 12; %plot waveform at single channel or the average across multiple channels (specify channel index - you can find the channel labels in 'chanlabels'); set 0 to plot all channels
% S.condition_n = 7; %plot single condition (e.g. S.condition=3) or average waveform across conditions (e.g. S.condition=2:7)
% % S.condition_label = {'Standard','Pitch','Slide','Intensity', 'Localization', 'Timbre', 'Rhythm'};
% S.groups = {'mus1','nonm1','mus0','nonm0'}; %Group labels
% S.colorline = {'k', 'r', 'b', 'g'}; %Select colorline for each group
% S.sensors = -1; % -1 for magnetometer, 0 for gradiometers (gradiometers have already been combined)
% S.x_lim = []; % Set x limits
% S.y_lim = [-120 120]; %Set y limits
% S.signtp = [0.1 0.3]; % A gray patch will highlight the area between the two time limits
% S.gsubj = gsubj;


% Make sure that all needed fields are filled
if isempty(S.chans_index)
    warning('Define the channels to be plotted!!');
end

if isempty(S.condition_n)
    warning('Define the conditions to be plotted!!');
end

if isempty(S.colorline)
    S.colorline = {'k', 'r', 'b', 'g', 'm', 'c'}; % If no colors are defined, set default colors
elseif length(S.colorline)<length(S.groups)
    warning('The n. of colors does not coincide with the n. of groups: the same color might be used for two or more groups.');
end

if S.chans_index == 0
    if S.sensors<-1 || S.sensors>0
        warning('Define sensors to be plotted!! Magnetometers = -1; Gradiometers = 0');
    end
end

% if isempty(S.signtp)
%     warning('Significant time points were not defined!');
% end


% channels = S.data(:,1,1,1);


%%% THE FOLLOWING LINES SHOULD BE BETTER WORKED OUT
% Assign sensortype to a new variable
countwave = S.sensors; %THIS IS FOR OSL STUFF (E.G. MAG-GRAD-MAG-GRAD)

%THIS IS FOR GRAD-GRAD-GRAD-MAG-MAG-MAG
% if S.sensors == -1
%     countwave = 102;
% else
%     countwave = 0;
% end
%%% UNTIL HERE




% left_mag = S.chans_index; %channels index. The confusing "left_mag" name is simply due to the laziness (also referred to as "historical reasons") in changing the name here and downstream this line
% mat_left_mag = zeros(length(left_mag),Sz(2),Sz(3),length(S.condition_n)); %Initialize a matrix with the following dimensions: channels x time x subjects x conditions




%Plot ALL channels
if sum(double(S.chans_index)) == 0 %if S.chans_index = 0, plot the waveform at all channels
    for jjk = 1:(round((length(S.chanlabels)/2)/9))+1 %N. of figures (depends on the number of channels and the n. of subplots in each figure)
        clf(figure(jjk)); %Clear current figure
        figure(jjk); %create new figure        
        for iik = 1:9
            
            %%% THE FOLLOWING LINES SHOULD BE BETTER WORKED OUT
            countwave = countwave + 2; % %THIS IS FOR OSL STUFF (E.G. MAG-GRAD-MAG-GRAD)
%             countwave = countwave + 1; %THIS IS FOR GRAD-GRAD-GRAD-MAG-MAG-MAG
            %%% UNTIL HERE

            if countwave <= length(S.chanlabels)
                subplot(3,3,iik); %3x3 plots within each window = 9 subplots (iik = tot. n. subplots)
                for aa = 1:length(S.groups)
                    % Plot average across conditions (if only one condition specified, only that condition will be plotted)
                    plot(S.time_real,squeeze(nanmean(nanmean(S.data(countwave,:,S.gsubj{aa},S.condition_n),3),4)),S.colorline{aa},'LineWidth',1,'DisplayName',S.groups{aa}); %average waveform across all subjects and all conditions (nanmean(var,3)) = nanmean across third dimension; (nanmean(var,4))=nanmean across fourth dimension)
                    % xlim(S.x_lim);
                    if ~isempty(S.y_lim)
                        ylim(S.y_lim);
                    end
                    hold on

                    ch_id = S.chanlabels(countwave);
                    if length(S.condition_n) == length(S.conditions)
                        title(sprintf('%s', ch_id{:}, ' All conditions'), 'Interpreter', 'none');
                    else
                        con = cell(1,length(S.condition_n));
                        for bb = 1:length(S.condition_n)
                            con{1,bb} = S.conditions(S.condition_n(bb));
                            a = [];
                            for cc =1:length(con)
                                a = horzcat( a, con{cc}{1});
                                %                 title(sprintf('%s', ch_id{:}, con{{bb}}));
                            end
                            title([ch_id{1} ' ' a])
                        end
                    end
                    set(gcf,'Color','w')
                    xlabel('ms')
                    if ~isempty(S.x_lim)
                        xlim(S.x_lim);
                    end
                end
                grid minor
                if S.legendl == 1
                    legend('show');
                end
            else
%                 disp('countwave exceeds n. of labels')
            end

        end
    end
    
    %Plot SINGLE channel(s)
elseif S.chans_index>0 %if channel index is specified, plot waveform at that channel
%     figure('InvertHardcopy','off','Color',[1 1 1], 'Renderer','painters'); %create new figure and apply rendering
    figure
         %Add gray patch to highlight significant time points (if specified)
        if ~isempty(S.signtp)
            %plotting the significant time-points with gray shades
            patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
            sgf2 = S.signtp; %Find the nanmean value across the significant time points and create the patch around it (before and after)
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
                caz = zeros(2,length(S.groups));
                for ab = 1:length(S.groups)
                    caz(1,ab) = max(squeeze(nanmean(nanmean(nanmean(S.data(S.chans_index,:,S.gsubj{ab},S.condition_n),1),3),4)));
                    caz(2,ab) = min(squeeze(nanmean(nanmean(nanmean(S.data(S.chans_index,:,S.gsubj{ab},S.condition_n),1),3),4)));
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
    if countwave <= length(S.chanlabels)
%         countwave = countwave + 2;
        
        for aa = 1:length(S.groups)
            % if more than one channel is specified, compute nanmean across the selected channels
            %Calculating the standard error nanmean (SEM) for each experimental group
            data_SEM = nanstd(squeeze(nanmean(nanmean(S.data(S.chans_index,:,S.gsubj{aa},S.condition_n),4),1)),0,2)/sqrt(length(S.gsubj{aa}));
            ananmean = squeeze(nanmean(nanmean(nanmean(S.data(S.chans_index,:,S.gsubj{aa},S.condition_n),1),3),4));
%             hold on
            %plot standard errors (shades..)
%             fillOut = fill([S.time_real fliplr(S.time_real)],[ananmean+data_SEM fliplr(ananmean-data_SEM)],S.colorline{aa}, 'FaceAlpha', 0.005,'linestyle','none');
            if S.STE == 1
                fill([S.time_real fliplr(S.time_real)],[ananmean+data_SEM fliplr(ananmean-data_SEM)],S.colorline{aa}, 'FaceAlpha', 0.005,'linestyle','none');
                %plot actual nanmean
                hold on;
            elseif S.STE == 2      
                plot(S.time_real,ananmean + data_SEM',':','Color',S.colorline{aa},'LineWidth',0.5); %upper std error
%                 xlim(x_lim_temp);
%                 if ~isempty(y_lim_ampl)
%                     ylim(y_lim_ampl);
%                 end
                hold on
                plot(S.time_real,ananmean - data_SEM',':','Color',S.colorline{aa},'LineWidth',0.5); %lower std error
%                 xlim(x_lim_temp);
%                 if ~isempty(y_lim_ampl)
%                     ylim(y_lim_ampl);
%                 end
                hold on
            end
                %             lineOut = plot(S.time_real,ananmean,S.colorline{aa},'linewidth',1.5,'DisplayName',S.groups{aa}); %% change color or linewidth to adjust nanmean line
                p1(aa) = plot(S.time_real,ananmean,'Color',S.colorline{aa},'linewidth',1.5); %% change color or linewidth to adjust nanmean line
            hold on
        end
        if S.legendl == 1
            %sanguinoso trick to get the proper legend..
            legend([p1],S.groups)
        end
        
        % Set condition name(s) as title to the plot
        if length(S.condition_n)== length(S.conditions)
            %                 title(sprintf('%s', ch_id{:}, ' All conditions'), 'Interpreter', 'none');
            title('All conditions')
        else
%             con = cell(1,length(S.condition_n));
%             for bb = 1:length(S.condition_n)
%                 con{1,bb} = S.conditions(S.condition_n(bb));
%                 a = [];
%                 for cc =1:length(con)
%                     a = horzcat(a,con{cc}{1});
%                     
%                 end
%                 if length(S.chans_index)==1
%                     title([S.chanlabels{S.chans_index(1)} ' ' a])
%                 else
%                     title(a)
%                 end
%             end
        end
        
        set(gcf,'Color','w')
        if ~isempty(S.y_lim)
            ylim(S.y_lim);
        elseif exist('ylims','var') == 1
            ylim(ylims);
        end
        if ~isempty(S.x_lim)
            xlim(S.x_lim);
        end
        grid minor
        %Set axis labels
        xlabel('s')
%         ylabel('fT')
        %Set tick labels for the x axis
%         get(gca, 'XTickLabel');
%         set(gca, 'LabelFontSizeMultiplier',1.2,'XTick',[-0.1 0 0.1 0.2 0.3 0.4],'XTickLabel',...
%             {'-100','0','100','200','300','400'});
    else
%         disp('countwave exceeds n. of labels')
    end
    
end


end 












% 
% 
% %%
% 
% function [] = plot_sensors_wavebis2(S)
% %% Plot waveforms at sensors
% 
% % Plot waveform at MEG sensors (magnetometers or gradiometers) of interest. The
% % signal can be plotted for single conditions (to be specified) or as the average across the selected conditions.
% % Waveforms can be plotted for single channels or for all channels (if S.chans_index=0).
% % Topographies can also be plotted.
% 
% 
% % INPUT SETTINGS
% % The argument of the function is the structure S. Before launching the
% % function, the following fields of S need to be specified:
% 
% 
% % MANDATORY FIELDS:
% %                    -S.data:           Data that you want to plot.
% %                                       data_mat is a 4D structure containing the data extracted from the
% %                                       spm object.
% %                                       It has the following configuration :
% %
% %                                                sensors x time points x subjects x conditions
% %
% %                                       (info can be visualized with size(data_mat))
% %                                       Please provide data as MAG - GRAD - MAG - GRAD - MAG - GRAD.. 
% %
% %                    -S.condition_n:    Plot single condition (e.g. S.condition=3)
% %                                       or average waveform across conditions (e.g. S.condition=2:7)
% %                    -S.chans_index:    Plot waveform at single channel or the average across multiple
% %                                       channels (specify channel index - you can find the channel labels
% %                                       in 'chanlabels'); set 0 to plot all channels
% %                    -S.groups:         group labels
% %                    -S.sensors:        -1 for magnetometer, 0 for gradiometers (gradiometers have already
% %                                       been combined).
% %                                       This is relevant only if you want to plot all channels.. otherwise, just submit the specific
% %                                       channel indices that you want to plot..
% %                    -S.gsubj:          cell containing the indices of the subjects. Each vector contains the
% %                                       indices of the subjects belonging to one of the groups.
% %                                       !! The order of the vectors must correspond to the order of the group
% %                                       labels!!
% %                    -S.signtp:          Significant time points (in seconds). If the field is filled, the significant time
% %                                        points will be highlighted with a gray patch. Significant time points
% %                                        are stored in a 1 x time-points array. If no significant time points are specified, define it as [].
% %                                        (only in the case of single channel or average across selected channels)
% %                    -S.chanlabels:     chanlabels provided as a 1 x channels cell array with the label of each channel in each cell
% %                    -S.legendl:        set 1 for having the legend; 0 otherwise
% %                    -S.STE:            set 1 if you want to plot the standard error (shadows, heavier, suggested in case of very clear differences).
% %                                       set 2 if you want to plot the standard error (dot lines, lighter, suggested for more subtle differences)
% %                                       set 0 not to plot any standard error
% 
% %
% 
% % OPTIONAL FIELDS:
% %
% %                    -S.x_lim:           Set temporal (x) limits in s
% %                    -S.y_lim:           Set amplitudel (y) limits in fT
% %                    -S.colorline:       Select colorline for each group
% %
% %                    -S.time_real:       time points - needed to plot the
% %                                        standard error (only in the case of single channel or average across selected channels)
% 
% 
% % !! When single channels (or average across specified channels) are
% % plotted, fancy features such as standard error and a gray patch
% % highlighting are added to the waveform
% 
% 
% % -------------------------------------------------------------------------
% % SILVIA E.P. BRUZZONE, LEONARDO BONETTI
% % July 2020
% % Center for Music in the Brain (MIB), Aarhus, Denmark
% % silviaelisabettaportis.bruzzone@studenti.units.it
% % leonardo.bonetti@clin.au.dk
% % -------------------------------------------------------------------------
% 
% 
% % INPUT SETTINGS - example
% 
% % S.data = data_mat; %data to be plotted
% % S.chans_index = 12; %plot waveform at single channel or the average across multiple channels (specify channel index - you can find the channel labels in 'chanlabels'); set 0 to plot all channels
% % S.condition_n = 7; %plot single condition (e.g. S.condition=3) or average waveform across conditions (e.g. S.condition=2:7)
% % % S.condition_label = {'Standard','Pitch','Slide','Intensity', 'Localization', 'Timbre', 'Rhythm'};
% % S.groups = {'mus1','nonm1','mus0','nonm0'}; %Group labels
% % S.colorline = {'k', 'r', 'b', 'g'}; %Select colorline for each group
% % S.sensors = -1; % -1 for magnetometer, 0 for gradiometers (gradiometers have already been combined)
% % S.x_lim = []; % Set x limits
% % S.y_lim = [-120 120]; %Set y limits
% % S.signtp = [0.1 0.3]; % A gray patch will highlight the area between the two time limits
% % S.gsubj = gsubj;
% 
% 
% % Make sure that all needed fields are filled
% if isempty(S.chans_index)
%     warning('Define the channels to be plotted!!');
% end
% 
% if isempty(S.groups)
%     warning('Define the groups to be plotted!!');
% end
% 
% if isempty(S.condition_n)
%     warning('Define the conditions to be plotted!!');
% end
% 
% if isempty(S.colorline)
%     S.colorline = {'k', 'r', 'b', 'g', 'm', 'c'}; % If no colors are defined, set default colors
% elseif length(S.colorline)<length(S.groups)
%     warning('The n. of colors does not coincide with the n. of groups: the same color might be used for two or more groups.');
% end
% 
% if S.chans_index == 0
%     if S.sensors<-1 || S.sensors>0
%         warning('Define sensors to be plotted!! Magnetometers = -1; Gradiometers = 0');
%     end
% end
% 
% if isempty(S.signtp)
%     warning('Significant time points were not defined!');
% end
% 
% 
% % channels = S.data(:,1,1,1);
% 
% 
% %%% THE FOLLOWING LINES SHOULD BE BETTER WORKED OUT
% % Assign sensortype to a new variable
% countwave = S.sensors; %THIS IS FOR OSL STUFF (E.G. MAG-GRAD-MAG-GRAD)
% 
% %THIS IS FOR GRAD-GRAD-GRAD-MAG-MAG-MAG
% % if S.sensors == -1
% %     countwave = 102;
% % else
% %     countwave = 0;
% % end
% %%% UNTIL HERE
% 
% 
% 
% 
% % left_mag = S.chans_index; %channels index. The confusing "left_mag" name is simply due to the laziness (also referred to as "historical reasons") in changing the name here and downstream this line
% % mat_left_mag = zeros(length(left_mag),Sz(2),Sz(3),length(S.condition_n)); %Initialize a matrix with the following dimensions: channels x time x subjects x conditions
% 
% 
% 
% 
% %Plot ALL channels
% if sum(double(S.chans_index)) == 0 %if S.chans_index = 0, plot the waveform at all channels
%     
%     for jjk = 1:(round((length(S.chanlabels)/2)/9))+1 %N. of figures (depends on the number of channels and the n. of subplots in each figure)
%         
%         clf(figure(jjk)); %Clear current figure
%         figure(jjk); %create new figure
%         
%         %Add gray patch to highlight significant time points (if specified)
% %         if ~isempty(S.signtp)
% %             %plotting the significant time-points with gray shades
% %             patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
% %             sgf2 = S.signtp; %Find the nanmean value across the significant time points and create the patch around it (before and after)
% %             if ~isempty(S.y_lim)
% %                            patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[S.y_lim(1) S.y_lim(1) S.y_lim(2) S.y_lim(2)], ...
% %                 patch_color,'EdgeColor','none','FaceAlpha',.5, 'HandleVisibility', 'off')
% %                 hold on
% %                 %at the moment there is no implementation for 
% % %             else
% % %                 ylims = max(abs(squeeze(nanmean(nanmean(nanmean(S.data(S.chans_index,:,:,S.condition_n)),3),4)))); %get y limits from the max absolute value among your data
% % %                 ylims= [-ylims-40 ylims+40];
% % %                 patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)], ...
% % %                     patch_color,'EdgeColor','none','FaceAlpha',.5, 'HandleVisibility', 'off')
% % %                 hold on
% %             end
% %             
% %         end
%         
%         
%         for iik = 1:9
%             
%             
%             %%% THE FOLLOWING LINES SHOULD BE BETTER WORKED OUT
%             countwave = countwave + 2; % %THIS IS FOR OSL STUFF (E.G. MAG-GRAD-MAG-GRAD)
% %             countwave = countwave + 1; %THIS IS FOR GRAD-GRAD-GRAD-MAG-MAG-MAG
%             %%% UNTIL HERE
% 
% 
%             
%             if countwave <= length(S.chanlabels)
%                 subplot(3,3,iik); %3x3 plots within each window = 9 subplots (iik = tot. n. subplots)
%                 for aa = 1:length(S.groups)
%                     
%                     % Plot average across conditions (if only one condition specified, only that condition will be plotted)
%                     plot(S.time_real,squeeze(nanmean(nanmean(S.data(countwave,:,S.gsubj{aa},S.condition_n),3),4)),S.colorline{aa},'LineWidth',1,'DisplayName',S.groups{aa}); %average waveform across all subjects and all conditions (nanmean(var,3)) = nanmean across third dimension; (nanmean(var,4))=nanmean across fourth dimension)
%                     grid minor
%                     
%                     % xlim(S.x_lim);
%                     if ~isempty(S.y_lim)
%                         ylim(S.y_lim);
% %                     else
% %                         ylim(ylims);
%                     end
%                     hold on
%                     
%                     ch_id = S.chanlabels(countwave);
%                     
%                     if length(S.condition_n)== length(S.conditions)
%                         title(sprintf('%s', ch_id{:}, ' All conditions'), 'Interpreter', 'none');
%                     else
%                         con = cell(1,length(S.condition_n));
%                         for bb = 1:length(S.condition_n)
%                             con{1,bb} = S.conditions(S.condition_n(bb));
%                             a = [];
%                             for cc =1:length(con)
%                                 a = horzcat( a, con{cc}{1});
%                                 %                 title(sprintf('%s', ch_id{:}, con{{bb}}));
%                             end
%                             title([ch_id{1} ' ' a])
%                         end
%                     end
%                     set(gcf,'Color','w')
%                     grid minor
%                     xlabel('ms')
% %                     ylabel('fT')
%                     get(gca, 'XTickLabel');
%                     set(gca, 'LabelFontSizeMultiplier',1.2,'XTick',[-0.1 0 0.1 0.2 0.3 0.4],'XTickLabel',...
%                         {'-100','0','100','200','300','400'});
%                 end
%             else
%                 disp('countwave exceeds n. of labels')
%             end
%             if S.legendl == 1
%                 legend('show');
%             end
%         end
%     end
%     
%     %Plot SINGLE channel(s)
% elseif S.chans_index>0 %if channel index is specified, plot waveform at that channel
% %     figure('InvertHardcopy','off','Color',[1 1 1], 'Renderer','painters'); %create new figure and apply rendering
%     figure
%          %Add gray patch to highlight significant time points (if specified)
%         if ~isempty(S.signtp)
%             %plotting the significant time-points with gray shades
%             patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
%             sgf2 = S.signtp; %Find the nanmean value across the significant time points and create the patch around it (before and after)
%             if ~isempty(S.y_lim)
%                 if S.STE == 1 %highlighting more the gray shadow if the standard errors are also plotted (just for visualization purposes)
%                     patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[S.y_lim(1) S.y_lim(1) S.y_lim(2) S.y_lim(2)], ...
%                         patch_color,'EdgeColor','none','FaceAlpha',.8, 'HandleVisibility', 'off')
%                 else
%                     patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[S.y_lim(1) S.y_lim(1) S.y_lim(2) S.y_lim(2)], ...
%                         patch_color,'EdgeColor','none','FaceAlpha',.5, 'HandleVisibility', 'off')
%                 end
%                 hold on
%             else
%                 %getting maximum and minimum values of the timeseries that will then be plotted
%                 caz = zeros(2,length(S.groups));
%                 for ab = 1:length(S.groups)
%                     caz(1,ab) = max(squeeze(nanmean(nanmean(nanmean(S.data(S.chans_index,:,S.gsubj{ab},S.condition_n),1),3),4)));
%                     caz(2,ab) = min(squeeze(nanmean(nanmean(nanmean(S.data(S.chans_index,:,S.gsubj{ab},S.condition_n),1),3),4)));
%                 end
%                 ylimsmax = max(max(caz));
%                 ylimsmin = min(min(caz));
%                 %adding and subtracting a percentage (suggested 10-15%) of the range amplitude (maximum and minimum values of the timeseries) from the maximum and minimum values of the same timeseries
%                 ylims= [ylimsmin-(((ylimsmax - ylimsmin)/100)*15) ylimsmax+(((ylimsmax - ylimsmin)/100)*15)];
%                 if S.STE == 1 %highlighting more the gray shadow if the standard errors are also plotted (just for visualization purposes)
%                     patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)], ...
%                         patch_color,'EdgeColor','none','FaceAlpha',.8, 'HandleVisibility', 'off')
%                 else
%                     patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)], ...
%                         patch_color,'EdgeColor','none','FaceAlpha',.5, 'HandleVisibility', 'off')
%                 end
%                 hold on
%             end
%         end
%     if countwave <= length(S.chanlabels)
% %         countwave = countwave + 2;
%         
%         for aa = 1:length(S.groups)
%             % if more than one channel is specified, compute nanmean across the selected channels
%             %Calculating the standard error nanmean (SEM) for each experimental group
%             data_SEM = nanstd(squeeze(nanmean(nanmean(S.data(S.chans_index,:,S.gsubj{aa},S.condition_n),4),1)),0,2)/sqrt(length(S.gsubj{aa}));
%             ananmean = squeeze(nanmean(nanmean(nanmean(S.data(S.chans_index,:,S.gsubj{aa},S.condition_n),1),3),4));
% %             hold on
%             %plot standard errors (shades..)
% %             fillOut = fill([S.time_real fliplr(S.time_real)],[ananmean+data_SEM fliplr(ananmean-data_SEM)],S.colorline{aa}, 'FaceAlpha', 0.005,'linestyle','none');
%             if S.STE == 1
%                 fill([S.time_real fliplr(S.time_real)],[ananmean+data_SEM fliplr(ananmean-data_SEM)],S.colorline{aa}, 'FaceAlpha', 0.005,'linestyle','none');
%                 %plot actual nanmean
%                 hold on;
%             elseif S.STE == 2      
%                 plot(S.time_real,ananmean + data_SEM',':','Color',S.colorline{aa},'LineWidth',0.5); %upper std error
% %                 xlim(x_lim_temp);
% %                 if ~isempty(y_lim_ampl)
% %                     ylim(y_lim_ampl);
% %                 end
%                 hold on
%                 plot(S.time_real,ananmean - data_SEM',':','Color',S.colorline{aa},'LineWidth',0.5); %lower std error
% %                 xlim(x_lim_temp);
% %                 if ~isempty(y_lim_ampl)
% %                     ylim(y_lim_ampl);
% %                 end
%                 hold on
%             end
%                 %             lineOut = plot(S.time_real,ananmean,S.colorline{aa},'linewidth',1.5,'DisplayName',S.groups{aa}); %% change color or linewidth to adjust nanmean line
%                 p1(aa) = plot(S.time_real,ananmean,'Color',S.colorline{aa},'linewidth',1.5); %% change color or linewidth to adjust nanmean line
%             hold on
%         end
%         if S.legendl == 1
%             %sanguinoso trick to get the proper legend..
%             legend([p1],S.groups)
%         end
%         
%         % Set condition name(s) as title to the plot
%         if length(S.condition_n)== length(S.conditions)
%             %                 title(sprintf('%s', ch_id{:}, ' All conditions'), 'Interpreter', 'none');
%             title('All conditions')
%         else
% %             con = cell(1,length(S.condition_n));
% %             for bb = 1:length(S.condition_n)
% %                 con{1,bb} = S.conditions(S.condition_n(bb));
% %                 a = [];
% %                 for cc =1:length(con)
% %                     a = horzcat(a,con{cc}{1});
% %                     
% %                 end
% %                 if length(S.chans_index)==1
% %                     title([S.chanlabels{S.chans_index(1)} ' ' a])
% %                 else
% %                     title(a)
% %                 end
% %             end
%         end
%         
%         set(gcf,'Color','w')
%         if ~isempty(S.y_lim)
%             ylim(S.y_lim);
%         elseif exist('ylims','var') == 1
%             ylim(ylims);
%         end
%         grid minor
%         %Set axis labels
%         xlabel('ms')
% %         ylabel('fT')
%         %Set tick labels for the x axis
%         get(gca, 'XTickLabel');
%         set(gca, 'LabelFontSizeMultiplier',1.2,'XTick',[-0.1 0 0.1 0.2 0.3 0.4],'XTickLabel',...
%             {'-100','0','100','200','300','400'});
%     else
%         disp('countwave exceeds n. of labels')
%     end
%     
% end
% 
% 
% end 
% 





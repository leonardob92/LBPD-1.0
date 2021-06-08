function [  ] = MEG_sensors_MCS_plottingclusters_LBPD_D( S )

% It plots topoplots for gradiometers and magnetometers clusters (given the
% cluster ID that you want to plot).
% Magnetometers positive and negative are plotted together, assuming that
% user is providing cluster IDs that are compatible. In other words,
% magnetometers have both positive and negative clusters that are two parts
% of the same process. It is very likely (but not necessarily true)
% that positive and negative magnetometers have similar cluster sizes and
% occur at (nearly) the same time. Therefore, to have a meaningful
% plot, user should look into the time-ranges of the clusters and selecting
% positive and negative clusters with similar sizes and happening around
% the same time. This requires critical thinking and will be done automatically
% in the future.
% It uses topoplot functionalities of FieldTrip.


%   INPUT:  -S.fieldtrip_mask:  path to fieldtrip mask file used to plot topoplot
%                               (tipicaly provided in the external folder of
%                               this collection of functions).
%           -S.PPall:           1 x 3 cell array containing significant clusters
%                               outputted by MEG_sensors_MonteCarlosim_LBPD_D.m.
%                               The data order MUST be:
%                               -GRAD (1st column)
%                               -MAG pos (2nd col)
%                               -MAG neg (3th col)
%                               Leave the cell empty PPall(x) = {{[]}} if you
%                               do not want to plot x, with x that can be
%                               either GRAD, MAG pos or MAG neg.
%           -S.clustnum_grad:   set gradiometers cluster that you want to
%                               plot.
%           -S.clustnum_magp:   same but for magnetometers positive
%           -S.clustnum_magn:   same but for magnetometers negative
%           -S.zlim:            cluster size limits for plotting.
%                               E.g. [0 19] set the colorbar from 0 to 19
%                               (this number refers to the size of the
%                               clusters).
%                               Leave empty [] for default colorbar (based
%                               on submitted data).
%           -S.time:            vector with time (in seconds) corresponding
%                               to the time-samples used in
%                               MEG_sensors_MonteCarlosim_LBPD_D.m

%   OUTPUT: -plots..




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 28/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








    
%load an example from a previous fieldtrip file obtained after some pre-processing
fieldtrip_example = load([S.fieldtrip_mask '/fieldmask.mat']);
%takes new labels from the fieldtrip example
label = fieldtrip_example.M_timb_lock_gradplanarComb.label;
%set the labels according to the ones (presumibly) in the layout..
%this should work when following the whole pipeline and working with OSL,
%SPM objects and data from 306-channel Elekta Neuromag Triux system..
%it may be that with different systems and/or preprocessing packages you
%may need to update your labels to be able to fit the ones in the FieldTrip
%layout
for ii = 1:204
    label{ii,1} = [label{ii,1}(1:3) label{ii,1}(5:end)];
end
PPall = S.PPall;
labeldummy = label(103:end,1); %getting some labels (labels are for magnetometers only since the computation is done for magnetometers spatial position only.. a bit contort but under control..)
%used in possible warnings later..
dfg{1} = 'gradiometers';
dfg{2} = 'magnetometers positive';
dfg{3} = 'magnetometers negative';

%preparing data
disp('preparing data..')
AVG = cell(length(PPall),1);
for jjh = 1:length(PPall)
    PP = PPall{jjh}; %extracting data
    if ~isempty(PP{1,1})
        avg = zeros(204,1); %initialising data vector
        if jjh == 1 %getting proper cluster number that you want to plot
            clustnum = S.clustnum_grad;
        elseif jjh == 2
            clustnum = S.clustnum_magp;
        else
            clustnum = S.clustnum_magn;
        end
        ff = size(PP{clustnum,3}); %get the number of significant clusters
        %reshaping the data
        for ii = 1:ff(1) %over significant channels
            azz = PP{clustnum,3}{ii,1}; %get the channel
            if jjh == 1 %if gradiometers
%                 if isnan(str2double(azz(7))) %if the channel has 3 figures
%                     azz = str2double(azz(4:6)) - 1; %-1 since the codes below works with magnetometers layout
%                     azz = num2str(azz); %later we need azz in string
%                 else %otherwise (if it has 4 figures..)
                    azz = str2double(azz(4:7)) - 1;
                    if azz/1000 > 1 %this is to check whether the channel is < or > than 1000
                        azz = num2str(azz);
                    else
                        azz = ['0' num2str(azz)]; %if not, add a 0
                    end
%                 end
            else 
%                 if length(azz) == 6 %if the channel has 3 figures (here I need to have this different statement since I have different labels for magnetometers and combined gradiometers)
%                     azz = azz(4:6); %here not -1 since the codes below works with magnetometers layout
%                 else %otherwise (if it has 4 figures..)
                    azz = azz(4:7);
%                 end
            end
            cpp = 0;
            flag = 0;
            while flag == 0
                cpp = cpp + 1;

%                 if cpp < 35
%                     ldm = labeldummy{cpp,1}(5:end);
%                 else
%                     ldm = labeldummy{cpp,1}(4:end);
%                 end
                ldm = labeldummy{cpp,1}(4:end);
                if strcmp(azz,ldm) == 1
                    avg(cpp,1) = length(cell2mat(PP{clustnum,3}(ii,2)));
                    flag = 1;
                end
            end
        end
        AVG(jjh,1) = {avg}; %storing these values..
    else
        warning(['either there are no significant clusters or no data has been provided for.. ' dfg{jjh}])
        AVG(jjh,1) = {[]}; %else assigning empty vector
    end
end

%actual plotting
disp('preparing plotting..')
for jjh = 1:length(AVG)-1 %over gradiometers and magnetometers
    if sum(double(cellfun(@isempty,AVG))) ~= 3 %if data exist..
%     if ~isempty(AVG{jjh,1}) %if data exist..
        avg = AVG{jjh,1}; %extracting data 
        if jjh == 1 %if gradiometers
            if ~isempty(AVG{jjh,1}) %if grad data exist
                %extracting time-points (in seconds) corresponding to maximum
                %and minimum temporal extent of the cluster
                PP2 = PPall{jjh};
                MAX = S.time(max(cellfun(@max,PP2{S.clustnum_grad,3}(:,2))));
                MIN = S.time(min(cellfun(@min,PP2{S.clustnum_grad,3}(:,2))));
                %creates the fake cfg field for the mask
                cfgdummy = fieldtrip_example.M_timb_lock_gradplanarComb.cfg;
                cfgdummy.previous = [];
                data = [];
                %create the mask for the data
                data.cfg = cfgdummy;
                data.time(1,:) = 1;
                data.label = label;
                data.dimord = fieldtrip_example.M_timb_lock_gradplanarComb.dimord;
                data.grad = fieldtrip_example.M_timb_lock_gradplanarComb.grad;
                data.avg = avg;
                %creates the cfg for the actual plotting (gradiometers)
                cfg = [];
                cfg.layout = [S.fieldtrip_mask '/neuromag306cmb.lay'];
                cfg.colorbar = 'yes';
                cfg.xlim = [1 1];
                if ~isempty(S.zlim) %if limit for colorbar (corresponding to size of the clusters) is provided by user
                    cfg.zlim = S.zlim;
                end
                cfg.colormap = 'jet'; %consider to let user decide this..
                figure
                %calling Fieldtrip function for plotting
                ft_topoplotER(cfg,data);
                title(['cluster time-range: ' num2str(MIN) ' to ' num2str(MAX)])
            end
        else %plotting magnetometers
            avg = zeros(204,1);
            avgdum = zeros(204,1);
            plo = 1;
            if ~isempty(AVG{jjh,1}) %if mag pos data exist..
                avg(103:end,1) = AVG{jjh,1}(1:102,1); %positive clusters
                PP2 = PPall{jjh}; %extracting values for positive clusters
                MAX(1) = S.time(max(cellfun(@max,PP2{S.clustnum_magp,3}(:,2)))); %maximum for positive clusters
                MIN(1) = S.time(min(cellfun(@min,PP2{S.clustnum_magp,3}(:,2)))); %minimum for positive clusters
                plo = plo + 1; %not elegant trick for calculating proper min and max if mag pos is not requested
            end
            if ~isempty(AVG{jjh+1,1}) %if mag neg data exist..
                avgdum(103:end,1) = (AVG{jjh+1,1}(1:102,1))*(-1); %negative clusters
                PP2 = PPall{jjh+1}; %extracting values for negative clusters
                MAX(plo) = S.time(max(cellfun(@max,PP2{S.clustnum_magn,3}(:,2)))); %maximum for negative clusters
                MIN(plo) = S.time(min(cellfun(@min,PP2{S.clustnum_magn,3}(:,2)))); %minimum for negative clusters
                avg = avg + avgdum; %putting together positive and negative clusters
            end
            if length(find(avg==0)) ~= 204 %if avg contains actual data and not only 0s
                %creates the fake cfg field for the mask
                cfgdummy = fieldtrip_example.M_timb_lock_gradplanarComb.cfg;
                cfgdummy.previous = [];
                data = [];
                %create the mask for the data (I may not recompute the first
                %lines below.. anyway..)
                data.cfg = cfgdummy;
                data.time(1,:) = 1;
                data.label = label;
                data.dimord = fieldtrip_example.M_timb_lock_gradplanarComb.dimord;
                data.grad = fieldtrip_example.M_timb_lock_gradplanarComb.grad;
                data.avg = avg;
                %creates the cfg for the actual plotting (gradiometers)
                cfg = [];
                cfg.layout = [S.fieldtrip_mask '/neuromag306mag.lay'];
                cfg.colorbar = 'yes';
                cfg.xlim = [1 1];
                if ~isempty(S.zlim) %if limit for colorbar (corresponding to size of the clusters) is provided by user
                    cfg.zlim = S.zlim;
                end
                cfg.colormap = 'jet'; %consider to let user decide this..
                figure
                %calling Fieldtrip function for plotting
                ft_topoplotER(cfg,data);
                title(['cluster time-range: ' num2str(min(MIN)) ' to ' num2str(max(MAX))])    
            end
        end
        set(gcf,'Color','w')
    end
end

end


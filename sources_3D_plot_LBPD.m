function [  ] = sources_3D_plot_LBPD( S )

% Plotting in 3D space neural signal over MEG sensors and in brain sources.
% If you do not provide information about the head model and MEG sensors
% position, previously computed structural/topologic data will be shown.
% For now, the function works in 8mm MNI152-T1 space, but this could be
% easily changed, if needed.
% The function uses FieldTrip and in-house-built codes.
% This function should be used in connection with LBPD.
% Please, contact me for more information (Leonardo Bonetti, leonardo.bonetti@clin.au.dk).



% INPUT:    -S.sensor_l:                    1 for magnetometers; 2 for combined planar gradiometers
%           -S.sources:                     vector with neural signal strenght for each source
%           -S.MEGsensordata:               vector with neural signal strenght over MEG sensors.
%                                           The order must be: first 102 combined planar gradiometer, then 102 magnetometers (as done in FieldTrip topoplot)
%           -S.time:                        time in seconds
%           -S.time_idxs:                   time indices for plotting purposes in seconds (e.g. [0.12 0.18]).
%                                           Interpreted as mean over the requested time-window.
%                                           Provide only 1 value (e.g. [0.12]) if you are interested in 1 specific time-point only.
%           -S.plot3D_l:                    1 for 3D plotting; 0 otherwise
%           -S.size_l:                      1 for plotting different sizes of brain sources with different stregth.
%                                           0 for the same concept, but with color
%           -S.topoplot_label:              1 for plotting also topoplot; 0 for not plotting topoplot
%           -S.scaleextr:                   xtreme values for scaling MEG sensors in plot (e.g.[0.0001 30])
%           -S.perce:                       percentage of active brain voxels to be 3d-plotted
%           -S.vol:                         description of the head volume conduction model as outputted by FieldTrip function ft_prepare_headmodel.
%                                           Leave empty [] to load a previously computed solution (sort of template).
%           -S.poschan:                     MEG channels position (their space must be consistent with the space of vol).
%                                           Leave empty [] to load a previously computed solution (sort of template).
%           -S.possource:                   brain sources position (their space must be consistent with the space of vol and poschan).
%                                           Leave empty [] to load a previously computed solution (sort of template).
%                                           Leaving empty makes sense only if you computed previous steps following LBPD suggestions..
%           -S.scatter_l:                   1 for scatter plot of brain sources


% OUTPUT:   -..:                            3D plotting solution showing neural activity over MEG sensors and brain sources







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 26/02/2021


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%finding correspondance between time-points in seconds and closest time-point in "time-samples"
[~,im1] = min(abs(S.time-S.time_idxs(1))); %time extreme 1
if length(S.time_idxs) == 2 %computing time-sample for time extreme 2, if existent
    [~,im2] = min(abs(S.time-S.time_idxs(2))); %time extreme 2
else %otherwise assigning im1 (and then computing the "mean" later only on 1 value
    im2 = im1;
end

%scatter plot of sources
if S.scatter_l == 1
    figure
    scatter(1:size(S.sources,1),mean(S.sources(:,im1:im2),2))
end

%position of MEG channels
if isempty(S.poschan)
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/plot3D_LF.mat','poschan')
else
    %extracting channels position
    poschan = S.poschan;
end
poschan = poschan(1:3:306,:); %getting only 1 set of coordinates for each of the MEG channels triplets

%preparing data
lf = mean(S.MEGsensordata(:,im1:im2),2);
if S.topoplot_label == 1
    avg = lf; %storing MEGsensor data for topoplot already here, if topoplot is required
end
%preparing MEG sensors
if S.sensor_l == 1 %magnetometers
    lf = lf(103:end); %magnetometers are located after gradiometers, as asked by this function and as done by FieldTrip when plotting topoplots
elseif S.sensor_l == 2 %gradiometers
    lf = lf(1:102); %gradiometers are located first
end
%scaling on the basis of MEG channels strength
rmin = min(abs(lf)); %getting minimum value of leadfield
rmax = max(abs(lf)); %getting maximum value of leadfield
%scaling (absolute value of) leadfield with extreme values (scaleextr) requested by user
LF = (abs(lf)-rmin)./(rmax-rmin).*(S.scaleextr(2)-S.scaleextr(1)) + S.scaleextr(1);
%getting values and indices in different variables for later plotting purposes
LFP = LF(lf>0); %positive values (on MEG sensors)
LFN = LF(lf<0); %negative values (on MEG sensors)
LFPi = find(lf>0); %positive indices (on MEG sensors)
LFNi = find(lf<0); %negative indices (on MEG sensors)

%scaling color or size for plotting the sources in the brain
sources = abs(mean(S.sources(:,im1:im2),2)); %absolute value since we do not care about the sign of the signal for this visualization system
[sortsour,idxf] = sort(sources,'descend');
sources3 = sortsour(1:round(length(sources)/100*S.perce),1); %getting most powerful sources
rmin = min(sources3); %getting minimum value of leadfield
rmax = max(sources3); %getting maximum value of leadfield
%scaling sources
if S.size_l == 1
    sources_scale = (sources3-rmin)./(rmax-rmin).*(S.scaleextr(2)-S.scaleextr(1)) + S.scaleextr(1);
    if isnan(sources_scale)
        sources_scale = ones(length(S.sources_scale),1)*S.scaleextr(2);
    end
else
    sources_scale = (sources3-rmin)./(rmax-rmin).*(1-0) + 0;
    if isnan(sources_scale)
        sources_scale = ones(length(sources_scale),1);
    end
end
cmp = zeros(length(sources),1); %producing a colormap with proressive darker shades of color
cmp(idxf(1:round(length(sources)/100*S.perce)),1) = 1;


%inspecting leadfield model using topoplot (if requested)
%topoplot must be done first for reasons related to Matlab handling figures..
if S.topoplot_label == 1
%     avg = lf;
    fieldtrip_example = load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/fieldmask.mat');
%     label = S.label;
    cfgdummy = fieldtrip_example.M_timb_lock_gradplanarComb.cfg;
    cfgdummy.previous = [];
    data = [];
    data.cfg = cfgdummy;
    data.time = 1;
    %taking new labels from the fieldtrip example
    label = fieldtrip_example.M_timb_lock_gradplanarComb.label;
    %setting labels according to the ones (presumibly) in the layout
    for ii = 1:204
        label{ii,1} = [label{ii,1}(1:3) label{ii,1}(5:end)];
    end
    data.label = label;
    data.dimord = fieldtrip_example.M_timb_lock_gradplanarComb.dimord;
    data.grad = fieldtrip_example.M_timb_lock_gradplanarComb.grad;
    data.avg = avg;
    %creating the cfg for the actual plotting (magnetometers)
    cfg = [];
    if S.sensor_l == 1 %magnetometers
        cfg.layout = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/neuromag306mag.lay';
    else %combined planar gradiometers
        cfg.layout = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/neuromag306cmb.lay';
    end
    cfg.colorbar = 'yes';
    % cfg.xlim = signt; %set temporal limits (in seconds)
    % cfg.zlim = [];
    cfg.colormap = 'jet';
    figure
    ft_topoplotER(cfg,data);
    set(gcf,'Color','w')
    %colormap with white for 0 values
    x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
    colormap(bluewhitered_PD(0,x))
end

%actual 3D plotting
if S.plot3D_l == 1
    figure
    if isempty(S.vol)
        load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/plot3D_LF.mat','vol')
    else
        vol = S.vol;
    end
    %brain with triangles
    ft_plot_vol(vol, 'edgecolor', [0 0 0], 'facealpha', 0);
    hold on
    %postion (coordinates) of sources from leadfield computation
    if isempty(S.possource)
        load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/plot3D_LF.mat','possource')
    else
        possource = S.possource;
    end
    %sources
    for ii = 1:length(S.sources)
        if cmp(ii,1) ~= 0
            if S.size_l == 1
                %size
                plot3(possource(ii, 1), possource(ii, 2), possource(ii, 3),'Color','g','Marker','.', 'MarkerSize', sources_scale(ii==idxf,1));
            else
                %color
                plot3(possource(ii, 1), possource(ii, 2), possource(ii, 3),'Color',[0 sources_scale(ii==idxf,1) 0],'Marker','.', 'MarkerSize', 20);
            end
        end
        hold on
    end
%     plot3(possource(S.source, 1), possource(S.source, 2), possource(S.source, 3), '.g', 'MarkerSize', 20);
    %positive MEG channels (strength recorded by MEG channels of neural signal from source of leadfield)
    for ii = 1:length(LFP)
        hold on
        plot3(poschan(LFPi(ii), 1), poschan(LFPi(ii), 2), poschan(LFPi(ii), 3), '.r', 'MarkerSize', LFP(ii));
    end
    %negative MEG channels (strength recorded by MEG channels of neural signal from source of leadfield)
    for ii = 1:length(LFN)
        hold on
        plot3(poschan(LFNi(ii), 1), poschan(LFNi(ii), 2), poschan(LFNi(ii), 3), '.b', 'MarkerSize', LFN(ii));
    end
    rotate3d on;
    axis off
    axis vis3d
    axis equal
    set(gcf,'Color','w')
end


end


function [  ] = LF_3D_plot_LBPD( S )

% Plotting in 3D space leadfield over the MEG sensors given a brain source
% and orientation.
% This function should work for each individual subject. If you do not have
% any subject information, just for visualization and demonstration purposes,
% a previously computed subject data will be shown.
% For now, the function works in 8mm MNI152-T1 space, but this could be
% easily changed, if needed.
% The function uses FieldTrip and in-house-built codes.
% This function should be used in connection with LBPD.
% Please, contact me for more information (Leonardo Bonetti, leonardo.bonetti@clin.au.dk).



% INPUT:    -S.sensor_l:                    1 for magnetometers; 2 for gradiometers; 3 for both
%           -S.source:                      source from the leadfield to be plotted
%                                           Must be both magnetomeres and gradiometers, sorted as:
%                                           Mag0111; Grad0112; Grad0113; Mag0121; Grad0122; Grad0123, etc.
%           -S.orientation:                 indicates 1 of the 3 orientations of the source (dipole) from the leadfield to be plotted
%           -S.topoplot_label:              1 for plotting also topoplot; 0 for not plotting topoplot
%           -S.L:                           leadfield model as outputted by FieldTrip function ft_prepare_leadfield
%           -S.scaleextr:                   xtreme values for scaling MEG sensors in plot (e.g.[0.0001 30])
%           -S.vol:                         description of the head volume conduction model as outputted by FieldTrip function ft_prepare_headmodel.
%                                           Leave empty [] to load a previously computed solution (sort of template).
%           -S.poschan:                     MEG channels position (their space must be consistent with the space of vol).
%                                           Leave empty [] to load a previously computed solution (sort of template).
%           -S.label:                       MEG channels labels (must refer to all 306 MEG channels, neuromag-elekta system)
%                                           Leave empty [] to load a previously computed solution which assumes standard label system for neuromag-elekta.
%           -S.possource:                   brain sources position (their space must be consistent with the space of vol and poschan).
%                                           Leave empty [] to load a previously computed solution (sort of template).


% OUTPUT:   -..:                            3D plotting solution showing source and leadfield over the MEG channels, with given orientation






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 26/02/2021


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%leadfield model
if isempty(S.L)
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/plot3D_LF.mat','L')
else
    %extracting channels position
    L = S.L;
end
%position of MEG channels
if isempty(S.poschan)
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/plot3D_LF.mat','poschan')
else
    %extracting channels position
    poschan = S.poschan;
end
%MEG sensor triplets are defined with overlapping positions, so here I rebuild an approximation of such triplets
poschan(2:3:306,:) = poschan(1:3:306,:) + 0.003; %2nd channel of the triplet
poschan(3:3:306,:) = poschan(1:3:306,:) - 0.003; %3th channel of the triple (the 1st channel of the triplet maintain the original position)
%preparing data from lead field model
lf = L{S.source}(:,S.orientation); %extracting source (and orientation) from leadfield
%preparing MEG sensors
if S.sensor_l == 1 %only magnetometers
    lf(2:3:306) = 0; %zeroing gradiometers
    lf(3:3:306) = 0;
elseif S.sensor_l == 2 %only gradiometers
    lf(1:3:306) = 0; %zeroing magnetometers
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
%inspecting leadfield model using topoplot (if requested)
%topoplot must be done first for reasons related to Matlab handling figures..
if S.topoplot_label == 1
    avg = lf;
    if ~exist('M_timb_lock_gradplanarComb','var') %if not already loaded..
        fieldtrip_example = load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/fieldmask.mat');
    end
    if ~isempty(S.label)
        label = S.label;
    else
        load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/plot3D_LF.mat','label')
    end
    cfgdummy = fieldtrip_example.M_timb_lock_gradplanarComb.cfg;
    cfgdummy.previous = [];
    data = [];
    data.cfg = cfgdummy;
    data.time = 1;
    data.label = label;
    data.dimord = fieldtrip_example.M_timb_lock_gradplanarComb.dimord;
    data.grad = fieldtrip_example.M_timb_lock_gradplanarComb.grad;
    data.avg = avg;
    %creating the cfg for the actual plotting (magnetometers)
    cfg = [];
    cfg.layout = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/neuromag306all.lay';
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
figure
if isempty(S.vol)
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/plot3D_LF.mat','vol')
else
    vol = S.vol;
end
%brain with triangles
ft_plot_vol(vol, 'edgecolor', [0 0 0], 'facealpha', 0);
% openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
hold on
%source from leadfield
if isempty(S.vol)
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/plot3D_LF.mat','possource')
else
    possource = S.possource;
end

%to switch back to original sources position
% load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/source_LBPD/MNI_back.mat')
% ab = find(S.source==MNI_sorted_index);
% plot3(possource(ab, 1), possource(ab, 2), possource(ab, 3), '.g', 'MarkerSize', 20);
%until here

%previous line
plot3(possource(S.source, 1), possource(S.source, 2), possource(S.source, 3), '.g', 'MarkerSize', 20);
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


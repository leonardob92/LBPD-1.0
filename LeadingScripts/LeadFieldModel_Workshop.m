%%

%%% LEADFIELD MODEL - WORKSHOP %%%


%instructions
%1)run the first section to start some functions up
%2)run the second section playing around with the three provided settings
%   -source_label: nunmbers corresponding to different brain sources (voxels)
%   -orientation: 1 of the 3 orientation of the brain sources (modelled as dipoles)
%   -topoplot_label: plotting also the topoplot, if you like (I think it is nice to show that as well to better understand the difference between the 2D and 3D representation of the MEG channels)



%% LBPD_startup_D

%starting up some functions
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
addpath(pathl);
LBPD_startup_D(pathl);

%% PLOTTING LEADFIELD MODEL IN 3D
%to be used after computing the beamforming (obviously it requires computation of leadfield model)

%%% user settings to be used in the demo (for DAVID)
source_label = 209; %source from the leadfield to be plotted
orientation = 2; %orientation of the source (dipole) from the leadfield to be plotted
%%%


%NOT TO BE USED
%additional user settings (not to be changed for the demo, I guess)
topoplot_label = 1; %1 for plotting also topoplot; 0 for not plotting topoplot
subject_n = 8; %subject number (referring to LearningBach2017_minor); numbers correspond to original subjects ID; 0 for previously computed subject automatically outputted by the function
loadsubj = 1; %if you inspect different sources of the same subject you do not need to load the subject every time, so switch it to 0
%datadir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/source_LBPD/LearningBach/Beam_abs_1_sens_1_freq_broadband_time_1_466'; %directory where the leadfields are stored
% more robust data dir (scratch folders may change over time):
% datadir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/source_LBPD/LearningBach/Beam_abs_1_sens_1_freq_broadband_time_1_466'; %directory where the leadfields are stored
% datadir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_1_sens_1_freq_broadband';
datadir = '/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/source_LBPD/LearningBach/Beam_abs_1_sens_1_freq_broadband_time_1_466';
original_source = 1; % 1 for original source in leadfield; 0 for sources with ID number as depicted in the following image: "/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz" (this is useful to index more easily the sources with visual inspection of the image for example in FSLEYES)
%user labels (most important settings)
sensor_l = 2; %1 for magnetometers; 2 for gradiometers; 3 for both
scaleextr = [0.0001 30]; %extreme values for scaling MEG sensors in plot
if subject_n == 0
    L = []; %leadfield model; leave empty [] to load a previously computed solution (sort of template)
    vol = []; %description of the head volume conduction model as outputted by FieldTrip function ft_prepare_headmodel; leave empty [] to load a previously computed solution (sort of template).
    poschan = []; %description of the head volume conduction model as outputted by FieldTrip function ft_prepare_headmodel; leave empty [] to load a previously computed solution (sort of template).
    possource = []; %brain sources position (their space must be consistent with the space of vol and poschan); leave empty [] to load a previously computed solution (sort of template).
    label = []; %MEG channels labels (must refer to all 306 MEG channels, neuromag-elekta system)
else
    if ~exist('OUT','var')
        load([datadir '/SUBJ_' num2str(subject_n) '_norm0_abs_1.mat'])
    end
    L = OUT.Ltot; %leadfield model; leave empty [] to load a previously computed solution (sort of template)
    vol = OUT.vol; %description of the head volume conduction model as outputted by FieldTrip function ft_prepare_headmodel; leave empty [] to load a previously computed solution (sort of template).
    poschan = OUT.sens.chanpos; %description of the head volume conduction model as outputted by FieldTrip function ft_prepare_headmodel; leave empty [] to load a previously computed solution (sort of template).
    possource = OUT.pos_brainsources_MNI8; %brain sources position (their space must be consistent with the space of vol and poschan); leave empty [] to load a previously computed solution (sort of template).
    label = OUT.sens.label; %MEG channels labels (must refer to all 306 MEG channels, neuromag-elekta system)
end
%loading original MNI coordinates
if ~exist('mni_coords','var') %if not already loaded..
    [ mni_coords, ~ ] = osl_mnimask2mnicoords('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz');
end
%building input structure
S = [];
%simply passing the inputs previously chosen by the user
S.L = L; S.vol = vol; S.poschan = poschan; S.possource = possource; S.label = label; S.sensor_l = sensor_l; S.source_label = source_label; S.orientation = orientation; S.topoplot_label = topoplot_label; S.scaleextr = scaleextr;
S.source = source_label;
%actual function
LF_3D_plot_LBPD(S)


%% 















































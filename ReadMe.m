%%% LBPD %%%

% CITATION
% If you find LBPD functions useful, please cite the following paper(s):
%
% 1) Spatiotemporal whole-brain dynamics of auditory patterns recognition. BiorXiv.
%    Bonetti, L., Brattico, E., Carlomagno, F., Cabral, J., Stevner, A., Deco, G., Whybrow, P.C., Pearce, M., Pantazis, D., Vuust, P., Kringelbach, M.L. (2020)
%    Pre-print: https://www.biorxiv.org/content/10.1101/2020.06.23.165191v3.full
% 2) Bonetti, L., Brattico, E., Carlomagno, F., Donati, G., Cabral, J., Haumann, N. T., Deco, G., Vuust, P., Kringelbach, M. L. (2021).
%    Rapid encoding of musical tones discovered in whole-brain connectivity. NeuroImage, 118735.
%    Paper: https://www.sciencedirect.com/science/article/pii/S1053811921010077
%    Leading script(s): https://github.com/leonardob92/RapidEncodingSoundsBrainConnectivityMEG_NeuroImage
% 3) Brain recognition of previously learned versus novel temporal sequences: a differential simultaneous processing. Cerebral Cortex.
%    Bonetti, L., Brattico, E., Bruzzone, S.E.P., Donati, G., Deco, G., Pantazis, D., Vuust, P., Kringelbach, M.L. (2022)
%    Paper: https://academic.oup.com/cercor/advance-article/doi/10.1093/cercor/bhac439/6809237
%    Leading script(s): https://github.com/leonardob92/Brain_Recognition_Temporal_Sequences_Differential_Simultaneous_Processing_Cerebral_Cortex
% 4) Bruzzone, S.E.P., Lumaca, M., Brattico, E., Vuust, P., Kringelbach, M.L., Bonetti, L. (2022)
%    Dissociated brain functional connectivity of fast versus slow frequencies underlying individual differences in fluid intelligence: a DTI and MEG study. Scientific Reports, 12, 4746.
%    Paper: https://www.nature.com/articles/s41598-022-08521-5
%    Leading script(s): https://github.com/leonardob92/FluidIntelligenceBrainConnectivity_MEG_DTI_ScientificReports
% 5) Whole-brain computation of cognitive versus acoustic errors in music: A mismatch negativity study.
%    Bonetti, L., Carlomagno, F., Kliuchko, M., Gold, B.P., Palva, S., Haumann, N.T., Tervaniemi, M., Huotilainen, M., Vuust, P., Brattico, E. (2022)
%    NeuroImage: Reports
%    Paper: https://www.sciencedirect.com/science/article/pii/S2666956022000691
%    Leading script(s): https://github.com/leonardob92/Whole_brain_computation_of_cognitive_vs_acoustic_errors_in_music_NeuroImageReports
% 6) Associations between abstract working memory abilities and brain activity underlying long-term recognition of auditory sequences.
%    PNAS Nexus, Volume 1, Issue 4, pgac216 (2022).
%    Fernández Rubio, G., Carlomagno, F., Vuust, P., Kringelbach, M.L., Bonetti, L.
%    Paper: https://doi.org/10.1093/pnasnexus/pgac216
%    Leading script(s): https://github.com/gemmaferu/MEG-auditory-memory/blob/main/papers/Associations_WM_long-term_recognition_auditory_sequences.m
% 7) Magnetoencephalography recordings reveal the spatiotemporal dynamics of recognition memory for complex versus simple auditory sequences
%    Fernández Rubio, G., Brattico, E., Kotz, S.A., Kringelbach, M.L., Vuust, P., Bonetti, L. (2022)
%    Communications Biology
%    Paper: https://www.nature.com/articles/s42003-022-04217-8#code-availability
%    Leading script(s): https://github.com/gemmaferu/Communications-Biology-MEG-recognition-auditory-sequences/tree/v1.0
% 8) Revealing the spacetime hierarchical whole-brain dynamics of auditory predictive coding. BiorXiv.
%    Bonetti, L., Fernandez Rubio, G., Carlomagno, F., Pantazis, D., Vuust, P., Kringelbach, M.L. (2022)
%    Pre-print: https://www.biorxiv.org/content/10.1101/2022.11.19.517195v2

% The folder "LeadingScripts" contains the pipeline of the above referenced studies which have not been published on peer-reviewed journals yet, but are available on BiorXiv as pre-prints.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INTRODUCTION TO LBPD
% LBPD is a collection of functions for MEG brain activity and connectivity
% data analysis.
% Summarizing some of my work in LBPD helped me to develop, order and better understand
% my analyses. Further, LBPD aims to allow other scientists to control and replicate
% my work and provide a few tools for analysis and plotting, in case anybody in the community is interested.
% LBPD strongly supports the idea that science should be based on codes
% that are highly commented and (hopefully) easy to understand and use.
% LBPD relies on a set of own functions and on external functions,
% sometimes taken from or corresponding to well-known toolboxes such as OSL,
% FieldTrip and SPM12.
% I do not want to take any credit for those functions and I deeply thank
% the brilliant scientists that developed them and made them publicly
% available.
% If you notice any mistake, typo, etc. in the codes, please feel free to
% contact me and I will fix/clarify them. The codes should work, but they
% are constantly been updated and improved.
% If you notice that I did not properly acknowledge the contribution of
% some external functions and their authors, I apologise. In that case,
% please let me know and I will made the requested changes.
% Please, use the Leading Scripts provided in the "LeadingScripts" folder
% to better understand how the LBPD (and other) functions (e.g.functions from SPM/FiledTrip/OSL,etc.)
% have been used to perform the analyses reported in the manuscripts mentioned above.
% Please, note that LBPD is free software and therefore you can redistribute
% it and/or modify it under the terms of the GNU General Public License as
% published by the Free Software Foundation, either version 3 of the License,
% or (at your option) any later version.

% leonardo.bonetti@clin.au.dk, Leonardo Bonetti, Aarhus, DK, 27/05/2021
% leonardo.bonetti@psych.ox.ac.uk, Leonardo Bonetti, Oxford, UK, 29/05/2022


% EXTERNAL FUNCTIONS
% Among the reported external functions, a key role is played by SPM12, OSL
% and FieldTrip. Even if LBPD uses only few of their functions, most of the LBPD codes are
% mainly meaningful in connection to those toolboxes. Especially, LBPD starts
% from reading data from MEEG SPM objects that are obtained by using SPM
% and/or OSL. Therefore, if you want to follow a proposed 'pipeline' of
% LBPD functions, you should start by preprocessing the data using OSL/SPM.
% (However, please note that some of the LBPD functions may be used even without
% relying on SPM objects. Please read the specific documentation provided
% for each function).
% To use the functions that require the OSL toolbox, please donwload OSL at the
% following link: https://ohba-analysis.github.io/osl-docs/
% Once you download OSL, please copy it in the directory: 'yourpath/LBPD/External/',
% then the 'LBPD_startup_D' function will set everything up.
% Remember that you need to have iOS or Linux to use OSL.
% If you have Microsoft Windows you can only run the codes that do not explicitely
% require OSL functions.
% Further, in that case, you should donwload FieldTrip and SPM.
% You find FielTrip (the current LBPD codes rely on the version: fieldtrip-20171231) here:
% http://www.fieldtriptoolbox.org/
% You find SPM here:
% https://www.fil.ion.ucl.ac.uk/spm/
% Then, please copy FieldTrip and SPM in the directory: 'yourpath/LBPD/External/',
% To start LBPD up, you only need to run the 'LBPD_startup_D' function
% (and use only the functions that do not explicitely require OSL).
% A special acknowledgment must be made for the codes developed by
% Dimitrios Pantazis that, even if not directly connected to LBPD
% functions, have been employed in the attached Leading Script(s) and can be
% used to understand the results of the following study:
% https://www.biorxiv.org/content/10.1101/2020.06.23.165191v3.full



% LBPD FUNCTIONS RECAP
% Please, note that the Leading Scripts uses different subsets of the functions reported below.
% All these functions have been made publicly available already now, and they
% will be further documented, explained and complemented with examples in
% future works that will be published soon.
% 
% 1)LBPD_startup_D

% IFC - phase synchrony - SFC (BRAIN FUNCTIONAL CONNECTIVITY ANALYSIS)
% 2)ps_preparedata_spmobj_LBPD_D
% 3)extracting_data_LBPD_D
% 4)preparing_baseline_from_restingstate_LBPD_D
% 5)phasesynchrony_LBPD_D
% 6)ps_statistics_LBPD_D
% 7)degree_segregation_MCS_LBPD_D
% 8)diff_conditions_matr_degree_MCS_LBPD_D
% 9)diff_conditions_matr_couplingROIs_MCS_LBPD_D
% 10)IFC_plotting_LBPD
% 11)signROIs_degree_connotherROIs_LBPD_D
% 12)signROIdegree_otherROI_prepareplotting_LBPD_D
% 13)generalCoupling_prepareplotting_LBPD_D
% 14)FC_Brain_Schemaball_plotting_LBPD_D
% 15)static_FC_MEG_LBPD_D

% MEG - sensors (BRAIN ACTIVITY ANALYSIS)
% 16)islands3D_LBPD_D
% 17)MEG_sensors_plotting_ttest_LBPD_D2
% 18)MEG_sensors_MCS_reshapingdata_LBPD_D
% 19)MEG_sensors_MonteCarlosim_LBPD_D
% 20)MEG_sensors_combining_magclustsign_LBPD_D
% 21)MEG_sensors_MCS_plottingclusters_LBPD_D

% Additional Monte Carlo simulation functions
% 22)oneD_MCS_LBPD_D
% 23)twoD_MCS_LBPD_D

% Additional functions (for source reconstruction using Beamforming, plotting networks in the brain, MCS on graph theoretical measures, spatial-functional k-means clustering to estimate a functional parcellation of the brain, etc.)
% 24)BrainSources_MonteCarlosim_3D_LBPD_D
% 25)GT_modul_plot_LBPD
% 26)LF_3D_plot_LBPD
% 27)MEG_SR_Beam_LBPD
% 28)MEG_SR_Stats1_Fast_LBPD
% 29)plot_sensors_wavebis2
% 30)plot_wave_conditions
% 31)schemaball_modularity_LBPD
% 32)sources_3D_plot_LBPD
% 33)cluster_DTI_perm_2groups_1
% 34)cluster_DTI_perm_2groups_2
% 35)cluster_DTI_perm_2groups_3
% 36)DTI_cluster_perm_2groups_LBPD
% 37)DTI_GT_MCS
% 38)FunctionalSpatialClustering_voxels2ROIs_LBPD_D
% 39)STC_plottingtimeseries_LBPD
% 40)MEG_SR_Stats_twogroups_LBPD
% 41)Extract_BrainCluster_Information_3D_LBPD_D
% 42)InducedResponses_Morlet_ROIs_LBPD_D
% 43)InducedResponses_Morlet_Coords_AALROIs_LBPD_D
% 44)FromCoordMatrix_2_3DNifti_8mm_LBPD_D
% 45)From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D
% 46)InducedResponses_Morlet_WholeBrain_LBPD_D
% 47)Extract_MEGSensor_Information_LBPD_D
% 48)waveplot_groups_server_LBPD (waveplot_groups_local for local computers and not Aarhus server)

% Additional scripts:
% - Workbench_Codes_LocalMac_Example
% - MEGSourceReconstruction_LeadFieldModel_Workshop

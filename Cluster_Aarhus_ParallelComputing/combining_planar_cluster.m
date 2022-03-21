
function O = combining_planar_cluster(input)

O = [];

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup


D = spm_eeg_combineplanar(input);



% perform analysis for subject n using data directory and subject nr


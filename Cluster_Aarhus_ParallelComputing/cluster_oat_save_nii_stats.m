function O = cluster_oat_save_nii_stats(S2)
O = []; 

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup
% addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')


[statsdir,times,count] = oat_save_nii_stats(S2);

end


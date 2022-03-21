function O = cluster_beamgrouplevel(oat)
O = []; 

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/ohba-external/fmt')



oat = osl_check_oat(oat);
oat.to_do = [0 0 0 1];
oat = osl_run_oat(oat);
% oat = osl_run_oat_F(oat);

end


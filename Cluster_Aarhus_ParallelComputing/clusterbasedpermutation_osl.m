function O = clusterbasedpermutation_osl(S)
O = []; 

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup


[ gstats ] = oat_cluster_permutation_testing(S);


end


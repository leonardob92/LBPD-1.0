function O = cluster_africa(S)
O = []; 

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup


D = S.D;
D = osl_africa(D,'do_ica',true,'do_ident',false,'do_remove',false,'used_maxfilter',true); 
D.save();

end


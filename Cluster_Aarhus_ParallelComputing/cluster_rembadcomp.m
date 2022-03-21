function O = cluster_rembadcomp(S)
O = []; 

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup


D = S.D;
% S = rmfield(S,'D');
D = osl_africa(D,'do_ident',false,'do_remove',true);
D.save();


% D = remove_bad_components_l(D,S);
% D.save();

end


%%

function O = coregfuncp(S) % apparently function need a 0 and no end to run in the cluster (!?), set the name of the function and the input variable 

O = [];

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')

%actual function
D1 = spm_eeg_load(S.D1);
S = rmfield(S,'D1');
ii = S.ii;
S = rmfield(S,'ii');
D = osl_headmodel(S);
% rhino_display(D);
D.save();

inv_rhino = D.inv;
D1.inv = inv_rhino;
D1.save();
clear D1
disp(num2str(ii))




% for using the server cluster we need to generate n output files with
% different name that's why we need to use num2str((Range)

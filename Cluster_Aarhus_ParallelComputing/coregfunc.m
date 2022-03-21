
function O = coregfunc(input) % apparently function need a 0 and no end to run in the cluster (!?), set the name of the function and the input variable 

O = [];


addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup



S = input;
%actual function
D = osl_headmodel(S);
D.save();
rhino_display(D);






% for using the server cluster we need to generate n output files with
% different name that's why we need to use num2str((Range)

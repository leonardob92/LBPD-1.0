function O = cluster_epoch_osl(S2)
O = []; 

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup


D_continuous = S2.D_continuous;
S2 = rmfield(S2,'D_continuous');

[epochinfo.trl, epochinfo.conditionlabels, MT] = spm_eeg_definetrial(S2);


%do epoching
S3 = [];
S3.prefix = S2.prefix;
S3 = epochinfo;
S3.D = D_continuous;
D = osl_epoch(S3);

end


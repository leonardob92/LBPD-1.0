function O = red_dim(S)
O = []; 

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup


p = S.p;
D = spm_eeg_load(S.filepath);                
D = D.montage('switch',S.montage);
D = ROInets.get_node_tcs(D,p.parcelflag(true),'PCA');
D.save();
                

end


function O = decoding(S)
O = []; 

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab/external');
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab/external/libsvm-3.21')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab/external/libsvm-3.21/matlab')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/SpatialLocation')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/SpatialLocation/plotchannel')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/permutationlab')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/permutationlab/private1')
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/permutationlab/developer')


% make



data = S.data;
condid = S.condid;
numperm = S.numperm;
kfold = S.kfold;
savedir = S.savedir;
ii = S.ii;
pairl = S.pairl;

%pairwise decoding
if pairl == 1
    d = sll_decodesvm_pairwise_haufe(data,condid,'numpermutation',numperm,'verbose',2,'kfold',kfold); %method = pairwise
else
    d = sll_decodesvm_pairwise_haufe(data,condid,'numpermutation',numperm,'method','multiclass','verbose',2,'kfold',kfold); %to get confusion matrix
end
save([savedir '/PD_SUBJ' num2str(ii) '.mat'],'d');

%temporal generalization
% d = sll_decodesvm(data,condid,'numpermutation',numperm,'method','temporalgen','verbose',1,'kfold',kfold); %max correlation classifier
% save([savedir '/TG_SUBJ' num2str(ii) '.mat'],'d');


end




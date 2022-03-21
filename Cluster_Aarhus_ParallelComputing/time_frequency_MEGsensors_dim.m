function O = time_frequency_MEGsensors_dim(S)
O = []; 

%starting up some of the functions that I wrote for the following analysis
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/BrainstormDemo_CodesTimeFrequencyDecomposition/functions') %to Dimitrios Morlet transform function

%preparing data
SS = size(S.data);
%actual computation
for cc = 1:length(S.conds) %over experimental conditions
    P = zeros(SS(1),length(S.f),SS(2),length(S.conds{cc}));
    for ii = 1:length(S.conds{cc}) %over trials of condition cc
        P(:,:,:,ii) = morlet_transform(S.data(:,:,S.conds{cc}(ii)),S.timesel,S.f); %frequency decomposition
        disp(ii)
    end
    S2 = S;
    S2 = rmfield(S2,'data');
    %averaging over trials already now to save space on disk..
    P = mean(P,4);
    %saving data on disk
    save([S.outdir '/Subject_' num2str(S.ii) '_' S.conds_label{cc} '.mat'],'P','S2','-v7.3')    
end

end
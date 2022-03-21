function O = cluster_Dlabel(S2)
O = []; 

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/osl-core'); %this add the path (not necessary being in the osl-core directory if you have this one)
osl_startup


D = S2.D; %file to be updated
D2 = S2.D2; %file with original labels

D=D.montage('switch',0);
for ii = 1:length(D.chanlabels)
    D = D.chanlabels(ii,[D2.chanlabels{ii}]);
    disp(ii)
end

D=D.montage('switch',1);
for ii = 1:length(D.chanlabels)
    D = D.chanlabels(ii,[D2.chanlabels{ii}]);
    disp(ii)
end

% D.sensors('MEG').label

D.save();



end


function O = RSA_APR2020_All( S )
O = []; 

% Aggregating RSAs for all comparisons of tones (Old and New)

% INPUT:    S.l:    1 = all together (but not same tones)
%                   2 = same tones
%                   3 = properly all together



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@psych.ox.ac.uk
% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 30/07/2023


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%building the list
list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/RSA/AveragedTrials/*.mat'); %dir to epoched files (encoding)
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat') %loading time
dumlist = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/RSA/AveragedTrials/OtherComparisons/*.mat');
carb = zeros(length(dumlist),1);
for ii = 1:length(dumlist) %over files
    carb(ii) = (double(dumlist(ii).name(25)) - double(dumlist(ii).name(20))); %storing distance between tones for each file
end
if S.l == 1 %if all together (but not same tones)
    list2 = dumlist(carb~=0); %getting new comparisons
    list3 = list(2:2:end); %getting previous comparisons (Old vs NewT1)
    list = cat(1,list3,list2); %combined list
elseif S.l == 2 %same tones
    list2 = dumlist(carb==0);
    list3 = list(3:2:end); %getting previous comparisons (Old vs NewT1)
    list = cat(1,list3,list2); %combined list
elseif S.l == 3 %properly all together
    list3 = list(2:end); %getting previous comparisons (Old vs NewT1)
    list2 = dumlist;
    list = cat(1,list3,list2); %combined list
end
%actual computation
RSA = zeros(89,89,2,83,length(list)); %time x time x conditions x subjects x pairs
for ii = 1:length(list) %over subsequent pairs of tones
    load([list(ii).folder '/' list(ii).name]); %loading data
    RSA(:,:,:,:,ii) = SS; %storing data
end
RSAm = mean(RSA,5); %mean over pairs

save(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/RSA/AveragedTrials/AggregatedResults/RSAm_case_' num2str(S.l) '.mat'],'RSAm','time');


end


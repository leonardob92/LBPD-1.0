function O = averaging_TRF(S)
O = []; 

S = [];
%OLD
PP = zeros(306,74,700,70);
for ii = 1:71
    if ii ~= 39
        load(['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/Time_frequency_MEGsensors/Subject_' num2str(ii) '_Old_Correct.mat']);
        PP(:,:,:,ii) = P;
    end
    disp(ii)
end
P = mean(PP,4);
save(['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/Time_frequency_MEGsensors/Old_Correct_SubjectsAverage.mat'],'P','-v7.3')    

%NEW
PP = zeros(306,74,700,70);
for ii = 1:71
    if ii ~= 39
        load(['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/Time_frequency_MEGsensors/Subject_' num2str(ii) '_New_Correct.mat']);
        PP(:,:,:,ii) = P;
    end
    disp(ii)
end
P = mean(PP,4);
save(['/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/Time_frequency_MEGsensors/New_Correct_SubjectsAverage.mat'],'P','-v7.3')    


end
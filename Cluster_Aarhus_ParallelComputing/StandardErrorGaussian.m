function O = StandardErrorGaussian(S)
O = []; 

%0.1-1Hz
indir = dir('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/LearningBach/source/source_localsphere_01_1Hz_AllTrials_BC_2.oat/subj*dir');
DAT = zeros(23,27,23,400,length(indir)); %3D image, time-points and then subjects
for ii = 1:length(indir) %over subjects
    T = load_nii([indir(ii).folder '/' indir(ii).name '/tstat3_8mm.nii.gz']); %loading nifti image
    DAT(:,:,:,:,ii) = T.img(:,:,:,1:400); %storing matrix with statistics from nifti image
    disp(ii)
end
STDE = std(DAT,0,5)./sqrt(length(indir)); %computing standard errors
save('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/Papers/GaussianBach/SlowFrew_STDE.mat','STDE')


end


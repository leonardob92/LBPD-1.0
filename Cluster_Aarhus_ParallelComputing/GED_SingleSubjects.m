function O = GED_SingleSubjects(S)
O = []; 



addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia')
comps2see = S.comps2see;
list = S.list;
shr   = 0.01;     % value between 0 and 1 (0 small perturbation, 1 big perturbation) - for Mattia regularisation (until 0.001 it works for restoring the full rank of the matrix)
%loading data
load('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1/Test_GED/bomba.mat');
for ii = 1:size(covR) %over subjects
    covS_subj = covS(:,:,ii);
    covR_subj = covR(:,:,ii);
    % regularisation (to bring covariance matrices to full rank)
    %narrow data
    covS_subj = covS_subj  + 1e-6*eye(size(covS_subj )); %making matrix full rank again by adding a small perturbation/noise (matrix regularisation)
    evalsR = eig(covR_subj );
    covR_subj  = (1-shr)*covR_subj  + shr * mean(evalsR) * eye(size(covR_subj ));
    
    % Eigendecomposition
    [evecs,evals] = eig(covS_subj ,covR_subj );
    [evals,sidx]  = sort( diag(evals),'descend' );
    evecs = evecs(:,sidx);
    evals2 = (evals.*100)./sum(evals);

    % Import broadband source data
    disp(['loading data subject ' list(ii).name])
    load([list(ii).folder '/' list(ii).name])
    broadData = OUT.sources_ERFs(:,:);
    
    % compute filter forward model and flip sign
    ncomps = 100;
    GEDmap = zeros(ncomps,size(broadData,1));
    GEDts  = zeros(ncomps,size(broadData,2));

    % computing time series of components
    disp('computing time series')
    for compi = 1:ncomps
        GEDmap(compi,:) = evecs(:,compi)'*covS_subj ; % get component
        [~,idxmax] = max(abs(GEDmap(:,1)));     % find max magnitude
        GEDmap(compi,:)  = GEDmap(compi,:)*sign(GEDmap(compi,idxmax)); % possible sign flip
        % compute time series (projections)
        GEDts(compi,:) = evecs(:,compi)'*broadData;
        disp(compi)
    end
    subjID = list(ii).name;
    disp(['saving subject ' list(ii).name])
    save(['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1/Test_GED/' list(ii).name(1:9) '.mat'],'covS_subj','covS_subj','evecs','evals','GEDmap','GEDts','srate','subjID','shr','-v7.3')    

    %plotting weights in the brain (nifti images)
    maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
    for iii = 1:length(comps2see) %over significant PCs
        %using nifti toolbox
        SO = GEDmap(iii,:)';
        %building nifti image
        SS = size(maskk.img);
        dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
        for jj = 1:size(SO,1) %over brain sources
            dumm = find(maskk.img == jj); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
            [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
            dumimg(i1,i2,i3,:) = SO(jj,:); %storing values for all time-points in the image matrix
        end
        nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
        nii.img = dumimg; %storing matrix within image structure
        nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
        disp(['saving nifti images - subject ' list(ii).name(1:9)])
        save_nii(nii,['/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1/Test_GED/' list(ii).name(1:9) '_Comp_' num2str(iii) '_Var_' num2str(round(evals2(iii))) '.nii']); %printing image
    end

end


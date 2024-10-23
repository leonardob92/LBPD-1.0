function [ OUT ] = PCA_LBPD( S )

% Principal Component Analysis (PCA) calculation (eigenvector solution on
% data covariance matrix) with Monte-Carlo simulation (MCS) implementation
% designed considering the variance of PCs (eigenvalues) in real and
% randomized data.
% Please, note that here only the weights of the PCA are provided and not
% the activation patterns (W*COV(data)) because no relevant differences
% emerged between them. It is instead necessary for other multivariate
% methods (e.g. SVM or generalised PCA).


% INPUT:
%         -S structure with fields:
%           -H:             data matrix (i.e. brain sources x time-points) (optionally it can be brain sources x time-points x experimental conditions).
%                           If a third dimension is provided  (e.g. for experimental conditions)
%                           PCA will be done on the data after average over the third dimension.
%           -rand_l:        strategy for randomisation:
%                           1 = randomising only time-points (independently for each time series)
%                           2 = randomising only space (independently for each time-point) 
%                           3 = randomising both time and space
%           -permnum:       number of permutations for Monte-Carlo simulations (e.g. 100 or 1000).
%                           If negative values, you take a fixed number of components (e.g. -3 = getting the first 3 components and not computing the MCS).  
%           -fig_l:         1 for plotting some figures illustrating the results; 0 otherwise
%           -sign_eig:      method for normalizing eigenvectors sign (SUGGESTED EITHER 'max_abs' or 'average'):
%                            -'occurrences' = using mean of the negative/positive values occurrences
%                            -'max_abs' = on the basis of the sign of the maximum value in absolute terms
%                            -'average' = on the basis of the sign of the average of weights
%           -namenii:       path plus name for nifti images (e.g. '/scratch5/MINDLAB2017_MEG-LearningBach/Leonardo/Analysis1'
%                               then: 'PC_ii.nii' will be automatically added for the ii significant PCs)
%           -time:          time in seconds (for plotting)
%           -onefig:        1 = time series of all significant components in the same figure (with all conditions as well).
%                           2 = one figure for each significant component

% OUTPUT:   -OUT:           OUT structure with:
%                            -variance explained by the components
%                            -timeseries of the significant components
%                            -weights (brain sources)
%                            -index/number of significant components






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 20/02/2021
% Leonardo Bonetti, Aarhus, DK, 10/03/2024


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%% EXAMPLE OF SETTINGS %%%
% load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
% 
% S = [];
% S.H = dum(:,:,1:2);
% S.permnum = 3;
% S.fig_l = 1;
% S.sign_eig = '0';
% S.namenii = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Test_PCA/test';
% S.time = time;
% S.rand_l = 2;
% S.onefig = 0;
% 
% [ OUT ] = PCA_LBPD( S )
%%%   %%%   %%%


%%% PCA (STEPS) - MANUAL SOLUTION
% H2 = bsxfun(@minus,H,mean(H)); %demeaning
% % [coeffp,scorep,latentp] = pca(H2); %PCA
% %PCA in single steps
% Hcov = cov(H2'); %covariance matrix (note that H2 is transposed, so time-points x sources)
% [V,D] = eig(Hcov); %eigenvector problem
% %getting eigenvalues from diagonal matrix (correspond to variance explained by the different princiapl components of PCA)
% Dv = D(D~=0);
% Dv = Dv(end:-1:1); %sorting eigenvalues (since they are original from the smallest to the highest explained variance)


%extracting inputs
time = S.time;
permnum = S.permnum;
rand_l = S.rand_l;

%actual PCA computation
H = mean(S.H(:,:,:),3);
[wcoeff,~,~,~,vare] = pca(H');

% PCA on randomised data
if permnum > 0
    mfD = zeros(permnum,1);
    for pp = 1:permnum %over permutations
        if rand_l == 1     %randomizing only time
            r_resh = zeros(size(H,1),size(H,2));
            for ss = 1:size(H,1) %over brain sources
                idx_dummy = randperm(size(H,2)); %create a permuted array of the indices of the temporal dimension
                r_resh(ss,:) = H(ss,idx_dummy);
            end
        elseif rand_l == 2 %randomizing only space
            r_resh = zeros(size(H,1),size(H,2));
            for ss = 1:size(H,2) %over time-points
                idx_dummy = randperm(size(H,1)); %create a permuted array of the indices of the spatial dimension
                r_resh(:,ss) = H(idx_dummy,ss);
            end
        elseif rand_l == 3 %randomizing both space and time
            %randomizing data..
            idx_dummy = randperm(size(H,1)*size(H,2)); %create a permuted array from 1 to size of original data vector
            r_dummy = zeros(size(H,1)*size(H,2),1); %initialise new vector
            r_dummy(1:size(H,1)*size(H,2)) = H(idx_dummy); %taking element in matrix data with index idx_dummy
            r_resh = reshape(r_dummy,[size(H,1),size(H,2)]); %reshaping the vector into a matrix shaped as the original data
            %running PCA in single steps on randomised data
        end
        [~,~,~,~,varef] = pca(r_resh');
        mfD(pp,1) = max(max(varef)); %storing maximum variance (max eigenvalue) occurring for PCs in randomized data
        disp(pp)
    end
    MFD = max(mfD); %absolute maximum variance
end

%ALTERNATIVE - MANUAL SOLUTION
% mfD = zeros(permnum,1);
% for pp = 1:permnum %over permutations
%     %randomizing data..
%     idx_dummy = randperm(size(H,1)*size(H,2)); %create a permuted array from 1 to size of original data vector
%     r_dummy = zeros(size(H,1)*size(H,2),1); %initialise new vector
%     r_dummy(1:size(H,1)*size(H,2)) = H(idx_dummy); %taking element in matrix data with index idx_dummy
%     r_resh = reshape(r_dummy,[size(H,1),size(H,2)]); %reshaping the vector into a matrix shaped as the original data
%     %running PCA in single steps on randomised data
%     fak = bsxfun(@minus,r_resh,mean(r_resh)); %demeaning
%     fakcov = cov(fak'); %covariance matrix
%     [~,fD] = eig(fakcov); %eigenvector problem
%     mfD(pp,1) = max(max(fD)); %storing maximum variance (max eigenvalue) occurring for PCs in randomized data
% end
% MFD = max(mfD); %absolute maximum variance
% %getting eigenvalues from diagonal matrix
% fDv = fD(fD~=0);
% fDv = fDv(end:-1:1);

%%% PLOTTING AND STATISTIC SOLUTIONS FOR REAL AND RANDOMIZED DATA
%plotting variances explained by first 20 PCs for actual and randomised data
if S.fig_l == 1
    figure
    plot(vare(1:100),'DisplayName','data')
    if permnum > 0
        hold on
        plot(varef(1:100),'DisplayName','rand')
    end
    grid minor
    legend('show')
    set(gcf,'color','w')
    title('variance (eigenvalues) of PCs')
    %plotting data using imagesc (sort of raster plots)
    MAX = max(max(max(H(:,:,:))));
    MIN = min(min(min(H(:,:,:))));
    for ii = 1:size(S.H,3)
        figure
        imagesc(time,1:3559,S.H(:,:,ii))
        set(gcf,'color','w')
        colorbar
        caxis([MIN MAX])
        title(['brain sources x time-points - cond ' num2str(ii)])
    end
end

if permnum > 0
    PCs = find(vare>MFD); %significant PCs obtained by getting PCs with variance (eigenvalues) higher than the maximum variance (eigenvalues) obtained from randomized data
else
    PCs = 1:(permnum*(-1));
end
disp('percentage of variance explained by significant PCs')
disp(vare(PCs))

%output structure
OUT = [];
OUT.variance_PCS = vare; %storing variance of significant PCs
if permnum > 0
    OUT.sign_comps_idx = PCs;
else
    OUT.sign_comps_idx = 'MCS has not been run';
end

%normalizing eigenvectors signs
dumones = ones(size(wcoeff,1),size(wcoeff,2)); %vector of 1s with lenngth of significant PCs
switch S.sign_eig
    case 'occurrences'
        % 1) normalizing eigenvectors sign by using mean of the negative/positive values occurrences
        dumones(:,mean(wcoeff > 0) < 0.5) = -1; %assigning -1 to eigenvectors that have more positive than negative weights
    case 'max_abs'
        % 2) normalizing eigenvectors sign on the basis of the sign of the maximum value in absolute terms
        [~,mi] = max(abs(wcoeff)); %getting maximum values indices
        ab = zeros(1,length(mi));
        for ii = 1:length(mi) %over significant PCs
            ab(1,ii) = wcoeff(mi(ii),ii); %storing original values with signs (corresponding to maximum in absolute terms)
        end
        dumones(:,sign(ab)<0) = -1;
    case 'average'
        % 3) normalizing eigenvectors sign on the basis of the sign of the average of weights
        mVV = mean(wcoeff);
        dumones(:,sign(mVV)<0) = -1;
end
VV2 = wcoeff .* dumones; %getting proper sign for eigenvectors
OUT.W = VV2; %weights of PCA
J = zeros(size(wcoeff,2),PCs(end),size(S.H,3));
for ii = 1:size(S.H,3) %over experimental conditions
    J(:,:,ii) = S.H(:,1:size(wcoeff,2),ii)' * wcoeff(:,PCs); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
end
OUT.PCA_timeseries = J; %storing PCA time series
OUT.time = time;

%plotting eigenvectors in the brain
if S.fig_l == 1
    %plotting weights in the brain (nifti images)
    maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
    for ii = 1:PCs(end) %over significant PCs
        %using nifti toolbox
        SO = wcoeff(:,ii);
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
        disp(['saving nifti image - component ' num2str(ii)])
        save_nii(nii,[S.namenii 'PC_' num2str(ii) '_Var_' num2str(round(vare(ii))) '.nii']); %printing image
    end
    %plotting time series of the significant PCA components
    if S.onefig == 1
        figure
        for pp = 1:length(PCs) %over significant PCA components
            for ii = 1:size(S.H,3) %over experimental conditions
                plot(time(1:size(wcoeff,2)),J(:,pp,ii),'LineWidth',2,'DisplayName',['Cond ' num2str(ii) ' PCA ' num2str(pp) ' Var ' num2str(vare(pp))])
                hold on
            end
        end
        xlim([time(1) time(end)])
        set(gcf,'color','w')
        legend('show')
        grid minor
    else
        for pp = 1:length(PCs) %over significant PCA components
            figure
            for ii = 1:size(S.H,3) %over experimental conditions
                plot(time(1:size(wcoeff,2)),J(:,pp,ii),'LineWidth',2,'DisplayName',['Cond ' num2str(ii) ' PCA ' num2str(pp) ' Var ' num2str(vare(pp))])
                hold on
            end
            xlim([time(1) time(end)])
            set(gcf,'color','w')
            legend('show')
            grid minor
        end
    end
end

end


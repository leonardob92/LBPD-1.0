function O = GED_AvgSingleCovs(S)
O = []; 


S = [];

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia')

% Settings
srate    = 250;   % recorded at 1000 Hz, downsampled by 4
stimfrex = 2.439; % stimulation frequency: hard-coded 1 / 0.410 s (range is 0.400 and 0.420 but 0.410 is the mode) 
fwhm  = .3;       % filter width
shr   = 0.01;     % value between 0 and 1 (0 small perturbation, 1 big perturbation) - for Mattia regularisation (until 0.001 it works for restoring the full rank of the matrix)

% For data import
path_sources = '/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1/';
cd(path_sources);
list = dir ([path_sources 'SUBJ*']);

% datapoints = [1:10000]';

% Compute covariance atrices
for subi = 1:length(list)
  
    disp(['Loading data and computing covariance of ' list(subi).name])   
    
    % Import broadband source data
    load([path_sources list(subi).name])
    broadData = OUT.sources_ERFs(:,:);
    
    % Initialize cube of covariance matrices upon first iteration
    if subi == 1
        [covR,covS] = deal( zeros(size(broadData,1), size(broadData,1), length(list)) );
    end
    
    % Filter data
    narrowData = filterFGx(broadData,srate,stimfrex,fwhm,0); % turn vis_flag to 0 for suppressing the plot 

    % Compute covariance matrices
    covS(:,:,subi) = cov(narrowData');
    covR(:,:,subi) = cov(broadData');

end

% Compute average covariance matrices, over subjects
covS_avg = squeeze( mean(covS, 3 )); 
covR_avg = squeeze( mean(covR, 3 ));

% regularisation (to bring covariance matrices to full rank)
%narrow data
covS_avg = covS_avg  + 1e-6*eye(size(covS_avg )); %making matrix full rank again by adding a small perturbation/noise (matrix regularisation)
evalsR = eig(covR_avg );
covR_avg  = (1-shr)*covR_avg  + shr * mean(evalsR) * eye(size(covR_avg ));

% visualize covariance matrices
% figure, clf
% subplot(121)
% imagesc(covS_avg ), axis square
% title('Covariance S (Narrowband)')
% colorbar
% subplot(122)
% imagesc(covR_avg ), axis square
% title('Covariance R (Broadband)')
% colorbar

% Eigendecomposition 
[evecs,evals] = eig(covS_avg ,covR_avg );
[evals,sidx]  = sort( diag(evals),'descend' );
evecs = evecs(:,sidx);
% compute filter forward model and flip sign
ncomps = 100;
GEDmap = zeros(ncomps,size(broadData,1));
GEDts  = zeros(ncomps,size(broadData,2));
disp('computing time series')
for compi = 1:ncomps
    GEDmap(compi,:) = evecs(:,compi)'*covS_avg ; % get component
    [~,idxmax] = max(abs(GEDmap(:,1)));     % find max magnitude
    GEDmap(compi,:)  = GEDmap(compi,:)*sign(GEDmap(compi,idxmax)); % possible sign flip
    % compute time series (projections)
    GEDts(compi,:) = evecs(:,compi)'*broadData;
    disp(compi)
end

save('/scratch7/MINDLAB2017_MEG-LearningBach/Mattia/after_maxfilter_mc/Source_LBPD/Beam_abs_0_sens_1_freq_broadband_invers_1/Test_GED/bomba.mat','covS_avg','covS_avg','evecs','evals','GEDmap','GEDts','covS','covR','srate','-v7.3')


end


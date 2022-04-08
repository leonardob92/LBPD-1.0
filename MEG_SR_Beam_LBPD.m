function [ OUT ] = MEG_SR_Beam_LBPD( S )

% Function to run the several steps required for source reconstruction
% using a beamforming algorithm.
% The algorithms used here are based on work of other brilliant people and
% I do not want to take any credit for that. This function was created simply
% to handle several steps more easily, make small different choices, try various solutions,
% and also use such algorithms in connection to the cluster of computers of Aarhus University
% (for this particular purpose, set the correspondent label below).
% Furthermore, writing these codes allowed me to deeply study the beamforming as well as to prepare
% materials for teaching purposes. Specifically, this function is based on
% a combination of codes and functions from FieldTrip, OSL and SPM
% (beamforming toolbox), plus some in-house-built solutions.
% Thus, to run this function you should first download OSL (https://github.com/OHBA-analysis/osl-core),
% which comprises also the required functions from FieldTrip and SPM. Then,
% you should also start OSL up so that all functions paths are known by Matlab.
% If you work with Aarhus University (and in some cases Oxford University)
% facilities, I would suggest you to contact me to set everything up.. best way could be using the set of functions named: LBPD.
% Leonardo Bonetti, leonardo.bonetti@clin.au.dk

% The function handles the following steps:
% I) NORMALIZING MEG SENSORS, MAINLY FOR COVARIANCE ESTIMATION (OPTIONAL BUT STRONGLY RECOMMENDED)
% II) GRID (i), LEADFIELD MODEL (ii), MEG SENSORS DATA EXTRACTION (iii), SPATIAL FILTERS COMPUTATION (iv)
% III) FULL INVERSION
% IV) SPATIAL AND/OR TEMPORAL SMOOTHING (OPTIONAL) (TO BE UPDATED SINCE I HAVE TO CHECK WHETHER THE NEW MASK IS COMPATIBLE WITH THE ORIGINAL ONE IN TERMS OF MNI COORDINATES, ETC.)
% V) CREATION OF NIFTI IMAGES (OPTIONAL)




% INPUT:    -S.norm_megsensors (I)
%            S.norm_megsensors.zscorel_cov: 1 = z-score normalization; 0 = no normalization;
%                                           (SUGGESTED 1!)
%            S.norm_megsensors.workdir:     working directory (where computation will be made and output stored; highly recommended
%                                           to have a different workdir per subject.
%            S.norm_megsensors.MEGdata_e:   path to MEG data (epoched)
%            S.norm_megsensors.freq:        frequency range in Hertz to be used for the source reconstruction (e.g. [2 8]); leave empty [] to not apply any filter
%            S.norm_megsensors.MEGdata_c:   path to MEG data (continuous).
%                                           This must be passed only if you ask more filtering (so if S.norm_megsensors.freq is not empty).
%                                           This is done since you should do the filtering on the continuous (and not on the epoched) data.
%            S.norm_megsensors.forward:     forward solution for leadfield; at the moment should be: 'Single Shell'

%           -S.beamfilters (II)
%            S.beamfilters.sensl:           1 = magnetometers; 2 = gradiometers; 3 = both MEG sensors (mag and grad) (SUGGESTED 3!)
%            S.beamfilters.maskfname:       path to brain mask: (e.g. 8mm MNI152-T1: '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz')

%           -S.inversion (III)
%            S.inversion.znorml:            1 for inverting MEG data using the zscored normalized data.
%                                           0 to normalize the original data with respect to maximum and minimum of the experimental conditions if you have both magnetometers and gradiometers.
%                                           0 to use original data in the inversion if you have only mag or grad (while e.g. you may have used zscored-data for covariance'matrix)
%                                           (SUGGESTED 0 IN BOTH CASES!)
%                                           (PLEASE, NOTE THAT THE COVARIANCE MATRIX IS ALWAYS CALCULATED BEFORE WITH THE OPTION PROVIDED WITH: "S.norm_megsensors.zscorel_cov")
%            S.inversion.timef:             data-points to be extracted (e.g. 1:300); leave it empty [] for working on the full length of the epoch
%            S.inversion.conditions:        cell with characters for the labels of the experimental conditions (e.g. {'Old_Correct','New_Correct'})
%            S.inversion.bc:                specify the time-samples to be used for baseline correction (as usual, subtraction of the mean of the baseline from the whole timeseries; independently done for each brain source).
%                                           e.g. [1 15] means first 15 time-samples (e.g. 100ms if you have sampling rate = 150Hz).
%                                           Leave empty [] for non calculating it.
%            S.inversion.abs:               1 for absolute values of brain sources values; 0 otherwise.
%                                           Please, note that the sign of the reconstructed signal is not really meaningful or at least it may be normalized after this function if you have knowledge of the polarity of the signal you are reconstructing.. THIS LAST STATEMENT MUST BE CHECKED WITH SOME TESTS!! IT SEEMED FINE IN OAT, BUT I'M NOT SURE YET HERE..
%                                           (SUGGESTED 1!)

%            S.inversion.effects:           Defining solution for computing main effects and contrasts:
%                                           -Main effects of experimental conditions:
%                                            1 = simple mean over trials
%                                            2 = mean divided by standard deviation
%                                            3 = mean divided by square root of standard deviation
%                                            4 = mean divided by (standard deviation divided by square root of number of trials) (actual t-value..)
%                                            5 = single trials (i.e. all trials independently)
%                                           -Contrasts between two experimental conditions:
%                                            1 = simple difference of means over trials
%                                            2-3-4 = two-sample t-tests,

%           -S.smoothing (IV)
%            S.smoothing.spatsmootl:        1 for spatial smoothing; 0 otherwise
%            S.smoothing.spat_fwhm:         spatial smoothing fwhm (suggested = 100)
%            S.smoothing.tempsmootl:        1 for temporal smoothing; 0 otherwise
%            S.smoothing.temp_param:        temporal smoothing parameter (suggested = 0.01)
%            S.smoothing.tempplot:          vector with sources indices to be plotted (original vs temporally smoothed timeseries; e.g. [1 2030 3269]).
%                                           Leave empty [] for not having any plot.

%           -S.nifti (V):                   1 for plotting nifti images of the reconstructed sources of the experimental conditions
%                                           Please, note that no statistics is computed; that will be done by other functions/scripts.

%           -S.Aarhus_cluster:              1 for using the cluster of computers (Haydes) of Aarhus (CFIN-MIB).
%                                           Please, note that to use this you need to have access to the Aarhus' facilities.
%                                           Please, contact me, Leonardo Bonetti (leonardo.bonetti@clin.au.dk), for further information.
%           -S.out_name:                    name (character) for saving output results and nifti images (conditions name is automatically detected and added) on disk


% OUTPUT:   -OUT:                           struct (both outputted and saved in S.norm_megsensors.workdir) with:
%                                           -MEG sensor data (ERFs)
%                                           -source reconstructed data (data before and after smoothing is stored)
%                                           -S structure with information used for the computations






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 19/02/2021


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







if S.Aarhus_cluster == 1
    %starting LBPD up (and therefore also OSL and a subset of SPM and FieldTrip functions)
    pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
    addpath(pathl);
    LBPD_startup_D(pathl);
end

%preparing output structure
OUT = [];
OUT.S = S;
%checking input for baseline correction
if ~isempty(S.inversion.bc) && length(S.inversion.bc) ~= 2
    error('S.inversion.bc (baseline correction) must be either empty [] or comprise two values (e.g. [1 15])')
end

    
%%%%% I - MEG SENSORS NORMALIZATION %%%%%
disp('%%%%% I - MEG SENSORS NORMALIZATION %%%%%')
%loading SPM object (D)
workdir = S.norm_megsensors.workdir;
if ~exist(workdir,'dir') %creating working folder if it does not exist
    mkdir(workdir)
end
%filtering, if requested
if ~isempty(S.norm_megsensors.freq) %if you have provided a frequency range
    if ~isfield(S.norm_megsensors,'MEGdata_c') || isempty(S.norm_megsensors.MEGdata_c) %checking if you have provided continuous data for filtering purposes
        error('if you want to apply the bandpass filter you need to provide continuous data (in S.norm_megsensors.MEGdata_c)')
    end   
    disp('filtering MEG sensor data')
    S3 = [];
    S3.D = S.norm_megsensors.MEGdata_c; %continuous data
    S3.freq = S.norm_megsensors.freq; %frequency range
    S3.band = 'bandpass';
%     S.prefix = 'f_LBPD';
    D = spm_eeg_filter(S3); %actual function
    disp('epoching MEG filtered sensor data with epochinfo from your previously epoched data')
    D_e = spm_eeg_load(S.norm_megsensors.MEGdata_e); %loading previously computed epoched data to get epochinfo (info about trials specification, labels, etc.)
    epochinfo = D_e.epochinfo; %getting the epochinfo information
    %building structure for spm_eeg_epochs
    S3 = [];
    S3.D = D; %filtered continuous data
    S3.trl = epochinfo.trl; %trial information from epochinfo
    S3.conditionlabels = epochinfo.conditionlabels; %condition labels
    S3.prefix = 'e_LBPD'; %prefix to clearly identify that I am doing a further epoching
    D = spm_eeg_epochs(S3); %actual function
    %store the epochinfo structure inside the D object
    D.epochinfo = epochinfo;
    %switch the montage to 1 (ICA denoised data; before was 0 (non-ICA-denoised data) since OSL people did it.. however should not really matter.. what is important in my understanding is that here it is switched back to 1)
    D = D.montage('switch',1);
    D.save(); %saving on disk
    Dfilte = D; %copying file for potential later usage when extracting MEG sensors for inversion
    clear D_e
else
    D = spm_eeg_load(S.norm_megsensors.MEGdata_e); %otherwise loading non-filtered (presumably broad band) data
end

%cloning SPM object
[workdiroriginal, fname, ext] = fileparts(D.fullfile); %getting information about the path of the data
if S.norm_megsensors.zscorel_cov == 1
    Dnew = D.clone( fullfile(workdir,['Clone' fname ext]) );
    idxMEGc1 = find(~cellfun(@isempty,strfind(D.chanlabels,'MEG'))); %elaborated way to get indices of MEG channels (even though I  am mainly interested in the first MEG channel index)
%     dat = D(idxMEGc1,:,:); %extracting only MEG channels
%     for jj = 1:size(dat,3) %over trials
%         dat(:,:,jj) = zscore(dat(:,:,jj),[],2);
%         disp(['normalizing trial ' num2str(jj) ' / ' num2str(size(dat,3))])
%     end
%     dat = zeros(length(idxMEGc1),size(D,2),size(D,3));
    for jj = 1:size(D,3) %over trials
        Dnew(idxMEGc1,:,jj) = zscore(D(idxMEGc1,:,jj),[],2);
%         dat(:,:,jj) = zscore(D(idxMEGc1,:,jj),[],2);
        disp(['normalizing trial ' num2str(jj) ' / ' num2str(size(D,3))])
    end
%     Dnew(idxMEGc1,:,:) = dat; %storing normalized MEG channels
    Dnew.save();
    D = Dnew;
    %%% OBS!!! NOTE THAT AFTER CLONE, HERE WE TAKE DATA FROM AFTER
    %%% AFRICA (SO THE ICA-DENOISED DATA) BUT WE HAVE TO STORE IT IN
    %%% THE MONTAGE = 0 OF THE CLONED SPM OBJECT.. SO IN THE CLONE
    %%% MONTAGE = 0 SHOULD CONTAIN THE ICA-DENOISED DATA AFTER ZSCORE
    %%% NORMALIZATION..
% else %otherwise just passing the data
%     Dnew(:,:,:) = D(:,:,:);
%     Dnew.save();
end
% clear D Dnew
% if S.norm_megsensors.zscorel_cov == 1 %this is just to make sure to properly work on normalized or non-normalized MEG data
%     D = spm_eeg_load([workdir '/Clone' fname ext]);
% else
if S.norm_megsensors.zscorel_cov ~= 1 %this is just to make sure to properly work on normalized or non-normalized MEG data
    if exist('Dfilte','var')
        D = Dfilte; %getting filtered file, if you previously computed it
    else %otherwise original epoched data provided by the user
        D = spm_eeg_load(S.norm_megsensors.MEGdata_e);
    end
end
S2 = [];
S2.forward_meg = S.norm_megsensors.forward; %- Specify forward model for MEG data ('Single Shell' or 'MEG Local Spheres')
D = osl_forward_model(D,S2);
D.save();
%previous specifications from SPM beamforming toolbox
space = 'MNI-aligned';
val = 1;
gradsource = 'inv';
%retrieving data from SPM object
BF.data = spm_eeg_inv_get_vol_sens(D, val, space, gradsource);
idxMEGc1 = find(~cellfun(@isempty,strfind(D.chanlabels,'MEG'))); %elaborated way to get indices of MEG channels (even though I  am mainly interested in the first MEG channel index)
BF.data.MEG.sens.label = D.chanlabels(idxMEGc1)'; %to prevent possible discrepancies between MEG channel labels..
%emulating SPM beamforming toolbox
BF.data.D = D;
BF.data.mesh = D.inv{val}.mesh;
BF.data.zscorel_cov = S.norm_megsensors.zscorel_cov; %storing information about normalization of MEG data
%1 BF file for each subject
save([workdir '/BF_' fname '.mat'],'BF');



%%%%% II - GRID (i), LEADFIELD MODEL (ii), MEG SENSORS DATA EXTRACTION (iii), SPATIAL FILTERS COMPUTATION (iv) %%%%%
disp('%%%%% II - GRID (i), LEADFIELD MODEL (ii), MEG SENSORS DATA EXTRACTION (iii), SPATIAL FILTERS COMPUTATION (iv) %%%%%')
sensl = S.beamfilters.sensl; %getting information about MEG sensors to be used
% timef = S.beamfilters.timef; %getting time-points to be used

%%% 1) COMPUTING GRID (8-MM MNI152 STANDARD MASK TRANSFORMED INTO SPACE WHERE WE DO THE BEAMFORMING)
disp('computing grid')
%loading prepared data (preprocessed, epoched, coregistered (rhino), (maybe normalized), forward model done (vol in ft))
load([workdir '/BF_' fname '.mat']);
[ mni_coords, ~ ] = osl_mnimask2mnicoords(S.beamfilters.maskfname);
%emulating SPM beamforming toolbox
S2 = [];
S2.pos = mni_coords;
% transform MNI coords in MNI space into space where we are doing the beamforming
M = inv(BF.data.transforms.toMNI);
gridv.pos = S2.pos;
%actual transformation of MNI coordinates into the space where we are doing the beamforming
res = ft_transform_geometry(M, gridv);
% establish index of nearest bilateral grid point (for potential use in lateral beamformer); NOT REALLY USED HERE..
res.bilateral_index = zeros(size(S2.pos,1),1);
for jj = 1:size(S2.pos,1)
    mnic = S2.pos(jj,:);
    mnic(1) = -mnic(1);
    res.bilateral_index(jj) = nearest_vec(S2.pos,mnic);
end;
BF.sources.mni = res;

%%% 2) COMPUTING LEADFIELD MODEL
disp('computing leadfield model')
%simulating SPM beamforming toolbox
BF.sources.pos = BF.sources.mni.pos;
BF.sources.ori = [];
nvert = size(BF.sources.pos, 1); %number of dipoles
modalities = {'MEG', 'EEG'};
reduce_rank = 2;
%getting some information
chanind = indchantype(BF.data.D, {'MEG', 'MEGPLANAR'}, 'GOOD'); %channels index
% chanunits = units(BF.data.D, chanind); %channels unit
%preparing data for computation of lead field model
[vol, sens] = ft_prepare_vol_sens(BF.data.(modalities{1}).vol, BF.data.(modalities{1}).sens, 'channel', chanlabels(BF.data.D, chanind));
%getting locations of dipoles (derived from original MNI coordinates from standard space mask)
pos = BF.sources.pos;
%plotting together grid in the brain, vertices and triangles from segmentation and MEG channels
%here I think I am plotting the native space with the dipoles derived from standard mask in MNI coordinates transformed into native space
figure
ft_plot_vol(vol, 'edgecolor', [0 0 0], 'facealpha', 0);
hold on
plot3(pos(:, 1), pos(:, 2), pos(:, 3), '.r', 'MarkerSize', 10);
hold on
plot3(BF.data.MEG.sens.chanpos(:, 1), BF.data.MEG.sens.chanpos(:, 2), BF.data.MEG.sens.chanpos(:, 3), '.b', 'MarkerSize', 10);
rotate3d on;
axis off
axis vis3d
axis equal
%getting label indices (the order of the triplets is fine since the
%leadfield function reads it from sens.label and now I am  also reading that
%from sens.label for extracting proper channels (e.g. magnetometers only)
%from leadfield computation
cntmag = 0; cntgrad2 = 0; cntgrad3 = 0; IDXMAG = zeros(102,1); IDXGRAD2 = zeros(102,1); IDXGRAD3 = zeros(102,1); 
for gg = 1:length(sens.label) %over MEG channels
    if strcmp(sens.label{gg}(end),'1') %if it is a magnetometer
        cntmag = cntmag + 1;
        IDXMAG(cntmag) = gg; %storing index
    elseif strcmp(sens.label{gg}(end),'2') %same concept for gradiometers (first of the couple)
        cntgrad2 = cntgrad2 + 1;
        IDXGRAD2(cntgrad2) = gg;
    elseif strcmp(sens.label{gg}(end),'3') %second of the couple
        cntgrad3 = cntgrad3 + 1;
        IDXGRAD3(cntgrad3) = gg;
    end
end
IDXGRAD = sort(cat(1,IDXGRAD2,IDXGRAD3)); %combining indices of planar gradiometers (according to progressive order in sens.label)
%preallocating space for lead field
L = cell(1,nvert);
Ltot = cell(1,nvert);
%actual computation of lead field
for ii = 1:nvert
    dumL = ft_compute_leadfield(BF.sources.pos(ii, :), sens, vol, 'reducerank', reduce_rank(1)); %computation of the leadfield using FieldTrip function
    if sensl == 1 %magnetometers
        %previous
%         L{ii} = dumL(1:3:306,:);
        %new
        L{ii} = dumL(IDXMAG,:);
    elseif sensl == 2 %gradiometers
        %previous
%         dumdumL = zeros(204,size(dumL,2));
%         dumdumL(1:2:204,:) = dumL(2:3:306,:);
%         dumdumL(2:2:204,:) = dumL(3:3:306,:);
%         L{ii} = dumdumL;
        L{ii} = dumL(IDXGRAD,:);
    else %all MEG sensors
        L{ii} = dumL;
    end
    Ltot{ii} = dumL; %storing full leadfield for later plotting computations
    %     L{i} = ft_compute_leadfield(BF.sources.pos(ii, :), sens, vol, 'reducerank', reduce_rank(1),'normalize', 'yes'); %trying to use the normalize option to reduce the bias towards the center of the head
    disp(['computing leadfield model.. source: ' num2str(ii) ' / ' num2str(nvert)])
end



%%%%%% CONSIDER TO REMOVE THIS IF EVERYTHING WORKS FINE
% %%% 3) EXTRACTING MEG SENSORS DATA
% disp('extracting MEG sensors data')
% %loading and extracting data (you need to decide whether the normalization should be used for all processes of source reconstruction or only for computing the covariance matrix)
% %z-scored normalized data (if you previously asked for normalization, otherwise simply clone of the original MEG data)
% Dc = spm_eeg_load([workdir '/Clone' fname ext]);
% if ~exist('timef','var')
%     warning('timef was not specified, so working on full length of the epoch.. you may want to check whether you need a subset of the time-window or the full epoch..')
%     timef = 1:size(D,2);
% end
% %extracting data (only MEG sensors and first timef time-points)
% idxMEGc1 = find(strcmp(Dc.chanlabels,'MEG0111')); %getting extreme channel index (first MEG channel, then I add 305 to get them all)
% chans = Dc.chanlabels(idxMEGc1:idxMEGc1+305); %extracting labels for MEG channels
% if sensl == 1 %magnetometers
%     datameg = Dc(idxMEGc1:3:idxMEGc1+305,timef,:); %data
% elseif sensl == 2 %gradiometers
%     datameg = zeros(204,length(timef),size(Dc,3));
%     datameg(1:2:204,:,:) = Dc(idxMEGc1+1:3:idxMEGc1+305,timef,:); %data
%     datameg(2:2:204,:,:) = Dc(idxMEGc1+2:3:idxMEGc1+305,timef,:); %data
% else
%     datameg = Dc(idxMEGc1:idxMEGc1+305,timef,:); %data
% end
% time = Dc.time(timef);
% ERFmean = mean(datameg(:,:,1:80),3);
% ERFmedian = median(datameg(:,:,1:80),3);
%%%%%%% UNTIL HERE

%RESHAPING ORDER OF MNI COORDINATES (FROM ORIGINAL STANDARD BRAIN TO BRAIN LAYOUT BUILT FOR BETTER VISUALIZATION PURPOSES, you find it here: /projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_coord_dyi.mat)
%this is important for printing the nifti image with the results to be visualized
disp('reshaping order of MNI coordinates.. just for handling better some visualization purposes..')
% if ~exist('MNI8','var') %if not already loaded..
warning('loading MNI152-T1 8mm brain newly sorted set of coordinates.. remember that if you want a different spatial resolution (e.g. 2mm), you need to create a new mask and update this line of code!!')
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_coord_dyi.mat') %loading coordinates for the new IDs
% end
% if ~exist('mni_coords','var') %if not already loaded..
%     [ mni_coords, ~ ] = osl_mnimask2mnicoords('/scratch5/MINDLAB2017_MEG-LearningBach/MNI152_T1_8mm_brain_PROVA.nii.gz');
% end
MNI_sorted_index = zeros(length(mni_coords),1);
for ii = 1:length(MNI8) %over coordinates (ii = coordinates from file: MNI152_8mm_coord_dyi)
    coorddum  = MNI8(ii,:); %getting new IDs coordinates
    source = find((double(mni_coords(:,1)==coorddum(1)) + double(mni_coords(:,2)==coorddum(2)) + double(mni_coords(:,3)==coorddum(3))) == 3); %finding old ID coordinates (from leadfield) corresponding to new IDs coordinates (image MNI152_8mm_coord_dyi.nii.gz)
    MNI_sorted_index(ii,1) = source; %storing correspondance between MNI8 and mni_coords
end
L2 = L(MNI_sorted_index); %sorting leadfield model (selected MEG channels only)
Ltot2 = Ltot(MNI_sorted_index); %sorting leadfield model (all MEG channels)
pos2 = pos(MNI_sorted_index,:); %sorting coordinates of brain sources for later plotting purposes
%storing outputs
OUT.L = L2; %leadfield model only for requested MEG channels (either magnetometers or gradiometers)
OUT.Ltot = Ltot2; %all MEG channels leadfield model (useful for later plotting purposes)
OUT.pos_brainsources_MNI8 = pos2; %position of coordinates
OUT.vol = vol; %volume conduction model (vol) from FieldTrip
OUT.sens = sens; %information about MEG channels
OUT.transforms = BF.data.transforms; %transform matrix from and to MNI and native
OUT.IDXMAG = IDXMAG; %storing magnetometers indices
OUT.IDXGRAD = IDXGRAD; %storing gradiometers indices


%%% 3) COMPUTING SPATIAL FILTERS
disp('extracting MEG sensors data for covariance estimation')
conditions = S.inversion.conditions;
%%% EXTRACTING MEG SENSORS DATA FOR COVARIANCE ESTIMATION
% disp('extracting MEG sensors data')
%loading and extracting data (z-scored normalized data)
if S.norm_megsensors.zscorel_cov == 1 %this is just to make sure to properly work on normalized or non-normalized MEG data
    Dc = spm_eeg_load([workdir '/Clone' fname ext]);
else
    warning('you are computing covariance matrix on MEG sensors which are not zscore-normalized.. so you should expect much stronger contribution of magnetometers than gradiometers in the reconstructed signal')
    warning('please, note that even if you use only magnetometers or gradiometers, I still suggest you to zscore-normalize the data for covariance matrix computation')
    if exist('Dfilte','var')
        Dc = Dfilte; %getting filtered file, if you previously computed it
    else %otherwise original epoched data provided by the user
        Dc = spm_eeg_load(S.norm_megsensors.MEGdata_e);
    end
end
% Dc = spm_eeg_load([workdir '/Clone' fname ext]);
if isempty(S.inversion.timef)
    warning('timef was not specified, so working on full length of the epoch.. you may want to check whether you need a subset of the time-window or the full epoch..')
    timef = 1:size(D,2);
else
    timef = S.inversion.timef;
end
%getting indices of conditions
condcat = [];
for ii = 1:length(conditions) %over condition labels
    condcat = cat(2,condcat,find(strcmp((Dc.conditions),conditions{ii}))); %getting the condition numbers ordered according to the specification of the user (concatenated together since the covariance matrix is computed on all experimental conditions trials)
end
%extracting data (only MEG sensors and first timef time-points)
%previous
% idxMEGc1 = find(strcmp(Dc.chanlabels,'MEG0111')); %getting extreme channel index (first MEG channel, then I add 305 to get them all)
% if isempty(idxMEGc1)
%     warning('there is no channel named MEG0111, is your data from Elekta-Neuromag?')
%     warning('trying to load MEG 0111, in case it is simply a problem of labels..')
%     idxMEGc1 = find(strcmp(Dc.chanlabels,'MEG 0111')); %getting extreme channel index (first MEG channel, then I add 305 to get them all)
% end
%new
% idxMEGc1 = find(strcmp(Dc.chanlabels,'MEG0111')); %getting extreme channel index (first MEG channel, then I add 305 to get them all)
idxMEGc1 = find(~cellfun(@isempty,strfind(Dc.chanlabels,'MEG'))); %elaborated way to get indices of MEG channels (even though I  am mainly interested in the first MEG channel index)
% chans = Dc.chanlabels(idxMEGc1:idxMEGc1+305); %extracting labels for MEG channels
if S.beamfilters.sensl == 1 %magnetometers
    %previous
    %     datameg = Dc(idxMEGc1:3:idxMEGc1+305,timef,condcat); %data
    %new
    datameg = Dc(idxMEGc1(1)-1 + IDXMAG,timef,condcat); %mag indices previously computed, here adding idxMEGc1(1), in case you have other channels before MEG channels; here I assume that once you have the first MEG channel then you have all the other ones, without any interruption
elseif sensl == 2 %gradiometers
    %this solution would be better.. here and in general..
    %gradi = sort(cat(2,[2:3:306],[3:3:306])); magi = 1:3:306;
    %datameg =  Dc(idxMEGc1+gradi,timef,condcat);
    %primitive codes.. previous
    %     datameg = zeros(204,length(timef),length(condcat));
    %     datameg(1:2:204,:,1:length(condcat)) = Dc(idxMEGc1+1:3:idxMEGc1+305,timef,condcat); %data
    %     datameg(2:2:204,:,1:length(condcat)) = Dc(idxMEGc1+2:3:idxMEGc1+305,timef,condcat); %data
    %new
    datameg = Dc(idxMEGc1(1)-1 + IDXGRAD,timef,condcat); %mag indices previously computed, here adding idxMEGc1(1), in case you have other channels before MEG channels; here I assume that once you have the first MEG channel then you have all the other ones, without any interruption 
else
    datameg = Dc(idxMEGc1,timef,condcat); %data
end
time = Dc.time(timef);
OUT.time = time; %storing time for later purposes

disp('computing spatial filters')
%YOU MAY CONSIDER TO TRY ALSO WITH COVARIANCE OVER TRIALS INSTEAD THAT OVER TIME.. BUT THIS IS MORE EXPERIMENTAL AND I WOULD DO THAT LATER..)
%please, note that you compute the covariance matrix on the same (amount of) data that you want to invert later
%THE SOLUTION BELOW (meancovl = 1 or meancovl = 3) SEEMS TO PROVIDE ALMOST THE SAME RESULTS (COVARIANCE OF CONCATENATED TRIALS SHOULD BE BETTER IN PRINCIPLE AND IT IS ALSO SLIGHTLY FASTER TO COMPUTE)
%conversely, as expected, meancovl = 2 is not fine
%I do not delete this since it may be useful to read it later on
%     meancovl = 3; %THIS SHOULD BE = 3
%     if meancovl == 1 %mean of trials covariance
%         Ct = zeros(size(datameg,1),size(datameg,1),size(datameg,3));
%         Cinvt = zeros(size(datameg,1),size(datameg,1),size(datameg,3));
%         for cc = 1:size(datameg,3) %over trials
%             Ct(:,:,cc) = cov(datameg(:,:,cc)'); %computing covariance matrix (time-sample x MEG channels, so covariance of MEG channels according to Matlab cov function)
%             Cinvt(:,:,cc) = pinv_plus(Ct(:,:,cc),50); %Mark Woolrich inv; 50 is the pca_dim related to Maxfilter computation (which makes the matrix of the data with a rank = 50)
%         end
%         C = mean(Ct,3); %mean over trials covariance
%         Cinv = mean(Cinvt,3); %mean over inv of trials coavariance
%     elseif meancovl == 2 %covariance calculated over mean of the trials
%         C = cov(ERFmean');
%         Cinv = pinv_plus(C,50); %Mark Woolrich inv; 50 is the pca_dim related to Maxfilter computation (which makes the matrix of the data with a rank = 50)
%     elseif meancovl == 3


%computing covariance between MEG sensors (trials concatenated)
datacov = reshape(datameg,size(datameg,1),size(datameg,2)*size(datameg,3))'; %data was extracted from SPM object previously loaded (2 sections above)
C = cov(datacov); %computing covariance matrix (time-sample x MEG channels, so covariance of MEG channels according to Matlab cov function)
Cinv = pinv_plus(C,50); %Mark Woolrich inv; 50 is the pca_dim related to Maxfilter computation (which makes the matrix of the data with a rank = 50)
%this function uses at most X (in this case 50) PCs to estimate inverted matrix; in this case 50 should stabilize the computation of pinv 
%     end
W = cell(1,length(L2)); %initializing spatial filters
WW = zeros(length(L2),size(L2{1},1));
%computing the spatial filters (independently for each MEG source)
for ii = 1:length(L2) %over brain sources
    llf = L2{ii}; %extracting brain source ii
    %combining dipole orientations
    tmp = llf' * Cinv * llf; %temporary estimation of covariance between sources
    nn = size(llf,2); %number of orientations
    [u, ~] = svd(real(pinv_plus(tmp(1:nn,1:nn),2,0)),'econ'); %single value decomposition, code from Mark Woolrich and Robert Oostenveld; (2 would be equal to the "reduce_rank" variable..)
    eta = u(:,1);
    llf = llf * eta; %actual reduction of the 3 orientations to 1 global vector
    %compuation of weights
    W{ii} = pinv_plus(llf' * Cinv * llf,2,0) * llf' * Cinv; %actual computation of spatial filter (Mark Woolrich pinv); I guess equivalent to: W{ii} = pinv_plus(llf' * Cinv * llf) * llf' * Cinv; pinv_plus(data,2,0) should not use 2 as maximum number of PCs.. anyway, I reproduce  previous codes here..
    W{ii} = W{ii}./sqrt(W{ii}*W{ii}'); %applying a weight normalization to counterbalance the bias towards the centre of the head, as suggested by Woolrich, Huang, Van Veen
    %     W{ii} = pinv(llf' * Cinv * llf) * llf' * Cinv; %actual computation of spatial filter
    WW(ii,:) = W{ii}(1,:);
%     disp(['computing spatial filters.. source: ' num2str(ii) ' / ' num2str(length(L2))])
end
%plotting images of filters (with weight normalisation)
figure
imagesc(WW)
colorbar
%storing output
OUT.W = W;



%%%%% III - FULL INVERSION %%%%%
disp('%%%%% III - FULL INVERSION %%%%%')
%%% EXTRACTING MEG SENSORS DATA
disp('extracting MEG sensors data for actual inversion')
%loading and extracting data
if S.inversion.znorml == 1
    %z-scored normalized data
    Dc = spm_eeg_load([workdir '/Clone' fname ext]);
else
    %original data (non-normalized)
    if exist('Dfilte','var')
        Dc = Dfilte; %getting filtered file, if you previously computed it
    else %otherwise epoched data provided by the user
        Dc = spm_eeg_load(S.norm_megsensors.MEGdata_e);
    end
end
%previous codes
%extracting data (only MEG sensors and first timef time-points)
% idxMEGc1 = find(strcmp(Dc.chanlabels,'MEG0111')); %getting extreme channel index (first MEG channel, then I add 305 to get them all)
% if isempty(idxMEGc1)
%     warning('there is no channel named MEG0111, is your data from Elekta-Neuromag?')
%     warning('trying to load MEG 0111, in case it is simply a problem of labels..')
%     idxMEGc1 = find(strcmp(Dc.chanlabels,'MEG 0111')); %getting extreme channel index (first MEG channel, then I add 305 to get them all)
% end
% % chans = Dc.chanlabels(idxMEGc1:idxMEGc1+305); %extracting labels for MEG channels
% if S.beamfilters.sensl == 1 %magnetometers
%     datameg = Dc(idxMEGc1:3:idxMEGc1+305,timef,:); %data
% elseif sensl == 2 %gradiometers
%     datameg = zeros(204,length(timef),size(Dc,3));
%     datameg(1:2:204,:,:) = Dc(idxMEGc1+1:3:idxMEGc1+305,timef,:); %data
%     datameg(2:2:204,:,:) = Dc(idxMEGc1+2:3:idxMEGc1+305,timef,:); %data
% else
%     datameg = Dc(idxMEGc1:idxMEGc1+305,timef,:); %data
% end
%new codes
idxMEGc1 = find(~cellfun(@isempty,strfind(Dc.chanlabels,'MEG'))); %elaborated way to get indices of MEG channels (even though I  am mainly interested in the first MEG channel index)
if S.beamfilters.sensl == 1 %magnetometers
    datameg = Dc(idxMEGc1(1)-1 + IDXMAG,timef,:); %mag indices previously computed, here adding idxMEGc1(1), in case you have other channels before MEG channels; here I assume that once you have the first MEG channel then you have all the other ones, without any interruption
elseif sensl == 2 %gradiometers
    datameg = Dc(idxMEGc1(1)-1 + IDXGRAD,timef,:); %mag indices previously computed, here adding idxMEGc1(1), in case you have other channels before MEG channels; here I assume that once you have the first MEG channel then you have all the other ones, without any interruption 
else
    datameg = Dc(idxMEGc1,timef,:); %data
end
time = Dc.time(timef);
%getting indices of conditions (Old_Correct and New_Correct)
conds = cell(1,length(conditions));
for ii = 1:length(conditions) %over condition labels
    conds{ii} = find(strcmp((Dc.conditions),conditions{ii})); %getting the condition numbers ordered according to the specification of the user in S.conditions (sum to avoid the problem of zero-vector occurring if it does not find the match between the specified label and the data label.. of course if you specify two times the same label in the structure S it crashes as well.. but you would need to be quite peculiar to do that so I assume you will not do that..)
end
%meanll = 1 or = 0 seem to return very similar (almost identical) results; this is slightly surprising, I would expect similar but not almost identical results.. anyway.. the procedure seems correct, so I'll just live with it..
% meanll = 1; % (SUGGESTED 1!) 1 for source reconstruction of mean over trials; 0 for mean over single-trials source reconstruction
%I did not delete these lines since it may useful to read this information in the future
% if meanll == 1 %mean over trials and then source reconstruction
%computing mean of ERFs (consider to compute also median later..)
ERFs = cell(1,length(conditions));
for ii = 1:length(conditions) %over experimental conditions
    if ~isempty(conds{ii}) %it mis empty if there are is no match between condition requested by user and condition in MEG data (it may happen for a few subjects..)
        if S.inversion.effects == 1 %simple mean over trials
            ERFs{ii} = mean(datameg(:,:,conds{ii}),3);
        elseif S.inversion.effects == 2 %mean divided by st
            ERFs{ii} = mean(datameg(:,:,conds{ii}),3) ./ std(datameg(:,:,conds{ii}),0,3);
        elseif S.inversion.effects == 3 %mean divided by square root of st
            ERFs{ii} = mean(datameg(:,:,conds{ii}),3) ./ sqrt(std(datameg(:,:,conds{ii}),0,3)); %(variation of t-value, dividing mean over trials by their standard error
        elseif S.inversion.effects == 4 %mean divided by (st divided by square root of number of trials) (actual t-value..)
            ERFs{ii} = mean(datameg(:,:,conds{ii}),3) ./ (std(datameg(:,:,conds{ii}),0,3)./sqrt(length(conds{ii}))); %t-value (sort of), dividing mean over trials by their standard error
        elseif S.inversion.effects == 5 %all trials
            ERFs{ii} = datameg(:,:,conds{ii});
        end
    end
end
%this should be ok, but a further check could be suggested..
if S.beamfilters.sensl == 3 && S.inversion.znorml ~= 1 %if you have both MEG sensors and they should be normalized with reference to their overall maximum value
    %getting maximum and minimum values experimental conditions together
    maxx = zeros(length(conditions),1);
    minn = zeros(length(conditions),1);
    rmaxx = zeros(length(conditions),1);
    rminn = zeros(length(conditions),1);
    indg = sort([2:3:306 3:3:306]); %index of gradiometers
    for cc = 1:length(conditions) %over conditions
        if S.inversion.effects < 5
            maxx(cc) = max(max(abs(ERFs{cc}))); %maximum for each conditions (all MEG sensors, so magnetometers..)
            minn(cc) = min(min(abs(ERFs{cc}))); %minimum for each condition (all MEG sensors, so magnetometers..)
            rmaxx(cc) = max(max(abs(ERFs{cc}(indg,:)))); %maximum for each conditions (only gradiometers)
            rminn(cc) = min(min(abs(ERFs{cc}(indg,:)))); %minimum for each condition (only gradiometers)
        else
            maxx(cc) = max(max(max(abs(ERFs{cc})))); %maximum for each conditions (all MEG sensors, so magnetometers..)
            minn(cc) = min(min(min(abs(ERFs{cc})))); %minimum for each condition (all MEG sensors, so magnetometers..)
            rmaxx(cc) = max(max(max(abs(ERFs{cc}(indg,:,:))))); %maximum for each conditions (only gradiometers)
            rminn(cc) = min(min(min(abs(ERFs{cc}(indg,:,:))))); %minimum for each condition (only gradiometers)
        end
    end
    rmax = max(maxx); %global maximum over all conditions (all MEG sensors, so magnetometers..)
    rmin = min(minn); %global minimum over all conditions (all MEG sensors, so magnetometers..)
    rmaxg = max(rmaxx); %global maximum over all conditions (only gradiometers)
    rming = min(rminn); %global minimum over all conditions (only gradiometers)
    for cc = 1:length(conditions) %over conditions
        %this should be important since we want the experimental conditions
        %to maintain their different relative strengths (if you zscore
        %normalized (X-mean(X)./std(X), you should loose their relative
        %different strengths.. thus now:
        % 1) COVARIANCE MATRIX COMPUTED ON Z-SCORE NORMALIZED MEG DATA
        % 2) INVERSION COMPUTED ON MEG DATA SCALED WITH RESPECT TO RELATIVE DIFFERENCES BETWEEN EXPERIMENTAL CONDITIONS
        if S.inversion.effects < 5
            dumones = ones(size(ERFs{cc},1),size(ERFs{cc},2));
            dumones(ERFs{cc}<0) = -1; %getting a matrix with values -1 and 1 according to the original sign in ERFs{cc}
            ERFs{cc} = abs(ERFs{cc}); %absolute value
            ERFs{cc}(indg,:) = (ERFs{cc}(indg,:)-rming)./(rmaxg-rming).*(rmax-rmin) + rmin; %normalizing gradiometers on the basis of magnetometers (with respect to the conditions different strengths)
            ERFs{cc} = ERFs{cc} .* dumones; %assigning the original sign to the data
        else
            dumones = ones(size(ERFs{cc},1),size(ERFs{cc},2),size(ERFs{cc},3));
            dumones(ERFs{cc}<0) = -1; %getting a matrix with values -1 and 1 according to the original sign in ERFs{cc}
            ERFs{cc} = abs(ERFs{cc}); %absolute value
            ERFs{cc}(indg,:,:) = (ERFs{cc}(indg,:,:)-rming)./(rmaxg-rming).*(rmax-rmin) + rmin; %normalizing gradiometers on the basis of magnetometers (with respect to the conditions different strengths)
            ERFs{cc} = ERFs{cc} .* dumones; %assigning the original sign to the data
        end
    end
end
%until here


disp('computing inversion')
if S.inversion.effects < 5
    sources_ERFs = zeros(length(W),length(time),length(conditions));
    for cc = 1:length(conditions) %over conditions
        if ~isempty(ERFs{cc}) %if there was no condition cc in MEG data
            for ii = 1:length(W) %over spatial filters (brain sources)
                for tt = 1:length(time) %over time-points
                    sources_ERFs(ii,tt,cc) = W{ii} * ERFs{cc}(:,tt); %actual inversion
                end
            end
        end
    end
else %single trials independently
    %     sources_ERFs = zeros(length(W),length(time),length(conditions));
    s_tr_ERFs = cell(length(conditions),1); %there may be a different number of trials in different conditions (e.g. in learningbach-recogminor), so you need to sotre different conditions in cells
    for cc = 1:length(conditions) %over conditions
        if ~isempty(ERFs{cc}) %if there was no condition cc in MEG data
            sources_ERFs = zeros(length(W),length(time),size(ERFs{cc},3));
            for ii = 1:length(W) %over spatial filters (brain sources)
                for tt = 1:length(time) %over time-points
                    for rr = 1:size(ERFs{cc},3) %over trials
                        sources_ERFs(ii,tt,rr) = W{ii} * ERFs{cc}(:,tt,rr); %actual inversion
                    end
                end
                disp(['spatial filter ' num2str(ii) ' / ' num2str(length(W)) ' - condition ' num2str(cc) ' / ' num2str(length(conditions))])
            end
            s_tr_ERFs{cc} = sources_ERFs;
        end
    end
end
%baseline correction (subtracting mean activity in the baseline)
if ~isempty(S.inversion.bc)
    disp('computing baseline correction')
    if S.inversion.effects < 5
        for cc = 1:length(conditions) %over conditions
            if ~isempty(ERFs{cc}) %if there was no condition cc in MEG data
                for ii = 1:length(W) %over spatial filters (brain sources)
                    sources_ERFs(ii,:,cc) = sources_ERFs(ii,:,cc) - mean(sources_ERFs(ii,S.inversion.bc(1):S.inversion.bc(2),cc),2);
                end
            end
        end
    else %single trials independently
        for cc = 1:length(conditions) %over conditions
            if ~isempty(ERFs{cc}) %if there was no condition cc in MEG data
                sources_ERFs = s_tr_ERFs{cc};
                for ii = 1:length(W) %over spatial filters (brain sources)
                    for rr = 1:size(s_tr_ERFs{cc},3) %over trials
                        sources_ERFs(ii,:,rr) = sources_ERFs(ii,:,rr) - mean(sources_ERFs(ii,S.inversion.bc(1):S.inversion.bc(2),rr),2);
                    end
                    disp(['baseline correction - condition ' num2str(cc) ' / ' num2str(length(conditions))])
                end
                s_tr_ERFs{cc} = sources_ERFs;
            end
        end
    end
end
% else %source reconstruction of single-trials and then mean
%     %extracting single-trial data
%     ERF_Old_st = datameg(:,:,conds{1}); %condition OLD
%     ERF_New_st = datameg(:,:,conds{2}); %condition NEW
%     sources_OLD = zeros(length(W),length(time));
%     sources_NEW = zeros(length(W),length(time));
%     for ii = 1:length(W) %over spatial filters (brain sources)
%         for tt = 1:length(time) %over time-points
%             %OLD condition
%             temptr = zeros(size(ERF_Old_st,3),1);
%             for rr = 1:size(ERF_Old_st,3) %over trials (OLD)
%                 temptr(rr,1) = W{ii} * ERF_Old_st(:,tt,rr); %actual inversion
%             end
%             sources_OLD(ii,tt) = mean(temptr); %mean over trials for source ii and time-point tt
%             %NEW condition
%             temptr = zeros(size(ERF_New_st,3),1);
%             for rr = 1:size(ERF_New_st,3) %over trials (NEW)
%                 temptr(rr,1) = W{ii} * ERF_New_st(:,tt,rr); %actual inversion
%             end
%             sources_NEW(ii,tt) = mean(temptr); %mean over trials for source ii and time-point tt
%         end
%         disp(['source ' num2str(ii) ' / ' num2str(length(W))])
%     end
% end
%scatter plotting (to check the distribution of the sources at a specific time-point)
% figure
% scatter(1:length(W),sources_OLD(:,timeselected),'o','r')
% hold on
% scatter(1:length(W),sources_NEW(:,timeselected),'o','b')

%storing outputs
OUT.data_MEG_sensors = ERFs;
if S.inversion.effects < 5
    if S.inversion.abs == 1 %absolute values, if requested
        sources_ERFs = abs(sources_ERFs);
    end
    OUT.sources_ERFs = sources_ERFs; %no absolute values
else
    if S.inversion.abs == 1 %absolute values, if requested
        for cc = 1:length(conditions) %over conditions
            s_tr_ERFs{cc} = abs(s_tr_ERFs{cc});
        end
    end
    OUT.sources_ERFs = s_tr_ERFs; %no absolute values
end


%%% MUST BE UPDATED SINCE NOW I GUESS IT WORKS WITH PREVIOUS MNI COORDINATES..
%%%%% IV - SMOOTHING %%%%%
%actual computations
%spatial smoothing
% if S.smoothing.spatsmootl == 1
%     disp('%%%%% IV - SMOOTHING %%%%%')
%     fwhm = S.smoothing.spat_fwhm; %spatial smoothing parameter
%     x1 = S.smoothing.temp_param;
%     %spatial smoothing
%     %preparing inputs
%     disp('computing spatial smoothing')
%     ERFs_spat_smooth = zeros(size(sources_ERFs,1),size(sources_ERFs,2),size(sources_ERFs,3));
%     for cc = 1:size(sources_ERFs,3) %over conditions
%         vol_as_matrix = sources_ERFs(:,:,cc); %data
%         %plotting original image
%         figure; imagesc(vol_as_matrix); title('original')
%         mask_fname = S.beamfilters.maskfname; %path to standard brain mask (MAYBE HERE THE NEW MASK??)
%         %actual function
%         ERFs_spat_smooth(:,:,cc) = smooth_vol_osl(vol_as_matrix, mask_fname, fwhm);
%         %plotting smoothed image
%         figure; imagesc(ERFs_spat_smooth(:,:,cc)); title(['spatially smoothed - condition ' S.inversion.conditions{cc}])
%     end
%     %storing outputs
%     OUT.spatial_smoothing = ERFs_spat_smooth;
%     sources_ERFs = ERFs_spat_smooth; %assigning spatially smoothed matrix to previously computed "full inversion" matrix since you want the spatially smoothed data to be passed to the next steps (temporally smoothing and printing nifti images) 
% else
%     OUT.spatial_smoothing = [];
% end
%temporal smoothing

if S.smoothing.tempsmootl == 1 && S.inversion.effects < 5
    disp('computing temporal smoothing')
    %preparing inputs (same for all experimental conditions)
    x2 = 1; %dimension 1 of matrix requested by osl_gauss
    x3 = size(sources_ERFs,2); %dimension 2 of matrix requested by osl_gauss
    f = fftshift(osl_gauss(x1/(time(2)-time(1)),x2,x3)'); %current_level.time_smooth_std/tres,1,length(datstd))');
    for cc = 1:size(sources_ERFs,3) %over conditions
        dat = sources_ERFs(:,:,cc); %data (condition cc)
        %plotting original image
        figure; imagesc(dat); title('original')
        ERFs_temp_smooth = zeros(size(sources_ERFs,1),size(sources_ERFs,2),size(sources_ERFs,3));
        for ii = 1:size(sources_ERFs,1) %over brain voxels
            ERFs_temp_smooth(ii,:,cc) = fftconv(dat(ii,:)',f);
        end
        %plotting temporally smoothed image
        figure; imagesc(ERFs_temp_smooth(:,:,cc)); title(['temporally smoothed - condition ' S.inversion.conditions{cc}])
        if ~isempty(S.smoothing.tempplot) %plotting difference between original and temporal smoothed data (only for requested brain voxels)
            for pp = 1:length(S.smoothing.tempplot)
                figure
                plot(dat(S.smoothing.tempplot(pp),:))
                hold on
                plot(fftconv(dat(S.smoothing.tempplot(pp),:)',f))
            end
        end
    end
    %storing output
    OUT.ERFs_temp_smooth = ERFs_temp_smooth;
    sources_ERFs = ERFs_temp_smooth; %assigning temporally smoothed matrix to previously computed "full inversion" matrix since you want the printing nifti images to be run on temporally (and maybe also spatially) smoothed data
else
    OUT.ERFs_temp_smooth = [];
    if S.smoothing.tempsmootl == 1 && S.inversion.effects == 5
        warning('spatial and temporal smoothing are implemented for averaged ERFs only')
    end
end




%%%%% V - PRINTING NIFTI IMAGE(s) %%%%%
if S.nifti == 1
    if S.inversion.effects < 5
        disp('%%%%% V - PRINTING NIFTI IMAGE(s) %%%%%')
        % %getting maximum and minimum values experimental conditions together (used for scaling images)
        % maxx = zeros(length(conditions),1);
        % minn = zeros(length(conditions),1);
        % for cc = 1:length(conditions)
        %     maxx(cc) = max(max(abs(sources_ERFs(:,:,cc)))); %maximum for each conditions
        %     minn(cc) = min(min(abs(sources_ERFs(:,:,cc)))); %minimum for each condition
        % end
        % rmax = max(maxx); %global maximum over all conditions
        % rmin = min(minn); %global minimum over all conditions
        warning('loading MNI152-T1 8mm brain newly sorted set of coordinates.. remember that if you want a different spatial resolution (e.g. 2mm), you need to create a new mask and update this line of code!!')
        maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
        for iii = 1:length(conditions)
            fnamenii = [workdir '/' S.out_name '_cond_' conditions{iii} '_norm' num2str(S.inversion.znorml) '_spatsmoot_' num2str(S.smoothing.spatsmootl) '_tempsmoot_' num2str(S.smoothing.tempsmootl) '_abs_' num2str(S.inversion.abs) '.nii.gz']; %path and name of the image to be saved
            
            %     %normalization of nifti image
            %     SO = (abs(sources_ERFs(:,:,iii))-rmin)./(rmax-rmin).*(100-0) + 0; %normalizing according to some scale (e.g. 0-100)
            %     if S.inversion.abs ~= 1 %absolute values, if requested
            %         dumones = ones(size(sources_ERFs,1),size(sources_ERFs,2)); %matrix of ones with dimensions of sources_ERFs(:,:,iii)
            %         dumones(sources_ERFs(:,:,iii)<0) = -1; %getting a matrix with values -1 and 1 according to the original sign in sources_ERFs(:,:,iii)
            %         SO = SO .* dumones; %rebuilding the original sign of the matrix with scaled values
            %     end
            %without normalization
            SO = sources_ERFs(:,:,iii);
            %building nifti image
            SS = size(maskk.img);
            dumimg = zeros(SS(1),SS(2),SS(3),length(time));
            for ii = 1:size(sources_ERFs,1) %over brain sources
                dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
                [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
                dumimg(i1,i2,i3,:) = SO(ii,:); %storing values for all time-points in the image matrix
            end
            nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
            nii.img = dumimg; %storing matrix within image structure
            nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
            disp(['saving nifti image - condition ' num2str(iii)])
            save_nii(nii,fnamenii); %printing image
        end
    else
        warning('printing nifti images is implemented only for averaged ERFs')
    end
end

%saving output on disk
disp('saving results')
if S.inversion.effects < 5 %normally small(ish) file
    save([S.norm_megsensors.workdir '/' S.out_name '_norm' num2str(S.inversion.znorml) '_abs_' num2str(S.inversion.abs) '.mat'],'OUT')
else %single trial data can be way bigger..
    save([S.norm_megsensors.workdir '/' S.out_name '_norm' num2str(S.inversion.znorml) '_abs_' num2str(S.inversion.abs) '.mat'],'OUT','-v7.3')
end

%deleting Clone data and variable BF
disp('deleting Clone data, BF and, if any, filtered MEG sensor data')
delete([workdir '/BF_' fname '.mat']) %BF
delete([workdir '/Clone' fname '.mat']) %Clone .mat
delete([workdir '/Clone' fname '.dat']) %Clone .dat
if ~isempty(S.norm_megsensors.freq) %if you have provided a frequency range
    delete([workdiroriginal '/' fname(7:end) '.mat']) %deleting the filtered MEG sensor data for saving space on disk (this does not affect the source reconstruction) .mat file
    delete([workdiroriginal '/' fname(7:end) '.dat']) %deleting the filtered MEG sensor data for saving space on disk (this does not affect the source reconstruction) .dat file
    delete([workdiroriginal '/' fname '.mat']) %epoched data by LBPD .mat
    delete([workdiroriginal '/' fname '.dat']) %epoched data by LBPD .dat
end

if S.Aarhus_cluster == 1
    OUT = []; %Aarhus cluster usually do not have output like normal Matlab functions, but saves results on disk
end


end


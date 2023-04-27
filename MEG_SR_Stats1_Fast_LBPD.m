function [ O ] = MEG_SR_Stats1_Fast_LBPD( S )
O = [];

% It computes statistics after beamforming (computed using MEG_SR_Beam_LBPD.m).
% Main effects here are simply mean(X)./sqrt(std(X)) (or just mean(X), while contrast
% between 2 conditions are actual t-tests.
% This function assumes that data is approximately normally distributed
% which is often the case of this kind of data.
% The function works for MEG sources, but also for MEG sensor results
% (loaded from the computation in MEG_SR_Beam_LBPD.m). In this case, it
% also combines planar gradiometers.
% Currently, it works only with source reconstructed data on the aggregated
% trials (i.e. averaged or similarly aggregated trials).



%  INPUT:   -S.workingdir:         path where the data from MEG_SR_Beam_LBPD.m is stored
%           -S.sensl:              1 = magnetometers only;
%                                  2 = gradiometers only;
%                                  3 = both magnetometers and gradiometers.
%           -S.plot_nifti:         1 to plot nifti images; 0 otherwise
%           -S.plot_nifti_name:    character with name for nifti files (and saved mat files) (it may be useful if you run separate analysis)
%                                  Leave empty [] to not  specify any name
%           -S.contrast:           one row per contrast (e.g. having 3 conditions, [1 -1 0; 1 -1 -1; 0 1 -1];
%                                  two or more -1 or 1 are interpreted as the mean over them first and then the contrast.
%                                  Leave empty [] for no contrasts.
%           -S.list:               list of subjects (mat files).
%                                  Leave empty [] (or not define the field) for automatic indexing of mat files of all subjects in S.workingdir
%           -S.effects:            Defining solution for computing main effects and contrasts:
%                                  -Main effects of experimental conditions:
%                                   1 = simple mean over subjects
%                                   2 = mean divided by standard deviation
%                                   3 = mean divided by square root of standard deviation
%                                   4 = mean divided by (standard deviation divided by square root of number of subjects) (actual t-value..)
%                                  -Contrasts between two experimental conditions:
%                                   1 = simple difference of means over subjects
%                                   2-3-4 = two-sample t-tests,
%                                   %%% MAYBE ALWASYS COMPUTE T-VALUES %%%
%           -S.subave:             cell matrix with requested subaverages in the first row (e.g. S.subave{1,1} = [1 3]; S.subave{1,2} = [2 4 5];
%                                  S.subave{1,3} = [6,8]; subaverages conditions 1 and 3, 2, 4 and 5, etc.).
%                                  Names in the second row (e.g. S.subave{2,1} = 'name1'; S.subave{2,2} = 'name2'; etc.).                                 
%                                  Leave empty S.subave = [] or do not provide field for no subaveraging.
%           -S.Aarhus_clust:       1 to use paralle computing (Aarhus University, contact me, Leonardo Bonetti, for more information)
%                                  leonardo.bonetti@clin.au.dk

%  OUTPUT:  -statistics (and possibly nifti images) saved on disk (in S.workingdir)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 06/03/2021
% Leonardo Bonetti, Oxford, UK, 27/04/2023


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    









if S.Aarhus_clust == 1
    %LBPD_startup_D
    pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
    addpath(pathl);
    LBPD_startup_D(pathl);
end

%getting a few inputs
sensl = S.sensl;
workingdir = S.workingdir;
mkdir(workingdir)
%MAIN EFFECTS CONDITIONS
if ~isfield(S,'list') || isempty(S.list)
    list = dir([workingdir '/SUBJ*mat']);
else
    list = S.list;
end
load([list(1).folder '/' list(1).name]) %assuming that all subjects have same dimensions and information.. loading first one to get some information about the data
SS_s = size(OUT.sources_ERFs); %getting size of sources
SR = zeros(SS_s(1),SS_s(2),SS_s(3),length(list)); %sources
SS_c = size(OUT.data_MEG_sensors{1}); %getting size of MEG channels data (only first condition here)
SR_c = zeros(SS_c(1),SS_c(2),length(OUT.data_MEG_sensors),length(list)); %sensors
for ii = 1:length(list) %over subjects
    load([list(ii).folder '/' list(ii).name]);
    %sources
    SR(:,:,:,ii) = OUT.sources_ERFs;
    %sensors
    for cc = 1:length(OUT.S.inversion.conditions) %over conditions
        if ~isempty(OUT.data_MEG_sensors{cc}) %it may be empty if condition cc for subject ii was not recorded during MEG data acquisition
            if sensl == 3 %magnetometers and gradiometers
                %magnetometers
                SR_c(103:end,:,cc,ii) = OUT.data_MEG_sensors{cc}(1:3:306,:); %condition cc (OLD first then NEW), magnetometers
                %gradiometers
                planar1 = OUT.data_MEG_sensors{cc}(2:3:306,:); %extracting planar 1
                planar2 = OUT.data_MEG_sensors{cc}(3:3:306,:); %extracting planar 2
                SR_c(1:102,:,cc,ii) = sqrt(planar1.^2 + planar2.^2); %combining planar gradiometers
            elseif sensl == 2 %gradiometers
                planar1 = OUT.data_MEG_sensors{cc}(1:2:204,:); %extracting planar 1
                planar2 = OUT.data_MEG_sensors{cc}(2:2:204,:); %extracting planar 2
                SR_c(:,:,cc,ii) = sqrt(planar1.^2 + planar2.^2); %combining planar gradiometers
            elseif sensl == 1 %magnetometers
                SR_c(:,:,cc,ii) = OUT.data_MEG_sensors{cc}; %condition cc (OLD first then NEW), magnetometers
            end
        end
    end
    disp(['loading subjects - ' num2str(ii)])
end
disp('computing and saving main effects - MEG sources and sensors')
SR(SR==0) = NaN; %if not all subjects had all conditions, replacing 0s with NaNs


%%% OBS!! %%%
%%% HERE I ALLOW TO DO SUBAVERAGES OF CONDITONS FOR MAIN EFFECT; THIS IS FINE BUT NOT GREATLY OPTIMIZED SINCE THE CODES ALREADY WORKED FOR DOING SUBAVERAGE BEFORE THE CONTRASTS..
%%% YOU SHOULD CONSIDER MAKE IT BETTER.. AND MAYBE REMOVE THE MEG SENSOR DATA FROM HERE.. %%%
if isfield(S,'subave') %this is to allow previous scripts to still work even if S.subave is not specified
    if ~isempty(S.subave) %if user requested subaveraging
        DUM = zeros(SS_s(1),SS_s(2),length(S.subave),length(list));
        for ii = 1:size(S.subave,2) %over requested subaverages
            DUM(:,:,ii,:) = mean(SR(:,:,S.subave{1,ii},:),3); %subaveraging set of conditions ii
        end
        SR = DUM;
        clear DUM
    end
end
%%% UNTIL HERE %%%

SRm = mean(SR,4,'omitnan'); %average over subjects
SRst = std(SR,0,4,'omitnan'); %computing standard deviations

%%% MAYBE HERE SMOOTHING OF VARIANCES AS DONE IN OSL, OPTION WITH ORIGINAL MASK IS PROPOSED BELOW
% for cc = 1:size(SRst,3) %over conditions
%     SRst(:,:,cc) = smooth_vol_osl(SRst(:,:,cc), '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz', 100);
% end
%%% UNTIL HERE.. I DO NOT SEE A PARTICULAR BENEFIT FROM THIS..

%computing main effects (several options) - MEG sources
if S.effects == 1 %simple mean over subjects
    t_val_s = SRm; 
elseif S.effects == 2 %mean divided by st
    t_val_s = SRm./SRst; 
elseif S.effects == 3 %mean divided by square root of st
    t_val_s = SRm./sqrt(SRst); 
elseif S.effects == 4 %mean divided by (st divided by square root of number of subjects) (actual t-value..)
    t_val_s = SRm./(SRst./sqrt(length(list)));
end

S_struct = OUT.S;
S_struct.S_stats = S; %storing information for the contrasts
if ~isfield(S,'plot_nifti_name') || isempty(S.plot_nifti_name)
    save([workingdir '/sources_main_effects.mat'],'t_val_s','S_struct')
else
    save([workingdir '/' S.plot_nifti_name '_sources_main_effects.mat'],'t_val_s','S_struct')
end
SR_c(SR_c==0) = NaN; %if not all subjects had all conditions, replacing 0s with NaNs
SR_cm = nanmean(SR_c,4); %average over subjects
SR_cst = nanstd(SR_c,0,4); %computing standard deviations
%%% MAYBE HERE SMOOTHING OF VARIANCES AS DONE IN OSL

%computing main effects (several options) - MEG sensors
if S.effects == 1 %simple mean over subjects
    t_val_c = SR_cm; 
elseif S.effects == 2 %mean divided by st
    t_val_c = SR_cm./SR_cst;
elseif S.effects == 3 %mean divided by square root of st
    t_val_c = SR_cm./sqrt(SR_cst);
elseif S.effects == 4 %mean divided by (st divided by square root of number of subjects) (actual t-value..)
    t_val_c = SR_cm./(SR_cst./sqrt(length(list)));
end

% t_val = SR_cm./sqrt(SR_cst); %computing t-values of main effect of the experimental conditions (OLD and NEW)
% t_val = SR_cm./(SR_cst./sqrt(length(list))); %computing t-values of main effect of the experimental conditions (OLD and NEW)

if ~isfield(S,'plot_nifti_name') || isempty(S.plot_nifti_name)
    save([workingdir '/sensors_plancomb_main_effects.mat'],'t_val_c','S_struct')
else
    save([workingdir '/' S.plot_nifti_name '_sensors_plancomb_main_effects.mat'],'t_val_c','S_struct')
end

maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
%nifti images, main effects (if requested)
if S.plot_nifti == 1
    %%%%% GROUP MAIN EFFECTS - PRINTING NIFTI IMAGE(s) %%%%%
    if isempty(S.plot_nifti_name) %if empty, I simply do not add any name to the automatic name of the nifti images
        pnn = 'Main';
    else
        pnn = [S.plot_nifti_name '_Main'];  
    end
    disp('%%%%% V - PRINTING NIFTI IMAGE(s) %%%%%')
    warning('loading MNI152-T1 8mm brain newly sorted set of coordinates.. remember that if you want a different spatial resolution (e.g. 2mm), you need to create a new mask and update this line of code!!')
    %     maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
    for iii = 1:size(SR,3) %length(OUT.S.inversion.conditions) %over experimental conditions
        if isfield(S,'subave') %this is to allow previous scripts to still work even if S.subave is not specified
            if ~isempty(S.subave) %if user requested subaveraging
                fnamenii = [workingdir '/' pnn '_cond_' S.subave{2,iii} '_abs_' num2str(OUT.S.inversion.abs) '.nii.gz']; %path and name of the image to be saved
            else
                fnamenii = [workingdir '/' pnn '_cond_' OUT.S.inversion.conditions{iii} '_abs_' num2str(OUT.S.inversion.abs) '.nii.gz']; %path and name of the image to be saved
            end
        end
        SO = t_val_s(:,:,iii);
        %building nifti image
        SS = size(maskk.img);
        dumimg = zeros(SS(1),SS(2),SS(3),size(SO,2));
        for ii = 1:size(SO,1) %over brain sources
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
end

%computing contrasts, if requested
if ~isempty(S.contrast)
    %computing actual statistical tests (t-tests)
    for aa = 1:size(S.contrast,1) %over the requested contrasts
        %extracting contrasts
        id1 = find(S.contrast(aa,:)==1); %getting index(es) of first element(s) of the contrast
        id2 = find(S.contrast(aa,:)==-1); %getting index(es) of second element(s) of the contrast
        %MEG sources
        if S.effects == 1 %simple mean over subjects
%             t_val_s = (mean(SRm(:,:,id1),3)-mean(SRm(:,:,id2),3)); %simple difference of means
        %%% MAYBE BETTER TO JUST DO REAL T-TESTS HERE.. CONSIDER THIS AND THEN UPDATE ALSO DOCUMENTATION IN THE BEGINNING OF THE FUNCTION %%%
            t_val_s = (mean(SRm(:,:,id1),3)-mean(SRm(:,:,id2),3))./sqrt(((mean(SRst(:,:,id1),3).^2)./length(list))+((mean(SRst(:,:,id2),3).^2)./length(list)));
        %%% UNTIL HERE %%%
        else
            %weighting difference of means by their standard deviations (independent samples t-values..)
            t_val_s = (mean(SRm(:,:,id1),3)-mean(SRm(:,:,id2),3))./sqrt(((mean(SRst(:,:,id1),3).^2)./length(list))+((mean(SRst(:,:,id2),3).^2)./length(list)));
        end
        
        if ~isfield(S,'name') || isempty(S.name)
            save([workingdir '/sources_contrast_' num2str(aa) '.mat'],'t_val_s','S','aa')
        else
            save([workingdir '/' S.name '_sources_contrast_' num2str(aa) '.mat'],'t_val_s','S','aa')
        end
        %MEG sensors
        if S.effects == 1 %simple mean over subjects
%             t_val_c = (mean(SR_cm(:,:,id1),3)-mean(SR_cm(:,:,id2),3)); %simple difference of means

        %%% MAYBE BETTER TO JUST DO REAL T-TESTS HERE.. CONSIDER THIS AND THEN UPDATE ALSO DOCUMENTATION IN THE BEGINNING OF THE FUNCTION %%%
            t_val_c = (mean(SR_cm(:,:,id1),3)-mean(SR_cm(:,:,id2),3))./sqrt(((mean(SR_cst(:,:,id1),3).^2)./length(list))+((mean(SR_cst(:,:,id2),3).^2)./length(list)));
        %%% UNTIL HERE %%%
        
        else
            %weighting difference of means by their standard deviations (independent samples t-values..)
            t_val_c = (mean(SR_cm(:,:,id1),3)-mean(SR_cm(:,:,id2),3))./sqrt(((mean(SR_cst(:,:,id1),3).^2)./length(list))+((mean(SR_cst(:,:,id2),3).^2)./length(list)));
        end
        if ~isfield(S,'plot_nifti_name') || isempty(S.plot_nifti_name)
            save([workingdir '/sensors_plancomb_contrast_' num2str(aa) '.mat'],'t_val_c','S','aa')
        else
            save([workingdir '/' S.plot_nifti_name '_sensors_plancomb_contrast_' num2str(aa) '.mat'],'t_val_c','S','aa')
        end
        if S.plot_nifti == 1
            if isempty(S.plot_nifti_name) %if empty, I simply do not add any name to the automatic name of the nifti images
                pnn = 'Contr_';
            else
                pnn = [S.plot_nifti_name '_Contr_'];
            end
            fnamenii = [workingdir '/' pnn  num2str(aa) '_abs_' num2str(OUT.S.inversion.abs) '.nii.gz']; %path and name of the image to be saved
            %building nifti image
            SS = size(maskk.img);
            dumimg = zeros(SS(1),SS(2),SS(3),size(t_val_s,2));
            for ii = 1:size(t_val_s,1) %over brain sources
                dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
                [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
                dumimg(i1,i2,i3,:) = t_val_s(ii,:); %storing values for all time-points in the image matrix
            end
            nii = make_nii(dumimg,[8 8 8]);
            nii.img = dumimg; %storing matrix within image structure
            nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
            disp(['saving nifti image - contrast ' num2str(aa) ' / ' num2str(size(S.contrast,1))])
            save_nii(nii,fnamenii); %printing image
        end
    end
end


end


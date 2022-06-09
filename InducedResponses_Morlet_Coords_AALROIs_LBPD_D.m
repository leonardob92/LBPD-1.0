function O = InducedResponses_Morlet_Coords_AALROIs_LBPD_D( S )
O = []; 

% It computes time-frequency analysis (induced responses) using Morlet
% wavelet.
% It computes Morlet wavelet transform independently on each voxel and each
% trial of the provided data. Then, it averages the results over voxels and trials.
% It does it on the voxels belonging to a series of region of interests (ROIs), as
% provided by the user.
% This should be used in connection with LBPD functions, after computing
% the beamforming (without using absolute values for the source
% reconstructed data).
% This function can work locally on any computer, but it was conceived for
% the cluster of computer of Aarhus University.



% INPUT:    -S.average_trials:    1: evoked responses (i.e. TF on averaged trials)
%                                 0: induced responses (i.e. TR on each voxel and trial and then average)
%           -S.f:                 frequencies to be used by Morlet
%           -S.conds:             indices of experimental conditions to be used (in connection to LBPD source reconstruction
%           -S.coordd:            vector (n of voxels x 3 (x,y,z in MNI space)) with coordinates to be used.
%                                 Leave empty [] if you want to use a full ALL ROI(s).
%           -S.subjlist:          list of subjects (output of LBPD source reconstruction) provided as the output of the "dir" Matlab function
%           -S.time:              vector of time in seconds (it assumes the source reconstructed data has the same time.
%           -S.outdir:            path and name (without ".mat") to store the results

% OUTPUT:   -file saved on disk with:
%                                   -P2:    time-frequency decomposition: voxels x frequencies x time-samples x conditions x subjects
%                                   -f:     frequencies used
%                                   -time:  time-points used (correspondent values in seconds)
%                                   -S:     structure to summarize inputs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Oxford, UK, 03/06/2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%adding path (meaningful if the function is used on the Aarhus University cluster of computers)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/BrainstormDemo_CodesTimeFrequencyDecomposition/functions') %path to Dimitrios function for Morlet transform
%extracting some information
average_trials = S.average_trials;
time = S.time;
conds = S.conds;
f = S.f;
list = S.subjlist;
coordd = S.coordd;
%getting indices of provided coordinates or AAL ROIs
if ~isempty(coordd)
    %matching coordinates and getting source indices (step by step)
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_coord_dyi.mat');
    idx = zeros(size(coordd,1),1);
    for ii = 1:size(coordd,1) %over voxels
        %previous
%         dum1 = double(coordd(ii,1) == MNI8(:,1)); %x
%         dum2 = double(coordd(ii,2) == MNI8(:,2)); %y
%         dum3 = double(coordd(ii,3) == MNI8(:,3)); %z
%         a = find((dum1+dum2+dum3)==3); %index
%         idx(ii) = a; %storing index   
        %updated
        %x
        dum1 = double(coordd(ii,1) == MNI8(:,1)); %x
        if isempty(find(dum1==1)) %if there is no pefect correspondance between given coordinate and mask coordinate, you look for an approximation
            dum = coordd(ii,1) - MNI8(:,1); %differnce between given coordinate and mask coordinate
            [m1] = min(dum(dum>0)); %getting closest value to 0 (positive)
            [m2] = max(dum(dum<0)); %getting closest value to 0 (negative)
            asd = [m1 m2]; %concatenating the values
            [~,i] = min(abs([m1 m2])); %getting the minimum value n absolute terms
            dum1 = zeros(size(MNI8,1),1);
            dum1(asd(i)==dum) = 1; %assigning 1 to the requested voxels (the ones that are closest to the given coordinate)
        end
        %y
        dum2 = double(coordd(ii,2) == MNI8(:,2)); %y
        if isempty(find(dum2==1)) %if there is no pefect correspondance between given coordinate and mask coordinate, you look for an approximation
            dum = coordd(ii,2) - MNI8(:,2); %differnce between given coordinate and mask coordinate
            [m1] = min(dum(dum>0)); %getting closest value to 0 (positive)
            [m2] = max(dum(dum<0)); %getting closest value to 0 (negative)
            asd = [m1 m2]; %concatenating the values
            [~,i] = min(abs([m1 m2])); %getting the minimum value n absolute terms
            dum2 = zeros(size(MNI8,1),1);
            dum2(asd(i)==dum) = 1; %assigning 1 to the requested voxels (the ones that are closest to the given coordinate)
        end
        %z
        dum3 = double(coordd(ii,3) == MNI8(:,3)); %z
        if isempty(find(dum3==1)) %if there is no pefect correspondance between given coordinate and mask coordinate, you look for an approximation
            dum = coordd(ii,3) - MNI8(:,3); %differnce between given coordinate and mask coordinate
            [m1] = min(dum(dum>0)); %getting closest value to 0 (positive)
            [m2] = max(dum(dum<0)); %getting closest value to 0 (negative)
            asd = [m1 m2]; %concatenating the values
            [~,i] = min(abs([m1 m2])); %getting the minimum value n absolute terms
            dum3 = zeros(size(MNI8,1),1);
            dum3(asd(i)==dum) = 1; %assigning 1 to the requested voxels (the ones that are closest to the given coordinate)
        end
        a = find((dum1+dum2+dum3)==3); %index
        idx(ii) = a; %storing index
    end
else
    AAL_ROIs = S.AAL_ROIs;
    aal = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'); %loading AAL nifti image
    aalt = aal.img; %extracting matrix with the values shown in the image
    mask = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %loading AAL nifti image
    maskt = mask.img; %extracting matrix with the values shown in the image
    idx = [];
    for ii = 1:length(AAL_ROIs) %over AAL ROIs
        dum = find(aalt == AAL_ROIs(ii)); %getting the indices of the selected AAL ROIs
        idx = cat(1,idx,maskt(dum)); %getting the correspondent indices of the mask upon which is based the source reconstruction (and therefore the data)
    end
    idx(idx==0) = []; %removing the few voxels that couldn't be indexed since the AAL ROIs and the mask used slightly different approximations and therefore a few voxels do not perfectly correspond (with MEG spatial resolution it does not really matter)
end
%loading data and computing time-frequency analysis
load([list(1).folder '/' list(1).name]);
ttt = length(OUT.S.inversion.timef);
P2 = zeros(length(idx),length(f),ttt,length(conds),length(list));
for ii = 1:length(list)
    disp(['loading data for subject ' num2str(ii)])
    if ii ~= 1
        load([list(ii).folder '/' list(ii).name]);
    end
    for jj = 1:length(conds) %over conditions
        disp(['subj ' num2str(ii) ' - cond ' num2str(jj)])
        if average_trials ~= 1
            dum = OUT.sources_ERFs{conds(jj)}; %getting single trial data for Old correct
            P = zeros(length(idx),length(f),ttt,size(dum,3));
            for xx = 1:size(dum,3) %over trials
                data = dum(idx,:,xx);
                P(:,:,:,xx) = morlet_transform(data,time(1:size(data,2)),f); %frequency decomposition
            end
            P2(:,:,:,jj,ii) = mean(P,4); %average over trials
            
        else
            data = OUT.sources_ERFs(idx,:,conds(jj));
            P2(:,:,:,jj,ii) = morlet_transform(data,time(1:size(data,2)),f); %frequency decomposition
        end
    end
    clear OUT dum P
end
time = time(1:ttt);
if average_trials ~= 1
    save([S.outdir '_avetr0.mat'],'P2','time','f','S');
else
    save([S.outdir '_avetr1.mat'],'P2','time','f','S');
end


end


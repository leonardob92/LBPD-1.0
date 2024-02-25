function O = InducedResponses_Morlet_ROIs_AarhusClust( S )
O = []; 

% It sends one job for each subject to Aarhus University cluster (for parallel
% computing purposes).
% It makes sense (and therefore it is implemented) only for single_subject
% option.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Oxford, UK, 12/06/2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%adding path.. nmeaningul only if you use this in the Aarhus Universiy cluster of computers setting
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/BrainstormDemo_CodesTimeFrequencyDecomposition/functions') %path to Dimitrios function for Morlet transform
%extracting information
ROI_n = S.ROI_n;
f = S.f;
maskk = S.mask;
list = S.subjlist;
time = S.time;
dop = length(find(maskk(:,ROI_n)==1));
ii = S.ii;
%actual computation
disp(['loading data for subject ' num2str(ii)])
load([list(ii).folder '/' list(ii).name]);
Psubj = zeros(length(f),length(time),length(OUT.sources_ERFs)); %frequencies x time-points x conditions
for cc = 1:length(OUT.sources_ERFs) %over conditions
    dum = OUT.sources_ERFs{cc}; %getting single trial data for Old correct
    data = dum(maskk(:,ROI_n)==1,:,:); %getting single voxels for VMPFC
    P2 = zeros(length(f),length(time),size(dum,3));
    for xx = 1:size(dum,3) %over trials
        Pdum = zeros(dop,length(f),length(time));
        if isfield(S,'baselcorr') %if baseline correction
            if ~isempty(S.baselcorr)
                dumbo = morlet_transform(data(:,:,xx),time,f); %frequency decomposition (for each voxel of ROI and each trial independently)
                Pdum(:,:,:) = dumbo - mean(dumbo(:,:,S.baselcorr(1):S.baselcorr(2)),3); %removing the indicating baseline
            else
                Pdum(:,:,:) = morlet_transform(data(:,:,xx),time,f); %frequency decomposition (for each voxel of ROI and each trial independently)
            end
        else
            Pdum(:,:,:) = morlet_transform(data(:,:,xx),time,f); %frequency decomposition (for each voxel of ROI and each trial independently)
        end
        Pdum = permute(Pdum,[2 3 1]); %permuting order so to have voxels in the 3th dimension
        P2(:,:,xx) = mean(Pdum,3); %average over voxels
        disp(['subject ' num2str(ii) ' - trial ' num2str(xx)])
    end
    Psubj(:,:,cc) = mean(P2,3); %average over trials
end
save([S.outdir(1:end-4) '_' list(ii).name(1:9) '.mat'],'Psubj','time','f','S');




end


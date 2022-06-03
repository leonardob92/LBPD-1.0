function O = InducedResponses_Morlet_ROIs_LBPD_D( S )
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



% INPUT:    -S.ROI_n:    number of ROI to be used (coherently with the ROIs order reported in S.mask)
%           -S.f:        frequencies to be used by Morlet
%           -S.mask:     path (including file) to a binary matrix (sources x ROIs) where 1s indicate
%                        the voxels to be taken into account.
%                        Voxels order is connected to LBPD system of coordinates.
%           -S.subjlist: list of subjects (output of LBPD source reconstruction)
%           -S.time:     vector of time in seconds (it assumes the source reconstructed data has the same time.
%                        If time is not the same, it takes from time-sample 1 up till the end of S.time).
%           -S.outdir:   path and name (including ".mat") to store the results

% OUTPUT:   -file saved on disk with:
%                                   -P3:    time-frequency decomposition: frequencies x time-samples x subjects
%                                   -f:     frequencies used
%                                   -time:  time-points used (correspondent values in seconds)
%                                   -S:     structure to summarize inputs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Oxford, UK, 03/06/2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%adding path.. nmeaningul only if you use this in the Aarhus Universiy cluster of computers setting
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/BrainstormDemo_CodesTimeFrequencyDecomposition/functions') %path to Dimitrios function for Morlet transform

ROI_n = S.ROI_n;
f = S.f;
load(S.mask);
list = S.subjlist;
time = S.time;
P3 = zeros(length(f),length(time),5,length(list)); %frequencies x time-points x conditions x subjects
for ii = 1:length(list)
    disp(['loading data for subject ' num2str(ii)])
    load([list(ii).folder '/' list(ii).name]);   
    for cc = 1:length(OUT.sources_ERFs)
        dop = length(find(ROIs_DCM(:,ROI_n)==1));        
        dum = OUT.sources_ERFs{cc}; %getting single trial data for Old correct
        data = dum(ROIs_DCM(:,ROI_n)==1,:,:); %getting single voxels for VMPFC
        P2 = zeros(length(f),length(time),size(dum,3));
        for xx = 1:size(dum,3) %over trials
            P = zeros(length(f),length(time),dop);
            for ss = 1:dop %over voxels of VMPFC
                P(:,:,ss) = morlet_transform(data(ss,:,xx),time,f); %frequency decomposition (for each voxel of VMPFC and each trial independently)
            end
            P2(:,:,xx) = mean(P,3); %average over voxels
            disp(['subject ' num2str(ii) ' - trial ' num2str(xx)])
        end
        P3(:,:,cc,ii) = mean(P2,3); %average over trials
    end
    clear OUT P2 P data dum
end
save(S.outdir,'P3','time','f','S');





end


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



% INPUT:    -S.ROI_n:           number of ROI to be used (coherently with the ROIs order reported in S.mask)
%           -S.f:               frequencies to be used by Morlet
%           -S.mask:            binary (0 and 1) matrix (sources x ROIs) where 1s indicate the voxels to be taken into account.
%                               Voxels order is connected to LBPD system of coordinates.
%           -S.subjlist:        list of subjects (output of LBPD source reconstruction)
%           -S.time:            vector of time in seconds (it assumes the source reconstructed data has the same time.
%                               If time is not the same, it takes from time-sample 1 up till the end of S.time).
%           -S.Aarhus_clust:    0 = working locally
%                                   integer number (i.e. 1) = sending one job for each subject to the Aarhus cluster (the number you insert
%                                   here corresponds to the slots of memory that you allocate for each job.
%           -S.outdir:          path and name (including ".mat") to store the results
%           -S.single_subj:     1 for saving single subject data; 0 for not saving it


% OUTPUT:   -file saved on disk with:
%                                   -P3:    time-frequency decomposition: frequencies x time-samples x subjects
%                                   -f:     frequencies used
%                                   -time:  time-points used (correspondent values in seconds)
%                                   -S:     structure to summarize inputs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Oxford, UK, 05/06/2022


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
if S.single_subj == 0
    load([list(1).folder '/' list(1).name]);
    P3 = zeros(length(f),length(time),length(OUT.sources_ERFs),length(list)); %frequencies x time-points x conditions x subjects
end
%setting up the cluster
if S.Aarhus_clust == 1
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
    clusterconfig('scheduler', 'cluster'); % 'none' or 'cluster'
    clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit
    clusterconfig('slot', S.Aarhus_clust); %slot in the queu
end
%actual computation (or jobs submission)
for ii = 1:length(list)
    if S.Aarhus_clust == 1
        S.ii = ii;
        jobid = job2cluster(@InducedResponses_Morlet_ROIs_AarhusClust,S);
    else
        disp(['loading data for subject ' num2str(ii)])
        %         if ii ~= 1 %you alredy loaded subject 1
        load([list(ii).folder '/' list(ii).name]);
        %         end
        if S.single_subj == 1
            Psubj = zeros(length(f),length(time),length(OUT.sources_ERFs)); %frequencies x time-points x conditions
        end
        for cc = 1:length(OUT.sources_ERFs) %over conditions
            dum = OUT.sources_ERFs{cc}; %getting single trial data for Old correct
            data = dum(maskk(:,ROI_n)==1,:,:); %getting single voxels for VMPFC
            P2 = zeros(length(f),length(time),size(dum,3));
            for xx = 1:size(dum,3) %over trials
                %             P = zeros(length(f),length(time),dop);
                %             for ss = 1:dop %over voxels of ROI
                %                 P(:,:,ss) = morlet_transform(data(ss,:,xx),time,f); %frequency decomposition (for each voxel of ROI and each trial independently)
                %             end
                %other solution
                Pdum = zeros(dop,length(f),length(time));
                Pdum(:,:,:) = morlet_transform(data(:,:,xx),time,f); %frequency decomposition (for each voxel of ROI and each trial independently)
                Pdum = permute(Pdum,[2 3 1]); %permuting order so to have voxels in the 3th dimension
                P2(:,:,xx) = mean(Pdum,3); %average over voxels
                disp(['subject ' num2str(ii) ' - trial ' num2str(xx)])
            end
            if S.single_subj == 0
                P3(:,:,cc,ii) = mean(P2,3); %average over trials
            else
                Psubj(:,:,cc) = mean(P2,3); %average over trials
            end
        end
        if S.single_subj == 1
            save([S.outdir(1:end-4) '_' list(ii).name(1:9) '.mat'],'Psubj','time','f','S');
        end
        clear OUT P2 Pdum data dum Psubj
    end
end
if S.single_subj == 0
    save(S.outdir,'P3','time','f','S');
end




end


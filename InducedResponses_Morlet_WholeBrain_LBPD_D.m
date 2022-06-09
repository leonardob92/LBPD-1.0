function O = InducedResponses_Morlet_WholeBrain_LBPD_D( S )
O = []; 

% It computes time-frequency analysis (induced responses) using Morlet
% wavelet.
% It computes Morlet wavelet transform independently on each voxel of the brain
% and each trial.
% This should be used in connection with LBPD functions, after computing
% the beamforming (without using absolute values for the source
% reconstructed data).
% This function can work locally on any computer, but it was conceived for
% the cluster of computer of Aarhus University.



% INPUT:    -S.f:                 frequencies to be used by Morlet
%           -S.conds:             indices of experimental conditions to be used (in connection to LBPD source reconstruction
%           -S.subjlist:          list of subjects (output of LBPD source reconstruction) provided as the output of the "dir" Matlab function
%           -S.time:              vector of time in seconds (it assumes the source reconstructed data has the same time)
%           -S.Aarhus_clust:      1 = Aarhus cluster for parallel computing (i.e. running one job per subject);
%                                 0 = work locally
%           -S.outdir:            path and name (without ".mat") to store the results

% OUTPUT:   -file saved on disk with:
%                                   -P2:    time-frequency decomposition: voxels x frequencies x time-samples x conditions
%                                   -f:     frequencies used
%                                   -time:  time-points used (correspondent values in seconds)
%                                   -S:     structure to summarize inputs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Oxford, UK, 09/06/2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%adding path (meaningful if the function is used on the Aarhus University cluster of computers)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/BrainstormDemo_CodesTimeFrequencyDecomposition/functions') %path to Dimitrios function for Morlet transform
%extracting some information
time = S.time;
conds = S.conds;
f = S.f;
list = S.subjlist;
if S.Aarhus_clust == 1
    %setting up the cluster
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
    clusterconfig('scheduler', 'cluster'); % 'none' or 'cluster'
    clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit
    clusterconfig('slot', 3); %slot in the queu
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')
end
for ii = 1:length(list)
    if S.Aarhus_clust ~= 1
        disp(['loading data for subject ' num2str(ii)])
        load([list(ii).folder '/' list(ii).name]);
        ttt = length(OUT.S.inversion.timef); %getting number of time-points
        sss = size(OUT.sources_ERFs{conds(1)},1); %getting number of voxels
        P2 = zeros(sss,length(f),ttt,length(conds));
        for jj = 1:length(conds) %over conditions
            disp(['subj ' num2str(ii) ' - cond ' num2str(jj)])
            dum = OUT.sources_ERFs{conds(jj)}; %getting single trial data for Old correct
            P = zeros(sss,length(f),ttt,size(dum,3));
            for xx = 1:size(dum,3) %over trials
                data = dum(:,:,xx);
                P(:,:,:,xx) = morlet_transform(data,time(1:size(data,2)),f); %frequency decomposition
            end
            P2(:,:,:,jj) = mean(P,4); %average over trials
        end
        time = time(1:ttt);
        save([S.outdir(1:end-4) '_' list(ii).name(1:9) '.mat'],'P2','time','f','S');
    else
        S.ii = ii;
        jobid = job2cluster(@Induced_Resp_SingleSubj_AarhusClust,S);
    end
end


end


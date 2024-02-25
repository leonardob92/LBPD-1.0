function O = Induced_Resp_SingleSub_AAL_Coordsj_AarhusClust( S )
O = []; 

% It sends one job per subject to the Aarhus cluster (for parallel
% computing)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Oxford, UK, 12/06/2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%adding path (meaningful if the function is used on the Aarhus University cluster of computers)
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/BrainstormDemo_CodesTimeFrequencyDecomposition/functions') %path to Dimitrios function for Morlet transform
%extracting some information
ii = S.ii;
idx = S.idx;
average_trials = S.average_trials;
time = S.time;
conds = S.conds;
f = S.f;
list = S.subjlist;
%loading data and computing time-frequency analysis
load([list(ii).folder '/' list(ii).name]);
ttt = length(OUT.S.inversion.timef);
Psubj = zeros(length(idx),length(f),ttt,length(conds));
%actual computation
for jj = 1:length(conds) %over conditions
    disp(['subj ' num2str(ii) ' - cond ' num2str(jj)])
    if average_trials ~= 1
        dum = OUT.sources_ERFs{conds(jj)}; %getting single trial data for Old correct
        P = zeros(length(idx),length(f),ttt,size(dum,3));
        for xx = 1:size(dum,3) %over trials
            data = dum(idx,:,xx);
            if isfield(S,'baselcorr') %if baseline correction
                if ~isempty(S.baselcorr)
                    dumbo = morlet_transform(data,time(1:size(data,2)),f); %computing power
                    P(:,:,:,xx) = dumbo - mean(dumbo(:,:,S.baselcorr(1):S.baselcorr(2)),3); %removing the indicated baseline
                else
                    P(:,:,:,xx) = morlet_transform(data,time(1:size(data,2)),f); %frequency decomposition
                end
            else
                P(:,:,:,xx) = morlet_transform(data,time(1:size(data,2)),f); %frequency decomposition
            end
        end
        Psubj(:,:,:,jj) = mean(P,4); %average over trials
    else
        data = OUT.sources_ERFs(idx,:,conds(jj));
        if isfield(S,'baselcorr') %if baseline correction
            if ~isempty(S.baselcorr)
                dumbo = morlet_transform(data,time(1:size(data,2)),f);
                Psubj(:,:,:,jj) = dumbo - mean(dumbo(:,:,S.baselcorr(1):S.baselcorr(2)),3);
            else
                Psubj(:,:,:,jj) = morlet_transform(data,time(1:size(data,2)),f); %frequency decomposition
            end
        else
            Psubj(:,:,:,jj) = morlet_transform(data,time(1:size(data,2)),f); %frequency decomposition
        end
    end
end
time = time(1:ttt);
%saving on disk
mkdir(S.outdir)
save([S.outdir '/' list(ii).name(1:9) '.mat'],'Psubj','time','f','S');



end


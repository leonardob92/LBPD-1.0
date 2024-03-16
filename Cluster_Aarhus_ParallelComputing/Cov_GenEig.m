function O = Cov_GenEig( S )
O = []; 

% Tests for computing covariance over averaged subjects for generalised
% eigenvector solution.
% Many solutions here are hard coded and the reason to use the cluster is
% to overtake problems about memory.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Aarhus, DK, 05/03/2024


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%main effect computed on averaged trials (for source reconstruction); this
%is much faster and easier since I already had this from previous
%computations
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/sources_main_effects.mat');
%adjusting polarity
timex = 45:52;
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if squeeze(mean(t_val_s(jj,timex,1),2)) > 0 %if the data in voxel jj is positive during N100 time
        vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
    else
        vect(jj,1) = 1;
    end
end



list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/SUBJ*.mat');
clear dum
cnt = 0;
for ss = 1:length(list) %over subjects
    cnt = cnt + 1; %counting subjects
    load([list(ss).folder '/' list(ss).name]);
    sources = OUT.sources_ERFs;
    for cc = 1:2%length(sources) %over conditions (only first 2)
        soudum = sources{cc};
        if size(soudum,3) > 14 %at least 15 trials correct (out of 27)
            dumbo = zeros(size(soudum,1),size(soudum,2),size(soudum,3));
            for rr = 1:size(soudum,3) %over trials
                for jj = 1:size(soudum,1) %over brain sources
                    dumbo(jj,:,rr) = soudum(jj,:,rr) .* vect(jj,1); %reversing (or not)..
                end
            end
            bumba = floor(size(soudum,3)/4); %getting number of trials to be subaveraged
            cntbumba = (bumba*(-1))+1; %initialise bumba
            for dd = 1:4 %over classes of subaveraged trials
                cntbumba = cntbumba + bumba;
                if dd < 4 %proper bumba trials
                    dum(:,:,dd,cc,ss) = mean(dumbo(:,:,cntbumba:cntbumba+(bumba-1)),3);
                else %last class may have more trials than bumba
                    dum(:,:,dd,cc,ss) = mean(dumbo(:,:,cntbumba:end),3);
                end
            end
        end
    end
    disp(ss)
end

dumma = mean(dum,5); %mean over subjects
save([list(1).folder '/GenEig_AverageOverSubjects.mat'],'dumma','cnt','-v7.3')



end


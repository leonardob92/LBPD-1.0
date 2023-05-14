function O = Replay_TG_Plot_Stats(S)
O = []; 

%It computes plots (and/or statistics against 50% chance level (time
%consuming)) for temporal generalisation (TG) decoding.
%It is designed to work on the Aarhus Cluster

% INPUT:

% -S structure with folloring fields:
%   -list_TG:                            directory with paths of files data. Typically obtained using dir() function
%   -stat:                               1 = stats against 50% change level; 0 = only computing and saving on disk average TG matrix
%   -t_thresh:                           threhsold to binarise the t-values in the TG matrix and run a 2-D cluster based correction for multiple comparisons
%   -time:                               time (in seconds)
%   -name:                               name for the file which will store the averageg TG matrix and the statistics
%   -output:                             output directory



% OUTPUT:
% -statistics and average over participants for later plotting (saved on disk)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Oxford, UK, 14/05/2023


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%% TO BE TESTED!! %%

list_TG = S.list_TG;
load([list_TG(1).folder '/' list_TG(1).name]); %load TG file
SS = size(d.d);
dd1 = zeros(SS(1),SS(2),length(list_TG)); %preallocate variable dd1 (time samples and TG files)
for ff = 1:length(list_TG) %over subjects (TG files)
    load([list_TG(ff).folder '/' list_TG(ff).name]); %load TG file
    dd1(:,:,ff) = d.d; %save d structure
    disp(ff)
end
ddTGm = mean(dd1,3); %mean over subjects
if isfield(S,'name')
    name = S.name;
else
    name = 'TG_Average');
end

if S.stat == 1 %
    %functions for Dimitrios' codes (decoding and statistics)
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding')
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab')
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab/external');
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab/external/libsvm-3.21')
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/scilearnlab/external/libsvm-3.21/matlab')
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/SpatialLocation')
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/SpatialLocation/plotchannel')
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/permutationlab')
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/permutationlab/private1')
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/MITDecodingToServer/Decoding/permutationlab/developer')
    %LBPD_startup_D
    pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions
    addpath(pathl);
    LBPD_startup_D(pathl);
    
    %computing statistics and plotting results
    stat = pl_permtest(dat,'alpha',0.05); %actual function
    %         time_selj = time_sel(1:776);
    %         figure; imagesc(time_selj,time_selj,stat.statmap); xlabel('Train time (s)'); ylabel('Test time (s)'); set(gca,'YDir','normal');
    %         colorbar
    %         set(gcf,'Color','w')
    %         caxis([-5 10])
    
    %testing contrast results (correcting for multiple comparison) by performing Monte Carlo simulations
    P2 = zeros(size(stat.statmap,1),size(stat.statmap,2));
    statt = stat.statmap;
    P2(statt>S.t_thresh) = 1; %threshold t-val = 2.6 corresponding to p-val < 0.01 (obtained by dividing 0.05 by the 4 comparisons employed here)
    thresh = 0;
    permut = 1000;
    threshMC = 0.001;
    perm_max = 1;
    t1 = S.time; t2 = t1;
    
    [ OUT ] = twoD_MCS_LBPD_D( P2, thresh, permut, threshMC, perm_max, t1 , t2 );
    
    PDn = cell2table(OUT); %table
    writetable(PDn,[S.output '/' name '.xlsx'],'Sheet',1); %printing excel file
    save([S.output '/' name '.mat'],'ddTGm','stat','OUT'); %saving averaged TG matrix, statistics and significant clusters
else %otherwise saving only average TG matrix 
    save([S.output '/' name '.mat'],'ddTGm');
end


end



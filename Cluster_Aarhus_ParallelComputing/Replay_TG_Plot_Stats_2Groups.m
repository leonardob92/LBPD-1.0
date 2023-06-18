function O = Replay_TG_Plot_Stats_2Groups(S)
O = []; 

%It computes plots for temporal generalisation (TG) decoding for two independent groups.
%It is designed to work on the Aarhus Cluster and it solves the OUT OF MEMORY issue.


% INPUT:

% -S structure with folloring fields:
%   -list_TG:                            directory with paths of files data. Typically obtained using dir() function
%   -time:                               time (in seconds)
%   -groups:                             cell array with two (or more) groups, with regards to list_TG. E.g. {[2,3,5],[4,6,7]}
%   -output:                             output directory and name for the output file
%   -S.stats:                            vector with groups to be contrasted in an independent-sample t-test (e.g. [1 2]).
%                                        Leave it empty [] for not running any statistics.



% OUTPUT:
% -average over participants for later plotting (saved on disk)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Oxford, UK, 18/06/2023


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%extracting inputs
list_TG = S.list_TG; %list of subjects
time = S.time;
groups = S.groups;
load([list_TG(1).folder '/' list_TG(1).name]); %load TG file
SS = size(d.d);

%actual computation
Dmean = zeros(SS(1),SS(2),length(groups)); %averaged file over groups of participants to be saved on disk later
for gg = 1:length(groups) %over groups
    dumg = groups{gg}; %vector with subjects in group gg
    D = zeros(SS(1),SS(2),length(dumg));
    for ii = 1:length(dumg) %over subjects
        load([list_TG(dumg(ii)).folder '/' list_TG(dumg(ii)).name]); %loading TG matrix
        D(:,:,ii) = d.d; %storing it
        disp(ii)
    end
    Dmean(:,:,gg) = mean(D,3); %average over subjects
end
clear D

if ~isempty(S.stats)
    %group 1
    d1group = groups{S.stats(1)};
    D1 = zeros(SS(1),SS(2),length(d1group));
    for ii = 1:length(d1group) %over subjects
        load([list_TG(d1group(ii)).folder '/' list_TG(d1group(ii)).name]); %loading TG matrix
        D1(:,:,ii) = d.d; %storing it
        disp(ii)
    end
    %group 2
    d1group = groups{S.stats(2)};
    D2 = zeros(SS(1),SS(2),length(d1group));
    for ii = 1:length(d1group) %over subjects
        load([list_TG(d1group(ii)).folder '/' list_TG(d1group(ii)).name]); %loading TG matrix
        D2(:,:,ii) = d.d; %storing it
        disp(ii)
    end
    %computing t-tests
    P = zeros(SS(1),SS(2));
    T = zeros(SS(1),SS(2));
    for ii = 1:size(P,1) %over time
        for jj = 1:size(P,1) %over time
            [~,p,~,stats] = ttest2(squeeze(D1(ii,jj,:)),squeeze(D2(ii,jj,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj) = p;
            T(ii,jj) = stats.tstat;
        end
        disp(ii)
    end
    save([S.output '.mat'],'Dmean','time','P','T'); %saving TG matrices averaged over the different groups
else
    save([S.output '.mat'],'Dmean','time'); %saving TG matrices averaged over the different groups
end

end



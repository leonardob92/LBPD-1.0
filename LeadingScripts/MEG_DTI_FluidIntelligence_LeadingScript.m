%% On the brain networks organization of individuals with high versus average fluid intelligence: a combined DTI and MEG study


%%% MEG - RESTING STATE - FUNCTIONAL CONNECTIVITY %%%


%BEFORE PROCEEDING, PLEASE NOTE:

%As follows, every analysis that we made has been reported to ensure full disclosure.
%Please, note that in this Leading script, I use built-in functions,
%functions from well-known toolboxes (e.g. OSL, SPM, FieldTrip) and
%in-house-built functions, which are reported together with this script in the collection named LBPD.
%Data can be provided according to the Danish regulation, upon reasonable request.
%If interested, please contact Leonardo Bonetti, leonardo.bonetti@clin.au.dk
%More information is provided in the ReadMe.mat file that I strongly advise to read.


%Leonardo Bonetti, Silvia Elisabetta Portis Bruzzone
%leonardo.bonetti@psych.ox.ac.uk
%silviaepbruzzone@clin.au.dk


%% DTI Portis pt. II: functional connectivity on MEG resting state (MEG_fcrs)

%1) Add LBPD functions
% starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl); % add to path LBPD and other OSL functions

%% Compute functional connectivity on MEG resting state in AAL space

%%% OBS!! TO BE RUN ONLY THE FIRST TIME %%%

%Make list with the frequencies folders
freq_list=dir('/scratch7/MINDLAB2017_MEG-LearningBach/Francesco/*rest10.oat');
order=[1:2:90 90:-2:2];
for ii=1:length(freq_list) %loop over folders containing different frequencies
    concat_list=dir([freq_list(ii).folder '/' freq_list(ii).name  '/concat*.mat']);
    path_fcrs= (['/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/MEG_fcrs/' freq_list(ii).name(20:end-11)]);
    mkdir(path_fcrs);
    for gg=1:length(concat_list) %loop over subjects within each folder
        D=spm_eeg_load([concat_list(gg).folder '/' concat_list(gg).name]); %load beamformed MEG files: 90 AAL regions x n. timepoints
        corM=corr(abs(hilbert(D(:,:)'))); %Create correlation matrix based on the envelopes computed with the hilbert transform - corr needs timepoints in the 1st dimension, that's why it was transposed
        corM=corM(order,order); %reshape corM to have it symmetrical (LEFT hemisphere=upper left triangle, RIGHT hem.=lower right triangle)
        save([path_fcrs '/' concat_list(gg).name(17:18) '.mat'], 'corM'); %save correlation matrix for each subject
        disp(['folder n. ' num2str(ii) ' subj n. ' num2str(gg)])
    end
end

%% PARCELLATION: compute centroids and prepare brain template figure

%%% OBS!! TO BE RUN ONLY THE FIRST TIME %%%

%this is for giving input such as parcellation, ROI names and template (otherwise it can also find by itself the proper template, in this case 8mm)
parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'; %this and the following files can be found in the provided codes folder
ROIsfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_ROIs.txt';
templatefile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz';
p = parcellation(parcelfile,ROIsfile,templatefile); %given required inputs, creating a parcellation object for SPM MEEG objects
p = p.remove_parcels(91:116);
%whole-brain MNI coordinates (2mm)
% [ WB_mni, ~ ] = osl_mnimask2mnicoords('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_T1_2mm_brain.nii.gz');
%getting centroids of all brain areas (2mm)
MNI_RC = zeros(90,3);
for ii = 1:2:90 %over AAL ROIs
    %MNI space AAL (LRLRLR)
    [ ROI_mni, ~ ] = osl_mnimask2mnicoords(['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/AAL_2mm_90ROIs/' num2str(ii) '.nii.gz']);
    MNI_RC(ii,:) = mean(ROI_mni,1); %left hemisphere
    MNI_RC(ii+1,1) = -MNI_RC(ii,1); %right hemisphere (x)
    MNI_RC(ii+1,2:3) = MNI_RC(ii,2:3); %right hemisphere (y,z)
    disp(ii)
end
order=[1:2:90 90:-2:2];
MNI_RC = MNI_RC(order,:); %order LLLRRR
save('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/MNI_RC.mat','MNI_RC')

%plotting brain template
p.plot_network(ones(90))
colorbar off
saveas(gcf,'/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/jes.fig')

%% GRAPH THEORETICAL MEASURES AND PLOTTING SOLUTIONS - BRAIN FUNCTIONAL CONNECTIVITY FROM MEG RESTING STATE DATA 

clear
label_cogmeas = 4; %1=PR, 2=WM, 3=SP, 4=PR+WM+SP, 5=BDI
label_freq = 5; %1=0.1-2Hz, 2=2-8Hz, 3=8-12Hz, 4=12-32Hz, 5=32-74Hz
label_test = 1; %1=Degree, 2=modularity intra and extra connectivity (connector and provincial hubs); 0 = no
brain_plot = 0; %1 = brain plot for label_test (degree and connector/provincial hubs); 0 = no
violins = 0; %1 = violins plot for label_test (degree and connector/provincial hubs); 0 = no
% label_ttest = 0; %1=t-test, 0=no
label_additiona_graph_measures = 0; %additional measures: 1 = char path length; 2 = glob eff; 3 = nod eff; 4 = loc eff; 5 = modularity
label_density = 0; % 1 = density; 0 = no (MAYBE DO SIMPLY SUM OF ALL CONNECTIONS)
perc=10; %percentage of values we want to remove for thresholding (for density)
label_modularity_brainplot = 0; %1 to plot modularity in the brain
label_intra = 1; %1 for intra subnetwork connectivity; 0 for inter subnetwork connectivity
label_schem_mod = 0; %1 for schemaball plotting of modularity; 0 for not
mod_perm_label = 0; %1 for MCS (1000 permutations) on modularity (establishing whether the brain is mode "modulable" than a random graph with its elements; 0 for not
label_FC_matrix_plot = 0; %1 to plot FC connectivity matrix and schemaball; 0 = no

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/DTI_Portis'); %adding path to function
close all
if label_cogmeas<5
    shiet=2;
else
    %loading xlsx file with all behavioral data you may be interested in for MINLABD2017-MEG_LearningBach
    shiet = 4; %select the excel sheet (1 = background; 2 = WAIS-IV; 3 = MET; 4 = BDI; 5 = GOLDSMITH; 6 = MEG behavioral task)
end
[~,~,raw] = xlsread('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/tbss_bis/BehavioralMeasuresLearningBach.xlsx',shiet);
%getting index of participants whose FA was actually used
if label_freq==1
    MEG_subjs =dir('/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/MEG_fcrs/Hz01_2/*.mat');
elseif label_freq==2
    MEG_subjs =dir('/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/MEG_fcrs/Hz2_8/*.mat');
elseif label_freq==3
    MEG_subjs =dir('/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/MEG_fcrs/Hz8_12/*.mat');
elseif label_freq==4
    MEG_subjs =dir('/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/MEG_fcrs/Hz12_32/*.mat');
elseif label_freq==5
    MEG_subjs =dir('/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/MEG_fcrs/Hz32_74/*.mat');
end
% WM-PR-SP and MET settings (for separate behavioral measures)
if label_cogmeas==1
    %building a design matrix for a two-sample unpaired t-test
    behav_ind = 19; %selecting the behavioral index (check the raw variable to get which index you need; e.g. WM = 24)
elseif label_cogmeas==2
    behav_ind=24;
elseif label_cogmeas==3
    behav_ind=29;
elseif label_cogmeas==5
    behav_ind=2;
end
if label_cogmeas ~= 4
    desmat = zeros(length(MEG_subjs),2);
    if label_cogmeas < 4
        for ii = 1:size(desmat)
            if ~isnan(raw{str2double(MEG_subjs(ii).name(1:2)) + 1,behav_ind}) %if nan just ignore and leave 0 (that should not take that participant into account in any of the groups
                desmat(ii,1) = raw{str2double(MEG_subjs(ii).name(1:2)) + 1,behav_ind}; %24 = WM
                if raw{str2double(MEG_subjs(ii).name(1:2)) + 1,behav_ind} == 0 %getting 1s for the second column if the value for the first column is 1
                    desmat(ii,2) = 1;
                end
            end
        end
    elseif label_cogmeas == 5
        % BDI settings
        for ii = 1:size(desmat)
            if ~isnan(raw{str2double(MEG_subjs(ii).name(1:2)) + 1,behav_ind})
                if (raw{str2double(MEG_subjs(ii).name(1:2)) + 1,behav_ind}<=4)
                    desmat(ii,1) = 1; %1 in the 1st column for BDI<=2 (NO tendency to depression)
                    desmat(ii,2) = 0;
                elseif (raw{str2double(MEG_subjs(ii).name(1:2)) + 1,behav_ind}>=9)   %getting 1s for the second column if the value for the first column is higher than 7
                    desmat(ii,2) = 1; % 1 in the 2nd column for BDI>7 (tendency to depression)
                    desmat(ii,1) = 0;
                end
            end
        end
    end
else
    %Defining different groups of subjects based on the average values of cognitive measures
    behav_ind = [17 22 27]; %select the behavioral index of the measures (check the raw variable to get which index you need; e.g. PR,WM,SP=[17 22 27])
    behav_mean = zeros(length(MEG_subjs),1); %Initialize vector where you will store the mean values across the behavioural measures of interest
    %Compute mean across conditions, for every subject
    for ii = 1:size(behav_mean)
        if ~isnan(raw{str2double(MEG_subjs(ii).name(1:2)) + 1,behav_ind(:,1)}) %if nan just ignore and leave 0 (that should not take that participant into account in any of the groups
            behav_mean(ii,1)=(raw{str2double(MEG_subjs(ii).name(1:2)) + 1,behav_ind(:,1)}+raw{str2double(MEG_subjs(ii).name(1:2)) + 1,behav_ind(:,2)}+raw{str2double(MEG_subjs(ii).name(1:2)) + 1,behav_ind(:,3)})/length(behav_ind); %Compute the average across the three cognitive measures
        end
    end
    all_behav_mean=mean(behav_mean(behav_mean~=0)); %Get mean value of cognitive measures across subjects
    desmat = zeros(length(MEG_subjs),2);
    for mm=1:size(desmat)
        if behav_mean(mm,1)~=0
            if behav_mean(mm,1)>all_behav_mean
                desmat(mm,1)= 1; %get 1 in desmat if the value is above the average
            elseif behav_mean(mm,1)<all_behav_mean
                desmat(mm,1) = 0; %get 0 in desmat if the value is below the average
                if behav_mean(mm,1)<all_behav_mean
                    desmat(mm,2)= 1; %get 1 in desmat if the value is below the average
                elseif behav_mean(ii,1)>all_behav_mean
                    desmat(mm,2) = 0; %get 0 in desmat if the value is above the average
                end
            end
        end
    end
end
% Graph Theory - Degree and degree-based hubs, modularity
% label_cogmeas='BDI';
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/DTI_Portis/GraphTheory/BCT/2019_03_03_BCT'); %add BCT functions to path (FFM)
louv_l = 0; %1 for Louvain; 0 for modularity
COL2 = [0.8 0 0; 0 0 0.8; 1 0 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0]; %vector with colors
clear datt datt2
% list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Graph_Theory_AAL90Matrix/Subj*'); %list of subjects
list = MEG_subjs;
%initialising matrices for several measures of graph theory
degree = zeros(length(list),90);
hubs=zeros(length(list),90);
provincial=zeros(length(list),90);
connector=zeros(length(list),90);
optstruct= zeros(90,length(list));
maxmod=zeros(length(list),1);
permnum= 1000;%num permutations
nodesum=zeros(length(list),90);
Mdiag = diag(1:90);
idx = find(Mdiag~=0);
CPL = zeros(length(list),1);
GE = zeros(length(list),1);
NE = zeros(length(list),90);
LE = zeros(length(list),1);
MOD = zeros(length(list),1);
%loading AAL labels
if ~exist('lab','var')
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat');
    order = [1:2:90 90:-2:2]; %to reshape the matrix symmetricaly
    lab2 = lab(order,:);
end
%actual computation
MM = zeros(90,90,length(list));
for ii = 1:length(list) %over subjects
    %loading matrix with strenght of connections (AAL, same ROIs used for distance, of course..)
    load([list(ii).folder '/' list(ii).name]); %Load matrix (M3) containing the strength of the connections across regions of interest
    M3=corM;
    MM(:,:,ii) = M3;
    degree(ii,:) = sum(M3); %compute degree of connectivity (weighted)
    %according to previous DTI-graph theory papers
    meandegreen=mean(degree(ii,:))/mean(degree(ii,:)); %mean degree, normalized = 1 (done just for the sake of formality)
    degreestd=std(degree(ii,:)./mean(degree(ii,:))); %std of degree, normalized
    degreen=degree(ii,:)./mean(degree(ii,:)); %normalized degree
    hubs(ii,:)=double(degreen>(meandegreen+degreestd));%compute degree-based hubs
    if louv_l == 1
        optstructP= zeros(90,permnum);
        maxmodP=zeros(permnum,1);
        for pp=1:permnum
            %                 [optstructP(:,pp), maxmodP(pp,1)] = modularity_und(M3);
            [optstructP(:,pp), maxmodP(pp,1)] = community_louvain(M3); %Compute modularity with Louvain algorithm
            %         disp(pp)
        end
        [maxval, maxind]=max(maxmodP); %getting maximum value showing the best modularity over the 1000 permutations
        maxmod(ii,1) =maxval; %storing such value for each subject
        optstruct(:,ii)=optstructP(:,maxind); %storing modularity corresponding to that value
    else
        
        optstruct(:,ii) = modularity_und(M3); %Newmann modularity (deterministic)
    end
    %strength of subnetworks as outputted by Louvain modularity
    for hh=1:90 %loop over the nodes (initial AAL areas)
        %       nodesum(:,hh)=degree(ii,hh)+degree(:,hh); %sum every connection with all the other connections, for each area (areas=columns)
        subnethh=optstruct(hh,ii); %Subnetwork to which the node hh belongs
        subidx=find(optstruct(:,ii)==subnethh); %Index of the nodes within the subnetwork
        nodesum(ii,hh)=(sum(M3(hh,subidx)))./sum(M3(hh,:)); %sum of the connections between node (hh) and all the
        %other nodes of its community, divided by all connections hh has with all the other nodes
    end
    %provincial and connector hubs based on modularity
    meanode=mean(nodesum(ii,:)); %Compute mean of the hubs (provincial and connectors), for each subject
    stdnode=std(nodesum(ii,:)); %Compute std of the hubs (provincial and connectors)
    provincial(ii,:)=double(nodesum(ii,:)>(meanode+stdnode));%compute provincial hubs (=values above mean value + std)
    connector(ii,:)=double(nodesum(ii,:)<(meanode-stdnode));%compute connector hubs (=values below mean value - std)
    
    if label_additiona_graph_measures == 1
        mM3=M3/(max(max(M3))); %Divide every value of M3 by the greatest value within M3 -> get weighted values between 0 and 1
        iM3=1-mM3; %create the inverse of mM3 (strongest connections = shortest paths)
        [distM, edges]=distance_wei(iM3); %Compute distance matrix containing the shortest path length between each pair of nodes
        distM(idx) = NaN; %assigning NaNs to diagonal
        CPL(ii,1) = nanmean(nanmean(distM)); %characteristic path length
        GE(ii,1) = nanmean(nanmean(1./distM)); %global efficiency
        NE(ii,:) = nanmean(1./distM);%nodal efficiency
        LE(ii,1) = mean(NE(ii,:)); %local efficiency
        [~,modd] = modularity_und(M3); %Newmann modularity (deterministic)
        MOD(ii,1) = modd;
    end
    disp(ii)
end
%plotting
if label_additiona_graph_measures == 1 %CPL
    highPR = CPL(desmat(:,1)==1,1); %characteristic path length of subjects with high PR
    lowPR = CPL(desmat(:,2)==1,1); %characteristic path length of subjects with low PR
    figure
    scatter(1:length(highPR), highPR, 'o', 'r');
    hold on
    scatter(1:length(lowPR), lowPR, 'o', 'b');
    [~,p,~,stats] = ttest2(highPR,lowPR);
    title(['CPL - pval ' num2str(p) ' - tval ' num2str(stats.tstat)])
    % Global efficiency
    highPR = GE(desmat(:,1)==1,1); %Global efficiency of subjects with high PR
    lowPR = GE(desmat(:,2)==1,1); %Global efficiency of subjects with low PR
    figure
    scatter(1:length(highPR), highPR, 'o', 'r');
    hold on
    scatter(1:length(lowPR), lowPR, 'o', 'b');
    [~,p,~,stats] = ttest2(highPR,lowPR);
    title(['GE - pval ' num2str(p) ' - tval ' num2str(stats.tstat)])
    %Nodal efficiency
    highPR = NE(desmat(:,1)==1,:); %nodal efficiency of subjects with high PR
    lowPR = NE(desmat(:,2)==1,:); %nodal efficiency of subjects with low PR
    diff2=mean(highPR,1)-mean(lowPR,1); %Difference between highPR and lowPR
    figure
    scatter(1:90, diff2);
    [~,p,~,stats] = ttest2(highPR,lowPR);
    disp('NE')
    lab2(find(p<0.05),:)
    %     title(['NE - pval ' num2str(p)])
    title('NE')
    % Local efficiency
    highPR = LE(desmat(:,1)==1,1); %Local efficiency of subjects with high PR
    lowPR = LE(desmat(:,2)==1,1); %Local efficiency of subjects with low PR
    figure
    scatter(1:length(highPR), highPR, 'o', 'r');
    hold on
    scatter(1:length(lowPR), lowPR, 'o', 'b');
    [~,p,~,stats] = ttest2(highPR,lowPR);
    title(['LE - pval ' num2str(p) ' - tval ' num2str(stats.tstat)])
    %Modularity
    highPR = MOD(desmat(:,1)==1,1); %characteristic path length of subjects with high PR
    lowPR = MOD(desmat(:,2)==1,1); %characteristic path length of subjects with low PR
    figure
    scatter(1:length(highPR), highPR, 'o', 'r');
    hold on
    scatter(1:length(lowPR), lowPR, 'o', 'b');
    [~,p,~,stats] = ttest2(highPR,lowPR);
    title(['mod - pval ' num2str(p) ' - tval ' num2str(stats.tstat)])
end
if label_density == 1
    Mones=ones(90);
    Mcat=[]; %initialize matrix which will contain the concatenated matrices of every subject
    for ii = 1:length(list) %over subjects
        %loading matrix with strenght of connections (AAL, same ROIs used for distance, of course..)
        load([list(ii).folder '/' list(ii).name]); %Load matrix (M3) containing the strength of the connections across regions of interest
        M3=corM;
        Mtriu=M3(find(triu(Mones,1)~=0));
        Mcat=cat(1,Mcat, Mtriu);
    end
    thresh1=sort(Mcat); %sort Mcat
    indthresh1=round(length(thresh1)/100*perc); %find indeces corresponding to 1% of Mcat
    if indthresh1 == 0
        indthresh1 = 1;
    end
    thresh1val=thresh1(indthresh1); %values corresponding to 1% of Mcat
    density = zeros(length(list),1);
    for ii = 1:length(list) %over subjects
        %loading matrix with strenght of connections (AAL, same ROIs used for distance, of course..)
        load([list(ii).folder '/' list(ii).name]); %Load matrix (M3) containing the strength of the connections across regions of interest
        M3=corM;
        M4=zeros(90);
        M4(M3>thresh1val)=M3(M3>thresh1val);
        density(ii,1) = density_und(M4); %compute degree of connectivity
    end
    %t-test
    highPR = density(desmat(:,1)==1,1); %density of subjects with high PR
    lowPR = density(desmat(:,2)==1,1); %density of subjects with low PR
    figure
    scatter(1:length(highPR), highPR, 'o', 'r');
    hold on
    scatter(1:length(lowPR), lowPR, 'o', 'b');
    [~,p,~,stats] = ttest2(highPR,lowPR);
    title(['density - pval ' num2str(p) ' - tval ' num2str(stats.tstat)])
end
%STAT TESTS ON GRAPH THEORY MEASURES
if label_test == 1 %computing MCS for degree and connector/provincial hubs
    % PERMUTATION TEST (DEGREE)
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/DTI_Portis'); %adding path to function
    S = [];
    S.thr = [];
    S.data{1,1} = degree(desmat(:,1)==1,:); %degrees of subjects with high PR
    S.data{1,2} = degree(desmat(:,2)==1,:); %degrees of subjects with low PR
    S.thr = std(median(S.data{1,1},1) - median(S.data{1,2},1)); %keeping only extreme data (leave empty [] for using all data)
    S.permtype = 1; %choose the permutation type: 1 for permuting the degrees of each couple of regions (high vs low PR) independently; 2 for permuting all degrees together
    S.permnum = 10000;
    disp('degree results')
    [p_val_pos, p_val_neg, posidx, negidx, diff1] = DTI_GT_MCS(S);
    disp(['p_val_pos: ' num2str(p_val_pos) ' - p_val_neg: ' num2str(p_val_neg)])
    %loading AAL labels
    if ~exist('lab','var')
        load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat');
        order = [1:2:90 90:-2:2]; %to reshape the matrix symmetricaly
        lab2 = lab(order,:);
    end
    disp('pos nodes')
    ROispos = lab2(posidx',:)
    valuespos = diff1(posidx)';
    disp('neg nodes')
    ROisneg = lab2(negidx',:)
    valuesneg = diff1(negidx)';
    disp('ST for degree')
    S.thr
    save(['Degree_freq' num2str(label_freq) '.mat'],'ROispos','ROisneg','p_val_pos','p_val_neg','valuespos','valuesneg')
    if violins == 1
        %violins
        po1 = (median(S.data{1,1},1))'; %median over subjects (group1)
        po2 = (median(S.data{1,2},1))'; %median over subjects (group2)
        datt{1,1} = po1; datt{1,2} = po2;
        cl = COL2(1:2, :); %same..
        figure
        h = rm_raincloud2(datt',cl)
        grid minor
        set(gcf,'color','w')
        set(gcf,'Position',[200,200,400,550])
        xlim([0 27])
        saveas(gcf,['Degree_freq' num2str(label_freq) '.svg'])
        %difference (violin plot)
        poo = (po1 - po2); %+ (mean(po1)+mean(po2))/2;
        datt2{1,1} = poo;
        figure
        h = rm_raincloud2(datt2',COL2(3, :))
        grid minor
        set(gcf,'color','w')
        set(gcf,'Position',[200,200,400,550])
        xlim([-6 6])
        saveas(gcf,['Degree_diff_freq' num2str(label_freq) '.svg'])
    end
    %     hold on
    %     jesui = xlim();
    %     x = linspace(jesui(1),jesui(2)); y = zeros(1,length(x));
    %     line(y,x)
    %     export_fig(['Degree_freq' num2str(label_freq) '.eps'])
    %PLOTTING DEGREE IN THE BRAIN (IN THE CENTROID LOCATIONS)
    %     if p_val_pos < 0.05 || p_val_neg < 0.05
    if brain_plot == 1
        clear clodum
        coldum(1) = 'r'; coldum(2) = 'b';
        if length(find(posidx==1)) < length(find(negidx==1)) %taking the side of the distribution where the difference was larger
            coldum(3) = 'b';
            posidx = negidx;
        else
            coldum(3) = 'r';
        end
        if isempty(find(posidx==1))
            posidx = negidx;
        end
        lo = (squeeze(median(S.data{1,1},1)) - squeeze(median(S.data{1,2},1)))'; %storing difference between median of the degree of the two groups
        S.data{1,3} = lo(posidx');
        S.data{1,1} = (median(S.data{1,1},1))'; %median over subjects (group1)
        S.data{1,2} = (median(S.data{1,2},1))'; %median over subjects (group2)
        rmin12 = min(cat(1,S.data{1},S.data{2})); %for later scaling..
        rmax12 = max(cat(1,S.data{1},S.data{2}));
        posidx = find(posidx==1)';
        for cc = 1:3 %group1 - group2 - difference of their medians
            limitt = [2,30]; %size of centroids for the 2 groups
            openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
            load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI_RC_AAL_2mm.mat')
            %scaling on the basis of connections strength
            dataa = abs(S.data{cc}); %to not have problems in case of negative distributions..)
            if cc == 3 %scaling difference  between groups
                rmin = min(dataa);
                rmax = max(dataa);
                %creating nifti image for WorkBench
                vector = zeros(90,1);
                vector(posidx) = dataa;
                create_AALnifti(vector,['FC_degree_freq_' num2str(label_freq) '_contr_' num2str(cc) '.nii'],1);
            else %otherwise scaling together group1 and group2
                rmin = rmin12; rmax = rmax12;
                %creating nifti image for WorkBench
                vector = dataa;
                create_AALnifti(vector,['FC_degree_freq_' num2str(label_freq) '_contr_' num2str(cc) '.nii'],1);
            end
            %scaling connection strengths with extreme values (limitt) requested by
            size_con = (dataa-rmin)./(rmax-rmin).*(limitt(2)-limitt(1)) + limitt(1);
            for gg = 1:length(dataa) %over AAL ROIs
                hold on
                if cc == 3
                    plot3(MNI_RC(posidx(gg),1), MNI_RC(posidx(gg),2), MNI_RC(posidx(gg),3), '.', 'Color', coldum(cc), 'MarkerSize', size_con(gg)); %centroid of ROI ii
                else
                    plot3(MNI_RC(gg,1), MNI_RC(gg,2), MNI_RC(gg,3), '.', 'Color', coldum(cc), 'MarkerSize', size_con(gg)); %centroid of ROI ii
                end
            end
            rotate3d on; axis off; axis vis3d; axis equal
            set(gcf,'color','w')
            title(['degree - contrast ' num2str(cc)])
            view([0 90])
%             saveas(gcf,['Dots_top_FC_degree_freq_' num2str(label_freq) '_contr_' num2str(cc) '.svg'])
            view([-90 0])
%             saveas(gcf,['Dots_left_FC_degree_freq_' num2str(label_freq) '_contr_' num2str(cc) '.svg'])
        end
    else
        %             disp('no significant degree..')
    end
    %     end
elseif label_test == 2
    %PERMUTATION TEST (NODESUM - INTRASUBNETWORK CONNECTIVITY DIVIDED BY ALL CONNECTIVITY FOR EACH NODE AT A TIME)
    S = [];
    S.data{1,1} = nodesum(desmat(:,1)==1,:); %degrees of subjects with high PR
    S.data{1,2} = nodesum(desmat(:,2)==1,:); %degrees of subjects with low PR
    S.thr = []; %keeping only extreme data (leave empty [] for using all data)
    S.thr = std(median(S.data{1,1},1) - median(S.data{1,2},1)); %keeping only extreme data (leave empty [] for using all data)
    % S.thr = std(mean(S.data{1,1},1) - mean(S.data{1,2},1)); %keeping only extreme data (leave empty [] for using all data)
    S.permtype = 1; %choose the permutation type: 1 for permuting the degrees of each couple of regions (high vs low PR) independently; 2 for permuting all degrees together
    S.permnum = 10000;
    disp('connector/provincial results')
    [p_val_pos, p_val_neg, posidx2, negidx2, diff2] = DTI_GT_MCS(S);
    disp(['p_val_pos: ' num2str(p_val_pos) ' - p_val_neg: ' num2str(p_val_neg)])
    %loading AAL labels
    if ~exist('lab','var')
        load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat');
        order = [1:2:90 90:-2:2]; %to reshape the matrix symmetricaly
        lab2 = lab(order,:);
    end
    disp('pos nodes')
    ROispos = lab2(posidx2',:)
    valuespos = diff2(posidx2)';
    disp('neg nodes')
    ROisneg = lab2(negidx2',:)
    valuesneg = diff2(negidx2)';
    disp('ST for connectors')
    S.thr
    save(['Connectors_freq' num2str(label_freq) '.mat'],'ROispos','ROisneg','p_val_pos','p_val_neg','valuespos','valuesneg')
    if violins == 1
        %violins
        po1 = (median(S.data{1,1},1))'; %median over subjects (group1)
        po2 = (median(S.data{1,2},1))'; %median over subjects (group2)
        datt{1,1} = po1; datt{1,2} = po2;
        cl = COL2(1:2, :); %same..
        figure
        h = rm_raincloud2(datt',cl)
        grid minor
        set(gcf,'color','w')
        set(gcf,'Position',[200,200,400,550])
        xlim([0.3 0.8])
        saveas(gcf,['Connector_freq' num2str(label_freq) '.svg'])
        %difference (violin plot)
        poo = (po1 - po2); %+ (mean(po1)+mean(po2))/2;
        datt2{1,1} = poo;
        figure
        h = rm_raincloud2(datt2',COL2(3, :))
        grid minor
        set(gcf,'color','w')
        set(gcf,'Position',[200,200,400,550])
        xlim([-0.15 0.15])
        saveas(gcf,['Connector_diff_freq' num2str(label_freq) '.svg'])
        %     export_fig(['Connector_freq' num2str(label_freq) '.eps'])
    end
    %PLOTTING CONNECTOR HUBS IN THE BRAIN (IN THE CENTROID LOCATIONS)
    %     if p_val_pos < 0.05 || p_val_neg < 0.05
    if brain_plot == 1
        clear clodum
        coldum(1) = 'r'; coldum(2) = 'b';
        if length(find(negidx2==1)) < length(find(posidx2==1)) %taking the side of the distribution where the difference was larger
            negidx2 = posidx2;
            coldum(3) = 'b';
        else
            coldum(3) = 'r';
        end
        lo = (squeeze(median(S.data{1,1},1)) - squeeze(median(S.data{1,2},1)))'; %storing difference between median of the degree of the two groups
        S.data{1,3} = lo(negidx2');
        S.data{1,1} = (median(S.data{1,1},1))'; %median over subjects (group1)
        S.data{1,2} = (median(S.data{1,2},1))'; %median over subjects (group2)
        rmin12 = min(cat(1,S.data{1},S.data{2})); %for later scaling..
        rmax12 = max(cat(1,S.data{1},S.data{2}));
        negidx2 = find(negidx2==1)';
        for cc = 1:3 %group1 - group2 - difference of their medians
            limitt = [2,30]; %size of centroids for the 2 groups
            openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
            load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI_RC_AAL_2mm.mat')
            %scaling on the basis of connections strength
            dataa = abs(S.data{cc});
            if cc == 3 %scaling difference between groups
                %             dataa = abs(dataa); %since we are interested in connector hubs (which are negative)
                rmin = min(dataa);
                rmax = max(dataa);
                %scaling connection strengths with extreme values (limitt) requested by
                size_con = (dataa-rmin)./(rmax-rmin).*(limitt(2)-limitt(1)) + limitt(1);
            else %otherwise scaling together group1 and group2
                rmin = rmin12; rmax = rmax12;
                %scaling connection strengths with extreme values (limitt) requested by
                size_con = (dataa-rmin)./(rmax-rmin).*(limitt(2)-limitt(1)) + limitt(1);
                size_con2 = limitt(1) + (limitt(2) - size_con);
            end
            if cc == 3
                %creating nifti image for WorkBench
                vector = zeros(90,1);
                vector(negidx2) = size_con;
                create_AALnifti(vector,['FC_connectorhubs_freq_' num2str(label_freq) '_contr_' num2str(cc) '.nii'],1);
            else
                %creating nifti image for WorkBench
                vector = size_con2;
                create_AALnifti(vector,['FC_connectorhubs_freq_' num2str(label_freq) '_contr_' num2str(cc) '.nii'],1);
            end
            for gg = 1:length(dataa) %over AAL ROIs
                hold on
                if cc == 3
                    plot3(MNI_RC(negidx2(gg),1), MNI_RC(negidx2(gg),2), MNI_RC(negidx2(gg),3), '.', 'Color', coldum(cc), 'MarkerSize', size_con(gg)); %centroid of ROI ii
                else
                    plot3(MNI_RC(gg,1), MNI_RC(gg,2), MNI_RC(gg,3), '.', 'Color', coldum(cc), 'MarkerSize', size_con2(gg)); %centroid of ROI ii
                end
            end
            rotate3d on; axis off; axis vis3d; axis equal
            set(gcf,'color','w')
            title(['connector hubs - contrast ' num2str(cc)])
            view([0 90])
            saveas(gcf,['Dots_top_FC_connector_freq_' num2str(label_freq) '_contr_' num2str(cc) '.svg'])
            view([-90 0])
            saveas(gcf,['Dots_left_FC_connector_freq_' num2str(label_freq) '_contr_' num2str(cc) '.svg'])
        end
    else
        %             disp('no significant connector hubs..')
    end
end
% end
% if label_ttest==1
%     %t-tests
%     
%     P = zeros(90,1);
%     T = zeros(90,1);
%     for ii = 1:90 %over AAL brain regions
%         a = degree(desmat(:,1)==1,ii); %degrees of subjects with high PR
%         b = degree(desmat(:,2)==1,ii); %degrees of subjects with low PR
%         [~,p,~,stats] = ttest2(a,b);
%         P(ii,1) = p;
%         T(ii,1) = stats.tstat;
%     end
%     %loading AAL labels
%     if ~exist('lab','var')
%         load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat');
%         order = [1:2:90 90:-2:2]; %to reshape the matrix symmetricaly
%         lab2 = lab(order,:);
%     end
%     disp('pos nodes')
%     lab2((P<0.05 & T>0),:)
%     disp('neg nodes')
%     lab2((P<0.05 & T<0),:)
% end
if label_modularity_brainplot == 1
    %loading MNI coordinates of AAL 2-mm centroids
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/MNI_RC.mat')
    %getting mean over groups of subjects and then computing modularity
    MM1 = mean(MM(:,:,desmat(:,1)==1),3); %group 1
    [com1,modl1] = modularity_und(MM1); %Newmann modularity (deterministic)
    COMM = zeros(length(com1),2);
    COMM(:,1) = com1; %storing ROIs communities ID for group 1
    MM2 = mean(MM(:,:,desmat(:,2)==1),3); %group 2
    [com1,modl2] = modularity_und(MM2); %Newmann modularity (deterministic)
    COMM(:,2) = com1; %storing ROIs communities ID for group 2
    %permutations to establish whether the brain is more "modulable" than a random graph with its same elements
    if mod_perm_label == 1
        %group1
        [~,m] = size(MM1);
        permut = 1000;
        MODP1 = zeros(permut,1);
        disp('computing permutations ofr MCS on modularity - group 1')
        for pp = 1:permut
            f = find(triu(ones(m),1) == 1); %find indexes of elements of the upper triangle of the square matrix
            idx_dummy = randperm(length(f)); %create a permuted array from 1 to length(f)
            idx_2 = f(idx_dummy); %shuffle the indexes in f according to the order in idx_dummy
            r_dummy = zeros(m*m,1); %initialise vector
            r_dummy(f) = MM1(idx_2); %taking element in matrix data with index idx_2 (shuffled indexes of the upper triangle)
            %rebuild 2D matrix with actual data
            r3 = reshape(r_dummy,[m,m]);
            %ricreate the simmetric matrix for calculating more easily the sums
            r3 = (r3 + r3') - eye(size(r3,1)).*diag(r3);
            %calculating segregation of the randomly rearranged matrix
            [~,mod_glob_r3] = modularity_und(r3); %only global value of segregation
            MODP1(pp) = mod_glob_r3;
%             disp(pp)
        end
        pval1 = length(find(modl1<MODP1))/permut;
        disp(['modularity p-value group 1 = ' num2str(pval1)])
        %group1
        [~,m] = size(MM2);
        permut = 1000;
        MODP2 = zeros(permut,1);
        disp('computing permutations ofr MCS on modularity - group 2')
        for pp = 1:permut
            f = find(triu(ones(m),1) == 1); %find indexes of elements of the upper triangle of the square matrix
            idx_dummy = randperm(length(f)); %create a permuted array from 1 to length(f)
            idx_2 = f(idx_dummy); %shuffle the indexes in f according to the order in idx_dummy
            r_dummy = zeros(m*m,1); %initialise vector
            r_dummy(f) = MM2(idx_2); %taking element in matrix data with index idx_2 (shuffled indexes of the upper triangle)
            %rebuild 2D matrix with actual data
            r3 = reshape(r_dummy,[m,m]);
            %ricreate the simmetric matrix for calculating more easily the sums
            r3 = (r3 + r3') - eye(size(r3,1)).*diag(r3);
            %calculating segregation of the randomly rearranged matrix
            [~,mod_glob_r3] = modularity_und(r3); %only global value of segregation
            MODP2(pp) = mod_glob_r3;
%             disp(pp)
        end
        pval2 = length(find(modl2<MODP2))/permut;
        disp(['modularity p-value group 2 = ' num2str(pval2)])
    end
    %elaborated way to match the color of most similar communities between groups
    if max(COMM(:,2)) > max(COMM(:,1)) %working finding communities of group with less communities in group with more communities
        labcom1 = 1;
        labcom2 = 2;
    else
        labcom1 = 2;
        labcom2 = 1;
    end
    matchh = zeros(max(COMM(:,labcom1)),2);
    for jj = 1:max(COMM(:,labcom1))
        dum = COMM(find(COMM(:,labcom1)==jj),labcom2); %getting elements of group2 with indices of group1 == community jj
        a = unique(dum); %getting unique values
        PP = zeros(1,length(a));
        for pp = 1:length(a) %over communities in dum
            PP(pp) = length(find(dum==a(pp))); %finding occurrencies of unique values (communities)
        end
        [~,ipp] = max(PP); %best match between communities in group1 and group2
        matchh(jj,1) = jj; %storing community for group1
        matchh(jj,2) = a(ipp); %with best corresponding community of group2
    end
    a4 = 1:max(COMM(:,labcom1));
    if labcom1 == 2 %matching communities color here
        %color codes using letters
%         COL{1} = ['r','m','k','g','y','c']; %assigning colors to communities of group1
%         COL{2} = COL{1}(matchh(:,2)); %colors to community of group2 from match between the two groups
%         COL{2}  = cat(2,COL{2},COL{1}(find(sum(double(a4==matchh(:,2)),1)==0)));
        %color codes using numbers (triplets)
        COL{1} = [0 0.6 1; 1 0.4 0; 1 0 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0]; %vector with colors
        COL{2} = COL{1}(matchh(:,2),:); %colors to community of group2 from match between the two groups
        COL{2}  = cat(1,COL{2},COL{1}(find(sum(double(a4==matchh(:,2)),1)==0),:));
    else
        %color codes using letters
%         COL{2} = ['r','m','k','g','y','c']; %assigning colors to communities of group2
%         COL{1} = COL{2}(matchh(:,2)); %colors to community of group1 from match between the two groups
%         COL{1}  = cat(2,COL{1},COL{2}(find(sum(double(a4==matchh(:,2)),1)==0)));
        %color codes using numbers (triplets)
        COL{2} = [0 0.6 1; 1 0.4 0; 1 0 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0]; %vector with colors
        COL{1} = COL{2}(matchh(:,2),:); %colors to community of group2 from match between the two groups
        COL{1}  = cat(1,COL{1},COL{2}(find(sum(double(a4==matchh(:,2)),1)==0),:));
    end
    for ii = 1:2 %over experimental groups
        S = [];
        if label_intra == 1
            S.intra_con = 1; %1 for intra-connectivity plots
            S.inter_con = 0; %1 for inter-connectivity plots
        else
            S.intra_con = 0; %1 for intra-connectivity plots
            S.inter_con = 1; %1 for inter-connectivity plots
        end
        S.perc_intra = 30; %percentage of the connections that you want to plot in each community (intra)
        S.perc_inter = 3; %percentage of the connections that you want to plot between each couple of communities (inter)
        S.MNI_RC = MNI_RC;
        S.conn_matrix = mean(MM(:,:,desmat(:,ii)==1),3); %connectivity matrix
        S.communities = COMM(:,ii);
        S.colors_mod = COL{ii};
        S.col_dot = 'b';
        S.limits = [0.1 5];
        S.title = ['group ' num2str(ii)];
        %actual plotting function
        GT_modul_plot_LBPD(S)
        
        %printing optimal community structure in excel files (high IQ)
        PD = cell(91,max(COMM(:,1)));
        for ll = 1:max(COMM(:,1))
            PD{1,ll} = ['Module ' num2str(ll)];
            PD(2:size(lab2((COMM(:,1)==ll),:),1)+1,ll) = cellstr(lab2((COMM(:,1)==ll),:));
        end
        PDnh = cell2table(PD); %remove the possible empty cell
        writetable(PDnh,['Freq_' num2str(label_freq) '_High_IQ.xlsx'],'Sheet',ll) %printing excel file
        %printing optimal community structure in excel files (average IQ)
        PD = cell(91,max(COMM(:,2)));
        for ll = 1:max(COMM(:,2))
            PD{1,ll} = ['Module ' num2str(ll)];
            PD(2:size(lab2((COMM(:,2)==ll),:),1)+1,ll) = cellstr(lab2((COMM(:,2)==ll),:));
        end
        PDna = cell2table(PD); %remove the possible empty cell
        writetable(PDna,['Freq_' num2str(label_freq) '_Average_IQ.xlsx'],'Sheet',ll) %printing excel file
        %         if label_intra == 1
        %             view([0 90])
        %             saveas(gcf,['IntraMod_FC_group_' num2str(ii) '_top_freq_' num2str(label_freq) '.svg'])
        %             view([-90 0])
        %             saveas(gcf,['IntraMod_FC_group_' num2str(ii) '_left_freq_' num2str(label_freq) '.svg'])
        %         else
        %             view([0 90])
        %             saveas(gcf,['InterMod_FC_group_' num2str(ii) '_top_freq_' num2str(label_freq) '.svg'])
        %             view([-90 0])
        %             saveas(gcf,['InterMod_FC_group_' num2str(ii) '_left_freq_' num2str(label_freq) '.svg'])
        %         end
        if label_schem_mod == 1
            %testing schemaball with modularity
            figure
            order = [1:2:90 90:-2:2];
            load('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat'); %loading labels
            lab2 = lab(order,:);
            r = S.conn_matrix./max(max(S.conn_matrix)); %scaling matrix
            perc_intra = 30; %TO BE DEFINED
            perc_inter = 3;
            limitt = [0.01 3]; %TO BE DEFINED
            %actual function
            schemaball_modularity_LBPD(r, lab2, S.communities, COL{ii}, perc_intra, perc_inter, limitt)
            title(['group ' num2str(ii)])
            set(gcf,'color','w')
            export_fig(['FC_schemball_Freq_' num2str(label_freq) '_group_' num2str(ii) '.eps'])
        end
    end
%     %CENTROIDS OF SUBNETWORKS (FROM MODULARITY)
%     %%% PROBABLY RELEVANT - TO BE REPORTED, maybe in supplementary materials
%     centroids_l = 0;
%     if centroids_l == 1
%         MS = [30,10]; %size of centroids for the 2 groups
%         openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
%         load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI_RC_AAL_2mm.mat')
%         for gg = 1:2 %over experimental groups
%             for cc = 1:max(COMM(:,gg)) %over communities of group gg
%                 ind = find(COMM(:,gg)==cc); %getting indices of community cc of group gg
%                 centroid = mean(MNI_RC(ind,:),1); %MNI 3D coordinates of centroid of community cc of group gg
%                 hold on
%                 plot3(centroid(1), centroid(2), centroid(3), '.', 'Color', COL{gg}(cc,:), 'MarkerSize', MS(gg)); %centroid of ROI ii
%             end
%         end
%         rotate3d on; axis off; axis vis3d; axis equal
%         set(gcf,'color','w')
%         title(['centroids of communities'])
%         view([0 90])
%         saveas(gcf,['/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021/Modularity/Centroids_FC_AllSubjs_top_freq_' num2str(label_freq) '.svg'])
%         view([-90 0])
%         saveas(gcf,['/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021/Modularity/Centroids_FC_AllSubjs_left_freq_' num2str(label_freq) '.svg'])
%     end
    %displaying here for convenience of the user..
    if exist('pval1','var')
        disp(['modularity p-value group 1 = ' num2str(pval1)])
    end
    if exist('pval2','var')
        disp(['modularity p-value group 2 = ' num2str(pval2)])
    end
end
%plotting average FC matrix
if label_FC_matrix_plot == 1
    Mdiag = diag(1:90);
    idx = find(Mdiag~=0);
    MM2 = mean(MM,3);
    MM2(idx) = 0; %zeros to diagonal.. maybe not necessary, but let's go with this.. (before I assigned NaNs to diagonal..)
    figure
    imagesc(MM2)
    colorbar
%     clims = [0 1];
%     caxis(clims)
    caxis([-1 1])
    set(gcf,'color','w')
    %changing colormap of figures (red-bue with white for 0 values)
    x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
    colormap(bluewhitered_PD(0,x))
    export_fig(['/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021/Matrix_FC_freq_' num2str(label_freq) '_AllSubjs.eps'])
    %schemaball plotting    
    S = [];
    S.MATF = MM2; %connectivity matrix ROIs x ROIs (can be both binary or non-binary) (double)
    S.symmetric_l = 1; %1 if MATF is LLLRRR; 0 if MATF is LRLRLR
    S.outpath = '/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021'; %path where saving images (characters)
    S.p = []; %parcellation object (created by using OSL)
    S.thr_cortex = []; %  strongest % of connectivity to be plotted (it is usually a quite high number, e.g. 0.993). If empty [] default is 0.993.
    S.name_gif = []; %name of gif/video (character). Leave empty [] if you do not want to plot the connectivity in the brain (e.g. you want only the schemaball).
    S.frame_vect = []; %vector with frames you want. If empty [] default is [15,90], another good option could be [1,15,59,90,130]
    S.fr_spec = []; %cell array with characters specifing the name of the requested frames. E.g. [15,90] corresponds to fr_spec = {'Posterior_L_R';'Frontal_R_L'} [1,15,59,90] corresponds to fr_spec = {'Orig_Angle';'Posterior_L_R';'Hem_R';'Frontal_R_L';'Hem_L'}
    S.schball_l = 1; %1 for having the schemaball; 0 for not having it.
    S.extr = []; %minimum and maximum values to be used for normalizing the matrix to be submitted to schemaball function.
        %Leave empty [] for scaling the matrix on its max(abs value).
    S.lab_parc_ph = '/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat'; %path to labels. It must be a .mat file with a character array named 'lab' (characters).
    S.colbacksch = []; %background color for shcemaball. If empty [], default is 'w', white.
    S.name_sch = [];
    S.schemplth = 100; %percentage of the strongest connections to be plotted in schemaball (5% to be submitted as '5').
%     S.name_sch  = [];
    S.name_sch = ['SchBal_FC_freq_' num2str(label_freq) '_allSubjs']; %name for saved schemaball image. Leave empty [] for not saving it.
    %actual function
    FC_Brain_Schemaball_plotting_LBPD_D(S)
    %plotting connectivity in the brain (using the modularity plotting with all ROIs belonging to 1 community)
    %loading MNI coordinates of AAL 2-mm centroids
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/MNI_RC.mat')
    S = [];
    S.intra_con = 1; %1 for intra-connectivity plots
    S.inter_con = 0; %1 for inter-connectivity plots
    S.perc_intra = 10; %percentage of the connections that you want to plot in each community (intra)
    S.perc_inter = 0; %percentage of the connections that you want to plot between each couple of communities (inter)
    S.MNI_RC = MNI_RC;
    S.conn_matrix = MM2; %connectivity matrix
    S.communities = ones(90,1); %assigning all ROIs to just 1 community
    S.colors_mod = [0.7 0 0]; %; 1 0.4 0; 1 0 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0];
    S.col_dot = [0 0 0.5];
    S.limits = [0.1 5];
    S.title = ['group ' num2str(ii)];
    %actual plotting function
    GT_modul_plot_LBPD(S)
    view([0 90])
    saveas(gcf,['/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021/Full_FC_AllSubjs_top_freq_' num2str(label_freq) '.svg'])
    view([-90 0])
    saveas(gcf,['/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021/Full_FC_AllSubjs_left_freq_' num2str(label_freq) '.svg'])
end

%%
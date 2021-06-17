function stat = pl_permtestcluster2(data1, data2, varargin)
%
% Performs 2-sample cluster-size statistical inference on 'data1' and 'data2' for difference against 0. 
% For a difference against any mean value m, use input data - m.
%
%   
% stat = pl_permclustertest2(data1,data2) uses default values
%
% stat = pl_permclustertest2(data,'param1',val1,'param2',val2,...) specifies one or
% more of the following name/value pairs
%
%       Parameter           Value
%       'alpha'             A value between 0-1 (default 0.05), specifying the significance level.
%       'statistic'         'tstat' (default), 'mean'
%                           Defines the test statistic, either a t-statistic or the mean  
%       'vartype'           'equal' (default), 'unequal', determines if the t-statistic assumes equal or unequal variance 
%       'numpermutation'    An integer value, the number of permutations (default 1000)
%       'tail'              'right' (default),'left' 
%                           Determines the alternative hypothesis mean > 0, mean < 0
%       'clusteralpha'      A value between 0-1 (default 0.05), specifying the cluster defining threshold
%       'clusterstatistic'  'maxsize', 'maxsum', 'wcm', determines the
%                           cluster statistic which can be equal to the number of voxels (maxsize),
%                           the sum of voxel values (maxsum), 
%                           or a weighted cluster mass (wcm, see 'Hayasaka & Nichols (NeuroImage, 2004)
%       'theta'             A scalar used when 'wcm' is requested (default 0.5). Determines whether the 
%                           cluster statistic is dominated by the cluster size (theta -> 0)
%                           or the largest voxel intensity in the cluster (theta -> 1)
%       'verbose'           Display iteration progress  
%
% stat  .statmap                    Statistical map (variable1 x variable2 x ...)
%       .criticalmap                Critical map (binary map, 1 significant, 0 non-significant
%       .clustermassdist            Distribution of the maximum cluster statistic
%       .clusters                   Indices of significant clusters
%       .clustermass                Mass of significant clusters
%       .clusteralphamap            Statistic values at cluster defining threshold (e.g. when clusteralpha = 0.05)
%       .tail                       Same as input or default values (for bookkeeping)
%       .alpha                      Same
%       .statistic                  Same
%       .numpermutation             Same
%       .clusteralpha               Same
%       .clusterstatistic           Same


% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is currently under development, please do not share


%% parse inputs 

numpermutation      = pl_inputparser(varargin,'numpermutation',1000,@(x) isscalar(x) && x>0 && x == round(x) );
tail                = pl_inputparser(varargin,'tail','right',{'right','left'}); %no 'both' allowed for clusters
alpha               = pl_inputparser(varargin,'alpha',0.05,@(x) isscalar(x) && (x > 0) && (x < 1));
statistic           = pl_inputparser(varargin,'statistic','tstat',{'tstat','mean'});
vartype             = pl_inputparser(varargin,'vartype','equal',{'equal','unequal'});
verbose             = pl_inputparser(varargin,'verbose',false);
clusteralpha        = pl_inputparser(varargin,'clusteralpha',0.05,@(x) isscalar(x) && (x > 0) && (x < 1));
clusterstatistic    = pl_inputparser(varargin,'clusterstatistic','maxsize',{'maxsum','maxsize','wcm'});
theta               = pl_inputparser(varargin,'theta',0.5,@(x) isscalar(x) && (x > 0) && (x < 1));


%% convert input data to cell arrays (if not already)

if ~iscell(data1) %convert data to cell (faster)
    data1 = pl_mat2cell(data1);
end
if ~iscell(data2) %convert data to cell (faster)
    data2 = pl_mat2cell(data2);
end


%% assign statistic function

switch statistic
    case 'mean'
        statisticfun = @pl_cellmean2;
    case 'tstat'
        if strcmp(vartype,'equal')
            statisticfun = @pl_celltstat2_equalvar;
        else
            statisticfun = @pl_celltstat2_unequalvar;
        end
end


%% assign tail function
switch tail
    case 'right'
        tailfun = @deal; %function handle of deal; it simply passes input to output
    case 'left'
        tailfun = @(x)(-x); %use negative to invert distributions
end


%% initialize variables

numobservation1 = length(data1); 
numobservation2 = length(data2); 
dn1 = size(data1{1});
dn2 = size(data2{1});
if any(dn1 - dn2) %check data1 data2 consistency
     error('The variables data1 and data2 have different variable dimensions.');
end
rng('shuffle'); %seed the random number generator based on the current time
data12 = [data1(:); data2(:)]; %concatenate data
clear data1 data2; %not needed any more


%% warning for processing speed

if verbose & ~isvector(data{1}) & ~exist('bwconncomp') & ~exist('spm_bwlabel')
    if ndims(data{1})<=3
        warning('Installing a) the image processing toolbox, or b) the spm software can significantly increase speed for multidimensional data!');
    else
        warning('Installing the image processing toolbox can significantly increase speed for multidimensional data!');
    end
end


%% compute statistic map values at clusteralpha level

if verbose
    disp(['Computing statistic values at cluster-defining threshold (clusteralpha = ' num2str(clusteralpha) ')...']);
end
slicesize = 50000; %slice data to reduce memory demand
[clusteralphamap] = pl_permalphamap2(data12(1:numobservation1), data12(numobservation1+1:end), 'numpermutation',numpermutation,'tail',tail,'alpha',clusteralpha,'statistic',statistic,'verbose',verbose,'slicesize',slicesize);


%% compute clustermass of original data

statmap = tailfun(statisticfun(data12(1:numobservation1),data12(numobservation1+1:end))); %statistic map of original data
statmap(statmap<=clusteralphamap) = 0; %thresholded statistic map 
[clustersize,clusters,numcluster] = pl_conncomp(statmap); %compute clusters
if numcluster == 0 %if no clusters found, return
    stat.statmap = statisticfun(data);
    stat.criticalmap = zeros(size(stat.statmap),'single');
    stat.clustermassdist = []; %clustermass max distribution
    stat.clusters = [];
    stat.clustermass = 0;
    stat.clusteralphamap = clusteralphamap;
    stat.tail = tail;
    stat.alpha = alpha;
    stat.statistic = statistic;
    stat.numpermutation = numpermutation;
    stat.clusteralpha = clusteralpha;
    stat.clusterstatistic = clusterstatistic;
    if strcmp(clusterstatistic,'wcm')
        stat.theta = theta;
    end
    return
end
switch clusterstatistic
    case 'maxsize'
        clustermass = clustersize;
    case 'maxsum'
        clustermass = zeros(numcluster,1,'single');
        for c = 1:numcluster
            clustermass(c) = sum(statmap(clusters{c}));
        end
    case 'wcm'
        clustermass = zeros(numcluster,1,'single');
        for c = 1:numcluster
            clustermass(c) = sum( abs(statmap(clusters{c}) - clusteralphamap(clusters{c})).^(theta / (1-theta)) );
        end
end


%% compute clustermass max distribution using permutations

if verbose
    disp('Computing clustermass distribution using permutations...');
end
clustermassdist = zeros(numpermutation,1,'single');
clustermassdist(1) = max(clustermass); %first permutation is original data
for i = 2:numpermutation %leave one permutation for original data
    %verbose output
    if verbose & ~rem(i,verbose)
        disp(['Permutation: ' num2str(i) ' out of ' num2str(numpermutation)]);
    end
    permndx = randperm(numobservation1 + numobservation2); %create permutation index
    statmapperm = tailfun( statisticfun( data12(permndx(1:numobservation1)) , data12(permndx(numobservation1+1:end)) )  ); %create permutation sample
    statmapperm(statmapperm<=clusteralphamap) = 0; %threshold permutation sample
    switch clusterstatistic
        case 'maxsize'
            [clustermassperm] = pl_conncomp(statmapperm);
        case 'maxsum'
            [~,clustersperm,numclusterperm] = pl_conncomp(statmapperm);
            clustermassperm = zeros(numcluster,1,'single');
            for c = 1:numclusterperm
                clustermassperm(c) = sum(statmapperm(clustersperm{c}));
            end
        case 'wcm'
            [~,clustersperm,numclusterperm] = pl_conncomp(statmapperm);
            clustermassperm = zeros(numclusterperm,1,'single');
            for c = 1:numclusterperm
                clustermassperm(c) = sum( abs(statmapperm(clustersperm{c}) - clusteralphamap(clustersperm{c})).^(theta / (1-theta)) );
            end
    end
    if isempty(clustermassperm) %if no clusters found
        clustermassdist(i) = 0;
    else
        clustermassdist(i) = max(clustermassperm); %distribution of maximum clustermass
    end
end


%% compute clustermass threshold and significant clusters

clustermassdist = sort(clustermassdist); 
clustermassth = clustermassdist(ceil(numpermutation*(1-alpha)));
criticalclusters = clustermass>clustermassth;
clustermass = clustermass(criticalclusters); %keep only critical clusters
clusters = clusters(criticalclusters); %keep only critical clusters
criticalndx = cat(1,clusters{:}); %find critical locations


%% parse output

%clusters
stat.statmap = statisticfun(data12(1:numobservation1),data12(numobservation1+1:end)); %recompute to avoid tailfun
stat.criticalmap = zeros(size(stat.statmap),'single'); stat.criticalmap(criticalndx) = 1;
stat.clustermassdist = clustermassdist; %clustermass max distribution
stat.clusters = clusters;
stat.clustermass = clustermass;
stat.clusteralphamap = clusteralphamap;

% store input parameters
stat.tail = tail;
stat.alpha = alpha;
stat.statistic = statistic;
stat.numpermutation = numpermutation;
stat.clusteralpha = clusteralpha;
stat.clusterstatistic = clusterstatistic;
if strcmp(clusterstatistic,'wcm')
    stat.theta = theta;
end














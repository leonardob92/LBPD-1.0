function stat = pl_permtest(data,varargin)

% Performs a permutation test on 'data' for difference against 0. 
% For a difference against any mean value m (e.g. baseline decoding accuracy 50), use input: data - m.
% The test returns the statistical map (using the requested 'statistic'),
% the pvalue statistical map, and significance using false discovery rate
% and familywise error rate
%
% The data array can have any number of dimensions N in the order [samples x variable1 x variable2 x ...]
% or be a cell in the form data{samples} = [variable1 x variable2 x ...]
% (the data will be converted internally in this form for faster execution)
%
%
% stat = pl_permtest(data) uses default values
%   
% stat = pl_permtest(data,'param1',val1,'param2',val2,...) specifies one or
% more of the following name/value pairs
%
%       Parameter           Value
%       'alpha'             A value between 0-1 (default 0.05), specifying the significance level.
%       'statistic'         'tstat' (default), 'mean'
%                           Defines the test statistic, either a t-statistic or the mean  
%       'numpermutation'    An integer value, the number of permutations (default 1000)
%       'tail'              'both' (default),'right','left' 
%                           Determines the alternative hypothesis mean is not 0, mean > 0, mean < 0
%       'verbose'           Display iteration progress  
%
% stat  .statmap                    Statistical map (variable1 x variable2 x ...)
%       .statmappv                  The pvalue statistical map
%       .FDR.criticalmap            False discovery rate (FDR) critical map (binary map, 1 significant, 0 non-significant
%       .FDR.statmappvadjusted      FDR adjusted pvalue map
%       .pvthreshold                FDR threshold for the requested 'alpha' level
%       .FWER.criticalmap           Familywise error rate critical map (binary map, 1 significant, 0 non-significant
%       .FWER.maxdist               Distribution of the maximal statistic
%       .FWER.statthreshold         FWER threshold of the requested 'alpha' level
%       .tail                       Same as input or default values (for bookkeeping)
%       .alpha                      Same
%       .statistic                  Same
%       .numpermutation             Same


% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is currently under development, please do not share


%% parse inputs 

alpha           = pl_inputparser(varargin,'alpha',0.05,@(x) isscalar(x) && (x > 0) && (x < 1));
statistic       = pl_inputparser(varargin,'statistic','tstat',{'tstat','mean'});
numpermutation  = pl_inputparser(varargin,'numpermutation',1000,@(x) isscalar(x) && x>0 && x == round(x) );
tail            = pl_inputparser(varargin,'tail','both',{'both','right','left'});
verbose         = pl_inputparser(varargin,'verbose',false);


%% convert input data to cell arrays (if not already)

if ~iscell(data) %convert data to cell (faster execution)
    data = pl_mat2cell(data);
end


%% assign statistic function

switch statistic
    case 'mean'
        statisticfun = @pl_cellmean;
    case 'tstat'
        statisticfun = @pl_celltstat;
end


%% assign tail function
switch tail
    case 'right'
        tailfun = @deal; %function handle of deal; it simply passes input to output
    case 'both'
        tailfun = @abs; %function handle of abs; useful for a two-sided test
    case 'left'
        tailfun = @(x)(-x); %use negative to invert distributions
end


%% initialize variables

numobservation = length(data); 
dn = size(data{1});
rng('shuffle'); %seed the random number generator based on the current time


%% perform permutations

statmap = tailfun(statisticfun(data)); %statistic map of original data
statmappv = ones(dn,'single'); %original sample is included in the pvalues 
maxdist = zeros(numpermutation,1,'single'); %distribution of max statistic for FWER
maxdist(1) = max(statmap(:)); %first permutation sample is original data
for i = 2:numpermutation
    %verbose output
    if verbose & ~rem(i,verbose)
        disp(['Permutation: ' num2str(i) ' out of ' num2str(numpermutation)]);
    end
    permndx = sign(randn(numobservation,1,'single')); %create permutation index
    statmapperm = tailfun(statisticfun(data,permndx)); %create permutation sample
    statmappv = statmappv + (statmapperm >= statmap); %update pvalue
    maxdist(i) = max(statmapperm(:)); % update maximum statistic distribution
end
statmappv = statmappv / numpermutation; %fix pvalues


%% parse output

%maps
stat.statmap = statisticfun(data);
stat.statmappv = statmappv;

% false discovery rate statistics
[pthr,pcor,padj] = fdr(stat.statmappv(:));
stat.FDR.criticalmap = zeros(size(statmap),'single'); stat.FDR.criticalmap(padj<alpha) = 1; %adjusted pvalues less than alpha
stat.FDR.statmappvadjusted = reshape(padj,size(stat.statmap));
stat.FDR.pvthreshold = pthr;

% familywise error rate statistics
maxdist = sort(maxdist);
statthreshold = maxdist(ceil(numpermutation*(1-alpha)));
stat.FWER.criticalmap = tailfun(stat.statmap) > statthreshold; 
stat.FWER.maxdist = maxdist;
stat.FWER.statthreshold = statthreshold;

% store input parameters
stat.tail = tail;
stat.alpha = alpha;
stat.statistic = statistic;
stat.numpermutation = numpermutation;








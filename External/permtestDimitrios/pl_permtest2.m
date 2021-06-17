function [stat] = pl_permtest2(data1,data2,varargin)

% Performs a 2-sample permutation test data1 versus data2 for difference against means. 
%
% stat = pl_permtest2(data1,data2) uses default values
%   
% stat = pl_permtest2(data1,data2,'param1',val1,'param2',val2,...) specifies one or
% more of the following name/value pairs
%
%       Parameter           Value
%       'alpha'             A value between 0-1 (default 0.05), specifying the significance level.
%       'statistic'         'tstat' (default), 'mean'
%                           Defines the test statistic, either a t-statistic or the mean  
%       'numpermutation'    An integer value, the number of permutations (default 1000)
%       'tail'              'both' (default),'right','left' 
%                           Determines the alternative hypothesis mean is not 0, mean > 0, mean < 0
%       'vartype'           'equal' (default), 'unequal', determines if the t-statistic assumes equal or unequal variance 
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
vartype         = pl_inputparser(varargin,'vartype','equal',{'equal','unequal'});
verbose         = pl_inputparser(varargin,'verbose',false);


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
    case 'both'
        tailfun = @abs; %function handle of abs; useful for a two-sided test
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


%% perform permutations
statmap = tailfun(statisticfun(data12(1:numobservation1),data12(numobservation1+1:end))); %statistic map of original data
statmappv = ones(dn1,'single'); %original sample is included in the pvalues 
maxdist = zeros(numpermutation,1,'single'); %distribution of max statistic for FWER
maxdist(1) = max(statmap(:)); %first permutation sample is original data
for i = 2:numpermutation
    %verbose output
    if verbose & ~rem(i,verbose)
        disp(['Permutation: ' num2str(i) ' out of ' num2str(numpermutation)]);
    end
    permndx = randperm(numobservation1 + numobservation2); %create permutation index
    statmapperm = tailfun( statisticfun( data12(permndx(1:numobservation1)) , data12(permndx(numobservation1+1:end)) )  ); %create permutation sample
    statmappv = statmappv + (statmapperm >= statmap); %update pvalue
    maxdist(i) = max(statmapperm(:)); % update maximum statistic distribution
end
statmappv = statmappv / numpermutation; %fix pvalues


%% parse output

%maps
stat.statmap = statisticfun(data12(1:numobservation1),data12(numobservation1+1:end));
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








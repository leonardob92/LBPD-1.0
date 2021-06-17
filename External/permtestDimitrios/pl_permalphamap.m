function [alphamap] = pl_permalphavalue(data, varargin)
%
% Returns the location-specific alpha values of a statistical map created
% from a statistic applied on 'data'

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is currently under development, please do not share


%% parse inputs 

numpermutation      = pl_inputparser(varargin,'numpermutation',1000,@(x) isscalar(x) && x>0 && x == round(x) );
tail                = pl_inputparser(varargin,'tail','right',{'right','left','both'}); 
alpha               = pl_inputparser(varargin,'alpha',0.05,@(x) isscalar(x) && (x > 0) && (x < 1));
statistic           = pl_inputparser(varargin,'statistic','tstat',{'tstat','mean'});
verbose             = pl_inputparser(varargin,'verbose',true);
slicesize           = pl_inputparser(varargin,'slicesize',50000);


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
    case 'left'
        tailfun = @(x)(-x); %use negative to invert distributions
end


%% initialize variables

numobservation = length(data); 
dn = size(data{1});
rng('shuffle'); %seed the random number generator based on the current time


%% calculate alpha map

alphamap = zeros(size(data{1}));
n = prod(dn);
[sndx,ssize,nslices] = pl_sliceblocks(n,slicesize);
for s = 1:nslices
    
    %verbose output
    if verbose & ~rem(s,verbose)
        disp(['Slice: ' num2str(s) ' out of ' num2str(nslices)]);
    end
    
    dataslice = pl_cellslice(data,sndx{s});
    statmap = tailfun(statisticfun(dataslice)); %first permutation sample is original data
    statmapperm = zeros(ssize(s),numpermutation,'single'); 
    statmapperm(:,1) = statmap'; %original sample is included in the permutations 
    for i = 2:numpermutation
        permndx = sign(randn(numobservation,1,'single')); %create permutation index
        statmapperm(:,i) = tailfun(statisticfun(dataslice,permndx)); %create permutation sample
    end
    if exist('prctile') %if image processing toolbox installed
        alphamap(sndx{s}) = prctile(statmapperm',100*(1-alpha));
    else
        statmapperm = sort(statmapperm,2);
        alphamap(sndx{s}) = statmapperm(:,ceil(numpermutation*(1-alpha)));
    end
end


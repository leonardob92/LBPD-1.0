function [alphamap] = pl_permalphavalue2(data1, data2, varargin)
%
% Returns the location-specific alpha values of a statistical map created
% from a 2-sample statistic between data1 and data2

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is currently under development, please do not share


%% parse inputs 

numpermutation      = pl_inputparser(varargin,'numpermutation',1000,@(x) isscalar(x) && x>0 && x == round(x) );
tail                = pl_inputparser(varargin,'tail','right',{'right','left','both'}); 
alpha               = pl_inputparser(varargin,'alpha',0.05,@(x) isscalar(x) && (x > 0) && (x < 1));
statistic           = pl_inputparser(varargin,'statistic','tstat',{'tstat','mean'});
vartype             = pl_inputparser(varargin,'vartype','equal',{'equal','unequal'});
verbose             = pl_inputparser(varargin,'verbose',true);
slicesize           = pl_inputparser(varargin,'slicesize',50000);


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


%% calculate alphamap

alphamap = zeros(size(data12{1}));
n = prod(dn1);
[sndx,ssize,nslices] = pl_sliceblocks(n,slicesize);
for s = 1:nslices
    
    %verbose output
    if verbose & ~rem(s,verbose)
        disp(['Slice: ' num2str(s) ' out of ' num2str(nslices)]);
    end
    
    data12slice = pl_cellslice(data12,sndx{s});
    statmap = tailfun(statisticfun(data12slice(1:numobservation1),data12slice(numobservation1+1:end))); %statistic map of original data
    statmapperm = zeros(ssize(s),numpermutation,'single'); 
    statmapperm(:,1) = statmap'; %original sample is included in the permutations 
    for i = 2:numpermutation
        permndx = randperm(numobservation1 + numobservation2); %create permutation index
        statmapperm = tailfun( statisticfun( data12slice(permndx(1:numobservation1)) , data12slice(permndx(numobservation1+1:end)) )  ); %create permutation sample
    end
    if exist('prctile') %if image processing toolbox installed
        alphamap(sndx{s}) = prctile(statmapperm',100*(1-alpha));
    else
        statmapperm = sort(statmapperm,2);
        alphamap(sndx{s}) = statmapperm(:,ceil(numpermutation*(1-alpha)));
    end
end


function P = morlet_transform(x,t,f,fc,FWHM_tc,squared)
% function P = morlet_transform(x,t,f,fc,FWHM_tc,squared)
%
% Applies complex Morlet wavelet transform to the timeseries stored in the 
% matrix x with size (ntimeseries x ntimes). It returns a wavelet coefficient map (by default
% squared)
%
% INPUTS:
%   x: (ntimeseries x ntimes) a vector of the timeseries
%   t: (1 x ntimes) a vector of times (in secs)
%   f: (1 x nfreqs) a vector of frequencies in which to estimate the
%   wavelet transform (in Hz). Default is 1:60Hz
%   fc: (default is 1) central frequency of complex Morlet wavelet in Hz
%   FWHM_tc: (default is 3) FWHM of complex Morlet wavelet in
%   time. Also see morlet_design.m. 
%   squared: 'y' (default) or 'n'. Flag that decided whether the function returns the
%           squared coefficients (y) or not (n). Squaring represents neural power in
%           the corresponding frequency.
%
% OUTPUT:
%   P: (ntimeseries x nfreqs x ntimes), or
%   P: (nFreqs x nTimes), if ntimeseries = 1
%      A matrix of wavelet coefficients (by default squared)
%
% Example
%   t = 0:.01:1;
%   f = 5:20;
%   P = morlet_transform(sin(2*pi*10*t),t,f);
%   Coefs = morlet_transform(sin(2*pi*10*t),t,f,[],[],'n');
%
% See also:
%   MORLET_PLOT_COEF
%
% Author: Dimitrios Pantazis, December 2008, USC

%test inputs
if nargin == 0 %if no input
    help morlet_transform %display help
    return
end
if(exist('fc') & ~isempty(fc))
else
    %default complex Morlet resolution:
    fc = 1; 
end
if(exist('FWHM_tc') & ~isempty(FWHM_tc))
else
    %default complex Morlet resolution:
    FWHM_tc = 3; 
end
%If fc = 1 and sigma_tc = 1.5, at 1Hz central frequency, it has 2*sigma = 3s temporal resolution and 0.212Hz frequency resolution
if(exist('squared') & ~isempty(squared))
else
    squared = 'y'; %by default return squared coefficients, ie neural power
end
if(exist('f') & ~isempty(f))
else
    f = 1:60; %default frequencies in Hz
end
if size(x,2)==1 %vector given as (ntimes x ntimeseries)
    x = x';
end


%signal parameters
Ts = t(2)-t(1); %sampling period of signal
Fs = 1/Ts; %sampling frequency of signal

%complex morlet wavelet parameters
scales = f ./ fc; %scales for wavelet
sigma_tc = FWHM_tc / sqrt(8*log(2));
sigma_t = fc*sigma_tc./f;
nscales = length(scales);

%compute wavelet kernels for each scale
precision = 3;
for s = 1:nscales
    xval = -precision*sigma_t(s): 1/Fs : precision*sigma_t(s);
    W{s} = sqrt(scales(s)) * morlet_wavelet(scales(s)*xval,fc,sigma_tc);
    %figure;plot(xval,real(W{s}));
end
    
%compute wavelet coefficients
nx = size(x,1); %number of timeseries
ntimes = size(x,2); %number of timepoints
P = zeros(nx,nscales,ntimes);
for s = 1:nscales
    %progress(s,1,nscales,1,'scales');
    %Convolution flips W{s}. So, I need to unflip it to become inner product.
    %Convolution does not conjugate, but wavelet transform does: <f,w> = f*conj(w)
    %So, in the convolution I need to use conj(fliplr(w)) instead of w
    %But for the complex Morlet wavelet they are equal!
    for ch = 1:size(x,1) %it is slower to use directly  P(:,s,:) = conv2(x,W{s},'same') * Ts; 
        P(ch,s,:) = conv2(x(ch,:),W{s},'same') * Ts; 
    end
end

%if only one timeseries, compress first dimension
if size(P,1)==1
    P = squeeze(P);
end

%if return squared coefficients
if strcmp(squared,'y')
    P = abs(P).^2; %return neural power
end


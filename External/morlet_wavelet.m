function W = morlet_wavelet(t,fc,sigma_tc);
% function W = morlet_wavelet(t,fc,sigma_tc)
%
% Returns the complex Morlet wavelet for a specified central frequency and standard deviation
%
% INPUTS:
%   t: timepoints where the wavelet will be calculated
%   fc: central frequency
%   sigma_tc: standard deviation of Gaussian kernel in time at the central
%   frequency. Decreasing sigma_tc improves temporal resolution (because
%   the wavelet has smaller temporal extent) at the expense of frequency
%   resolution
%
% Example values:
%  fc = 1;
%  sigma_tc = 1.5;
%
% Author: Dimitrios Pantazis, December 2008, USC

%complex Morlet wavelet
W = (sigma_tc*sqrt(pi))^(-0.5) * exp( -(t.^2)/(2*sigma_tc^2) ) .* exp(i*2*pi*fc*t);

function [  ] = LBPD_startup_D( pathl )

%It adds some required paths and startup OSL..


% INPUT:    -pathl:             path to functions (character)
%
% OUTPUT:   -starting up.. 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 24/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%adding paths..
addpath([pathl '/External']); %to some external functions..
addpath([pathl '/External/creatingAALnifti']);
addpath([pathl '/External/altmany-export_fig-412662f']);
addpath([pathl '/External/giftools']);
addpath([pathl '/External/schemaball-master']);
addpath([pathl '/External/permtestDimitrios']); %path to a bunch of functions written by Dimitrios Pantazis for statistical analysis
a = dir([pathl '/External/*osl']);
if ~isempty(a)
    addpath([pathl '/External/osl/osl-core']); %to a (very) sligthly edited version (by Leonardo) of OSL
    osl_startup %starting up OSL
    %if OSL exists, than remove SPM12 and FieldTrip to use the SPM12 and FieldTrip codes redistributed by OSL
    %rmdir([pathl '/External/spm12'],'s')
    %rmdir([pathl '/External/fieldtrip-20171231'],'s')
else %otherwise starting up SPM12 and FielTrip (assuming that the user donwloaded them and placed them in "LBPD/External"
    addpath([pathl '/External/spm12']); %path to spm12
    % starting up FieldTrip
    addpath([pathl '/External/fieldtrip-20171231']);
    ft_defaults
end


end


function [ Dofl ] = extracting_data_LBPD_D( listfl )

% Extracts data from each subject (SPM objects) and store it (with a
% corresponding label) in a single file.
% The SPM object should be a 2D or 3D data (channels x time-points x
% trials).
% This function is thought to be used before preparing_baseline_from_restingstate_LBPD_D.mat.



%   INPUT:  -listfl:        -list of paths to subject files.
%                            It is thought to be the output of the 'dir'
%                            function in matlab.

%   OUTPUT: -Dofl:          -cell with data (1st column) and corresponding
%                            label (2nd column)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% calumaxs@gmail.com
% leonardo.bonetti@clin.au.dk
% Francesco Carlomagno, Aarhus, DK, 20/03/2019
% Leonardo Bonetti, Aarhus, DK, 01/10/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








Dofl = cell(length(listfl),2); %preallocating space
for p = 1:length(listfl) %over subjects
    D = spm_eeg_load([listfl(p).folder '/' listfl(p).name]); %loading data
    Dofl(p,1) = {D(:,:,:)}; %extracting and storing data
    Dofl(p,2) = {listfl(p).name}; %storing labels
    disp(['you have just computed the subj num ' num2str(p)]) %showing progress..
end


end


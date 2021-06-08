function [ OUT ] = preparing_baseline_from_restingstate_LBPD_D( x, y, datapath )

% Given a cell matrix with resting state data for all subjects, it creates
% simulated trials, averages them, corrects for source leakage and computes
% their envelopes.
% Future updates may give user more options (e.g. calculating or not mean,
% source leakage correction, envelope..).



%   INPUT:  -x:         length of the simulated trial (in time-samples)
%           -y:         number of simulated trials
%           -datapath:  path to the .mat data file (character) containing a 
%                       variable named: Drest.
%                       Drest is a number of subjects x 2 cell matrix with:
%                        1st column: actual data (ROIs x time-points double in each cell)
%                        2nd column: subject IDs

%   OUTPUT: -OUT:       number of subjects x 2 cell matrix with:
%                        1st column: processed data (ROIs x time-points double in each cell)
%                        Data has been submitted to:
%                          -simulated trials creation
%                          -mean across simulated trials
%                          -source leakage correction (orthogonalisation)
%                          -calculation of envelope
%                        2nd column: subject IDs





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 24/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%preparing onsets of simulated trials
d = zeros(y,2);
countx = x*(-1) + 1;
for ii = 1:y
    countx = countx + x;
    d(ii,1) = countx;
    d(ii,2) = countx + x - 1;
end
%loading data
load(datapath); %this loads a .mat file containing a cell matrix called Drest
if ~exist('Drest','var')
    error('the resting state file you gave must be named "Drest".. ')
end

%computation
OUT = cell(length(Drest),2); %preallocating space
for ii = 1:length(Drest) %over subjects (assuming that you have at least two subjects..) 
    a = zeros(90,x,y); %preallocating space
    D2 =  Drest{ii,1}; %extracting data for subject ii
    for jj = 1:y %over simulated trials
        a(:,:,jj) = D2(:,d(jj,1):d(jj,2)); %using simulated onsets for creating simulated trials
    end
    if find(a == 0) %check for 0s.. this warning controls several common and dramatic mistakes, but if it does not sayig anything it does not necessarily mean that everything is fine..
        warning(['data for subject ' num2str(ii) ' contains 0s.. you may have made some mistakes in previous computations..']); %warning message
    end
    dum = mean(a,3); %averaging simulated trials
    dum2 = remove_source_leakage_b(dum,'symmetric'); %correcting for source leakage by computing orthogonalisation
    OUT(ii,1) = {hilbenv(dum2)}; %calculating envelope and storing data
    OUT(ii,2) = Drest(ii,2); %storing IDs
    disp(['subject ' num2str(ii) ' has just been computed']) %showing progress of computation
end



end


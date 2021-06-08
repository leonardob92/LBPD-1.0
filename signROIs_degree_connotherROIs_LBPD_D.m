function [ POS, NEG ] = signROIs_degree_connotherROIs_LBPD_D( INT, COUP, lab_ab, r, X )

% It identifies the connections between the ROIs significantly central within
% the brain network and the other ROIs.
% Currently, it is implemented for the contrast between two conditions or for
% working with only condition 1 or condition 2 (for having both cond 1 and
% 2 just run the function two times).
% This is not a statistical test, but simply a function that retrievs
% information after that statistics has been computed in previous functions
% such as: 'diff_conditions_matr_degree_MCS_LBPD_D.m' and
% 'ps_statistics_LBPD_D.m'.
% This function should be used after 'diff_conditions_matr_degree_MCS_LBPD_D.m'.


%   INPUT:  -INT:       ROIs x time x conditions, describing when ROIs were
%                       significantly central within the brain network
%                       (binary double)
%           -COUP:      ROIs x ROIs x time x conditions, describing the
%                       significant coupling between ROIs over time (binary
%                       double)
%           -lab_ab:    ordered label of ROIs for conditions 1 and 2
%                       (character).
%                       These labels are in the order of INT and COUP
%                       data..
%           -r:         statistics outputted by
%                       'diff_conditions_matr_degree_MCS_LBPD_D.m'.
%                       It can be either the stats for cond1 or cond2 or
%                       their contrast.
%           -X:         time-points to be used (they must match the ones used in
%                       'degree_segregation_MCS_LBPD_D.m'.
%                       Set = 0 to have all time-points.

%   OUTPUT: -POS:       for condition 1; cell describing the ROIs that were significantly
%                       central within the brain network, followed by:
%                       ROIs that were significantly coupled with original ROIs and how
%                       many times.
%                       Coupled ROIs are sorted in descendent order.
%           -NEG:       same as POS but for condition 2






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 16/01/2019
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%setting some parameters and checking if you submitted a contrast or simply
%cond1 or cond2 only
S = size(INT);
R = size(r);
if R(2) == 3 %this is a trick for using this function also for condition 1 or 2 without the contrast..
    rdu = r(:,3);
elseif R(2) == 2
    rdu = r(:,2);
end
if X == 0
   X = 1:S(2);
end

%cond1 vs cond2 (if you have a contrast, otherwise this works as cond1 or cond2)
Ipos = find(double(abs(double(cellfun('isempty',rdu))-1) > 0) + double(cell2mat(r(:,2)) > 0) == 2);
labs = string(lab_ab);
if ~isempty(Ipos)
    POS = cell(length(Ipos),S(1)*2+1);
    for ii = 1:length(Ipos)
        Id = find(strcmp(r{Ipos(ii),1},labs)); %finding index of ROI in the previous Order
        I2pos = find(INT(Id,X,1) == 1); %finding time-point where ROI Id was significantly central
        coup_sum = sum(COUP(Id,:,I2pos),3);
        [V,I] = sort(coup_sum,'descend');
        POS(ii,1) = {string(labs(Id))};
        count2 = 0;
        for pp = 1:length(I)
            count2 = count2 + 2;
            POS(ii,count2) = {labs(I(pp))};
            POS(ii,count2+1) = {V(pp)};
        end
    end
else
    POS = [];
end

%cond2 vs cond1 (if you have a contrast)
Ineg = find(double(abs(double(cellfun('isempty',rdu))-1) > 0) + double(cell2mat(r(:,2)) < 0) == 2); %this is not empty only in case of significance after contrast
labs = string(lab_ab);
if ~isempty(Ineg)
    NEG = cell(length(Ineg),S(1)*2+1);
    for ii = 1:length(Ineg)
        Id = find(strcmp(r{Ineg(ii),1},labs)); %finding index of ROI in the previous Order
        I2neg = find(INT(Id,X,2) == 1); %finding time-point where ROI Id was significanlty central
        coup_sum = sum(COUP(Id,:,I2neg),3);
        [V,I] = sort(coup_sum,'descend');
        NEG(ii,1) = {string(labs(Id))};
        count2 = 0;
        for pp = 1:length(I)
            count2 = count2 + 2;
            NEG(ii,count2) = {labs(I(pp))};
            NEG(ii,count2+1) = {V(pp)};
        end
    end
else
    NEG = [];
end




end


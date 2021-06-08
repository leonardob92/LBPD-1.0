function [ MATF ] = generalCoupling_prepareplotting_LBPD_D( POS, NEG, p_thresh, label_path, pos_lab_only )

% It creates an ROIs x ROIs matrix, ordering ROIs according to labels
% provided (example that I provide consists in AAL ordered like LRLRLR),
% showing the sum over time of significant connections between ROIs.
% This should be used after 'diff_conditions_matr_couplingROIs_MCS_LBPD_D.m'
% and prepares data for subsequent plotting.


%   INPUT:  -POS:           matrix outputted by 'diff_conditions_matr_couplingROIs_MCS_LBPD_D.m'
%                           with in 1st column ROIs significantly central
%                           within the brain network and then ROIs and
%                           corresponding number of how many times those ROIs
%                           were significantly coupled with original ROI in
%                           1st column (POS refers to cond1 > cond2)
%           -NEG:           same but for cond2 > cond1
%           -p_thresh:      option for changing alpha level.
%                           Leave empty [] for default = 0.05.
%           -label_path:    path for ROIs labels (must be a .mat file with
%                           character variable named 'lab')
%           -pos_lab_only:  1 for using POS only
%                           0 for both POS and NEG

%   OUTPUT: -MATF:          connectivity matrix between ROIs, ordered
%                           according to labels.
%                           In my proposed example it is LRLRLR.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, MIT, Boston, USA, 13/08/2019
% Leonardo Bonetti, Aarhus, DK, 26/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%preparing data
dum = POS; %input (P_pos or P_neg);
if isempty(p_thresh)
    p_thresh = 0.05;
end
matf = zeros(90,90);
load(label_path);
labs = string(lab);
%pos (cond1 vs cond2)
W = find(cell2mat(dum(:,3)) < p_thresh); %getting the p_values after possible threshold if required (it may be possible that you want to use a Bonferroni correction, maybe if you run the analysis in different frequency bands)
for ii = 1:length(W)
    Idi1 = find(strcmp(dum{ii,1},labs)); %get original index of ROI1 in the couple
    Idi2 = find(strcmp(dum{ii,2},labs)); %get original index of ROI2 in the couple
    matf(Idi1,Idi2) = dum{ii,4};
    matf(Idi2,Idi1) = dum{ii,4};
end

%if you want also cond2 > cond1
if pos_lab_only == 0
    %neg
    %same
    dum = NEG;
    matfn = zeros(90,90);
    W = find(cell2mat(dum(:,3)) < p_thresh); 
    for ii = 1:length(W)
        Idi1 = find(strcmp(dum{ii,1},labs));
        Idi2 = find(strcmp(dum{ii,2},labs));
        matfn(Idi1,Idi2) = dum{ii,4};
        matfn(Idi2,Idi1) = dum{ii,4};
    end
    %this is for putting together positive and negative contrasts (in the same connectivity matrix)
    MATF = zeros(90,90);
    for ii = 1:90
        for jj = 1:90    
            if matf(ii,jj) ~= 0 %the idea is if one matrix (pos) has some values and the other (neg) not take them (pos); precedence given to biggest element in case of overlaps..
                if matfn(ii,jj) == 0
                    MATF(ii,jj) = matf(ii,jj);
                else
                    if matf(ii,jj) > matfn(ii,jj)
                        MATF(ii,jj) = matf(ii,jj);
                    else
                        MATF(ii,jj) = matfn(ii,jj);
                    end
                end
            else
                if matfn(ii,jj) ~= 0
                    MATF(ii,jj) = matfn(ii,jj);
                end
            end
        end
    end
else
    MATF = matf;
end




end


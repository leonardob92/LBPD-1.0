function [ MATF ] = signROIdegree_otherROI_prepareplotting_LBPD_D( POS, NEG, label_path )

% It creates an ROIs x ROIs matrix, ordering ROIs according to labels
% provided (example that I provide consists in AAL ordered like LRLRLR),
% showing the number of significant connections between couple of ROIs 
% over time.
% This should be used after 'diff_conditions_matr_couplingROIs_MCS_LBPD_D.m' and
% prepares data for subsequent plotting.



%   INPUT:  -POS:           matrix outputted by 'diff_conditions_matr_couplingROIs_MCS_LBPD_D.m'
%                           with in 1st column ROIs significantly central
%                           within the brain network and then ROIs and
%                           corresponding number of how many times those ROIs
%                           were significantly coupled with original ROI in
%                           1st column (POS refers to cond1 > cond2)
%           -NEG:           same but for cond2 > cond1
%           -label_path:    path for ROIs labels (must be a .mat file with
%                           character variable named 'lab')

%   OUTPUT: -MATF:          connectivity matrix between ROIs, ordered
%                           according to labels (LRLRLR in my provided
%                           example). I provide this since for later
%                           plotting functions it is useful.






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, MIT, Boston, USA, 13/08/2019
% Leonardo Bonetti, Aarhus, DK, 26/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%cond1 > cond2
dum = POS;
matf = zeros(90,90);
SC = size(dum);
load(label_path)
labs = string(lab);
for ii = 1:SC(1)
    Idii = find(strcmp(dum{ii,1},labs)); %get original index of ROI significanlty central within the brain network
    for jj = 1:SC(2)-2
        Idjj = find(strcmp(dum{ii,jj+1},labs)); %get original lab index of the ROIs coupled with with ROI Idii
        if matf(Idjj,Idii) ~= 0 %dealing with the problem that ROI 1 might be coupled x time with ROI2, while ROI2 might be coupled y time with ROI1.. to deal with this if this condition exists..
            matf(Idii,Idjj) = (matf(Idjj,Idii) + dum{ii,jj+2})/2; %do the average
            matf(Idjj,Idii) = (matf(Idjj,Idii) + dum{ii,jj+2})/2; %also simmetrically (better to do like that instead of recreating the symmetric matrix otherwise you could lose some values.. since here I am storing some values in the upper triangle and some others in the lower triangle..
        else
            matf(Idii,Idjj) = dum{ii,jj+2}; %otherwise get normally the single value..
            matf(Idjj,Idii) = dum{ii,jj+2}; %also simmetrically
        end %not ideal but good enough..
    end
end

%cond2 > cond1
%same
dum = NEG;
matfn = zeros(90,90);
SC = size(dum);
labs = string(lab);
for ii = 1:SC(1)
    Idii = find(strcmp(dum{ii,1},labs));
    for jj = 1:SC(2)-2
        Idjj = find(strcmp(dum{ii,jj+1},labs)); 
        if matfn(Idjj,Idii) ~= 0 
            matfn(Idii,Idjj) = (matfn(Idjj,Idii) + dum{ii,jj+2})/2;
            matfn(Idjj,Idii) = (matfn(Idjj,Idii) + dum{ii,jj+2})/2;
        else
            matfn(Idii,Idjj) = dum{ii,jj+2};
            matfn(Idjj,Idii) = dum{ii,jj+2};
        end 
    end
end

%this is for putting together positive and negative contrasts (in the same connectivity matrix)
matfn = matfn*(-1); %changing the sing
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


end


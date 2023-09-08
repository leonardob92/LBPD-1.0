function [ PDn ] = Extract_BrainCluster_Information_3D_LBPD_D( S )

% It extracts information about the clusters of the provided brain image.




%  INPUT:   -S.fname:       path and name of the figure (nifti format) depicting the brain clusters
%           -S.stats:       label for the statistics (or other category of values) which form the clusters in the brain
%                           (leave empy [] or do not provide the field to get "T-val" by default)


%  OUTPUT:  -PDn:           table with information on the clusters (e.g. hemisphere, t-value, coordinates, AAL label)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Oxford, UK, 29/05/2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    








%getting path and name of the brain figure (with clusters)
fname = S.fname;
if ~isfield(S,'stats') || isempty(S.stats)
    stats = 'T-val';
else
    stats = S.stats;
end
[ mni_coords, xform ] = osl_mnimask2mnicoords(fname); %getting MNI coordinates of significant voxels within the provided image
V = nii.load(fname); %loading the image
VV = V(V~=0); %extracting statistics
VI = find(V~=0); %indices of non-zero values of nifti image
parcelfile = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'; %path to AAL template %load this from the provided codes folder
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat'); %loading AAL labels %load this from the provided codes folder
K = nii.load(parcelfile); %extracting AAL coordinates information
%sorting results in order to have strongest voxels (in absolute terms) at the top
if abs(max(VV)) > abs(min(VV)) %if positive t-values are bigger in absolute terms than negative t-values
    [VV2, II] = sort(VV,'descend'); %positive t-values on top
else
    [VV2, II] = sort(VV,'ascend'); %otherwise negative t-values on top
end
VI = VI(II);
mni_coords = mni_coords(II,:);
PD = cell(length(VV2)+2,6);  %final cell
ROI = zeros(length(VI),1); %getting AAL indices
PD{1,1} = 'ROI'; PD{1,2} = 'Hemisphere'; PD{1,3} = stats; PD{1,4} = 'MNI Coordinates';
PD{2,4} = 'x'; PD{2,5} = 'y'; PD{2,6} = 'z'; %coordinate labels are not perfectly aligned with the numbers.. think if and how to fix this..
cnt = 2;
for ii = 1:length(VI)
    ROI(ii) = K(VI(ii));
    if ROI(ii) > 0 && ROI(ii) < 91
        cnt = cnt + 1;
        PD(cnt,1) = {lab(ROI(ii),3:end)}; %storing ROI
        PD(cnt,4) = {mni_coords(ii,:)}; %storing MNI coordinates
        if mni_coords(ii,1) > 0 %storing hemisphere
            PD(cnt,2) = {'R'};
        else
            PD(cnt,2) = {'L'};
        end
        PD(cnt,3) = {round(VV2(ii),2)}; %storing t-statistics
    end
end
% PDn = cell2table(PD(~any(cellfun('isempty',PD),2),:)); %remove the possible empty cell
PDn = cell2table(PD); %table


end


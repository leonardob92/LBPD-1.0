function idx = From3DNifti_OrMNICoords_2_CoordMatrix_8mm_LBPD_D( S )
 
% It translates the voxels of a 3D brain nifti image or sets of MNI coordinates
% (8mm) in LBPD coordinates (matrix with 3559 voxels).




% INPUT:    -S.input:       1 = MNI coordinates, e.g. [-45 -32 6];
%                           2 = AAL ROIs e.g. [45 50 89];
%                           3 = general image with non-zero values (provided as a nifti image)
%           -S.coordd:      MNI coordinates x y z, e.g. [-45 -32 6]
%           -S.AAL_ROIs:    AAL ROIs numbers you want to use, e.g. [45 50 89] (the order is LRLRLR..);
%           -S.image:       path to nifti image with non-zero values

% OUTPUT:   -idx:           vector with LBPD coordinates (in case of
%                           S.input == 3 you get the indices in the 1st
%                           column and the original values in the 2nd columns)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Oxford, UK, 08/06/2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%getting indices of provided coordinates or AAL ROIs
if S.input == 1
    coordd = S.coordd;
    %matching coordinates and getting source indices (step by step)
    load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_coord_dyi.mat');
    idx = zeros(size(coordd,1),1);
    for ii = 1:size(coordd,1) %over voxels
        %x
        dum1 = double(coordd(ii,1) == MNI8(:,1)); %x
        if isempty(find(dum1==1)) %if there is no pefect correspondance between given coordinate and mask coordinate, you look for an approximation
            dum = coordd(ii,1) - MNI8(:,1); %differnce between given coordinate and mask coordinate
            [m1] = min(dum(dum>0)); %getting closest value to 0 (positive)
            [m2] = max(dum(dum<0)); %getting closest value to 0 (negative)
            asd = [m1 m2]; %concatenating the values
            [~,i] = min(abs([m1 m2])); %getting the minimum value n absolute terms
            dum1 = zeros(size(MNI8,1),1);
            dum1(asd(i)==dum) = 1; %assigning 1 to the requested voxels (the ones that are closest to the given coordinate)
        end
        %y
        dum2 = double(coordd(ii,2) == MNI8(:,2)); %y
        if isempty(find(dum2==1)) %if there is no pefect correspondance between given coordinate and mask coordinate, you look for an approximation
            dum = coordd(ii,2) - MNI8(:,2); %differnce between given coordinate and mask coordinate
            [m1] = min(dum(dum>0)); %getting closest value to 0 (positive)
            [m2] = max(dum(dum<0)); %getting closest value to 0 (negative)
            asd = [m1 m2]; %concatenating the values
            [~,i] = min(abs([m1 m2])); %getting the minimum value n absolute terms
            dum2 = zeros(size(MNI8,1),1);
            dum2(asd(i)==dum) = 1; %assigning 1 to the requested voxels (the ones that are closest to the given coordinate)
        end
        %z
        dum3 = double(coordd(ii,3) == MNI8(:,3)); %z
        if isempty(find(dum3==1)) %if there is no pefect correspondance between given coordinate and mask coordinate, you look for an approximation
            dum = coordd(ii,3) - MNI8(:,3); %differnce between given coordinate and mask coordinate
            [m1] = min(dum(dum>0)); %getting closest value to 0 (positive)
            [m2] = max(dum(dum<0)); %getting closest value to 0 (negative)
            asd = [m1 m2]; %concatenating the values
            [~,i] = min(abs([m1 m2])); %getting the minimum value n absolute terms
            dum3 = zeros(size(MNI8,1),1);
            dum3(asd(i)==dum) = 1; %assigning 1 to the requested voxels (the ones that are closest to the given coordinate)
        end
        a = find((dum1+dum2+dum3)==3); %index
        idx(ii) = a; %storing index
    end
elseif S.input == 2
    %AAL case
    AAL_ROIs = S.AAL_ROIs;
    aal = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/aal_8mm_try5.nii.gz'); %loading AAL nifti image
    aalt = aal.img; %extracting matrix with the values shown in the image
    mask = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %loading AAL nifti image
    maskt = mask.img; %extracting matrix with the values shown in the image
    idx = [];
    for ii = 1:length(AAL_ROIs) %over AAL ROIs
        dum = find(aalt == AAL_ROIs(ii)); %getting the indices of the selected AAL ROIs
        idx = cat(1,idx,maskt(dum)); %getting the correspondent indices of the mask upon which is based the source reconstruction (and therefore the data)
    end
    if isempty(idx) || sum(idx) == 0 %if no correspondance is found (it can happen if the ROI is very smal and has e.g. only 1 voxel which does not perfectly align with my template)
        cnt = 0;
        while isempty(idx) || sum(idx) == 0  %in case it happens, I look for the nearest voxel which matches my template (this small imprecision is negligible with MEG)
            cnt = cnt + 1;
            idx = cat(1,idx,maskt(dum+cnt));
        end
    end
    idx(idx==0) = []; %removing the few voxels that couldn't be indexed since the AAL ROIs and the mask used slightly different approximations and therefore a few voxels do not perfectly correspond (with MEG spatial resolution it does not really matter)
elseif S.input == 3 %case with general image
    IMG = load_nii(S.image);
    imgdum = IMG.img; %extracting matrix with the values shown in the image
    mask = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %loading nifti image
    maskt = mask.img; %extracting matrix with the values shown in the image
    dum = find(imgdum ~= 0); %getting the indices of the non-zero voxels
    idx = zeros(length(dum),2);
    idx(:,1) = maskt(dum); %proper indices
    idx(:,2) = imgdum(dum); %associated values (it may be useful to have them as well)
    idx(idx(:,1)==0,:) = []; %deleting 0s from indices.. this may happen for different approximations in the images but it does not represent a concern
end

end


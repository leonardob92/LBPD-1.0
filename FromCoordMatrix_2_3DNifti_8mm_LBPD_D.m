function O = FromCoordMatrix_2_3DNifti_8mm_LBPD_D( S )
O = []; 

% It translates LBPD coordinates (8mm) from matrix (3559 voxels) into a 3D
% configuration normally used to generate nifti images with brain template.




% INPUT:    -S.data:        data (voxels x ROIs (time-points))
%           -S.fname:       path and name for the image to be saved (no '.nii.gz')
%           -S.names:       cell array with names for the different images to be saved (e.g. time-points, ROIs, etc.).
%                           If you do not provide the field, it creates images using progressive numbers (i.e. 1-2-3).
%           -S.singleimage: 1 if you want to save the images in a single image (i.e. putting together different time-points or ROIs, etc.).
%                           0 or do not provide field for saving one image per time-point (or ROI, etc.).

% OUTPUT:   -3D nifti image



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@psych.ox.ac.uk
% Leonardo Bonetti, Oxford, UK, 08/06/2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







data = S.data;
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
%the following lines could be made more compact, but for historical reasons it is not faster to do like that..
if isfield(S,'singleimage')
    if S.singleimage == 1 %if you want to have all time=points (or ROIs) in one single image
        if isfield(S,'names')
            fnamenii = [S.fname '/' S.names{1} '.nii.gz']; %path and name of the image to be saved
        else
            fnamenii = [S.fname '/SingleImage.nii.gz']; %path and name of the image to be saved
        end
        SO = data(:,1);
        if ~isempty(SO)
            %building nifti image
            SS = size(maskk.img);
            dumimg = zeros(SS(1),SS(2),SS(3),size(data,2)); %no time here so size of 4D = 1
            for ii = 1:size(SO,1) %over brain sources
                dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
                [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
                dumimg(i1,i2,i3,:) = data(ii,:); %storing values
            end
            nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
            nii.img = dumimg; %storing matrix within image structure
            nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
            disp(['saving nifti single image'])
            save_nii(nii,fnamenii); %printing image
        end
    else
        for iii = 1:size(data,2) %over time-points (or ROIs, etc.)
            if isfield(S,'names')
                fnamenii = [S.fname '/' S.names{iii} '.nii.gz']; %path and name of the image to be saved
            else
                fnamenii = [S.fname '/Image_' num2str(iii) '.nii.gz']; %path and name of the image to be saved
            end
            SO = data(:,iii);
            if ~isempty(SO)
                %building nifti image
                SS = size(maskk.img);
                dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
                for ii = 1:size(SO,1) %over brain sources
                    dum = find(maskk.img == ii); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
                    [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dum); %getting subscript in 3D from index
                    dumimg(i1,i2,i3,:) = SO(ii,:); %storing values
                end
                nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
                nii.img = dumimg; %storing matrix within image structure
                nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
                disp(['saving nifti image - condition ' num2str(iii)])
                save_nii(nii,fnamenii); %printing image
            end
        end
    end
end

end


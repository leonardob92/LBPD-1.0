function create_HCP_MMP1_nifti_LBPD(vector,name)

% script to create HCP-MMP1 nifti file for rendering

% INPUTS:       vector:     a 1x90 vector with a value for each AAL region
%               name:       the name of the output file





if length(vector) == 180
    MMP1 = load_nii('HCP-MMP1_on_MNI152_ICBM2009a_nlin.nii.gz');
elseif length(vector) == 44
    MMP1 = load_nii('HCP-MMP1_22_on_MNI152_ICBM2009a_nlin.nii.gz');
else
    error('wrong vector length..')
end

% create new image
newnii = MMP1;
% set datatype to correct single type (rather than uint8)
newnii.img = zeros(size(MMP1.img),'single');
%     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32 
newnii.hdr.dime.datatype = 16;
newnii.hdr.dime.bitpix = 32;


vals = vector;


for ii = 1:length(vector)
    % find the indices of voxels for a given aal region
    vox = find(MMP1.img == ii);
    % set value for each voxel in aal region 
    newnii.img(vox) = vals(ii);
end;


%save new nifti file
save_nii (newnii, [name '.nii.gz']);


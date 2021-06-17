function create_AALnifti(aalvector,name,symm)
% script to create AAL nifti file for rendering
%  aalvector is a 1x90 vector with a value for each AAL region
%  name is the name of the output file
%  symm is a flag that if >0 means that ordering is symmetrical

aalnii=load_nii('aal_2mm.nii.gz')

% create new image
newnii=aalnii;
% set datatype to correct single type (rather than uint8)
newnii.img=zeros(size(aalnii.img),'single');
%     8 Signed integer                (int32, bitpix=32) % DT_INT32, NIFTI_TYPE_INT32 
newnii.hdr.dime.datatype=16;
newnii.hdr.dime.bitpix=32;

if (symm>0) 
    % convert symmetrical to standard
    aalvals(1:2:90)=aalvector(1:45);
    aalvals(2:2:90)=aalvector(90:-1:46);
else
    aalvals=aalvector;
end;

for i=1:90
    % find the indices of voxels for a given aal region
    vox=find(aalnii.img==i);
    % set value for each voxel in aal region 
    newnii.img(vox)=aalvals(i);
end;


%save new nifti file
save_nii (newnii, name);


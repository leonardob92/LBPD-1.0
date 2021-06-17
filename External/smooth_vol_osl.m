function v2 = smooth_vol_osl(vol_as_matrix, mask_fname, fwhm)

%computing spatial smoothing
% vol_as_matrix = data (voxels x time-points)
% mask_fname = character with path to a standard brain mask
% fwhm = paramter for spatial smoothing (e.g. 100)


    [mask,res,xform] = nii.load(mask_fname);
    mask(mask>0) = 1; % Turn it into a mask

    sd = fwhm/2.3; % std spatial smoothing, in mm
    sd = sd/res(1); % standard deviation in voxels

    smooth_mask = smooth3(mask,'gaussian',5,sd);
    smooth_mask(~mask) = 0;

    v2 = nan(size(vol_as_matrix));

    for j = 1:size(vol_as_matrix,2)
        v = matrix2vols(vol_as_matrix(:,j),mask);
        v = smooth3(v,'gaussian',5,sd);
        v(~mask) = 0;
        v = v./smooth_mask;
        v2(:,j)=vols2matrix(v,mask);
    end
end
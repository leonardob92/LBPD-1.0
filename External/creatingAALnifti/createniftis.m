load empiricalAwake.mat
load empiricalSleep.mat

create_AALnifti(mean(Hierarchy),'aw_hierarchy.nii.gz',1);
create_AALnifti(mean(Hierarchy2),'sl_hierarchy.nii.gz',1);
create_AALnifti(mean(Hierarchyn2),'sln2_hierarchy.nii.gz',1);

create_AALnifti(mean(SpatialTemporalBroadness_sink),'aw_sink.nii.gz',1);
create_AALnifti(mean(SpatialTemporalBroadness_sink2),'sl_sink.nii.gz',1);
create_AALnifti(mean(SpatialTemporalBroadness_sinkn2),'sln2_sink.nii.gz',1);

create_AALnifti(mean(SpatialTemporalBroadness_source),'aw_source.nii.gz',1);
create_AALnifti(mean(SpatialTemporalBroadness_source2),'sl_source.nii.gz',1);
create_AALnifti(mean(SpatialTemporalBroadness_sourcen2),'sln2_source.nii.gz',1);

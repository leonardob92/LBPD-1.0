function O = cluster_DTI_perm_2groups_1(S_ori)
O = []; 

%starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl);
% addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')

V_bin_1 = S_ori.V_bin_1;
V_bin_2 = S_ori.V_bin_2;
V_t = S_ori.V_t;

%cond 1 > cond 2
%actual clustering (islands) algorithm
[~,k_info_t,Idx1t,Idx2t,Idx3t] = islands3D_LBPD_D(V_bin_1,1); %extracting binary patterns (islands)
cluster_ori_1 = zeros(size(k_info_t,1),2); %matrix to store cluster characteristics in (column 1 = cluster size; column 2 = cluster mass)
[cls,idx] = sort(k_info_t(:,2),'descend'); %getting sorted original cluster sizes (number of voxels forming the clusters)
cluster_ori_1(:,1) = cls; %storing cluster sizes
%getting global statistics of the networks (sum of voxels statistic); similar to what FSL people call cluster mass (instead of size which is the simple number of voxels forming the cluster)
for ii = 1:length(Idx1t) %over clusters
    for jj = 1:length(Idx1t{idx(ii)}) %over voxels forming the cluster ii
        cluster_ori_1(ii,2) = cluster_ori_1(ii,2) + V_t(Idx1t{idx(ii)}(jj),Idx2t{idx(ii)}(jj),Idx3t{idx(ii)}(jj)); %summing progressively the statistic of each voxel forming the cluster ii
    end
end
%sorting information before saving it
k_info_t = k_info_t(idx,:);
Idx1t = Idx1t(idx);
Idx2t = Idx2t(idx);
Idx3t = Idx3t(idx);
%creating directory if not existent and saving file
if ~exist(S_ori.outdir,'dir')
    mkdir(S_ori.outdir)
end
save([S_ori.outdir '/' S_ori.analysis_name '_origclust_cond1cond2.mat'],'cluster_ori_1','k_info_t','Idx1t','Idx2t','Idx3t')


%cond 2 > cond 1
%actual clustering (islands) algorithm
[~,k_info_t,Idx1t,Idx2t,Idx3t] = islands3D_LBPD_D(V_bin_2,1); %extracting binary patterns (islands)
cluster_ori_2 = zeros(size(k_info_t,1),2); %matrix to store cluster characteristics in (column 1 = cluster size; column 2 = cluster mass)
[cls,idx] = sort(k_info_t(:,2),'descend'); %getting sorted original cluster sizes (number of voxels forming the clusters)
cluster_ori_2(:,1) = cls; %storing cluster sizes
%getting global statistics of the networks (sum of voxels statistic); similar to what FSL people call cluster mass (instead of size which is the simple number of voxels forming the cluster)
for ii = 1:length(Idx1t) %over clusters
    for jj = 1:length(Idx1t{idx(ii)}) %over voxels forming the cluster ii
        cluster_ori_2(ii,2) = cluster_ori_2(ii,2) + V_t(Idx1t{idx(ii)}(jj),Idx2t{idx(ii)}(jj),Idx3t{idx(ii)}(jj)); %summing progressively the statistic of each voxel forming the cluster ii
    end
end
%sorting information before saving it
k_info_t = k_info_t(idx,:);
Idx1t = Idx1t(idx);
Idx2t = Idx2t(idx);
Idx3t = Idx3t(idx);

save([S_ori.outdir '/' S_ori.analysis_name '_origclust_cond2cond1.mat'],'cluster_ori_2','k_info_t','Idx1t','Idx2t','Idx3t')



end


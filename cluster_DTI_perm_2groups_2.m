function O = cluster_DTI_perm_2groups_2(S_perm)
O = []; 


%starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl);


% addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/Cluster')
perm = S_perm.perm;
V_t = S_perm.V_t;
SS = size(V_t);
VTOT = S_perm.VTOT;

%doing permutations
MAX = zeros(perm,2);
for pp = 1:perm
    %actual clustering (islands) algorithm
    ind = find(V_t ~= 0); %getting indices of non-0 voxels
    for ll = 1:S_perm.kk %desperate trick to try to prevent the apparently "same" randomization happening when splitting the jobs on the cluster
        idx_dummy = randperm(length(ind)); %create a permuted array from 1 to length(ind)
    end
    idx_2 = ind(idx_dummy); %shuffle the indexes in ind according to the order in idx_dummy
    r_dummy = zeros(SS(1)*SS(2)*SS(3),1); %initialise vector
    t_dummy = zeros(SS(1)*SS(2)*SS(3),1); %initialise vector
    r_dummy(ind) = VTOT(idx_2); %taking element in matrix data with index idx_2
    t_dummy(ind) = V_t(idx_2); %same but for t-values
    t_rand = reshape(t_dummy,[SS(1),SS(2),SS(3)]); %reshaping the vector into a matrix shaped as the original data
    V_rand = reshape(r_dummy,[SS(1),SS(2),SS(3)]); %reshaping the vector into a matrix shaped as the original data
    V_rand(V_rand==0) = NaN; %removing 0 values as done above for the actual data
    [~,k_info_r,Idx1t,Idx2t,Idx3t] = islands3D_LBPD_D(V_rand,1); %extracting binary patterns (islands)
    %cluster size
    MAX(pp,1) = max(k_info_r(:,2)); %getting maximum random cluster size
    %cluster mass
    dummy = zeros(length(Idx1t),1);
    for ii = 1:length(Idx1t) %over clusters
        for jj = 1:length(Idx1t{ii}) %over voxels forming the cluster ii
            dummy(ii,1) = dummy(ii,1) + t_rand(Idx1t{ii}(jj),Idx2t{ii}(jj),Idx3t{ii}(jj)); %summing progressively the statistic of each voxel forming the cluster ii
        end
    end
    MAX(pp,2) = max(dummy); %getting maximum random cluster mass
    disp(pp)
end
if ~exist(S_perm.outdir,'dir')
    mkdir(S_perm.outdir)
end
save([S_perm.outdir '/' S_perm.analysis_name '_clust_perm_' num2str(S_perm.kk) '.mat'],'MAX')


end


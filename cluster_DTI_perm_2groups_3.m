function OUT = cluster_DTI_perm_2groups_3(S_fin)
%combining original and permuted clusters
%this function works in connection with the previous one

% Leonardo Bonetti, Aarhus, 22/12/2020


%indexing and loading permuting clusters
a = dir([S_fin.outdir '/' S_fin.analysis_name '_clust_perm*']);
if ~isempty(a) %if there are no permuted clusters (i.e. permuting the data gave rise to no significant clusters)
    b = [];
    for ii = 1:length(a)
        load([a(ii).folder '/' a(ii).name],'MAX');
        b = cat(1,b,MAX); %concatenating maximum permuted cluster sizes and masses
    end
else %creating a vector with zeros
    b = zeros(1,S_fin.perm);
end
if length(b) ~= S_fin.perm %if some of the permutations did not gie rise to significant clusters (but some other permutations yes)
    b2 = zeros(S_fin.perm,2); %creating a dummy variable with length equal to the number of permutations
    b2(1:length(b),:) = b; %assigning vector b to it (so having a vector with values from permutations and zeros to complete it..)
    b = b2; %assigning that to b
end
%sorting cluster sizes and masses
perm_size_clust = sort(b(:,1),'descend');
perm_mass_clust = sort(b(:,2),'descend');
%if no clusters were detected during permutations you simply assign to the index threshold the end of the vector
if perm_size_clust == 0
    th_psc = 1;
else
    %getting threshold in permuted clusters to be used for original clusters
    th_psc = perm_size_clust(S_fin.thresh_mc*length(perm_size_clust));
end
if perm_mass_clust == 0
    th_pmc = 1;
else
    th_pmc = perm_mass_clust(S_fin.thresh_mc*length(perm_mass_clust));
end

%loading original clusters
a2 = dir([S_fin.outdir '/' S_fin.analysis_name '_origclust*']);
for ii = 1:length(a2)
    load([a2(ii).folder '/' a2(ii).name],['cluster_ori_' num2str(ii)]);
end

%comparing original vs permuted clusters
cluster_ori_2sorted_n = sort(cluster_ori_2(:,2)); %sorting original clusters according to cluster mass (cond2>cond1)
cluster_ori_2sorted_p = sort(cluster_ori_1(:,2),'descend'); %sorting original clusters according to cluster mass (cond1>cond2)

%cluster size (getting the indices of the closest values to thresholds
[~,ret] = min(abs(th_psc - cluster_ori_1(:,1)));
[~,retn] = min(abs(th_psc - cluster_ori_2(:,1)));
%cluster mass (getting the indices of the closest values to thresholds
[~,retm] = min(abs(th_pmc - cluster_ori_2sorted_p));
[~,retnm] = min(abs((th_pmc*-1) - cluster_ori_2sorted_n));

%cluster size cond1>cond2
if ret > 1
    clust_size_1 = cluster_ori_1(1:ret-1,1); %here ret-1 is to avoid to get also the significant clusters on the threshold.. of course if you have no clusters in the permutations, this should probably be only "ret" and not "ret-1".. anyway.. better to be a bit more strict probably so I leave it like this.. (same concept for the lines below!!)
else
    clust_size_1 = [];
end
%cluster size cond2>cond1
if retn > 1
    clust_size_2 = cluster_ori_2(1:retn-1,1);
else
    clust_size_2 = [];
end
%cluster mass cond1>cond2
if retm > 1
    clust_mass_1 = cluster_ori_2sorted_p(1:retm-1,1);
else
    clust_mass_1 = [];
end
%cluster mass cond2>cond1
if retnm > 1
    clust_mass_2 = cluster_ori_2sorted_n(1:retnm-1,1);
else
    clust_mass_2 = [];
end

%saving results
OUT = cell(3,2);
OUT(1,1) = {clust_size_1}; %first column cond1>cond2
OUT(2,1) = {clust_mass_1};
OUT(1,2) = {clust_size_2}; %second column cond2>cond1
OUT(2,2) = {clust_mass_2};
OUT(3,1) = {'1,1: cond1>cond2 size; 2,1: cond1>cond2 mass; 1,2: cond2>cond1 size; 2,2: cond2>cond1 mass'}; %sort of legend..


end


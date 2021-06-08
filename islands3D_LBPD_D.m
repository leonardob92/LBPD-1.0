function [b2b,k_info,Ind1,Ind2,Ind3] = islands3D_LBPD_D(a,infor)

% It finds all islands in a 3D binary matrix.
% An island is here defined as a cluster of same neighbouring values along the 3 dimensions.
% A single element is not considered a cluster.
% The function works ignoring NaN if any.


% INPUT:     -a:        3D matrix with data
%            -infor:    1 to show progressive cluster estimation
%                       0 not to show it
%                         this can be useful to check slow computations of very large datasets

% OUTPUT:    -b2b:      a 3D matrix of the same size of a where the found clusters are stored.
%                       Elements belonging to the same cluster are identified by the same number.
%            -k_info:   3 x n matrix withÆ
%                            -1st clm: cluster progressive number
%                            -2nd clm: cluster size
%                            -3th clm: corresponding value of a
%            -Ind1:     cells with indices for dimension 1 of each element of each cluster
%            -Ind2:     cells with indices for dimension 2 of each element of each cluster
%            -Ind3:     cells with indices for dimension 3 of each element of each cluster





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 06/09/2018
% Leonardo Bonetti, Aarhus, DK, 27/04/2020


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    








%random values for testing..
% a = rand(4,5,3);
% a(a>0.7) = 1;
% a(a<1) = 0;


S = size(a);
%creating a new matrix a2 with NaN values as a frame for the original a matrix
a2 = zeros(S(1)+2,S(2)+2,S(3)+2);
a2(1,:,:) = NaN;
a2(S(1)+2,:,:) = NaN;
a2(:,1,:) = NaN;
a2(:,S(2)+2,:) = NaN;
a2(:,:,1) = NaN;
a2(:,:,S(3)+2) = NaN;
b2 = a2; %in b2 I store the new cluster identifications
a2(2:S(1)+1,2:S(2)+1,2:S(3)+1) = a; %inscribing a in a2

%checking if information about clustering algorithm label exists
if ~exist('infor','var') %if it does not exist, assigning 0
    infor = 0;
end

%looking for islands
cnt = 0;
for ii = 2:S(1)+1
    for jj = 2:S(2)+1
        for pp = 2:S(3)+1
            %comparing one element with all the other ones (1st and 2nd dimensions)
            %checking up straight
            if a2(ii,jj,pp) == a2(ii-1,jj,pp) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii-1,jj,pp))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii-1,jj,pp) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii-1,jj,pp); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii-1,jj,pp) = b2(ii,jj,pp);
                    end
                elseif b2(ii-1,jj,pp) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii-1,jj,pp); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii-1,jj,pp) = cnt;
                end
            end
            %checking up right diagonal
            if a2(ii,jj,pp) == a2(ii-1,jj+1,pp) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii-1,jj+1,pp))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii-1,jj+1,pp) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii-1,jj+1,pp); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii-1,jj+1,pp) = b2(ii,jj,pp);
                    end
                elseif b2(ii-1,jj+1,pp) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii-1,jj+1,pp); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii-1,jj+1,pp) = cnt;
                end
            end
           %checking right straight
            if a2(ii,jj,pp) == a2(ii,jj+1,pp) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii,jj+1,pp))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii,jj+1,pp) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii,jj+1,pp); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii,jj+1,pp) = b2(ii,jj,pp);
                    end
                elseif b2(ii,jj+1,pp) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii,jj+1,pp); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii,jj+1,pp) = cnt;
                end
            end
            %checking bottom right diagonal
            if a2(ii,jj,pp) == a2(ii+1,jj+1,pp) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii+1,jj+1,pp))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii+1,jj+1,pp) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii+1,jj+1,pp); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii+1,jj+1,pp) = b2(ii,jj,pp);
                    end
                elseif b2(ii+1,jj+1,pp) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii+1,jj+1,pp); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii+1,jj+1,pp) = cnt;
                end
            end
            %checking bottom straight
            if a2(ii,jj,pp) == a2(ii+1,jj,pp) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii+1,jj,pp))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii+1,jj,pp) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii+1,jj,pp); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii+1,jj,pp) = b2(ii,jj,pp);
                    end
                elseif b2(ii+1,jj,pp) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii+1,jj,pp); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii+1,jj,pp) = cnt;
                end
            end
            %checking bottom left diagonal
            if a2(ii,jj,pp) == a2(ii+1,jj-1,pp) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii+1,jj-1,pp))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii+1,jj-1,pp) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii+1,jj-1,pp); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii+1,jj-1,pp) = b2(ii,jj,pp);
                    end
                elseif b2(ii+1,jj-1,pp) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii+1,jj-1,pp); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii+1,jj-1,pp) = cnt;
                end
            end
            %checking left straight
            if a2(ii,jj,pp) == a2(ii,jj-1,pp) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii,jj-1,pp))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii,jj-1,pp) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii,jj-1,pp); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii,jj-1,pp) = b2(ii,jj,pp);
                    end
                elseif b2(ii,jj-1,pp) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii,jj-1,pp); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii,jj-1,pp) = cnt;
                end
            end
            %checking up left diagonal
            if a2(ii,jj,pp) == a2(ii-1,jj-1,pp) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii-1,jj-1,pp))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii-1,jj-1,pp) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii-1,jj-1,pp); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii-1,jj-1,pp) = b2(ii,jj,pp);
                    end
                elseif b2(ii-1,jj-1,pp) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii-1,jj-1,pp); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii-1,jj-1,pp) = cnt;
                end
            end
            
            %comparing one element with all the other ones (1st, 2nd and 3d dimensions (pp+1))
            %checking up straight
            if a2(ii,jj,pp) == a2(ii-1,jj,pp+1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii-1,jj,pp+1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii-1,jj,pp+1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii-1,jj,pp+1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii-1,jj,pp+1) = b2(ii,jj,pp);
                    end
                elseif b2(ii-1,jj,pp+1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii-1,jj,pp+1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii-1,jj,pp+1) = cnt;
                end
            end
            %checking up right diagonal
            if a2(ii,jj,pp) == a2(ii-1,jj+1,pp+1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii-1,jj+1,pp+1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii-1,jj+1,pp+1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii-1,jj+1,pp+1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii-1,jj+1,pp+1) = b2(ii,jj,pp);
                    end
                elseif b2(ii-1,jj+1,pp+1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii-1,jj+1,pp+1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii-1,jj+1,pp+1) = cnt;
                end
            end
           %checking right straight
            if a2(ii,jj,pp) == a2(ii,jj+1,pp+1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii,jj+1,pp+1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii,jj+1,pp+1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii,jj+1,pp+1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii,jj+1,pp+1) = b2(ii,jj,pp);
                    end
                elseif b2(ii,jj+1,pp+1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii,jj+1,pp+1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii,jj+1,pp+1) = cnt;
                end
            end
            %checking bottom right diagonal
            if a2(ii,jj,pp) == a2(ii+1,jj+1,pp+1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii+1,jj+1,pp+1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii+1,jj+1,pp+1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii+1,jj+1,pp+1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii+1,jj+1,pp+1) = b2(ii,jj,pp);
                    end
                elseif b2(ii+1,jj+1,pp+1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii+1,jj+1,pp+1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii+1,jj+1,pp+1) = cnt;
                end
            end
            %checking bottom straight
            if a2(ii,jj,pp) == a2(ii+1,jj,pp+1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii+1,jj,pp+1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii+1,jj,pp+1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii+1,jj,pp+1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii+1,jj,pp+1) = b2(ii,jj,pp);
                    end
                elseif b2(ii+1,jj,pp+1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii+1,jj,pp+1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii+1,jj,pp+1) = cnt;
                end
            end
            %checking bottom left diagonal
            if a2(ii,jj,pp) == a2(ii+1,jj-1,pp+1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii+1,jj-1,pp+1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii+1,jj-1,pp+1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii+1,jj-1,pp+1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii+1,jj-1,pp+1) = b2(ii,jj,pp);
                    end
                elseif b2(ii+1,jj-1,pp+1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii+1,jj-1,pp+1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii+1,jj-1,pp+1) = cnt;
                end
            end
            %checking left straight
            if a2(ii,jj,pp) == a2(ii,jj-1,pp+1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii,jj-1,pp+1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii,jj-1,pp+1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii,jj-1,pp+1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii,jj-1,pp+1) = b2(ii,jj,pp);
                    end
                elseif b2(ii,jj-1,pp+1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii,jj-1,pp+1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii,jj-1,pp+1) = cnt;
                end
            end
            %checking up left diagonal
            if a2(ii,jj,pp) == a2(ii-1,jj-1,pp+1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii-1,jj-1,pp+1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii-1,jj-1,pp+1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii-1,jj-1,pp+1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii-1,jj-1,pp+1) = b2(ii,jj,pp);
                    end
                elseif b2(ii-1,jj-1,pp+1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii-1,jj-1,pp+1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii-1,jj-1,pp+1) = cnt;
                end
            end
            %checking usame element but in 3th dimension pp+1
            if a2(ii,jj,pp) == a2(ii,jj,pp+1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii,jj,pp+1)) 
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii,jj,pp+1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii,jj,pp+1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii,jj,pp+1) = b2(ii,jj,pp);
                    end
                elseif b2(ii,jj,pp+1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii,jj,pp+1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii,jj,pp+1) = cnt;
                end
            end
            
            
            %comparing one element with all the other ones (1st, 2nd and 3d dimensions (pp-1))
            %checking up straight
            if a2(ii,jj,pp) == a2(ii-1,jj,pp-1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii-1,jj,pp-1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii-1,jj,pp-1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii-1,jj,pp-1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii-1,jj,pp-1) = b2(ii,jj,pp);
                    end
                elseif b2(ii-1,jj,pp-1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii-1,jj,pp-1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii-1,jj,pp-1) = cnt;
                end
            end
            %checking up right diagonal
            if a2(ii,jj,pp) == a2(ii-1,jj+1,pp-1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii-1,jj+1,pp-1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii-1,jj+1,pp-1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii-1,jj+1,pp-1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii-1,jj+1,pp-1) = b2(ii,jj,pp);
                    end
                elseif b2(ii-1,jj+1,pp-1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii-1,jj+1,pp-1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii-1,jj+1,pp-1) = cnt;
                end
            end
           %checking right straight
            if a2(ii,jj,pp) == a2(ii,jj+1,pp-1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii,jj+1,pp-1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii,jj+1,pp-1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii,jj+1,pp-1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii,jj+1,pp-1) = b2(ii,jj,pp);
                    end
                elseif b2(ii,jj+1,pp-1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii,jj+1,pp-1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii,jj+1,pp-1) = cnt;
                end
            end
            %checking bottom right diagonal
            if a2(ii,jj,pp) == a2(ii+1,jj+1,pp-1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii+1,jj+1,pp-1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii+1,jj+1,pp-1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii+1,jj+1,pp-1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii+1,jj+1,pp-1) = b2(ii,jj,pp);
                    end
                elseif b2(ii+1,jj+1,pp-1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii+1,jj+1,pp-1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii+1,jj+1,pp-1) = cnt;
                end
            end
            %checking bottom straight
            if a2(ii,jj,pp) == a2(ii+1,jj,pp-1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii+1,jj,pp-1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii+1,jj,pp-1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii+1,jj,pp-1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii+1,jj,pp-1) = b2(ii,jj,pp);
                    end
                elseif b2(ii+1,jj,pp-1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii+1,jj,pp-1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii+1,jj,pp-1) = cnt;
                end
            end
            %checking bottom left diagonal
            if a2(ii,jj,pp) == a2(ii+1,jj-1,pp-1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii+1,jj-1,pp-1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii+1,jj-1,pp-1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii+1,jj-1,pp-1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii+1,jj-1,pp-1) = b2(ii,jj,pp);
                    end
                elseif b2(ii+1,jj-1,pp-1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii+1,jj-1,pp-1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii+1,jj-1,pp-1) = cnt;
                end
            end
            %checking left straight
            if a2(ii,jj,pp) == a2(ii,jj-1,pp-1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii,jj-1,pp-1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii,jj-1,pp-1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii,jj-1,pp-1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii,jj-1,pp-1) = b2(ii,jj,pp);
                    end
                elseif b2(ii,jj-1,pp-1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii,jj-1,pp-1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii,jj-1,pp-1) = cnt;
                end
            end
            %checking up left diagonal
            if a2(ii,jj,pp) == a2(ii-1,jj-1,pp-1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii-1,jj-1,pp-1))
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii-1,jj-1,pp-1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii-1,jj-1,pp-1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii-1,jj-1,pp-1) = b2(ii,jj,pp);
                    end
                elseif b2(ii-1,jj-1,pp-1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii-1,jj-1,pp-1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii-1,jj-1,pp-1) = cnt;
                end
            end
            %checking usame element but in 3th dimension pp+1
            if a2(ii,jj,pp) == a2(ii,jj,pp-1) && ~isnan(a2(ii,jj,pp)) && ~isnan(a2(ii,jj,pp-1)) 
                if b2(ii,jj,pp) ~= 0 %if 1st voxel ~= 0
                    if b2(ii,jj,pp-1) ~= 0 %if 2nd voxel ~= 0
                        dum = b2(ii,jj,pp-1); %get cluster ID of 2nd voxel
                        b2(b2==dum) = b2(ii,jj,pp); %and overwrite on those voxels cluster ID of 1st voxel
                    else %else simply write cluster ID of voxel 1 in voxel 2
                        b2(ii,jj,pp-1) = b2(ii,jj,pp);
                    end
                elseif b2(ii,jj,pp-1) ~= 0 %else, if 1st voxel == 0 and 2nd voxel ~= 0
                    b2(ii,jj,pp) = b2(ii,jj,pp-1); %1st voxel has no ID, so it gets 2nd voxel ID
                else %else, if both voxels == 0, simply create a new cluster ID
                    cnt = cnt + 1;
                    b2(ii,jj,pp) = cnt; %and assign it to the 2 voxels
                    b2(ii,jj,pp-1) = cnt;
                end
            end
        end
    end
    if infor == 1
        disp([num2str(ii) '/' num2str(S(1)+1)])
    end
end

%preparing output
b2b = b2(2:S(1)+1,2:S(2)+1,2:S(3)+1); %extracting real values outside NaN frame..
%finding cluster elements in the b2 matrix and store cluster information
countk = 0;
if max(max(max(b2b))) > 0 %if there was at least 1 cluster (maybe this check is not really necessary..)
    for kkk = 1:max(max(max(b2b)))
        [I1,I2,I3] = findND(b2b == kkk); %getting indices of cluster kkk
        if ~isempty(I1)
            countk = countk + 1;
            k_info(countk,1) = countk; %progressive number of the cluster
            k_info(countk,2) = length(I1); %size of the cluster
            k_info(countk,3) = a(I1(1),I2(1),I3(1)); %original number in matrix a associated to cluster kkk
            Ind1(countk,1) = {I1}; %stores the indexes of each cluster element
            Ind2(countk,1) = {I2};
            Ind3(countk,1) = {I3};
        end
    end
else
    k_info = 0;
    Ind1 = 0;
    Ind2 = 0;
    Ind3 = 0;
    
end

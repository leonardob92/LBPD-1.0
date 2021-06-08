function [ ] = GT_modul_plot_LBPD( S )

% It plots intra and inter subnetworks connectivity (after finding ideal
% subnetworks; e.g. using modularity algorithms).



%  INPUT:   -S.intra_con:          1 for connectivity intra subnetwork; 0 otherwise
%           -S.inter_con:          1 for connectivity inter subnetwork; 0 otherwise
%           -S.perc_intra:         percentage of intra subnetwork connections to be plotted (e.g. 10 means 10%)
%           -S.perc_inter:         percentage of inter subnetwork connections to be plotted (e.g. 10 means 10%)
%           -S.communities:        brain regions (ROIs) vector with community ID for each ROI
%           -S.MNI_RC:             ROIs x 3D coordinates matrix with coordinates of the centroid of all ROIs
%           -S.conn_matrix:        ROIs x ROIs undirected connectivity matrix
%           -S.colors_mod:         vector with colors for connections of subnetworks (e.g. ['r','k','g'] or RGB triplets, e.g. [0 0.6 1; 1 0.4 0; 1 0 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0])
%           -S.col_dot:            color for centroids of brain regions (e.g. 'b')
%           -S.title:              name to be added to the title of the  images (character)
%           -S.limits:             limits to be used for rescaling strength of the connections for plotting purposes (e.g.  (suggested) [0.1 5]).
%                                  If you provide 1 number only, there is no rescaling of the relative strength of the different connections.

%  OUTPUT:  -plotting solutions






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 01/04/2021


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    










%getting inputs
intra_con = S.intra_con;
inter_con = S.inter_con;
perc_intra = S.perc_intra;
perc_inter = S.perc_inter;
MNI_RC = S.MNI_RC;
com1 = S.communities;
colors = S.colors_mod;
col_dot = S.col_dot;
limitt = S.limits;
MM1 = S.conn_matrix;
SS = size(MM1);

%actual computations
thresh2 = round((((SS(1)*SS(1))-SS(1))/2/100)*perc_intra);
D12 = sort(reshape(MM1,[SS(1)*SS(1),1]),'descend');
rmin = D12(thresh2);
% rmin = min(min(r(r~=0)));
rmax = max(max(MM1));
%actual computation
%intra subnetwork connectivity
if intra_con == 1
    %plotting brain template and AAL nodes (centroids of areas)
    openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
    hold on
    for ii = 1:90 %over ROIs
        plot3(MNI_RC(ii, 1), MNI_RC(ii, 2), MNI_RC(ii, 3), ['.'], 'Color', col_dot, 'MarkerSize', 10); %centroid of ROI ii
        hold on
    end
    rotate3d on; axis off; axis vis3d; axis equal
    for ii = 1:max(com1) %over communities
        Mfak = zeros(90);
        Mfak(com1==ii,com1==ii) = MM1(com1==ii,com1==ii); %getting only connections intra-subnetwork ii
        Mfakt = triu(Mfak,1); %upper triangle
        ind = find(Mfakt~=0);
        a = sort(Mfakt(ind),'descend');
        if round((length(a)/100)*perc_intra) == 0 %the threshold must take at least the first element (if it is 0, it of course crashes..)
            thresh = a(round((length(a)/100)*perc_intra)+1);
        else
            thresh = a(round((length(a)/100)*perc_intra));
        end
        Mfakt(Mfakt<thresh) = 0; %zeroing elements below the threshold
        ind = find(Mfakt~=0); %getting indices of non-zero elements
        [row,col] = ind2sub([90 90],ind); %converting indices in 2D
        %if requested, scaling
        if length(limitt) == 2
            [a2,i2] = sort(Mfakt(ind),'descend'); %sorting only elements that survived the threshold
            row = row(i2); %sorting row from highest to lowest connection value
            col = col(i2); %sorting col from highest to lowest connection value
            %scaling on the basis of connections strength
%             rmin = a2(end);
%             rmax = a2(1);
            %scaling connection strengths with extreme values (limitt) requested by user
            size_con = (a2-rmin)./(rmax-rmin).*(limitt(2)-limitt(1)) + limitt(1);
        else
            size_con = ones(length(row),1).*limitt(1);
        end
        %the following line is a reasonable solution that, however, must be better fixed in the future
        size_con(size_con<0) = 0.0000000001; %barbaric trick to avoid the problem that the scaling may return negative values in some cases..
        for jj = 1:length(row) %over ROI jj of community ii
%             for zz = 1:length(col) %over ROI zz of community ii
                vdum = [MNI_RC(row(jj),:); MNI_RC(col(jj),:)]; %vector with the two selected nodes (jj and zz)
                hold on
                plot3(vdum(:,1),vdum(:,2),vdum(:,3),'Color',colors(ii,:),'LineWidth',size_con(jj))  %plotting line connecting the two points
%             end
            disp(['intra-subnetwork - community ' num2str(ii) ' - connection ' num2str(jj) ' / ' num2str(length(row))])
        end
        
        %previous wrong solution..
%         for jj = 1:length(row) %over ROI jj of community ii
%             for zz = 1:length(col) %over ROI zz of community ii
%                 vdum = [MNI_RC(row(jj),:); MNI_RC(col(zz),:)]; %vector with the two selected nodes (jj and zz)
%                 hold on
%                 plot3(vdum(:,1),vdum(:,2),vdum(:,3),colors(ii),'MarkerSize', size_con(zz))  %plotting line connecting the two points
%             end
%             disp(['intra-subnetwork - community ' num2str(ii) ' - connection ' num2str(jj) ' / ' num2str(length(row))])
%         end
    end
    set(gcf,'color','w')
    title(['intra subnetwork - ' S.title])
end

%inter subnetwork connectivity
% thresh = round((4005/100)*perc_inter);
% D12 = sort(reshape(MM1,[8100,1]),'descend');
% rmin = D12(thresh);
% rmin = min(min(r(r~=0)));
% rmax = max(max(MM1));
if inter_con == 1
    %plotting brain template and AAL nodes (centroids of areas)
    openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
    hold on
%     for ii = 1:90 %over ROIs
%         plot3(MNI_RC(ii, 1), MNI_RC(ii, 2), MNI_RC(ii, 3), '.b', 'MarkerSize', 10); %centroid of ROI ii
%         hold on
%     end
    rotate3d on; axis off; axis vis3d; axis equal
    for ii = 1:max(com1) %over communities
        vac = find(com1==ii);
        for iii = 1:length(vac) %over ROIs
            plot3(MNI_RC(vac(iii), 1), MNI_RC(vac(iii), 2), MNI_RC(vac(iii), 3),'.','Color',colors(ii,:), 'MarkerSize', 15); %centroid of ROI ii
            hold on
        end
        
        Mfak = zeros(90);
        Mfak(com1==ii,com1~=ii) = MM1(com1==ii,com1~=ii); %getting only connections intra-subnetwork ii
        %         Mfakt = triu(Mfak,1); %upper triangle
        Mfakt = Mfak; %we do not need the upper triangle here.. I already index only one connection (i.e. only a(2,3) and not also a(3,2))
        ind = find(Mfakt~=0);
        a = sort(Mfakt(ind),'descend');
        if round((length(a)/100)*perc_inter) == 0 %the threshold must take at least the first element (if it is 0, it of course crashes..)
            thresh = a(round((length(a)/100)*perc_inter)+1);
        else
            thresh = a(round((length(a)/100)*perc_inter));
        end
        Mfakt(Mfakt<thresh) = 0; %zeroing elements below the threshold
        ind = find(Mfakt~=0); %getting indices of non-zero elements
        [row,col] = ind2sub([90 90],ind); %converting indices in 2D
        %if requested, scaling
        if length(limitt) == 2
            [a2,i2] = sort(Mfakt(ind),'descend'); %sorting only elements that survived the threshold
            row = row(i2); %sorting row from highest to lowest connection value
            col = col(i2); %sorting col from highest to lowest connection value
            %scaling on the basis of connections strength
            rmin = a2(end);
            rmax = max(max(MM1));
%             rmax = a2(1);
            %scaling connection strengths with extreme values (limitt) requested by user
            size_con = (a2-rmin)./(rmax-rmin).*(limitt(2)-limitt(1)) + limitt(1);
        else
            size_con = ones(length(row),1).*limitt(1);
        end
        %the following line is a reasonable solution that, however, must be better fixed in the future
        size_con(size_con<0) = 0.0000000001; %barbaric trick to avoid the problem that the scaling may return negative values in some cases..
        for jj = 1:length(row) %over ROI jj of community ii
%             for zz = 1:length(col) %over ROI zz of community ii
                vdum = [MNI_RC(row(jj),:); MNI_RC(col(jj),:)]; %vector with the two selected nodes (jj and zz)
                hold on
                plot3(vdum(:,1),vdum(:,2),vdum(:,3),'Color',[0.4 0.4 0.4],'LineWidth',size_con(jj))  %plotting line connecting the two points
%                 plot3(vdum(:,1),vdum(:,2),vdum(:,3),'Color',colors(ii,:),'LineWidth',size_con(jj))  %plotting line connecting the two points
%             end
            disp(['inter-subnetwork - community ' num2str(ii) ' - connection ' num2str(jj) ' / ' num2str(length(row))])
        end
        
        %previous wrong solution..
%         for jj = 1:length(row) %over ROI jj of community ii
%             for zz = 1:length(col) %over ROI zz of community ii
%                 vdum = [MNI_RC(row(jj),:); MNI_RC(col(zz),:)]; %vector with the two selected nodes (jj and zz)
%                 hold on
%                 plot3(vdum(:,1),vdum(:,2),vdum(:,3),colors(ii),'MarkerSize', size_con(zz))  %plotting line connecting the two points
%             end
%             disp(['intra-subnetwork - community ' num2str(ii) ' - connection ' num2str(jj) ' / ' num2str(length(row))])
%         end
    end
    set(gcf,'color','w')
    title(['inter subnetwork - ' S.title])
end


end


function h = schemaball_modularity_LBPD(r, lbls, modul, COL, perc_intra, perc_inter, limitt)

% Plots modularity of graphs (such as brain networks) in schemaball.
% This work is based on the brilliant function "schemaball" created by Oleg
% Komarov (oleg.komarov@hotmail.it), even though the codes here have been
% considerably modified.
% I suggest to plot only the strongest connections (e.g. 50% of
% intra-subnetworks connectivity and 3% of inter-subnetworks connectivity).
% The codes could be better optimized, however they are reasonably fast if
% you do not plot all connections.



%  INPUT:   -r:                    connectivity square matrix (undirected)
%           -lbls:                 M x N character matrix with custom labels M labels for each node of the network
%           -modul:                M vector with communities (e.g. after modularity) for each node
%           -COL:                  n. communities x triplet RGB color, to plot intra-subnetworks connectivity
%           -perc_intra:           percentage of intra-subnetwork connectivity to be plotted (e.g. 10 = 10%)
%           -perc_inter:           percentage of inter-subnetwork connectivity to be plotted (e.g. 10 = 10%)
%           -limitt:               limit for line width in schemaball (e.g. [1 3]).
%                                  A single value (e.g. [1]) plots lines with the same width.

%  OUTPUT:  -schamaball plot





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 09/04/2021


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 






% Setting a few parameters to be used later
%Oleg's codes
%Nodes edge color
ecolor = [.25 .103922 .012745];
%Text color
tcolor = [0 0 0]; %leonardo (black color)
%preparing values for later Bezier curve
t = (0.025: 0.05 :1)';
t2 = [1-t, t].^2;
ncolor2 = COL; %getting colors
if ischar(lbls)
    lbls = cellstr(lbls);
end
sz = size(r);
% Use tau http://tauday.com/tau-manifesto
tau = 2*pi;
% Positions of nodes on the circle starting from (0,-1), useful later for label orientation
step = tau/sz(1);
theta = -.25*tau : step : .75*tau - step;
% Get cartesian x-y coordinates of the nodes
x = cos(theta);
y = sin(theta);


%leonardo's codes
%sorting nodes (and correspondent colors) according to their community (i.e. nodes belonging to community 1 are now placed close to each other in the circle, etc.)
ROIs_sorted = [];
NCUL = [];
% rdum = triu(r,1); %upper triangle
for ii = 1:max(modul) %over communities
    ok = find(modul==ii); %getting indices of community ii
    [Z,idumdeg] = sort(sum(r(modul==ii,modul==ii)),'descend'); %computing degree in the community ii and sorting the ROIs accordingly
    ok = ok(idumdeg); %sorting indices of ROIs in community ii according to their degree (first ROIs in vector 'ok' have stronger degree in community ii)
    ROIs_sorted = cat(1,ROIs_sorted,ok); %concatenating indices to be plotted later
    Z = (Z-min(Z)+0.01)/(max(Z)-min(Z)+0.01); Z = Z'; %scaling degree of ROIs internal to community ii
%         ncul = hsv2rgb([repmat(ncolor2(ii,1:2), length(Z),1) Z*ncolor2(ii,3)]); %scaling and sorting colors
    ncul = repmat(ncolor2(ii,1:3), length(Z),1); %simply passing colors..
    NCUL = cat(1,NCUL,ncul); %concatenating colors to be plotted later
    %     jk = find(ncolor2(ii,:)~=0);
end
%generating figure
figure
s.s = scatter(x,y,[], NCUL,'fill','MarkerEdgeColor',ecolor,'LineWidth',1); %plotting nodes (brain ROIs) (as dots)

%Oleg's codes
% PLACE TEXT LABELS such that you always read 'left to right'
lbls2 = lbls(ROIs_sorted); %sorting labels
ipos       = x > 0;
s.t        = zeros(sz(1),1);
s.t( ipos) = text(x( ipos)*1.1, y( ipos)*1.1, lbls2( ipos),'Color',tcolor,'fontsize',6); %leonardo added the 'fontsize' setting
set(s.t( ipos),{'Rotation'}, num2cell(theta(ipos)'/tau*360))
s.t(~ipos) = text(x(~ipos)*1.1, y(~ipos)*1.1, lbls2(~ipos),'Color',tcolor,'fontsize',6); %leonardo added the 'fontsize' setting
set(s.t(~ipos),{'Rotation'}, num2cell(theta(~ipos)'/tau*360 - 180),'Horiz','right')
% ADJUST FIGURE height width to fit text labels (same codes of Oleg)
xtn        = cell2mat(get(s.t,'extent'));
post       = cell2mat(get(s.t,'pos'));
sg         = sign(post(:,2));
posfa      = cell2mat(get([gcf gca],'pos'));
% Calculate xlim and ylim in data units as x (y) position + extension along x (y)
ylims      = post(:,2) + xtn(:,4).*sg;
ylims      = [min(ylims), max(ylims)];
xlims      = post(:,1) + xtn(:,3).*sg;
xlims      = [min(xlims), max(xlims)];
% Stretch figure
posfa(1,3) = (( diff(xlims)/2 - 1)*posfa(2,3) + 1) * posfa(1,3);
posfa(1,4) = (( diff(ylims)/2 - 1)*posfa(2,4) + 1) * posfa(1,4);
% Position it a bit lower (movegui slow)
posfa(1,2) = 100;
% Axis settings
set(gca, 'Xlim',xlims,'Ylim',ylims, 'color', 'w','XColor','none','YColor','none',...
    'clim',[-1,1])
set(gcf, 'pos' ,posfa(1,:),'Visible','on')
axis equal

%leonardo's codes
%sorting matrix with connectivity data and communities indices
R = r(ROIs_sorted,ROIs_sorted); %connectivity matrix
modul2 = modul(ROIs_sorted); %vector with communities
%maximum and minimum of original matrix, for later scaling purposes (this will provide a decent estimation of the values to be used for scaling the plots)
% thresh = a(round((length(a)/100)*perc_intra))
thresh = round((((sz(1)*sz(1))-sz(1))/2/100)*perc_intra); %getting threshold index (from percentage)
D12 = sort(reshape(r,[sz(1)*sz(1),1]),'descend'); %sorting connectivity matrix reshaped as vector
rmin = D12(thresh); %getting minimum value (corresponding to threshold)
% rmin = min(min(r(r~=0)));
rmax = max(max(r)); %getting maximum value
%computing and plotting Bezier curve
for ii = 1:max(modul2) %over communities
    %intra-subnetworks connectivity
    rdum = zeros(90);
    rdum(modul2==ii,modul2==ii) = R(modul2==ii,modul2==ii); %getting square connectivity matrix of community ii
    Mfakt = triu(rdum,1); %upper triangle
    ind = find(Mfakt~=0); %only value in the upper triangle
    %computing the threshold for zeroing smaller connections that the user does not want to be plotted 
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
%         rmin = a2(end);
%         rmax = a2(1);
%         AZZ(ii) = rmax;
        %scaling connection strengths with extreme values (limitt) requested by user
        size_con = (a2-rmin)./(rmax-rmin).*(limitt(2)-limitt(1)) + limitt(1);
    else
        size_con = ones(length(row),1).*limitt(1);
    end
    %the following line is a reasonable solution that, however, must be better fixed in the future
    size_con(size_con<0) = 0.0000000001; %barbaric trick to avoid the problem that the scaling may return negative values in some cases..
    for jj = 1:length(row) %over ROI jj of community ii        
        %computing and plotting Bezier curve
        Bx = t2*[x(col(jj)); x(row(jj))];
        By = t2*[y(col(jj)); y(row(jj))];
        hold on
        plot(Bx(:),By(:),'Color',ncolor2(ii,:),'LineWidth',size_con(jj))
        disp(['intra-subnetwork - community ' num2str(ii) ' - connection ' num2str(jj) ' / ' num2str(length(row))])
    end
end
%INTER-SUBNETWORK CONNECTIVITY
%it could probably be better optimized.. anyway..
for ii = 1:max(modul2) %over communities
    %here same concept used for intra-subnetwork connectivity
    rdum = zeros(90);
    rdum(modul2==ii,modul2~=ii) = R(modul2==ii,modul2~=ii); %getting inter-subnetworks connectivity
%     Mfakt = triu(rdum,1); %upper triangle
    Mfakt = rdum;
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
        rmin = a2(end); %using minimum of inter-subnetworks values, since the intra-subnetworks values are too bigger
        %Thus, inter and intra connectivity should be visually compared with caution..
        rmax = max(max(r));
%         rmax = a2(1);
        %scaling connection strengths with extreme values (limitt) requested by user
        size_con = (a2-rmin)./(rmax-rmin).*(limitt(2)-limitt(1)) + limitt(1);
    else
        size_con = ones(length(row),1).*limitt(1);
    end
    %the following line is a reasonable solution that, however, must be better fixed in the future
    size_con(size_con<0) = 0.0000000001; %barbaric trick to avoid the problem that the scaling may return negative values in some cases..
    for jj = 1:length(row) %over ROI jj of community ii
        %computing and plotting Bezier curve
        Bx = t2*[x(col(jj)); x(row(jj))];
        By = t2*[y(col(jj)); y(row(jj))];
        hold on
        plot(Bx(:),By(:),'Color',[0.4 0.4 0.4],'LineWidth',size_con(jj))
        disp(['inter-subnetwork - community ' num2str(ii) ' - connection ' num2str(jj) ' / ' num2str(length(row))])
    end
end

end
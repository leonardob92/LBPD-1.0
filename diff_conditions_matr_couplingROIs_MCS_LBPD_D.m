function [ P_pos_fin, P_neg_fin, PM ] = diff_conditions_matr_couplingROIs_MCS_LBPD_D( a, b, ROI_lab, X, permut, threshMC, plot_label )

% It calculates the difference between two conditions of strength of the significant
% coupling between each couple of ROIs. Then, it calculates Monte Carlo
% simulations to assess the probability that a specific difference was
% randomly occurred. Here I think that the p-values obtained from MCS
% should be divided according to how many differences you are looking at.
% This is because the single p-value is showing the probability after MCS of only one
% difference given the input data. Therefore, a good idea would be dviding
% your alpha level by the number of differences between ROIs that you are
% looking at. E.g. if you look at whole-brain level in 90-ROI AAL
% non-cerebellar parcellation, assuming an aplha level = 0.05, you should
% probably adjust it as follows: 0.05/4005 = 1.2e-05. In this case it would
% be recommended to run something like 10000 permutations (in other words a number
% of permutations that is consistent with the adjusted alpha level). This looks very
% strict but from my tests it may be a more reliable approach.



%   INPUT:  -a:          matrix n x n x t, where n are ROIs (square matrix)
%                        and t are time-points of condition 1 (double)
%           -b:          matrix n x n x t, where n are ROIs (square matrix)
%                        and t are time-points of condition 2 (double)
%                        Currently, a and b must have the same amount of
%                        time-points.
%           -ROI_lab:    labels of the ROIs corresponding to a and b
%                        values (cell)
%           -X:          -set 0 for working on all time-points (double)
%                        -set the time-window you want to work on (e.g. X = 1:180)
%           -permut:     number of permutations (double)
%           -threshMC:   threshold for MCS (5% must be expressed as 0.05)
%           -plot_label: -0 for computation only
%                        -1 for plotting only
%                        -2 for both computation and plotting
%                        Remember that if you want to use this plotting you
%                        should give the matrix reshaped in the order
%                        LLLRRR and not LRLRLR.

%   OUTPUT: -P_pos_fin:  cell with the two ROIs significantly coupled
%                        and the corresponding p-value (positive means that
%                        these specific couplings occurred more in the condition
%                        1 than in the condition 2)
%           -P_neg_fin:  basically the same, but with coupling occurring
%                        more in condition 2 than in condition 1
%           -PM:         cell with matrices of sums for condition 1, 2 and their difference





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 10/01/2019
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 25/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%getting which time-points the user requested
if X == 0
    P1 = a;
    P2 = b;
else
    P1 = a(:,:,X);
    P2 = b(:,:,X);
end
%calculating sum over time
P1sum = sum(P1,3);
P2sum = sum(P2,3);
%working only with the upper triangle
P1sum = triu(P1sum,1);
P2sum = triu(P2sum,1);
%difference between conditions
P_diff = P1sum - P2sum;
%here I assign the number of total permutations to every element of the matrix
[row,col,prof] = size(P1);
P_pos_d = ones(row)*permut; %here I assign the number of total permutations to every element of the matrix
P_neg_d = ones(row)*permut;

%computing permutations
if plot_label == 0 || plot_label == 2
    for pp = 1:permut
        disp(['computing permutation number ' num2str(pp)])
        %working randomizing only the upper triangle of each matrix over time
        %cond1
        f = find(triu(ones(row),1) == 1); %find indexes of elements of the upper triangle of the square matrix
        countf = -1*(row*row);
        F = zeros(length(f)*prof,1);
        countlenf = -length(f);
        for ii = 1:prof %calculating indices of the upper triangle of the square matrix over time (in the 3th dimension)
            countf = countf + row*row;
            countlenf = countlenf + length(f);
            F(countlenf+1:countlenf+length(f)) = f + countf; %in F are stored the indexes of the upper triangle of matrices corresponding to each time-point
        end
        idx_dummy = randperm(length(F)); %create a permuted array from 1 to length(F)
        idx_2 = F(idx_dummy); %shuffle the indexes in F according to the order in idx_dummy
        r1_dummy = zeros(row*row*prof,1); %preallocation of space
        r1_dummy(F) = P1(idx_2); %putting values of P1 in r1_dummy according to the randomized order if idx_2 (randomized only taking into account the indexes of the upper triangle of each matrix)
        [I_l,J_l,K_l] = ind2sub([row,col,prof],1:length(r1_dummy)); %convert indexes from 1D array to 3D matrix
        %rebuild 3D matrix with actual data (only upper triangle)
        r1 = zeros(row,col,prof);
        count = 0;
        for ii = 1:length(I_l)
            count = count + 1;
            r1(I_l(ii),J_l(ii),K_l(ii)) = r1_dummy(count);
        end
        r1sum = sum(r1,3); %sum over time
        %cond2
        f = find(triu(ones(row),1) == 1); %find indexes of elements of the upper triangle of the square matrix
        countf = -1*(row*row);
        F = zeros(length(f)*prof,1);
        countlenf = -length(f);
        for ii = 1:prof %calculating indices of the upper triangle of the square matrix over time (in the 3th dimension)
            countf = countf + row*row;
            countlenf = countlenf + length(f);
            F(countlenf+1:countlenf+length(f)) = f + countf; %in F are stored the indexes of the upper triangle of matrices corresponding to each time-point
        end
        idx_dummy = randperm(length(F)); %create a permuted array from 1 to length(F)
        idx_2 = F(idx_dummy); %shuffle the indexes in F according to the order in idx_dummy
        r2_dummy = zeros(row*row*prof,1); %preallocation of space
        r2_dummy(F) = P2(idx_2); %putting values of P1 in r1_dummy according to the randomized order if idx_2 (randomized only taking into account the indexes of the upper triangle of each matrix)
        [I_l,J_l,K_l] = ind2sub([row,col,prof],1:length(r2_dummy)); %convert indexes from 1D array to 3D matrix
        %rebuild 3D matrix with actual data (only upper triangle)
        r2 = zeros(row,col,prof);
        count = 0;
        for ii = 1:length(I_l)
            count = count + 1;
            r2(I_l(ii),J_l(ii),K_l(ii)) = r2_dummy(count);
        end
        r2sum = sum(r2,3); %sum over time
        r_diff_max = max(max(r1sum - r2sum)); %here I take the max and the min for a problem of memory.. otherwise I would need to save too many numbers.. I think that would be better to generate the reference distribution using all numbers, however also this is a good enough solution I think, also because this should be even more stringent and conservative than using the global distribution
        r_diff_min = min(min(r1sum - r2sum)); %note that the distributions of maxs and mins do not look as a smooth normal distribution, but they should work fine
        [I,J] = find(P_diff > r_diff_max); %finding when the original data was higher than the randomized one
        for ii = 1:length(I) %this loop may be removed, but I left it for clarity..
            P_pos_d(I(ii),J(ii)) = P_pos_d(I(ii),J(ii)) - 1; %here the idea is subtracting the significant elements to themselves.. and then dividing them by the total number of permutations
        end
        [I,J] = find(P_diff < r_diff_min);
        for ii = 1:length(I) %the same here..
            P_neg_d(I(ii),J(ii)) = P_neg_d(I(ii),J(ii)) - 1;
        end
    end
    %testing and storing significance
    %pos
    [I,J] = find(P_pos_d/permut < threshMC/2); %dividing by the total number of permutations and therefore obtaining final P-values (in this case getting already the indexes for writing a more readable output); here threshMC is divided by 2 since I compute this operation 2 times (for positive and negative differences..)
    P_pos_fin = cell(length(I),4);
    %sorting significant results
    dum = zeros(length(I),1);
    for ii = 1:length(I)
        dum(ii) = P_pos_d(I(ii),J(ii))/permut;
    end
    [~,Js] = sort(dum);
    I = I(Js);
    J = J(Js);
    %storing results in cell
    for ii = 1:length(I)
        P_pos_fin(ii,1) = ROI_lab(I(ii)); %label ROI1
        P_pos_fin(ii,2) = ROI_lab(J(ii)); %label ROI2
        P_pos_fin(ii,3) = {P_pos_d(I(ii),J(ii))/permut}; %p-values
        P_pos_fin(ii,4) = {P_diff(I(ii),J(ii))}; %number of couplings between ROI1 and ROI2 over time
    end
    %neg
    [I,J] = find(P_neg_d/permut < threshMC/2); %dividing by the total number of permutations and therefore obtaining final P-values (in this case getting already the indexes for writing a more readable output); here threshMC is divided by 2 since I compute this operation 2 times (for positive and negative differences..)
    P_neg_fin = cell(length(I),4);
    %sorting significant results
    dum = zeros(length(I),1);
    for ii = 1:length(I)
        dum(ii) = P_neg_d(I(ii),J(ii))/permut;
    end
    [~,Js] = sort(dum);
    I = I(Js);
    J = J(Js);
    %storing results in cell
    for ii = 1:length(I)
        P_neg_fin(ii,1) = ROI_lab(I(ii)); %label ROI1
        P_neg_fin(ii,2) = ROI_lab(J(ii)); %label ROI2
        P_neg_fin(ii,3) = {P_neg_d(I(ii),J(ii))/permut}; %p-values
        P_neg_fin(ii,4) = {P_diff(I(ii),J(ii))}; %number of couplings between ROI1 and ROI2 over time
    end
end

%plotting
if plot_label == 1 || plot_label == 2
    %ricreate the simmetric matrices for better plots
    P1sum = (P1sum + P1sum') - eye(size(P1sum,1)).*diag(P1sum);
    P2sum = (P2sum + P2sum') - eye(size(P2sum,1)).*diag(P2sum);
    max1 = zeros(2,1);
    max1(1) = max(max(P1sum));
    max1(2) = max(max(P2sum));
    MAX = max(max1);
    clim = [0 MAX];
    PPsum = P1sum - P2sum;
    maxc = zeros(2,1);
    maxc(1) = max(max(PPsum));
    maxc(2) = abs(min(min(PPsum)));
    MAXC = max(maxc);
    climc = [MAXC*(-1) MAXC];
    %plotting
    %old
    figure
    colormap(jet)
    imagesc(P1sum,clim)
    colorbar
    if col == 90 %this is meaningful (and correct) only for a 90-ROI parcellation (such as AAL) with ROIs sorted as LLLRRR (and not LRLRLR)
        labelY = {'19','39','59','79','82','62','42','22','2'};
    end
    labelX = labelY;
    ax = gca;
    ax.YTickLabel = labelY;
    ax.XTickLabel = labelX;
    %new
    figure
    colormap(jet)
    imagesc(P2sum,clim)
    colorbar
    if col == 90
        labelY = {'19','39','59','79','82','62','42','22','2'};
    end
    labelX = labelY;
    ax = gca;
    ax.YTickLabel = labelY;
    ax.XTickLabel = labelX;
    %difference
    figure
    colormap(jet)
    imagesc(PPsum,climc)
    colorbar
    if col == 90
        labelY = {'19','39','59','79','82','62','42','22','2'};
    end
    labelX = labelY;
    ax = gca;
    ax.YTickLabel = labelY;
    ax.XTickLabel = labelX;
end
%checking if nothing is significant.. giving empty [] output
if ~exist('P_pos_fin','var')
    P_pos_fin = [];
end
if ~exist('P_neg_fin','var')
    P_neg_fin = [];
end
%further output
PM{1} = P1sum;
PM{2} = P2sum;
PM{3} = PPsum;

end







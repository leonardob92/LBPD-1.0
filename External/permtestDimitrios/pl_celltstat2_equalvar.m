function [t] = pl_cell_tstat_2sample(Cell1,Cell2)
% Computes a 2-sample t-statistic assuming equal variance. 
% 
% The variables are expressed in cells for speed

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is currently under development, please do not share

N1 = length(Cell1);
N2 = length(Cell2);

%mean of cell1
m1 = Cell1{1};
for i = 2:N1
    m1 = m1 + Cell1{i};
end
m1 = m1/N1; %element mean

%sd of cell 1
var1 = (Cell1{1}-m1).^2;
for i = 2:N1
    var1 = var1 + (Cell1{i}-m1).^2;
end
var1 = var1/(N1-1); %element variance

%mean of cell2
m2 = Cell2{1};
for i = 2:N2
    m2 = m2 + Cell2{i};
end
m2 = m2/N2; %element mean

%sd of cell 2
var2 = (Cell2{1}-m2).^2;
for i = 2:N2
    var2 = var2 + (Cell2{i}-m2).^2;
end
var2 = var2/(N2-1); %element variance

%t-test, equal variance
t = 1 / sqrt(1/N1 + 1/N2) *  (m1 - m2) ./ sqrt( ((N1-1)*var1 + (N2-1)*var2 ) / (N1 + N2 - 2) );






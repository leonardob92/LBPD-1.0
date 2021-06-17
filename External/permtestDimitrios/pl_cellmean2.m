function m = fl_cell_mean2(Cell1,Cell2)
% Computes the mean difference between two cell elements
% 
% The cells must be 1-dimensional
% (The function does not make an internal copy of the cell variable, which
% is critical for large data sets)

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is currently under development, please do not share

N1 = length(Cell1);
N2 = length(Cell2);

%compute mean
m1 = Cell1{1};
for i = 2:N1
    m1 = m1 + Cell1{i};
end
m1 = m1/N1;

%compute mean
m2 = Cell2{1};
for i = 2:N2
    m2 = m2 + Cell2{i};
end
m2 = m2/N2;

m = m1-m2;


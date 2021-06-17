function Cell_slice = pl_cell_slice(Cell,ndx)
%
% Slices cell elements keeping the indices 'ndx'

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is currently under development, please do not share

N = length(Cell);

%compute mean
for i = 1:N
    Cell_slice{i} = Cell{i}(ndx);
end



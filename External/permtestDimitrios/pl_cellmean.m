function m = fl_cell_mean(Cell,w)
% Computes the mean of a cell
% 
% The cell must be 1-dimensional
% (The function does not make an internal copy of the cell variable, which
% is critical for large data sets)

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is currently under development, please do not share

wflag = exist('w');

N = length(Cell);

if wflag %if weighted (useful for permutations)
    
    %compute mean
    if w(1)==1
        m = Cell{1};
    else
        m = -Cell{1};
    end
    for i = 2:N
        if w(i) == 1
            m = m + Cell{i};
        else
            m = m - Cell{i};
        end
    end
    m = m/N;
    
else
    
    %compute mean
    m = Cell{1};
    for i = 2:N
        m = m + Cell{i};
    end
    m = m/N;
    
end



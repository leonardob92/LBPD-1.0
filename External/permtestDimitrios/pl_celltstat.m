function [t] = fl_cell_tstat(Cell,w)
% Computes a 1-sample t-statistic using cells
% 
% The variables are expressed in cells for speed

% This function is part of the permutationlab software:
% Author: Dimitrios Pantazis
% The code is currently under development, please do not share

wflag = exist('w'); %w is a vector of 1s and 0s.

N = length(Cell);

if wflag %if weighted
    
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
    
    %compute std
    if w(1)==1
        sd = (Cell{1}-m).^2;
    else
        sd = (-Cell{1}-m).^2;
    end
    for i = 2:N
        if w(i)==1
            sd = sd + (Cell{i}-m).^2;
        else
            sd = sd + (-Cell{i}-m).^2;
        end
    end
    sd = sqrt(sd)/sqrt(N-1); %element standard deviation
       
    t = m ./ (sd/sqrt(N)); %t-statistic (mean / standard deviation of mean)
  
else
    
    m = Cell{1};
    for i = 2:N
        m = m + Cell{i};
    end
    m = m/N; %element mean
    
    sd = (Cell{1}-m).^2;
    for i = 2:N
        sd = sd + (Cell{i}-m).^2;
    end
    sd = sqrt(sd)/sqrt(N-1); %element standard deviation
    
    t = m ./ (sd/sqrt(N)); %t-statistic (mean / standard deviation of mean)
    
end 
    

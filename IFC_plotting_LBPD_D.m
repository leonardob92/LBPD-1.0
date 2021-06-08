function [  ] = IFC_plotting_LBPD_D( S )

% It plots IFC matrices, either difference between conditions or matrices after
% statistics, or matrices of a single condition, allowing to do some
% subaveraging.
% PLease note that if you have very few time-points t (like 2-10) you should
% try to ask the function to plot a number of matrices x that is a divisor
% of t, otherwise the function makes an approximation that could discard
% relevant information. If you have a large number of time-points t, that
% does not matter very much anymore (the larger the ts the less this is a
% problem..).


%   INPUT:  -either:                    S.contr_1dum: data condition 1 (n x n x t)
%                                       S.contr_2dum: data condition 2 (n x n x t)
%                or:                    S.STAT: tstats or signrank stats
%                                               between cond1 and cond2 (matrix n x n x t)
%           -S.ordresh:                 1 for reshaping ROIs order from
%                                        LRLRLR to LLLRRR.
%                                       0 for keeping order LRLRLR.
%           -S.cond_contr_meanbaseline: 1 for calculating average across
%                                       baseline (if you have S.STAT the
%                                       function assumes that you
%                                       calculated the statistics and
%                                       already had averaged baseline over
%                                       time)
%           -S.Num_window:              number of time-windows to be used for averaging
%                                       successive chunks of data
%           -S.imag_contr_lim:          limits for plotting (e.g. [0 4])
%                                       Leave empty [] for automatic
%                                       adjustment of limits.
%           -S.indimagesnumber:         set which image(s) you want to have plotted
%                                       in independent figures.
%                                       Leave empty [] for not having it(them).

%   OUTPUT: -plots




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 10/12/2018
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 26/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%checking some inputs
limm = S.imag_contr_lim;    
if isfield(S,'contr_1dum') && isfield(S,'contr_2dum')
    contr_1dum = S.contr_1dum;
    contr_2dum = S.contr_2dum;
    contr_stat = 0;
    SE = size(contr_1dum);
elseif isfield(S,'STAT')
    STAT = S.STAT;
    contr_stat = 1;
    SE = size(STAT);
    for kk = 1:SE(3) %this is useful for better visualization of later plots (in some cases..)
        STAT(:,:,kk) = STAT(:,:,kk).*~(eye(SE(1))); %assigning 0s to diagonal..
    end
else
    error('wrong inputs..')
end
if S.ordresh == 1 %if you want to reshape the order of the matrices
    Order = [1:2:SE(1) SE(1):-2:2];
else
    Order = [1:1:SE(1)];
end

%plotting phase synchrony matrices for contrasts or differences between conditions..
figure 
colormap(jet)
Num_window = S.Num_window;
bb = ceil(sqrt(Num_window)); %this better optimise the space for plotting the matrices
wind_lenght = floor((SE(3))/Num_window);
if S.cond_contr_meanbaseline == 0
    if contr_stat == 0 %if comparing same time-points simply subtracting one matrix from the other one
        IFC = contr_1dum(:,:,:) - contr_2dum(:,:,:);
    else %the code here is redundant (because of its historical issues..) and could have been simplified in its tructure.. however it works and for now it is fine like this..
        IFC = STAT(:,:,:); 
    end
    for window = 1:Num_window
        subplot(bb,bb,window)
        if isempty(limm)
            imagesc(mean(IFC(Order,Order,(window - 1) * wind_lenght + 1:window * wind_lenght),3));
        else
            imagesc(mean(IFC(Order,Order,(window - 1) * wind_lenght + 1:window * wind_lenght),3),limm);
        end
        colorbar
        xlabel([num2str(window)])
        set(gcf,'color','w')
    end
else %otherwise
    if contr_stat == 0 %here the same but with averaged baseline over time-points..
        IFC_1 = contr_1dum(:,:,:); %extracting data from cond1
        IFC_2 = mean(contr_2dum,3); %against baseline averaged across the whole time-window
        for window = 1:Num_window
            subplot(bb,bb,window + (kk - 1) * Num_window)
            IFC_dumm = mean(IFC_1(Order,Order,(window - 1) * wind_lenght + 1:window * wind_lenght),3) - IFC_2;
            if isempty(limm)
                imagesc(IFC_dumm);
            else
            	imagesc(IFC_dumm,limm);
            end
            xlabel([num2str(window)])
            set(gcf,'color','w')
        end
    else
        IFC = STAT(:,:,:);
        for window = 1:Num_window
            subplot(bb,bb,window)
            IFC_dumm = mean(IFC(Order,Order,(window - 1) * wind_lenght + 1:window * wind_lenght),3);
            if isempty(limm)
                imagesc(IFC_dumm);
            else
            	imagesc(IFC_dumm,limm);
            end
            colorbar
            xlabel([num2str(window)])
            set(gcf,'color','w')
        end
    end
end

%if user requires single images
if ~isempty(S.indimagesnumber)
    cnt = 0;
    for ii = S.indimagesnumber
        cnt = cnt + 1;
        figure
        colormap(jet)
        DUM = mean(IFC(Order,Order,(ii - 1) * wind_lenght + 1:ii * wind_lenght),3);
        if isempty(limm)
            imagesc(DUM);
        else
            imagesc(DUM,limm);
        end
        labelY = {'19','39','59','79','82','62','42','22','2'};
        labelX = labelY;
        ax = gca;
        ax.YTickLabel = labelY;
        ax.XTickLabel = labelX;
        colorbar
        title([num2str(S.indimagesnumber(cnt))])
        set(gcf,'color','w')
    end
end



end


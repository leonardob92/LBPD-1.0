function [  ] = plotting_ROIs_horzbars_LBPD_D( R )


% It plots horizontal bars describing the frequence of a specific
% parameter (originally thought for verteces degree) for each ROI.
% This relies on previous codes by Joana Cabral.


%   INPUT:  -R: cell array with for each cell a 3 x number of ROIs cell
%               matrix (C)
%               C: first column: ROIs (characters..)
%                  second: numbers to be plotted (double)
%                  third: (only for contrast..) number if p-value was
%                  significant, empty [] if it is not.
%                  This is a manual setting that I do not like very much,
%                  but it can be useful.

%   OUTPUT: -plots..





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Oxford, UK, 16/01/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 26/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%getting maximum value for setting later proper limits for the plots of the different conditions
dum = zeros(length(R),1);
for ii = 1:length(R)
    dum(ii) = max(cell2mat(R{ii}(:,2)));
end
MAX = max(dum);

figure
for ii = 1:length(R)
    subplot(1,length(R),ii)        
    V2 = cell2mat(R{ii}(:,2)); %extracting numbers to be plotted
    V2 = sort(V2,'ascend'); %sorting them for having the most positive up in the plots
    V1 = R{ii}(:,1); %getting corresponding labels
    i_red=find(V2>-1); %finding indeces
    i_blue=find(V2<1);
    barh(i_red,V2(i_red),'FaceColor','r','EdgeColor','none')
    hold on
    barh(i_blue,V2(i_blue),'FaceColor','b','EdgeColor','none')
    if ~isempty(find(V2<0)) %if it is the contrast
        %plotting a star for significance
        i_pos = find(double(abs(double(cellfun('isempty',(R{ii}(end:-1:1,3))))-1) > 0) + double((V2(90:-1:1) > 0)) == 2); %finding indeces of significant ROIs (cond1 vs cond2)
        i_neg = find(double(abs(double(cellfun('isempty',(R{ii}(end:-1:1,3))))-1) > 0) + double((V2(90:-1:1) < 0)) == 2); %cond2 vs cond1

%             original.. it is not working because of the random width.. boh
%             hold on
%             barh(i_yel,V2(i_yel),'FaceColor','k','EdgeColor','none') %and change their color from red or blue to a new colour

        title(['Contrast'])
        xlim([(min(V2)-1) max(V2)+1])
        plot(V2(i_pos) - 0.5,i_pos,'*k') %here it is confusing but since the indeces are reversed the i_pos in the data corresponds to 'pos' but in the barh corresponds to 'neg'..
        plot(V2(i_neg) + 0.5,i_neg,'*k')
    else
        title(['Condition ' num2str(ii-1)])
        xlim([0 MAX+1])
    end
    set(gca,'YTick',1:length(V1),'Fontsize',5)  %original fontsize 6  %original y-tick 1:91
    set(gca,'YTickLabel',V1(end:-1:1))
    ylim([0 (length(V1)+1)])
    grid on
end

end









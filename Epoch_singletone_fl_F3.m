function [OUT ] = Epoch_singletone_fl_F3 ( C,x,y,z,d,task,por)


% Epoching, concatenating trials, reshaping, averaing and removing source leakage correction
% of the memory encoding task data.
% This function was made for didactic reasons and it is not to be
% considered a general function that is very likely to be applied to many more
% general cases..


%   INPUT: C    = actual data (ROIs x time-samples x musical piece repetition
%          x    = duration of every single note
%          y    = number of total notes
%          z    = number of trials to be concatenated
%          d    = onsets of the notes
%          task = 1 for task; 0 for resting state (used as baseline)
%          por  = repetitions of the musical piece you want to consider
%
%
%   OUTPUT:A = data epoched, averaged, source leakage corrected and with
%              envelope computed (with trials concatenation)



% Francesco Carlomagno, Leonardo Bonetti
% Last modified: 12/10/2019, Aarhus






OUT = cell(length(C),2);
for ii = 1:length(C)
    %epoching and creating trials on the basis of onsets stored in d matrix
    D2 = C{ii,1};
%     disp('epoching..')
    if task == 1 %for task
        PD = size(D2); 
%         D3 = reshape(D2,PD(1),PD(2)*PD(3)); 
        if PD(3) == 4 %hard coding since 1 subject has no fourth repetition of the musical piece
            por2 = por;
        else
            por2 = por;
            por2(por2==4) = [];
        end
        d1o = zeros(90,x,z,length(por2));
        for pp = por2%1:PD(3)
            uu = 0;
            a = zeros(90,x,y);
            for jj = 1:y 
                uu = uu + 1; 
                a(:,:,uu) = D2(:,d(jj,2):d(jj,3),pp);
            end
            if find(a == 0)
                warning(['subject ' num2str(ii) ' has 0s']);
            end
            A = reshape(a(:,:,1:(floor(y/z)*z)),90,x,floor(y/z),z);
            A(A == 0) = NaN;% set the 0s to NaN (just in case..)
            d1o(:,:,:,pp) = squeeze(mean(A,3,'omitnan'));
        end
        d1o2 = mean(d1o,4);
    else %for resting state (used as baseline)
        por2 = por;
        d1o = zeros(90,x,z,length(por2));
        for pp = por2 %loop until 4 epochs
            uu = 0;
            a = zeros(90,x,y);
            for jj = 1:y %loop until y
                uu = uu + 1; %update the counter
                a(:,:,uu) = D2(:,(x*y)*(pp-1)+d(jj,2):(x*y)*(pp-1)+d(jj,3)); %getting a mixture of note onsets and randomised starting points for creating fake trials (to be used as baseline)
            end
            if find(a == 0)
                warning(['subject ' num2str(ii) ' has 0s']);
            end
            %doing same procedure as done before for task data
            A = reshape(a(:,:,1:(floor(y/z)*z)),90,x,floor(y/z),z); % reshaping process with floor for made the average possible
            A(A == 0) = NaN;% set the 0s to NaN
            d1o(:,:,:,pp) = squeeze(mean(A,3,'omitnan'));
        end
        d1o2 = mean(d1o,4);
    end
%     if find(a == 0) 
%         warning(['subject ' num2str(ii) ' has 0s']);
%     end
%     %reshaping matrices for doing subaveraging on some trials and then concatenating the resulting subaveraged trials
%     if task == 1 %hard coding for specific situation in task database..
% %         if PD(3) == 3 %hard coding because 1 subj had only 3 repetitions of c minor Bach
% %             a2 = a(:,:,1:y*3); %creating another variable for making reshaping possible
% %             A = reshape(a2(:,:,1:(floor(y*3/z)*z)),90,x,floor(y*3/z),z); % reshaping process with flor for made the average possible
% %         else
% %             A = reshape(a(:,:,1:(floor(y*4/z)*z)),90,x,floor(y*4/z),z); % reshaping process with floor for made the average possible
%             A = reshape(a(:,:,1:(floor(y*length(porcodio)/z)*z)),90,x,floor(y*length(porcodio)/z),z); % reshaping process with floor for made the average possible
% 
% %         end
%     else
%         A = reshape(a(:,:,1:(floor(y*4/z)*z)),90,x,floor(y*4/z),z); % reshaping process with floor for made the average possible
%     end
    %subaveraging
%     A(A == 0) = NaN;% set the 0s to NaN
%     dio = squeeze(mean(A,3,'omitnan'));
    %actual concatenation of subaveraged trials
    diodo = reshape(d1o2,90,x*z);
    %removing source leakage and calculating envelope
    disp('removing source leakage correction..')
    tso_lk = ROInets.remove_source_leakage(diodo,'symmetric'); %leakage correction
    Ab = hilbenv(tso_lk); %envelope
    %storing output
    OUT(ii,1) = {Ab}; % create OUT variable with epoching data in the first column
    OUT(ii,2) = C(ii,2); % create OUT variable with ID data in the second column
    disp(['You have just done the subj num ' num2str(ii)]) % display confirm process
end

end


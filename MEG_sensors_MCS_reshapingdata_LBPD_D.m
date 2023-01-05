function [ MAG_data_pos, MAG_data_neg, GRAD_data, MAG_GRAD_tval ] = MEG_sensors_MCS_reshapingdata_LBPD_D( S )

% It reshapes the order of data according to labels that are
% submitted, providing 2D approximation of MEG sensors layout.
% This function works in combination with a following function that finds
% spatio-temporal patterns in the data.
% The 2D approximation provided by this function is not ideal, but
% simulations and testing showed it to work well enough.
% Of course, you need to have the same order between labels and
% data. Furthermore, right now it assumes that magnetometers and
% gradiometers are ordered in the same way (e.g. order of MAG 1211 == order
% of GRAD 1212 + 1213). Therefore, input data needs to be submitted
% following this logic. The actual reshaping is done on the basis of
% magnetometers. Therefore you need to submit magnetometer labels, even if
% you work with gradiometers.
% OBS!! Moreover, you need to provide separately the p-values related to
% POSITIVE or NEGATIVE t-values for magnetometers (because of the double
% polarity of the magnetometers).
% For clarity, you MUST provide binarized data (1 corresponds to
% significant elements, 0 to non-significant ones).
% E.g. TSTAT_mag_pos will have 1s on the channels where the t-values
% were positive and significant and 0s in the other channels.



%  INPUT:   -S.label:         magnetometers label sorted in the same order as the
%                             actual data
%           -S.TSTAT_mag_pos: magnetometers data sorted in the same order
%                             as the label (p-values of the POSITIVE t-values).
%           -S.TSTAT_mag_neg: magnetometers data sorted in the same order
%                             as the label (p-values of the NEGATIVE t-values)
%           -S.TSTAT_grad:    same but for gradiometers data
%           -S.raw_channels:  2D approximation of magnetometers MEG layout
%           -S.TVal:          matrix with actual t-values (102 channels x time-points x sensor-type (first mag second combined planar grad)).
%                             If you do not want to provide it, simply do not write the field.

%  OUTPUT:  -MAG_data_pos:    reshaped POSITIVE magnetometers data
%           -MAG_data_neg:    reshaped NEGATIVE magnetometers data
%           -GRAD_data:       reshaped gradiometers data
%           -MAG_GRAD_tval:   (if requested..), reshaped t-values for magnetometers (x,x,x,1) and gradiometers (x,x,x,2)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 10/04/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 28/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%inputs
label = S.label;
TSTAT_mag_pos = S.TSTAT_mag_pos;
TSTAT_mag_neg = S.TSTAT_mag_neg;
TSTAT_grad = S.TSTAT_grad;
raw_channels = S.raw_channels;

%preallocating variables
[~,coltstat] = size(TSTAT_mag_pos);
[rawmeg,colmeg] = size(raw_channels);
MAG_data_pos = zeros(rawmeg,colmeg,coltstat);
MAG_data_neg = zeros(rawmeg,colmeg,coltstat);
GRAD_data = zeros(rawmeg,colmeg,coltstat);
if isfield(S,'TVal')
    MAG_GRAD_tval = zeros(rawmeg,colmeg,coltstat,2);
    Tval = S.TVal;
else
    MAG_GRAD_tval = [];
end

%actual reshaping
for pp = 1:coltstat %iterates over time
    for ii = 1:rawmeg %rows of the MEG layout
        for jj = 1:colmeg %columns of the MEG layout
            if raw_channels{ii,jj} ~= 0 %if the MEG layout has an actual channel
                [Inx,~] = find(label(:,1) == cell2mat(raw_channels(ii,jj))); %find the index of the label corresponding to the index of the specific channel in the dataset
                MAG_data_pos(ii,jj,pp) = TSTAT_mag_pos(Inx,pp); %move the data (mag)
                MAG_data_neg(ii,jj,pp) = TSTAT_mag_neg(Inx,pp); %move the data (mag)
                GRAD_data(ii,jj,pp) = TSTAT_grad(Inx,pp); %move the data (grad)
                if isfield(S,'TVal')
                    MAG_GRAD_tval(ii,jj,pp,1) = Tval(Inx,pp,1); %move the data (mag tvals)
                    MAG_GRAD_tval(ii,jj,pp,2) = Tval(Inx,pp,2); %move the data (grad tvals)
                end
            end
        end
    end
end



end


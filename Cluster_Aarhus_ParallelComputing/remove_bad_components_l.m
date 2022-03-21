function D = remove_bad_components_l(D,S)
    % Take in a D object with D.ica.bad_components
    % Make a new montage where all the channels stay the same but 
    % they have had the bad ICA components subtracted from them

    D = D.montage('switch',0);

    if strcmp(S.modality,'EEG')
        chantype = 'EEG';
        modality = 'EEG';
    else
        chantype = {'MEG','MEGANY'};
        modality = 'MEG';
    end

    chan_inds = D.ica.chan_inds; % This is a record of which channels have had components subtracted from them
    bad_components = unique(D.ica.bad_components);
    megdata        = D(chan_inds,:,:);
    megdata        = reshape(megdata,size(megdata,1),[]);

    tc = (D(chan_inds,:,:)'*pinv(D.ica.sm(chan_inds,:))').';
    tra = eye(D.nchannels);
    dat_inv = pinv_plus(megdata', D.ica.params.num_ics);
    tra(chan_inds,chan_inds) = (eye(numel(chan_inds)) - dat_inv*(tc(bad_components,:)'*D.ica.sm(chan_inds,bad_components)'))';
    
    tmp = struct(D);
    D = add_montage(D,tra,S.montagename,D.chanlabels,tmp.channels);
end
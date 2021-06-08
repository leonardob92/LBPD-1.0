function [ OUT ] = MEG_sensors_plotting_ttest_LBPD_D_marina( S )

% It prepares data (loading data or extracting it from SPM objects (usually
% obtained from previous preprocessing in OSL)), optionally calculates t-tests
% between two conditions and then provides multiple plotting options:
% 1)single-channel waveforms
% 2)averaged channels waveforms (independently for left and right
%   hemispheres)
% 3)topoplots
% More details about input requirements are provided below.
% Note that the labels of conditions are graphically not very good-looking
% if they present underscores.. (e.g. Old_Correct). Consider to either
% provide labels without underscores or to manually change the labels in the
% plot.. or in the codes of the function.



% INPUT:

% -preparing data settings
%   -S.data:                             data to be used for plotting.
%                                        Order must be: mag 0111,grad 0112+0113, mag 0121, grad 0122+0123, etc.
%                                        Leave empty [] if you want to extract data
%                                        from the SPM objects (in that case
%                                        an independent path for subject
%                                        must be provided).
%                                        Gradiometers must be passed
%                                        combined (e.g. 0112+0113).
%   -S.outdir:                           character path for output directory. If not existent, it automatically creates it.
%   -S.chanlabels:                       channel labels to be provided coherently with the
%                                        data order.
%                                        If data = [] chanlabels are
%                                        automatically extracted from SPM
%                                        objects.
%   -S.time_real:                        time in seconds describing data time-series.
%                                        If data = [] time_real is automatically extracted
%                                        from SPM objects.
%   -S.timeextract:                      vector with time-points (samples)
%                                        to be extracted from SPM objects.
%                                        If empty [], default is all time-points.
%                                        This is meaningfull only if S.data = [].
%   -S.centerdata0:                      1 to make the data starting from 0
%   -S.spm_list:                         path to SPM object for each subject.
%                                        This is meaningfull only if S.data = [].
%   -S.conditions:                       conditions to be extracted (cell
%                                        array with characters).
%                                        E.g. {'Old';'New'}.
%                                        Here I assume that a condition A
%                                        is represented in the data only
%                                        once.
%   -S.save_data:                        1 for saving the data     
%   -S.save_name_data:                   name for data file to be saved (character)

% -1)individual waveform plotting
%   -S.waveform_singlechannels_label:    1 to plot single channel waveforms
%   -S.wave_plot_conditions_together:    1 for plotting the average of all
%                                        conditions
%   -S.mag_lab:                          1 for magnetometers
%                                        2 for gradiometers
%   -S.x_lim_temp_wave:                  limits for time (in secs)
%                                        E.g. [-0.1 3.4]
%   -S.y_lim_ampl_wave:                  limit for amplitude
%                                        E.g. [0 120] magnetometes, [0 6] gradiometers

% -2)averaged waveform plotting
%   -S.waveform_average_label:           1 to plot averaged waveforms
%                                        across specific channels
%   -S.left_mag:                         channels to be plotted together
%                                        (indeces from original channel labels).
%                                        The name 'left_mag' is confusing 
%                                        and due to historical reasons of
%                                        this function.
%   -S.avewave_contrast:                 1 to plot the contrast between
%                                        S.cond_ttests_tobeplotted_topoplot(1) and
%                                        S.cond_ttests_tobeplotted_topoplot(2).
%   -S.signtp:                           significant time-points to be
%                                        plotted (in seconds).
%                                        It must be a 1 x significant time-points cell array.
%                                        E.g. S.signtp(1,1) = {1:3}; S.signtp(1,2) = {6:10}; etc.
%                                        Leave empty [] if you want no
%                                        signficant time-points to be highlighted.
%   -S.legc:                             1 to plot legend.
%                                        0 not to plot legend.
%   -S.save_label_waveaverage:           1 to save averaged waveforms plots
%   -S.label_plot:                       name for saved plots (character)

% -ttests settings
%   -S.t_test_for_permutations:          1 to compute univariate t-tests for each channel and time-point.
%   -S.cond_ttests_tobeplotted_topoplot: conditions to be contrasted and then plotted in the
%                                        topoplot and averaged waveform difference between conditions.
%                                        E.g. [2 1] contrasts cond2 against cond1.
%                                        Maximum 2 conditions can be contrasted.
%                                        In case of one condition only (e.g. [1]) that condition is
%                                        contrasted with the mean value of its own baseline.

% -3)topoplotting settings:
%   -S.topoplot_label:                   1 to plot topoplot
%   -S.fieldtrip_mask:                   path to fieldtrip mask file used
%                                        to plot topoplot (tipicaly
%                                        provided in the external folder of
%                                        this collection of functions).
%   -S.topocontr:                        1 for plotting contrasts
%                                        0 for plotting mean across
%                                        conditions or single condition.
%   -S.topocondsing:                     empty [] for mean across conditions.
%                                        E.g. [1] for plotting only
%                                        S.conditions{1}.
%                                        Meaningful only if S.topocontr = 0. 
%   -S.xlim:                             set limit time for topoplot (in seconds)
%                                        E.g. [0.95 1.05]
%   -S.zlimmag:                          set limit amplitude for topoplot (magnetometers)
%                                        E.g. [0 6]
%                                        Leave empty [] for automatic determination.
%   -S.zlimgrad:                         same, but for gradiometers
%   -S.colormap_spec:                    set specification for Matlab colormap.
%                                        E.g. 'jet'.
%   -S.topoplot_save_label:              1 for saving the outputted topoplots



% OUTPUT:
% -preparing data and plotting..






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, Aarhus, DK, 25/06/2019
% Leonardo Bonetti, MIT, Boston, USA, 12/08/2019
% Leonardo Bonetti, Aarhus, DK, 08/10/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











%preparing data - computation
if isfield(S,'data') %it requires that S.data exists
    if isempty(S.data) %if you want to load data from disk leave it empty
        spm_list = S.spm_list;
        D = spm_eeg_load(spm_list{1}); %load data from disk for knowing how many time-points are in the data
        SS = size(D);
        if ~isempty(S.timeextract) %if user requires specific time-points to be extracted
            tex = S.timeextract; %extracting them
        else
            tex = 1:SS(2); %otherwise extracting all time-points
        end
        time_sel = D.time(tex); %assuming SPM objects have the same time extent..
        data_mat = zeros(204,length(tex),length(spm_list),length(S.conditions)); %preallocation of space
        for jj = 1:length(spm_list) %over subjects
            D = spm_eeg_load(spm_list{jj}); %loading data for each subject
            conds = zeros(1,length(S.conditions));
            for ii = 1:length(S.conditions)
                conds(ii) = sum(find(strcmp((D.conditions),S.conditions{ii}))); %getting the condition numbers ordered according to the specification of the user in S.conditions (sum to avoid the problem of zero-vector occurring if it does not find the match between the specified label and the data label.. of course if you specify two times the same label in the structure S it crashes as well.. but you would need to be quite peculiar to do that so I assume you will not do that..)
            end
            if sum(double(conds == 0)) > 0 %this used to be an error; now it is converted in warning since user may be interested in plotting (mainly for exploratory reasons) a condition that occurs very seldomly and not for all participants (for instance ERF associated to wrong answers in a simple behavioral task)
                warning(['for subject ' num2str(jj) ' you have specified a condition name that does not match with the condition labels in the data..'])
            end
            chans = D.chanlabels;
            idxMEGc1 = find(strcmp(chans,'MEG0111')); %getting extreme channels indexes
            idxMEGc2 = find(strcmp(chans,'MEG2642+2643'));          
            if isempty(idxMEGc1) || isempty(idxMEGc2)
                error('the label of the MEG channels do not correspond to the standard ones used in this function.. you need to change either the labels of the function or the ones of your data.. check for spaces between MEG and 0111.. sometimes that is the problem..')
            end
            chanlabels = D.chanlabels(idxMEGc1:idxMEGc2);
            for hh = 1:length(conds)
                if conds(hh) ~= 0
                    data_mat(:,:,jj,hh) = D(idxMEGc1:idxMEGc2,tex,conds(hh));
                end
            end
            disp(['subject ' num2str(jj) ' computed'])
        end
        %centering data on 0
        if S.centerdata0 == 1
            data2 = zeros(size(data_mat,1),size(data_mat,2),size(data_mat,3),size(data_mat,4)); %preallocating space
            for ccc = 1:size(data_mat,4) %over conditions
                for ppp = 1:size(data_mat,3) %over subjects
                    for sss = 1:size(data_mat,1) %over MEG channels
                        data2(sss,:,ppp,ccc) = data_mat(sss,:,ppp,ccc) - data_mat(sss,25,ppp,ccc); %subtracting the value of the first time-point to make the time-serie starting from 0. Here the data was actually already been baseline-corrected, so this passage may not be that crucial.. however, I think it is a good idea to do that
                    end
                end
            end
            data_mat = data2;
        end
        %if requested, saving data
        if S.save_data == 1
            %creating the directory if it does not exist already
            if ~exist(S.outdir,'dir')
                mkdir(S.outdir);
            end
            save([S.outdir '/' S.save_name_data '.mat'],'data_mat','chanlabels','time_sel');
            condsl = S.conditions;
            save([S.outdir '/conditions.mat'],'condsl');
        end
    else %otherwise taking the data and the condition labels passed in by the user
        data_mat = S.data;
        conds = S.conditions;
        chanlabels = S.chanlabels;
        time_sel = S.time_real;
    end
else
    error('you did not provide any data.. if you want to read data from SPM objects saved on disk, please give as input: S.data = []; ')
end
%making 0s of data_mat equal to NaN for preventing later averages and t-tests to be affected by 0s that are not actual data                        
data_mat(data_mat(:,2:end,:,:) == 0) = NaN;
Sz = size(data_mat); %used later for calculating standard errors

clcc = S.colorcode; %
%color_line = {'b','r','g','k','y','m','c','b','r','k'};%it provides 8 different colours.. this is for working around the problem of indexing different colours for different conditions.. here it provides 8 different options..
%if you have more then 8 the plot would anyway be not really understandle, so I did not really care here.. of course this could be easily improved if someone feels like it is actually useful..    
if S.waveform_singlechannels_label == 1
    %getting a few inputs
    kkk = S.mag_lab - 1;
    x_lim_temp = S.x_lim_temp_wave;
    y_lim_ampl = S.y_lim_ampl_wave;
    %single-channel waveform plotting
    if isempty(x_lim_temp) %if time-window for plotting is not specified, it takes the extreme of D.time (data)
        x_lim_temp = [time_sel(1) time_sel(end)];
    end
    data_mean = squeeze(nanmean(data_mat,3)); %averaging data over participants
    data_stde = squeeze(nanstd(data_mat,0,3))/sqrt(Sz(3)); %calcolating standard error
    cond_n = length(conds); %get the number of conditions
    if S.wave_plot_conditions_together == 1
        data_mean = nanmean(data_mean,3); %averaging over conditions
        data_stde = nanmean(data_stde,3);
        plotlab = 'avg conditions';
        cond_n = 1; %and therefore the number of conditions becomes one
    end
    countwave = -1;
    for jjk = 1:12
        clf(figure(jjk));
        figure(jjk);
        for iik = 1:9
            kkk = kkk + 1;
            countwave = countwave + 1;
            if countwave + kkk < 205
                legcell = cell(1,cond_n);
                for aa = 1:cond_n
                    subplot(3,3,iik);
                    if S.wave_plot_conditions_together == 1
                        %original h1 = plot(time_sel,data_mean((countwave + kkk),:,aa),color_line{aa},'LineWidth',2,'DisplayName',plotlab); %average waveform
                        h1 = plot(time_sel,data_mean((countwave + kkk),:,aa),'color',clcc{aa},'LineWidth',2,'DisplayName',plotlab); %average waveform %marina
                    else
                        %original h1 = plot(time_sel,data_mean((countwave + kkk),:,aa),color_line{aa},'LineWidth',2,'DisplayName',S.conditions{aa}); %average waveform
                        h1 = plot(time_sel,data_mean((countwave + kkk),:,aa),'color',clcc{aa},'LineWidth',2,'DisplayName',S.conditions{aa}); %average waveform %marina
                    end
                    xlim(x_lim_temp);
                    if ~isempty(y_lim_ampl)
                        ylim(y_lim_ampl);
                    end
                    hold on
                    %original plot(time_sel,data_mean((countwave + kkk),:,aa) + data_stde((countwave + kkk),:,aa),[':' color_line{aa}],'LineWidth',0.5); %upper std error
                    plot(time_sel,data_mean((countwave + kkk),:,aa) + data_stde((countwave + kkk),:,aa),[':' 'color',clcc{aa}],'LineWidth',0.5); %upper std error %marina
                    xlim(x_lim_temp);
                    if ~isempty(y_lim_ampl)
                        ylim(y_lim_ampl);
                    end
                    hold on
                    %original plot(time_sel,data_mean((countwave + kkk),:,aa) - data_stde((countwave + kkk),:,aa),[':' color_line{aa}],'LineWidth',0.5); %lower std error
                    plot(time_sel,data_mean((countwave + kkk),:,aa) - data_stde((countwave + kkk),:,aa),[':' 'color',clcc{aa}],'LineWidth',0.5); %lower std error %marina
                    xlim(x_lim_temp);
                    if ~isempty(y_lim_ampl)
                        ylim(y_lim_ampl);
                    end
                    hold on
                    legendv(1,aa) = h1; %this is for working around the problem of the multiple legends for stde that we do not want..
                    if S.wave_plot_conditions_together == 1
                        legcell(1,aa) = {plotlab};
                    else
                        legcell(1,aa) = {S.conditions{aa}}; %same..
                    end
                    grid on
                    grid minor
                end                   
                legend([legendv],legcell) %same..
                legend('show');
                title(chanlabels(countwave + kkk));
                set(gcf,'Color','w')
            end
        end
    end
    clear legendv
end
        
        
%computing t-tests..
if S.t_test_for_permutations == 1
    disp('computing t-tests..')
    if length(S.cond_ttests_tobeplotted_topoplot) == 2
        disp(['contrast ' S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} '_vs_' S.conditions{S.cond_ttests_tobeplotted_topoplot(2)}])
    else
        disp(['contrast ' S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} '_vs_ its own baseline'])
    end
    disp('this can take few seconds up to several minutes..')
    if length(S.cond_ttests_tobeplotted_topoplot) > 2
        error('you cannot compute t-test for more than 2 conditions..')
    end
    %computing contrasts
    Sz = size(data_mat);
    deftstat = zeros(Sz(1),Sz(2));
    deftstatp = zeros(Sz(1),Sz(2));
    for jsk = 1:Sz(1)
        disp(['calculating t-tests for channel ' num2str(jsk) ' / ' num2str(Sz(1))])
        for jsl = 1:Sz(2)
            if length(S.cond_ttests_tobeplotted_topoplot) == 2 %if you want to compare two conditions
                a = squeeze(data_mat(jsk,jsl,:,S.cond_ttests_tobeplotted_topoplot(1)));
                b = squeeze(data_mat(jsk,jsl,:,S.cond_ttests_tobeplotted_topoplot(2)));
                ad = find(isnan(a)); %this is for preventing t-tests to be calculated considering NaNs.. if there are any subject in a or in b corresponding to NaN both a and b are deprivated of that subject
                ab = find(isnan(b));
                abd = [ab ad];
                a(abd) = [];
                b(abd) = [];
                [~,pttest,~,tstats] = ttest(a,b);
                deftstat(jsk,jsl) = tstats.tstat;
                deftstatp(jsk,jsl) = pttest;
            else %if you want to compare one condition with its own averaged baseline (tha
                a = squeeze(data_mat(jsk,jsl,:,S.cond_ttests_tobeplotted_topoplot(1)));
                b = squeeze(mean(data_mat(jsk,1:find(time_sel==0),:),2)); %looking for time 0 seconds and taking time from time-sample 1 to time-sample corresponding to 0 seconds (so to the onset of the stimulus).. in other words averaging the baseline of each subject and each channel over time
                ad = find(isnan(a)); %this is for preventing t-tests to be calculated considering NaNs.. if there are any subject in a or in b corresponding to NaN both a and b are deprivated of that subject
                ab = find(isnan(b));
                abd = [ab ad];
                a(abd) = [];
                b(abd) = [];                
                [~,pttest,~,tstats] = ttest(a,b);
                deftstat(jsk,jsl) = tstats.tstat;
                deftstatp(jsk,jsl) = pttest;
            end
        end
    end
    warning('now the MAGMEG and GRADMEG are being moved assuming that their original order is: MEG0111 - MEG0112+0113 - MEG0121 - MEG0122+0123 etc. as it is quite often.. CHECK if this applies to your case, otherwise just compute the t-tests yourself!')
    TSTAT_grad = zeros(102,Sz(2));
    TSTAT_mag = zeros(102,Sz(2));
    TSTATP_grad = zeros(102,Sz(2));
    TSTATP_mag = zeros(102,Sz(2));
    disp('reshaping order of the data..')
    %moves the gradiometers
    count = 0;
    for ii = 1:102
        count = count + 2;
        for kk = 1:Sz(2) %1:length(D.time)
            TSTAT_grad(ii,kk) = deftstat(count,kk); %D((count),kk,5);
            TSTATP_grad(ii,kk) = deftstatp(count,kk);
        end
    end
    %moves the magnetometers
    count = -1;
    for ii = 1:102
        count = count + 2;
        for kk = 1:Sz(2) %1:length(D.time)
            TSTAT_mag(ii,kk) = deftstat(count,kk); %D((count),kk,5);
            TSTATP_mag(ii,kk) = deftstatp(count,kk);
        end
    end
    %preparing OUT structure
    OUT.TSTAT_mag = TSTAT_mag;
    OUT.TSTATP_mag = TSTATP_mag;
    OUT.TSTAT_grad = TSTAT_grad;
    OUT.TSTATP_grad = TSTATP_grad;
    OUT.chanlabels = chanlabels;
    OUT.time_sel = time_sel;
    disp('saving statistics..') %this is always done automatically
    if length(S.cond_ttests_tobeplotted_topoplot) == 2 %if you want to compare two conditions        
        save([S.outdir '/' S.save_name_data '_OUT_' S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} '_vs_' S.conditions{S.cond_ttests_tobeplotted_topoplot(2)} '.mat'],'OUT')
    else %if you want to compare one condition with its own baseline
        save([S.outdir '/' S.save_name_data '_OUT_' S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} '_vs_baseline.mat'],'OUT')        
    end
else
    OUT.TSTAT_mag = [];
    OUT.TSTATP_mag = [];
    OUT.TSTAT_grad = [];
    OUT.TSTATP_grad = [];
end
%averaged waveforms plotting
if S.waveform_average_label == 1
    %getting some inputs
    x_lim_temp = S.x_lim_temp_wave;
    y_lim_ampl = S.y_lim_ampl_wave;
    %single-channel waveform plotting
    if isempty(x_lim_temp) %if time-window for plotting is not specified, it takes the extreme of D.time (data)
        x_lim_temp = [time_sel(1) time_sel(end)];
    end
    left_mag = S.left_mag;
    
    %originally I wrote codes for plotting one picture per hemisphere since
    %due to the often-happpening symmetry of brain activity that sounded
    %like a good idea. However, I think it is now more reasonable to provide only one plotting option and run it again as many times as needed by user 
%     right_mag = S.right_mag;
%     if S.mag_lab == 2
%         left_mag = left_mag + 1;
% %         right_mag = right_mag + 1;
%     elseif S.mag_lab ~= 1
%         warning('are you sure you know what you want to plot..? magnetometers or gradiometers or what..?')
%     end
    %preparing left and right hemisphere data    
    Sz = size(data_mat);
    mat_left_mag = zeros(length(left_mag),Sz(2),Sz(3),length(conds));
%     mat_right_mag = zeros(length(right_mag),Sz(2),Sz(3),length(conds));
    for count_lx = 1:length(left_mag) %count for both left and right
        for gg = 1:length(conds)
            mat_left_mag(count_lx,:,:,gg) = data_mat(left_mag(count_lx),:,:,gg);
%             mat_right_mag(count_lx,:,:,gg) = data_mat(right_mag(count_lx),:,:,gg);
        end
    end
    %mean
    mean_left_mag = squeeze(nanmean(squeeze(nanmean(mat_left_mag,1)),2));
%     mean_right_mag = squeeze(mean(squeeze(mean(mat_right_mag,1)),2));
    %standard error
    stde_left_mag = squeeze(std(squeeze(nanmean(mat_left_mag,1)),0,2)/sqrt(Sz(3)));
%     stde_right_mag = squeeze(std(squeeze(mean(mat_right_mag,1)),0,2)/sqrt(Sz(3)));
    %if you want to plot the average over conditions
    if S.wave_plot_conditions_together == 1
        mean_left_mag = nanmean(mean_left_mag,2);
%         mean_right_mag = mean(mean_right_mag,2);
        stde_left_mag = nanmean(stde_left_mag,2);
%         stde_right_mag = mean(stde_right_mag,2);
        plotlab = 'avg conditions';
        cond_n = 1; %and therefore the number of conditions becomes one
    else
        cond_n = length(conds);
    end
    %actual plotting
    legcell = cell(1,cond_n);
    clf(figure(13))
    figure(13)
    thc = ones(cond_n)*2;
%     %%% hard coding.. just for one specific purpose.. not for more general
%     %%% ones..
%     color_line = {'b','r','b','r'};
%     thc = [3,3,2,2];
%     %%%
    
    for aa = 1:cond_n
        if S.wave_plot_conditions_together == 1
            %original h1 = plot(time_sel,mean_left_mag(:,aa),color_line{aa},'LineWidth',thc(aa),'DisplayName',plotlab); %average waveform 
            h1 = plot(time_sel,mean_left_mag(:,aa),'color',clcc{aa},'LineWidth',thc(aa),'DisplayName',plotlab); %average waveform %marina
        else
            %original h1 = plot(time_sel,mean_left_mag(:,aa),color_line{aa},'LineWidth',thc(aa),'DisplayName',S.conditions{aa}); %average waveform
            h1 = plot(time_sel,mean_left_mag(:,aa),'color',clcc{aa},'LineWidth',thc(aa),'DisplayName',S.conditions{aa}); %average waveform %marina
        end
        xlim(x_lim_temp);
        if ~isempty(y_lim_ampl)
            ylim(y_lim_ampl);
        end
        hold on
        %original plot(time_sel,mean_left_mag(:,aa) + stde_left_mag(:,aa),[':' color_line{aa}],'LineWidth',0.5); %upper std error
        plot(time_sel,mean_left_mag(:,aa) + stde_left_mag(:,aa),':', 'color', clcc{aa},'LineWidth',0.5); %upper std error %marina
        xlim(x_lim_temp);
        if ~isempty(y_lim_ampl)
            ylim(y_lim_ampl);
        end
        hold on
        %original plot(time_sel,mean_left_mag(:,aa) - stde_left_mag(:,aa),[':' color_line{aa}],'LineWidth',0.5); %lower std error
        plot(time_sel,mean_left_mag(:,aa) - stde_left_mag(:,aa),':', 'color', clcc{aa},'LineWidth',0.5); %lower std error %marina
        xlim(x_lim_temp);
        if ~isempty(y_lim_ampl)
            ylim(y_lim_ampl);
        end
        hold on
        legendv(1,aa) = h1; %this is for working around the problem of the multiple legends for stde that we do not want..
        if S.wave_plot_conditions_together == 1
            legcell(1,aa) = {plotlab};
        else
            legcell(1,aa) = {S.conditions{aa}}; %same..
        end
    end
%     grid on
%     grid minor
    if ~isempty(S.signtp{1})
        %plotting the significant time-points with gray shades
        patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
        ylims = get(gca,'YLim');
        for ii = 1:length(S.signtp)
%             sgf2 = S.signtp{ii}./S.sr;
            sgf2 = S.signtp{ii};
            patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.4)
            hold on
        end
        %plotting again in order to have a better visualization of the
        %time-series.. currently this implementation looks quite awful, but it
        %works and takes a very reasonable amount of time to compute..
        for aa = 1:cond_n
            if S.wave_plot_conditions_together == 1
                %original h1 = plot(time_sel,mean_left_mag(:,aa),color_line{aa},'LineWidth',thc(aa),'DisplayName',plotlab); %average waveform
                h1 = plot(time_sel,mean_left_mag(:,aa),'color',clcc{aa},'LineWidth',thc(aa),'DisplayName',plotlab); %average waveform %marina
            else
                %original h1 = plot(time_sel,mean_left_mag(:,aa),color_line{aa},'LineWidth',thc(aa),'DisplayName',S.conditions{aa}); %average waveform
                h1 = plot(time_sel,mean_left_mag(:,aa),'color',clcc{aa},'LineWidth',thc(aa),'DisplayName',S.conditions{aa}); %average waveform %marina
            end
            xlim(x_lim_temp);
            if ~isempty(y_lim_ampl)
                ylim(y_lim_ampl);
            end
            hold on
            %original plot(time_sel,mean_left_mag(:,aa) + stde_left_mag(:,aa),[':' color_line{aa}],'LineWidth',0.5); %upper std error
            plot(time_sel,mean_left_mag(:,aa) + stde_left_mag(:,aa),[':' 'color',clcc{aa}],'LineWidth',0.5); %upper std error %marina
            xlim(x_lim_temp);
            if ~isempty(y_lim_ampl)
                ylim(y_lim_ampl);
            end
            hold on
            %original plot(time_sel,mean_left_mag(:,aa) - stde_left_mag(:,aa),[':' color_line{aa}],'LineWidth',0.5); %lower std error
            plot(time_sel,mean_left_mag(:,aa) - stde_left_mag(:,aa),[':' 'color',clcc{aa}],'LineWidth',0.5); %lower std error %marina
            xlim(x_lim_temp);
            if ~isempty(y_lim_ampl)
                ylim(y_lim_ampl);
            end
            hold on
            legendv(1,aa) = h1; %this is for working around the problem of the multiple legends for stde that we do not want..
            if S.wave_plot_conditions_together == 1
                legcell(1,aa) = {plotlab};
            else
                legcell(1,aa) = {S.conditions{aa}}; %same..
            end
        end
    end
    grid on
    grid minor
    set(gcf,'Color','w')
    if S.legc == 1
        legend([legendv],legcell) %same..
        legend('show');
    end
%     title('left averaged channels');
    %right MEG channels
%     legcell = cell(1,cond_n);
%     clf(figure(14))
%     figure(14)
%     for aa = 1:cond_n     
%         if S.wave_plot_conditions_together == 1
%             h1 = plot(time_sel,mean_right_mag(:,aa),color_line{aa},'LineWidth',2,'DisplayName',plotlab); %average waveform            
%         else
%             h1 = plot(time_sel,mean_right_mag(:,aa),color_line{aa},'LineWidth',2,'DisplayName',S.conditions{aa}); %average waveform
%         end  
%         xlim(x_lim_temp);
%         if ~isempty(y_lim_ampl)
%             ylim(y_lim_ampl);
%         end
%         hold on
%         plot(time_sel,mean_right_mag(:,aa) + stde_right_mag(:,aa),[':' color_line{aa}],'LineWidth',0.5); %upper std error
%         xlim(x_lim_temp);
%         if ~isempty(y_lim_ampl)
%             ylim(y_lim_ampl);
%         end
%         hold on
%         plot(time_sel,mean_right_mag(:,aa) - stde_right_mag(:,aa),[':' color_line{aa}],'LineWidth',0.5); %lower std error
%         xlim(x_lim_temp);
%         if ~isempty(y_lim_ampl)
%             ylim(y_lim_ampl);
%         end
%         hold on
%         legendv(1,aa) = h1; %this is for working around the problem of the multiple legends for stde that we do not want..
%         if S.wave_plot_conditions_together == 1
%             legcell(1,aa) = {plotlab};
%         else
%             legcell(1,aa) = {S.conditions{aa}}; %same..
%         end
%         grid on
%     end                   
%     legend([legendv],legcell) %same..
%     legend('show');
%     title('right averaged channels');
    if S.avewave_contrast == 1 %checking if the contrasts have been previously calculated
        if S.wave_plot_conditions_together == 0 %if you do not want conditions averaged together, it is assumed that you may want not only cond1 and cond2 independently but also their contrast
            if isempty(OUT.TSTAT_mag) %loading contrasts if they are not in the workspace (very likely)
                if length(S.cond_ttests_tobeplotted_topoplot) == 2 %if you want to compare two conditions
                    load([S.outdir '/' S.save_name_data '_OUT_' S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} '_vs_' S.conditions{S.cond_ttests_tobeplotted_topoplot(2)} '.mat'],'OUT');
                    disp(['loading contrast ' S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} '_vs_' S.conditions{S.cond_ttests_tobeplotted_topoplot(2)}])
                else %if you want to compare one condition with its own baseline
                    load([S.outdir '/' S.save_name_data '_OUT_' S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} '_vs_baseline.mat'],'OUT')
                    disp(['loading contrast ' S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} '_vs_ its own baseline'])
                end
                TSTAT_grad = OUT.TSTAT_grad;
                TSTAT_mag = OUT.TSTAT_mag;
                %extracting requested channels
                if mod(left_mag(1),2) == 1 %if indices are odd we have magnetometers
                    datdiff_l = TSTAT_mag(floor(left_mag/2 + 1),:);
%                     datdiff_r = TSTAT_mag(floor(right_mag/2 + 1),:);
                else %otherwise we have gradiometers
                    datdiff_l = TSTAT_grad(left_mag/2,:);
%                     datdiff_r = TSTAT_grad(right_mag/2);
                end
                %mean of t-values
                mddl = squeeze(nanmean(datdiff_l,1));
                mddl(isnan(mddl)) = 0; %here the first column may be NaN because it may have occurred to calculate t-tests with 0s against 0s.. therefore here the NaN values are converted to 0 for plotting purposes
%                 mddr = squeeze(mean(datdiff_r,1));
                %actual plotting
                clf(figure(15))
                figure(15)
                plot(time_sel,mddl,'k','LineWidth',2);
                xlim([time_sel(1) time_sel(end)])
                %significant time-points
                if ~isempty(S.signtp{1})
                    %plotting the significant time-points with gray shades
                    patch_color = [.85 .85 .85]; % Color of grey box marking time-ranges
                    ylims = get(gca,'YLim');
                    for ii = 1:length(S.signtp)
            %             sgf2 = S.signtp{ii}./S.sr;
                        sgf2 = S.signtp{ii};
                        patch([sgf2(1) sgf2(end) sgf2(end) sgf2(1)],[ylims(1) ylims(1) ylims(2) ylims(2)],patch_color,'EdgeColor','none','FaceAlpha',.5)
                        hold on
                    end
                end
                %replotting for better visualization
                plot(time_sel,mddl,'k','LineWidth',2);
                xlim([time_sel(1) time_sel(end)])
                %additional plottin details                
                grid on
                grid minor
                if S.legc == 1
                    legend('show');
                end
                if length(S.cond_ttests_tobeplotted_topoplot) == 2 %if you want to compare two conditions
                                    
                    title([S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} ' vs ' S.conditions{S.cond_ttests_tobeplotted_topoplot(2)}])
                else
                    title([S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} ' vs its own baseline'])
                end
                set(gcf,'Color','w')
%                 clf(figure(16))
%                 figure(16)
%                 plot(time_sel,mddr,'k','LineWidth',2,'DisplayName','right hem');
%                 grid on
%                 legend('show');
%                 title([S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} ' vs ' S.conditions{S.cond_ttests_tobeplotted_topoplot(2)}])
            end
        end
    end
    %saving if requested
    if S.save_label_waveaverage == 1
        if S.mag_lab == 1
            mgh = 'mag';
        else
            mgh = 'grad';
        end
        saveas(figure(13),[S.outdir '/wave_' mgh '_' S.label_plot '.jpg' ],'jpg');
%         saveas(figure(14),[S.outdir '/wave_' mgh '_right_' S.label_plot '.jpg' ],'jpg');
        if ~isempty(S.cond_ttests_tobeplotted_topoplot)
            if S.wave_plot_conditions_together == 0
                saveas(figure(15),[S.outdir '/wave_diff_' mgh '_' S.label_plot '.jpg' ],'jpg');
%                 saveas(figure(16),[S.outdir '/wave_diff_' mgh '_right_' S.label_plot '.jpg' ],'jpg');
            end
        end
    end
end


%topoplot
if S.topoplot_label == 1
    if S.topocontr == 1 %loading contrasts if you want contrasts (otherwise later it will be calculated a simple mean difference..)
        if isempty(OUT.TSTAT_mag)     
           if length(S.cond_ttests_tobeplotted_topoplot) == 2 %if you want to compare two conditions
               load([S.outdir '/' S.save_name_data '_OUT_' S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} '_vs_' S.conditions{S.cond_ttests_tobeplotted_topoplot(2)} '.mat'],'OUT');
               disp(['loading contrast ' S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} '_vs_' S.conditions{S.cond_ttests_tobeplotted_topoplot(2)}])
           else %if you want to compare one condition with its own baseline
               load([S.outdir '/' S.save_name_data '_OUT_' S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} '_vs_baseline.mat'],'OUT')
               disp(['loading contrast ' S.conditions{S.cond_ttests_tobeplotted_topoplot(1)} '_vs_ its own baseline'])
           end
            TSTAT_grad = OUT.TSTAT_grad;
            TSTAT_mag = OUT.TSTAT_mag;
        end     
    end
    Sz = size(data_mat);
    %loading an example from a previous fieldtrip file obtained after some pre-processing
    fieldtrip_example = load([S.fieldtrip_mask '/fieldmask.mat']);
    %taking new labels from the fieldtrip example
    label = fieldtrip_example.M_timb_lock_gradplanarComb.label;
    %setting labels according to the ones (presumibly) in the layout
    for ii = 1:204
        label{ii,1} = [label{ii,1}(1:3) label{ii,1}(5:end)];
    end
    dummy_time_window = length(time_sel);
    %defining if user wants to plot a single condition or the difference (t-values) between two specific conditions
    avg = zeros(Sz(1),Sz(2));
    if S.topocontr == 1 %here we do not move GRAD and MAG since we already calculated the t-tests earlier and separately stored the GRAD and MAG (if you want the contrasts)
        avg(1:102,:) = TSTAT_grad;
        avg(103:204,:) = TSTAT_mag;
    else
        diff_matrix_dum = squeeze(nanmean(data_mat,3));
        if isempty(S.topocondsing) %otherwise if you want the mean across all conditions
            diff_matrix = nanmean(diff_matrix_dum,3);
        else %or a specific condition only
            diff_matrix = diff_matrix_dum(:,:,S.topocondsing);
        end
        %moves the gradiometers
        count = 0;
        for ii = 1:length(label)/2
            count = count + 2;
            for kk = 1:dummy_time_window %1:length(D.time)
                avg(ii,kk) = diff_matrix(count,kk); %D((count),kk,5);
            end
        end
        %moves the magnetometers
        count = -1;
        for ii = 103:length(label)
            count = count + 2;
            for kk = 1:dummy_time_window %1:length(D.time)
                avg(ii,kk) = diff_matrix(count,kk); %D((count),kk,5);
            end
        end
    end
 
    %actual plotting
    %creating the mask for the data
    cfgdummy = fieldtrip_example.M_timb_lock_gradplanarComb.cfg;
    cfgdummy.previous = [];
    data.cfg = cfgdummy;
    data.time(1,:) = time_sel(:);
    data.label = label;
    data.dimord = fieldtrip_example.M_timb_lock_gradplanarComb.dimord;
    data.grad = fieldtrip_example.M_timb_lock_gradplanarComb.grad;
    data.avg = avg;
    %creating the cfg for the actual plotting (magnetometers)
    cfg = [];
    cfg.layout = [S.fieldtrip_mask '/neuromag306mag.lay'];
    cfg.colorbar = 'yes';
    cfg.xlim = S.xlim;
    if ~isempty(S.zlimmag)
        cfg.zlim = S.zlimmag;
    end
    cfg.colormap = S.colormap_spec;
    %clear and open the figure
    clf(figure(17));
    figure(17);
    ft_topoplotER(cfg,data);
    set(gcf,'Color','w')
    %saving if requested
    if S.topoplot_save_label == 1
        saveas(figure(17),[S.outdir '/' S.label_plot '_' num2str(S.xlim(1)) '_' num2str(S.xlim(2)) '_mag.jpg' ],'jpg');
    end
    %creating the cfg for the actual plotting (gradiometers)
    cfg = [];
    cfg.layout = [S.fieldtrip_mask '/neuromag306cmb.lay'];
    cfg.colorbar = 'yes';
    cfg.xlim = S.xlim;
    if ~isempty(S.zlimgrad)
        cfg.zlim = S.zlimgrad;
    end
    cfg.colormap = S.colormap_spec;
    %clear and open the figure
    clf(figure(18));
    figure(18);
    ft_topoplotER(cfg,data);
    set(gcf,'Color','w')
    %saving if requested
    if S.topoplot_save_label == 1
        saveas(figure(18),[S.outdir '/' S.label_plot '_' num2str(S.xlim(1)) '_' num2str(S.xlim(2)) '_cmb.jpg' ],'jpg');
    end
end




end


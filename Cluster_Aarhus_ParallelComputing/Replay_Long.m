function O = Replay_Long(labell)
O = []; 


if labell == 1
    load('/scratch7/MINDLAB2022_MEG-EncodingMusicSeq/gemma/decoding/Localizer_19_20time.mat'); %time in seconds
    list_TG = dir(['/scratch7/MINDLAB2022_MEG-EncodingMusicSeq/gemma/decoding/Localizer_1_2/TG*']); %list of TG files
    dummean = zeros(length(time_sel),length(time_sel),20); %preallocate variable dd1 (time samples and TG files)
    for tt = 1:20 %over tones
        cnt = 0;
        dd2 = zeros(length(time_sel),length(time_sel),length((tt+1):20),length(list_TG)); %preallocate variable dd1 (time samples and TG files)
        for nn = (tt+1):20 %over the other tones
            cnt = cnt + 1;
            list_TG = dir(['/scratch7/MINDLAB2022_MEG-EncodingMusicSeq/gemma/decoding/Localizer_' num2str(tt) '_' num2str(nn) '/TG*']); %list of TG files
            for ff = 1:length(list_TG) %over subjects
                load([list_TG(ff).folder '/' list_TG(ff).name]); %load TG file
                dd2(:,:,cnt,ff) = d.d;
            end
        end
        dd2(dd2==0) = NaN;
        dummean(:,:,tt) = nanmean(nanmean(dd2,4),3); %mean over participants and tones nn
        disp(tt)
    end
    dd2TGm = nanmean(dummean,3); %average participants
    figure; imagesc(time_sel,time_sel,dd2TGm); xlabel('Train time (s)'); ylabel('Test time (s)'); set(gca,'YDir','normal'); %settings for figure
    colorbar %add colorbar
    set(gcf,'Color','w')
    save('/scratch7/MINDLAB2022_MEG-EncodingMusicSeq/gemma/decoding/TG_Localizer_average.mat','dd2TGm');
%     export_fig('/scratch7/MINDLAB2022_MEG-EncodingMusicSeq/gemma/loc_alltones.png')
%     saveas(gcf,'/scratch7/MINDLAB2022_MEG-EncodingMusicSeq/gemma/loc_alltones2.png')
elseif labell == 2
    
    block = 1; %1 = auditory vs visual numbers; 2 = auditory numbers vs rest; 3 = visual numbers vs rest
    
    if block == 1
        list_TG = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/decoding_replay/Aud_vs_Vis/TG*.mat');
    elseif block == 2
        list_TG = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/decoding_replay/Aud_vs_Rest/TG*.mat');
    elseif block == 3
        list_TG = dir('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/decoding_replayVis_vs_Rest/TG*.mat');
    end
    load('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/decoding_replay/Aud_vs_Resttime.mat');
    dd1 = zeros(length(time_sel),length(time_sel),length(list_TG)); %preallocate variable dd1 (time samples and TG files)
    for ff = 1:length(list_TG) %over TG files
        load([list_TG(ff).folder '/' list_TG(ff).name]); %load TG file
        dd1(:,:,ff) = d.d; %save d structure
        disp(ff)
    end
    ddTGm = mean(dd1,3); %average participants
    figure; imagesc(time_sel,time_sel,ddTGm); xlabel('Train time (s)'); ylabel('Test time (s)'); set(gca,'YDir','normal'); %settings for figure
    colorbar %add colorbar
    set(gcf,'Color','w')
    if block == 1
            save('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/decoding_replay/Aud_vs_Vis/TG_MemList_average.mat','ddTGm');
%         export_fig('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/decoding_replay/Aud_vs_Vis.png')
%         saveas(gcf,'/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/decoding_replay/Aud_vs_Vis.png')
    elseif block == 2
        save('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/decoding_replay/Aud_vs_Rest/TG_MemList_average.mat','ddTGm');
%         export_fig('/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/decoding_replay/Aud_vs_Rest.png')
%         saveas(gcf,'/scratch7/MINDLAB2021_MEG-TempSeqAges/leonardo/after_maxfilter/decoding_replay/Aud_vs_Rest.png')
    end
elseif labell == 3
    list_TG = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/decoding_replay_samemel_visualpat/Aud_vs_Vis/TG*.mat');
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/decoding_replay_samemel_visualpat/Aud_vs_Vistime.mat'); %time in seconds
    dd1 = zeros(length(time_sel),length(time_sel),length(list_TG)); %preallocate variable dd1 (time samples and TG files)
    for ff = 1:length(list_TG) %over TG files
        load([list_TG(ff).folder '/' list_TG(ff).name]); %load TG file
        dd1(:,:,ff) = d.d; %save d structure
        disp(ff)
    end
    ddTGm = mean(dd1,3); %average participants
    save('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/decoding_replay_samemel_visualpat/TG_average.mat','ddTGm');
end


end




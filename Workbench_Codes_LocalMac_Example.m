%% *** Temporal pattern recognition in the human brain: a dual simultaneous processing *** %%

% EXAMPLE OF AUTOMATIC CREATION OF LEFT AND RIGHT HEMISPHERE IMAGES FOR WORKBENCH

% Leonardo Bonetti
%Aarhus, DK, 29/01/2021


%% CREATING ONE IMAGE PER HEMISPHERE USING WORKBENCH BRAIN TEMPLATES

%defining all paths where you have nifti images to be used by wb_command
paths{1} = '/Users/au550322/Documents/AarhusUniversitet/PhDThesis/PaperFinal3_PCAKmeans/Figures_Tables/Figures/Figure_1_S3';
paths{2} = '/Users/au550322/Documents/AarhusUniversitet/PhDThesis/PaperFinal3_PCAKmeans/Figures_Tables/Figures/Figure_1_S3/SpatioTemporalParcels_TempClust3';
paths{3} = '/Users/au550322/Documents/AarhusUniversitet/PhDThesis/PaperFinal3_PCAKmeans/Figures_Tables/Figures/Figure_2';
paths{4} = '/Users/au550322/Documents/AarhusUniversitet/PhDThesis/PaperFinal3_PCAKmeans/Figures_Tables/Figures/Figure_4/Brain_Parcels';
paths{5} = '/Users/au550322/Documents/AarhusUniversitet/PhDThesis/PaperFinal3_PCAKmeans/Figures_Tables/Figures/Figure_S2';
paths{6} = '/Users/au550322/Documents/AarhusUniversitet/PhDThesis/PaperFinal3_PCAKmeans/Figures_Tables/Figures/Figure_S4';
paths{7} = '/Users/au550322/Documents/AarhusUniversitet/PhDThesis/PaperFinal3_PCAKmeans/Figures_Tables/Figures/Figure_S5';
paths{8} = '/Users/au550322/Documents/AarhusUniversitet/PhDThesis/PaperFinal3_PCAKmeans/Figures_Tables/Figures/Figure_S6/Brain_Parcels';
paths{9} = '/Users/au550322/Documents/AarhusUniversitet/PhDThesis/PaperFinal3_PCAKmeans/Figures_Tables/Figures/Figure_S7/Parcels_NEWCondition__OLDConditionParcelsInFigure_S5';
paths{10} = '/Users/au550322/Downloads/MCS_Figure1F';

%masks to be used
mask_left = '/Users/au550322/Documents/osl/std_masks/ParcellationPilot.L.midthickness.32k_fs_LR.surf.gii';
mask_right = '/Users/au550322/Documents/osl/std_masks/ParcellationPilot.R.midthickness.32k_fs_LR.surf.gii';

%actual computation
for ii = 10 %over paths
    list = dir([paths{ii} '/*.gz']); %list of nifti file within path ii
    for jj = 1:length(list) %over images within path ii
        inima = [list(jj).folder '/' list(jj).name]; %input name
        %left hemisphere
        outimaleft = [list(jj).folder '/' list(jj).name(1:end-7) '_L.func.gii']; %output image
        %call to wb_command
        cmd = ['/Users/au550322/Applications/workbench/bin_macosx64/wb_command -volume-to-surface-mapping ' inima ' ' mask_left ' ' outimaleft ' -trilinear'];
        system(cmd) %running the job
        %right hemisphere
        outimaright = [list(jj).folder '/' list(jj).name(1:end-7) '_R.func.gii']; %output image
        %call to wb_command
        cmd = ['/Users/au550322/Applications/workbench/bin_macosx64/wb_command -volume-to-surface-mapping ' inima ' ' mask_right ' ' outimaright ' -trilinear'];
        system(cmd) %running the job
    end
    disp(ii)
end


%%

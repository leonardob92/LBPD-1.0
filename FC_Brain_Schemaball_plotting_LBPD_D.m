function [  ] = FC_Brain_Schemaball_plotting_LBPD_D( S )

% Given a connectivity matrix and the corresponding parcellation, it uses
% oslspinningbrain and schemaball for plotting connectivity within the
% brain (both video and specific frames) and as schemaball.


%   INPUT:  -S.MATF:        connectivity matrix ROIs x ROIs (can be both
%                           binary or non-binary) (double)
%           -S.symmetric_l: 1 if MATF is LLLRRR
%                           0 if MATF is LRLRLR
%           -S.outpath:     path where saving images (characters)
%           -S.p:           parcellation object (created by using OSL)
%           -S.thr_cortex:  strongest % of connectivity to be plotted (it
%                           is usually a quite high number, e.g. 0.993).
%                           If empty [] default is 0.993.
%           -S.name_gif:    name of gif/video (character)
%                           Leave empty [] if you do not want to plot the
%                           connectivity in the brain (e.g. you want only
%                           the schemaball).
%           -S.frame_vect:  vector with frames you want.
%                           If empty [] default is [15,90], another good
%                           option could be [1,15,59,90,130]
%           -S.fr_spec:     cell array with characters specifing the name
%                           of the requested frames.
%                           E.g. [15,90] corresponds to fr_spec = {'Posterior_L_R';'Frontal_R_L'}
%                           [1,15,59,90] corresponds to fr_spec = {'Orig_Angle';'Posterior_L_R';'Hem_R';'Frontal_R_L';'Hem_L'}
%           -S.spinlim:     limits for spinning brain (e.g. [-1 1].
%                           If you want automatic calculation of limits, do
%                           not provide this field or leave it empty.
%           -S.schball_l:   1 for having the schemaball
%                           0 for not having it.
%           -S.extr:        minimum and maximum values to be used for
%                           normalizing the matrix to be submitted to
%                           schemaball function.
%                           Leave empty [] for scaling the matrix on its
%                           max(abs value).
%           -S.lab_parc_ph: path to labels.
%                           It must be a .mat file with a character array
%                           named 'lab' (characters).
%           -S.colbacksch:  background color for shcemaball.
%                           If empty [], default is 'w', white.
%           -S.schemplth:   percentage of the strongest connections to be
%                           plotted in schemaball (5% to be submitted as
%                           '5').
%                           Percentage of the actual connections (0s (null connections between ROIs) are
%                           removed before calculating x%..
%           -S.name_sch:    name for saved schemaball image.
%                           Leave empty [] for not saving it.

%   OUTPUT: -plots..





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% leonardo.bonetti@clin.au.dk
% Leonardo Bonetti, MIT, Boston, USA, 13/08/2019
% Leonardo Bonetti, Aarhus, DK, 27/08/2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%preparing data
SS = size(S.MATF);
order = [1:2:SS(1) SS(1):-2:2];
if S.symmetric_l == 0
    MATFbrain = S.MATF;
    MATFschemb = S.MATF(order,order);
else
    MATFschemb = S.MATF;
    %the following two lines make the order LRLRLR from LLLRRR; they should be fine but an extra check may be a good idea
    [~,I] = sort(order);
    MATFbrain = S.MATF(I,I);     
end

%plotting connectivity in the brain
if ~isempty(S.name_gif) %if you want to plot (and presumably save) connectivity in the brain
    %assigning default if value is not specified by user
    if isempty(S.thr_cortex)
        S.thr_cortex = 0.993;
    end
    if isempty(S.frame_vect) %if it is empty assigning two specific frames by default
        S.frame_vect = [15,90];
    end
    %plotting network using an OSL function related to the object p
    %parcellation
    p = S.p;
    [~,h_scatter] = p.plot_network(MATFbrain,S.thr_cortex);
    set(h_scatter,'SizeData',15); %for now keeping fixed these parameters..
    if isfield(S,'spinlim') && ~isempty(S.spinlim)
        set(gca,'CLim',S.spinlim)
    end
    % colormap('jet')
    colormap(bluewhitered) %same here
    %creating the directory if it does not exist already
    if ~exist(S.outpath,'dir')
        mkdir(S.outpath);
        cd(S.outpath);
    else
        cd(S.outpath);
    end
    %creating gif (video) by using an OSL function
    osl_spinning_brain([S.outpath '/' S.name_gif '.gif']) %for now I assume you always want both video and frames of connectivity in the brain, even if you may want to simply plot the first frame of it only..
    
    %getting specific images corresponding to different frames of the gif
    obj = gifread([S.outpath '/' S.name_gif '.gif']);
    for ii = 1:length(S.frame_vect)
%         figure(ii)
%         clf(figure(ii))
        figure
        imshow(obj(:,:,:,S.frame_vect(ii)))
        export_fig([S.outpath '/' S.name_gif '_fig_' S.fr_spec{ii} '.eps']) %here you may want to either use this export_fig or the normal saveas function.. depending on your matlab/server, etc. it seems to produce higher quality images one or the other function.. I do not know why, so I just left export_fig with saveas as comment
    %     saveas(figure(ii),[S.outpath '/' S.name_gif '_fig_' fr_spec{ii} '.jpg'],'jpg');
    end
end

%schemaball
if S.schball_l == 1
    %if background is not specified, assigning white background by default
    if isempty(S.colbacksch)
        S.colbacksch = 'w';
    end  
    load(S.lab_parc_ph); %labels submitted in LRLRLR order
    if isempty(S.extr) 
        matrx_z_norm = MATFschemb/(max(max(abs(MATFschemb)))); %normalizing data matrix with its maximum value
    elseif isempty(S.extr)
        matrx_z_norm = MATFschemb/(max(max(abs(MATFschemb)))); %normalizing data matrix with its maximum value
    else
        if max(abs(S.extr)) < max(max(abs(MATFschemb)))
            error('the maximum extreme (in absolute value) must be >= than the maximum (absolute value) of the data matrix')
        end
        matrx_z_norm = MATFschemb/(max(abs(S.extr))); %otherwise with values given by the user
    end
    labsymm = lab(order,:); %ordering labels
    %plotting only the strongest connections
    dd = triu(matrx_z_norm,1); %getting only upper triangle
    dd2 = reshape(dd,[SS(1)*SS(1),1]); %reshaping the matrix into a vector
    dd21 = dd2(dd2~=0); %removing 0s
    dumabsd = sort(abs(dd21),'descend'); %sorting the absolute values for getting maximum connections independently by the sign
    dumthr = dumabsd(floor(length(dumabsd)*S.schemplth/100)); %getting the x% threshold (S.schemplth)        
    matf1 = zeros(SS(1),SS(1)); %creating variable
    matf1(matrx_z_norm >= dumthr) = matrx_z_norm(matrx_z_norm >= dumthr); %getting values higher than threshold (positive connections)
    matf1(matrx_z_norm <= ((-1) * dumthr)) = matrx_z_norm(matrx_z_norm <= ((-1) * dumthr)); %getting values lower then threshold * (-1) (negative connections)    
    matf1(matf1 == 0) = NaN; %assigning NaN to 0s..
    %         schemaball(matrx_z_norm,LBLS); %original
    %         schemaball(matf1,LBLS,[1 0 0; 0 0 1]);
    %trick to have left ROIs in the "left" part of the sphere.. (ball..)
    %labels
    labsymm2 = labsymm;
    labsymm2(1:45,:) = labsymm(46:90,:);
    labsymm2(46:90,:) = labsymm(1:45,:);
    %data
    matf2 = matf1;
    matf2(1:45,1:45) = matf1(46:90,46:90);
    matf2(46:90,46:90) = matf1(1:45,1:45);
    matf2(1:45,46:90) = matf1(46:90,1:45);
    matf2(46:90,1:45) = matf1(1:45,46:90);
    schemaball(matf2,labsymm2,[0 0 1; 1 0 0],[1 0 0]); %using external function 'schemaball' with a couple of small modifications of mine
    %original call with right hemisphere ROIs in the "right" part f the sphere..
    % schemaball(matf1,labsymm,[0 0 1; 1 0 0],[1 0 0]); %using external function 'schemaball' with a couple of small modifications of mine
    set(gca,'Color',S.colbacksch)
    colorbar
    if ~isempty(S.name_sch) %saving option for schemaball
        export_fig([S.outpath '/' S.name_sch '.eps'])
    end
end   



end


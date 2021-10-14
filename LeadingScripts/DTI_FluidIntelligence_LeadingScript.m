%% On the brain networks organization of individuals with high versus average fluid intelligence: a combined DTI and MEG study


%%% DTI - STRUCTURAL CONNECTIVITY %%%


%BEFORE PROCEEDING, PLEASE NOTE:

%As follows, every analysis that we made has been reported to ensure full disclosure.
%Please, note that in this Leading script, I use built-in functions,
%functions from well-known toolboxes (e.g. OSL, SPM, FieldTrip) and
%in-house-built functions, which are reported together with this script in the collection named LBPD.
%Data can be provided according to the Danish regulation, upon reasonable request.
%If interested, please contact Leonardo Bonetti, leonardo.bonetti@clin.au.dk
%More information is provided in the ReadMe.mat file that I strongly advise to read.


%Leonardo Bonetti, Silvia Elisabetta Portis Bruzzone
%leonardo.bonetti@psych.ox.ac.uk
%silviaepbruzzone@clin.au.dk


%% PIPELINE OF ALL SUBJECTS - DTI - FSL %%

%In the first part of the script, we have used a standard FSL pipeline
%adapted to the cluster of computers of MINDLAB, MIB-CFIN, Aarhus University.

%%OBS!! remember to open the terminal, write "use anaconda", then open matlab in the same terminal

%% INITIAL CONVERSION FROM DICOM (RAW) FILE TO NIFTI (dcm2nii) - fast

path='/projects/MINDLAB2017_MEG-LearningBach/raw'; %path to the folder containing all the raw data
path_out='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/'; %set the path to the folder where you will save your converted data
path1=dir([path '/00*']); %query the content of the folder "raw": it contains one folder for each subject
for ii=1:length(path1) %loop across all subjects
    path2=dir([path1(ii).folder '/' path1(ii).name '/20*']); %for each subject, search folders named '2018*'. There should be 1 or 2 folders per participant: one contains the folder 'MR' and 'SR' the other one contains the folder 'MEG'. We need the folder 'MR'.
    mkdir(path_out, path1(ii).name);  %for each subject (ii), create a new folder in (DTI_Portis) and call it with the same name as the subject's (path1(ii).name)
    for dd=1:length(path2)
        path3=dir([path2(dd).folder '/' path2(dd).name '/MR/014*AP']); %choose 'MR' rather than 'SR' - AP
        path4=dir([path2(dd).folder '/' path2(dd).name '/MR/023*PA']); %choose 'MR' rather than 'SR' - PA
        if ~isempty(path3)
            pathAP=[path3.folder '/' path3.name '/files']; %set path to diff data (AP)
            %Convert AP to .nii.gz
            cmd = ['dcm2nii -o ' path_out path1(ii).name ' ' pathAP]; %local
            system(cmd) %send command to the cluster
        end
        if ~isempty(path4)
            pathPA=[path4.folder '/' path4.name '/files'];%set path to PA
            %Convert PA
            cmd = ['dcm2nii -o ' path_out path1(ii).name ' ' pathPA];
            system(cmd)
        end
    end
end
% ! Subjects 7,33 and 41 do not have the MRI files: the corresponding folders were deleted

%% Getting ready for TOPUP

path_out='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/'; %current folder
cd(path_out);

% 1) NODIF (fslroi) - fast
%creating reference volume (nodif) based on the first image of the DTI data 
path5=dir([path_out '00*']);
for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder (specify 20* to avoid hidden files)
    cd([path5(gg).folder '/' path5(gg).name]); %change directory for each subject, so that nodif is saved in that of the corresponding subject
    cmd = ['fslroi ' path6(1).folder '/' path6(3).name ' nodif 0 1']; %apply fslroi to AP, which is the 3rd file in the folder
    system(cmd)
    disp(gg)
end 
%fslroi command
%input DTI data (AP)
%output name of the reference volume (nodif)
%0 (minimum index of the volumes within the nifti file (indexing as in python=starting from 0))
%1 (number of images to be taken (in this case only first))

% 2) NODIF_PA (copyfile) - fast
%renaming PA as follows: 
%1)copy the DTI PA file (20180123_102505cmrrmbep2ddiffmultidirPASubjectNo0001s023a1001.nii.gz)
%2)rename it as "nodif_PA.nii.gz" 
for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    cd([path5(gg).folder '/' path5(gg).name]); %change directory for each subject, so that a nodif_PA.nii.gz is created for each subject
    copyfile([path6(4).folder '/' path6(4).name],[path6(4).folder '/nodif_PA.nii.gz']); %copy PA (4th file in the folder) and rename it 
    disp(gg)
end

% 3) AP_PA_b0 (fslmerge) - fast
%merging AP with PA data with the reference volume (nodif) 
for mm=1:length(path5)
    path6=dir([path5(mm).folder '/' path5(mm).name '/20*']); %loop across each subject and get the content of each subject's folder
    cd([path5(mm).folder '/' path5(mm).name]); %change directory for each subject, so that a nodif_PA.nii.gz is created for each subject  
    cmd = ['fslmerge -t ' path6(4).folder '/AP_PA_b0 nodif nodif_PA'];
    system(cmd)
    disp(mm)
end
%fslmerge command
%-t (parameter of fslmerge)
%output name (AP_PA_b0)
%reference AP (nodif)
%reference PA (nodif_PA)

%% Creating acquisition parameters file (acqparams.txt) - To be run ONLY for the FIRST SUBJECT!!! - fast
%txt file with information about acquisition parameters (the same for all subjects)
%In this case: 2x4 matrix 
%row 1: 0 -1 (meaning AP) 0 (from the beginning of time) 0.104 (echo time in seconds (104ms))

matrixacqpar = [0 -1 0 0.104; 0 1 0 0.104]; %creating double matrix with values
t = table(matrixacqpar); %converting matrix to table
writetable(t,'acqparams.txt') %saving table as .txt file (acqparams.txt MUST be the name)

%% Copying acqparams.txt to each subject's folder - To be run for ALL THE OTHER SUBJECTS

path_acqparams='/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/DTI001/acqparams.txt';
path_out='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/'; %set the path to the folder where you will save your converted data
path5=dir([path_out '00*']);
for gg=1:length(path5) 
    copyfile(path_acqparams,[path5(gg).folder '/' path5(gg).name '/acqparams.txt']); %copy acqparams.txt from DTI001 in every subject's folder  
    disp(gg)
end

%% Call to TOPUP command (topup) - quite fast

path_out='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/'; %set the path to the folder where you will save your converted data
path5=dir([path_out '00*']);
for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    %cmd = 'topup --imain=AP_PA_b0 --datain=acqparams.txt --config=b02b0.cnf --out=topup_AP_PA_b0 --iout=topup_AP_PA_bo_iout --fout=topup_AP_PA_bo_fout'; %apply topup to every subject
    cmd = ['submit_to_cluster -q all.q -n 1 -p MINDLAB2017_MEG-LearningBach "/usr/local/fsl/bin/topup --imain=' path6(1).folder '/AP_PA_b0 --datain=' path6(1).folder '/acqparams.txt --config=b02b0.cnf --out=' path6(1).folder '/topup_AP_PA_b0 --iout=' path6(1).folder '/topup_AP_PA_b0_iout --fout=' path6(1).folder '/topup_AP_PA_b0_fout"']; %apply topup to every subject    
    system(cmd)
    disp(gg)
end
%cmd = 'topup --imain=AP_PA_b0 --datain=acqparams.txt --config=b02b0.cnf --out=topup_AP_PA_b0 --iout=topup_AP_PA_bo_iout --fout=topup_AP_PA_bo_fout';
%topup command
%DTI data after fslmerg (--imain=AP_PA_b0)
%acquisition parameters (--datain=acqparams.txt)
%configuration file specifying command line arguments (--datain=acqparams.txt); this file is already in the FSL directory and it contains some (default) specifications
%output file (--out=topup_AP_PA_b0)
%output file with unwarped images (--iout=topup_AP_PA_bo_iout)
%output file with field (Hz) (--fout=topup_AP_PA_bo_fout)

%% Generating a brain mask from the corrected b0 (fslmaths, bet) - fast

path_out='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/'; %set the path to the folder where you will save your converted data
path5=dir([path_out '00*']);

for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    cd([path5(gg).folder '/' path5(gg).name]); %change directory for each subject, so that a nodif_PA.nii.gz is created for each subject   
    cmd = 'fslmaths topup_AP_PA_b0_iout.nii.gz -Tmean hifi_nodif'; %apply fslmaths to create a mask from the corrected b0   
    system(cmd)
    disp(gg)
end
% cmd = 'fslmaths topup_AP_PA_b0_iout.nii.gz -Tmean hifi_nodif';
%fslmaths command
%input (topup_AP_PA_b0_iout.nii.gz)
%specification -Tmean
%output (hifi_nodif)
% system(cmd)

%extracting the brain from b0 - fast
for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    cd([path5(gg).folder '/' path5(gg).name]); %change directory for each subject, so that a nodif_PA.nii.gz is created for each subject
    cmd = 'bet hifi_nodif hifi_nodif_brain -m -f 0.2'; %apply BET to extract the brain from b0   
    system(cmd)
    disp(gg)
end
% cmd = 'bet hifi_nodif hifi_nodif_brain -m -f 0.2';
%bet command
%input brain mask (hifi_nodif)
%output (hifi_nodif_brain)
%parameters (-m -f 0.2); -f means "fraction intensity threshold"
% system(cmd)

%% Creating index file - to be run ONLY THE 1st TIME!
%.txt file with 1s on the first column
%as many 1s as number of volumes (images in the DICOM format) (in this case 211)

ind = ones(211,1); %creating the double vector
t = table(ind); %converting vector to table
writetable(t,'index.txt') %saving table as .txt file (acqparams.txt MUST be the name)
%OBS!! delete the first row which indicates the legend of the table ("ind") and save the file!

%% Copying index.txt to every subject's folder - to be run ALL THE OTHER TIMES

path_out='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/'; %set the path to the folder where you will save your converted data
path5=dir([path_out '00*']);
path_index='/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/DTI001/index.txt';

for gg=1:length(path5)  
    copyfile(path_index,[path5(gg).folder '/' path5(gg).name '/index.txt']); %copy index.txt from DTI001 into every subject's folder
    disp(gg)
end

%% Correcting for EDDY currents (eddy) - slow: about 10-15 hours/subject
%(currents generated in the MRI machine because of a rapid change of the magnetic field direction during the acquisition: echo planar images are acquired rapidly in different orientations)

path_out='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/';
path5=dir([path_out '00*']);

for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    cmd = ['submit_to_cluster -q long.q -n 1 -p MINDLAB2017_MEG-LearningBach "/usr/local/fsl/bin/eddy --imain=' path6(3).folder '/' path6(3).name ' --mask=' path6(1).folder '/hifi_nodif_brain_mask --index=' path6(1).folder '/index.txt --acqp=' path6(1).folder '/acqparams.txt --bvecs=' path6(2).folder '/' path6(2).name ' --bvals=' path6(1).folder '/' path6(1).name ' --fwhm=0 --topup=' path6(1).folder '/topup_AP_PA_b0 --flm=quadratic --out=' path6(3).folder '/eddy_unwarped_images"'];
    system(cmd)
    disp(gg)
end
% cmd = 'eddy --imain=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.nii.gz --mask=hifi_nodif_brain_mask --index=index.txt --acqp=acqparams.txt --bvecs=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bvec --bvals=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bval --fwhm=0 --topup=topup_AP_PA_b0 --flm=quadratic --out=eddy_unwarped_images';
%command eddy
%main input DTI data from initial nifti convertion (--imain=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.nii.gz)
%mask outputted from BET (hifi_nodif_brain_mask); NOTE that the name of the output that we use is "hifi_nodif_brain_mask" and not "hifi_nodif_brain"
%file with indices of the volumes (--index=index.txt)
%acquisition parameters (--acqp=acqparams.txt)
%bvecs outputted from the initial nifti convertion (--bvecs=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bvec); bvecs indicate the direction of diffusion
%bvals outputted from the initial nifti convertion (--bvals=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bval); bvals indicate the amount of the diffusion
%parameter telling that no smoothing should be applied (--fwhm=0)
%output from TOPUP algorithm (--topup=topup_AP_PA_b0)
%parameter that assumes a quadratic model for the EC-fields (--flm=quadratic)
%output name (--out=eddy_unwarped_images)
% system(cmd)

%% TENSOR FITTING - fast (for later TBSS)
%fitting the tensor into the DTI data (this is to get fractional anisotropy (FA) to see whether there are microstructural changes or differences between groups, etc. this is usually used in connection with TBSS)

path_out='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/';
path5=dir([path_out '00*']);
for gg=1:length(path5)
    path6=dir([path5(gg).folder '/' path5(gg).name '/20*']); %loop across each subject and get the content of each subject's folder
    cd([path5(gg).folder '/' path5(gg).name]); %change directory for each subject, so that a nodif_PA.nii.gz is created for each subject
    cmd = ['dtifit --data=' path6(1).folder '/eddy_unwarped_images.nii.gz --mask=' path6(1).folder '/hifi_nodif_brain_mask.nii.gz --bvecs=' path6(2).folder '/' path6(2).name ' --bvals=' path6(1).folder '/' path6(1).name ' --out=' path6(1).folder '/dti_fitted_tensors'];   
    system(cmd)
    disp(gg)
end
% cmd = ['dtifit --data=eddy_unwarped_images.nii.gz --mask=hifi_nodif_brain_mask.nii.gz --bvecs=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bvec --bvals=20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bval --out=dti_fitted_tensors'];
%dtifit command
%input data from eddy (--data=eddy_unwarped_images.nii.gz)
%mask that you got before (--mask=hifi_nodif_brain_mask.nii.gz)
%bvecs (same as for eddy)
%bvals (same as for eddy)
%output file (--out=dti_fitted_tensors)
% system(cmd)
%returns several outputs (starting with "dti_fitted_tensors..", NOTE in particular the ones ending with "_V1" "_V2" "_V3" and open them in fsleyes to conduct a further inspection

%% TBSS - Tract-Based Spatial Statistics

%To compare FA values between two groups of subjects
%ALL the FA files (one per subject) must be in the SAME directory
path_tbss='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/tbss_2mm'; %set the path to the folder where you will save your converted data
path_out='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/'; %path to the folder containing all the dti data that are being analyzed, divided by subject
path5=dir([path_out '00*']);
%0) Copy and rename *_FA.nii.gz files from each subject into a new folder, where tbss will be run - fast
for gg=1:length(path5)
    path_FA=dir([path5(gg).folder '/' path5(gg).name '/*FA.nii.gz']); %loop across each subject and get the content of each subject's FA file
    copyfile([path_FA(1).folder '/' path_FA(1).name],[path_tbss '/' path5(gg).name  '_FA.nii.gz']); %copy acqparams.txt from DTI001 in every subject's folder
    disp(gg)
end

%% 1) Removing likely outliers by removing brain-edge artefacts and the zero end slices (tbss_1_preproc) - FAST process

% tbss_1_preproc *nii.gz 
%It creates a new folder called FA, containing the sub-directory "origdata", with the original images
cd(path_tbss); %set tbss folder as current path
cmd = ['tbss_1_preproc *nii.gz']; %run tbss preprocessing
system(cmd)
% index.html shows the slices of every single subject

%% 2) Non-linear registration: aligning all the FA data across subjects (tbss_2_reg ) - moderately to very SLOW process (see options)
% tbss_2_reg 
% It estimates warping parameters to standardize space (later applied in step 3)
% Two ways to do it:
% a) align to FSL's "FMRIB58_FA" Template (about 10 min/participant) - RECOMMENDED by FSL guys. Use -T flag
% b) automatic search for the most representative image to be used as a template across subjects (days/weeks). Use -t flag to choose your own target image

% 2.1) Creating new FMRIB58_FA template with the wanted resolution (flirt)
% - not necessary
flag_2mm=1; %Select template resolution: 1=2mm or 0=1mm (WE CHOSE 1mm, as suggested by FSL)
if flag_2mm~=1
    varT='T';
else
    template_path='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/tbss_2mm/template';
    cmd=['flirt -in ' template_path '/FMRIB58_FA_1mm.nii.gz -ref ' template_path '/FMRIB58_FA_1mm.nii.gz -out ' template_path '/FMRIB58_FA_2mm.nii.gz -applyisoxfm 2']; %convert the 1mm template into a 2mm one (matching the DTI acquisition parameters)
    system(cmd);
    varT=['t ' template_path '/FMRIB58_FA_2mm.nii.gz']; %select the desired template by using -t and path to template
end
%2.2)Running non-linear registration (tbss_2_reg)
cd(path_tbss); %set tbss folder as current path
cmd=['submit_to_cluster -q short.q -n 1 -p MINDLAB2017_MEG-LearningBach "/usr/local/fsl/bin/tbss_2_reg -' varT '"']; %run non-linear registration
system(cmd);

%% 3) Post-registration: applying the previous registration to take all subjects into 1x1x1mm standard space (tbss_3_postreg) - relatively fast
%Transforming the original FA image into MNI152 (1x1x1mm) space and merges all of the subjects' images into a single 4D image called "all_FA" (saved in "stats" folder) and calculates the mean of all images (mean_FA), which is used to create the "mean_FA_skeleton"

%OBS! It's always better to VISUALLY INSPECT that the skeleton is well aligned with MNI152 image!! --> FSLview or FSLeyes
% Options:
% a) -S : use the mean across all subjects as a skeleton (recommended).
% b) -T : use the FMRIB58_FA mean FA image and its derived skeleton.
cd(path_tbss); %set tbss folder as current path
cmd=['tbss_3_postreg -S'];
system(cmd);
%  Check MATLAB command window for info about the status
%when visualizing the output on fsl, remember to add also the MNI152
%template and change the colourscale for the skeleton

%% 4) Project the pre-aligned data into the skeleton (tbss_4_prestats)

%Thresholding the mean FA skeleton image at the chosen threshold (0.3).
% the script takes the "all_FA" image (containing all subjects' aligned FA data) and, for each subject, projects the FA data onto the mean FA skeleton.
%The output is a binary skeleton mask (4D image) with the projected skeletonised FA data = set of voxels that will be used for the actual statistics
cd(path_tbss); %set tbss folder as current path
cmd=['tbss_4_prestats 0.2']; %0.2 is a value that generally works well
system(cmd);

%% 5) Voxelwise statistics

%% 5.1 Generating a design matrix (design.mat) and contrast files

%GENERAL SETTINGS
%loading xlsx file with all behavioral data you may be interested in for MINLABD2017-MEG_LearningBach
shiet = 2 %select the excel sheet (1 = background; 2 = WAIS-IV; 3 = MET; 4 = BDI; 5 = GOLDSMITH; 6 = MEG behavioral task)
[~,~,raw] = xlsread('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/tbss_bis/BehavioralMeasuresLearningBach.xlsx',shiet);
%getting index of participants whose FA was actually used
FA_subjs = dir('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/tbss_bis/origdata/0*');

%% WM-PR-SP and MET settings (for separate behavioral measures)

% elaborated.. and not ideally coded, but it works

% %building a design matrix for a two-sample unpaired t-test
% behav_ind = 19; %selecting the behavioral index (check the raw variable to get which index you need; e.g. WM = 24)
% desmat = zeros(length(FA_subjs),2);
% for ii = 1:size(desmat)
%     if ~isnan(raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind}) %if nan just ignore and leave 0 (that should not take that participant into account in any of the groups
%         desmat(ii,1) = raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind}; %24 = WM
%         if raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind} == 0 %getting 1s for the second column if the value for the first column is 1
%             desmat(ii,2) = 1;
%         end
%     end
% end
%Av. WM-PR-SP and Av. MET settings (Mean across behavioural measures) 
% WM+PR+SP
behav_ind = [17 22 27]; %selecting the behavioral index of the measures (check the raw variable to get which index you need; e.g. PR,WM,SP=[17 22 27])
behav_mean = zeros(length(FA_subjs),1); %Initialize vector where you will store the mean values across the behavioural measures of interest
%Compute mean across conditions, for every subject
for ii = 1:size(behav_mean)
    if ~isnan(raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind(:,1)}) %if nan just ignore and leave 0 (that should not take that participant into account in any of the groups
        behav_mean(ii,1)=(raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind(:,1)}+raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind(:,2)}+raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind(:,3)})/length(behav_ind); %Compute the average across the three cognitive measures
    end
    
end
all_behav_mean=mean(behav_mean(behav_mean~=0)); %Get mean value of cognitive measures across subjects
%Fill desmat based on the average values: 1st col:
%1=behav_mean(mm,1)>average, 0=behav_mean(mm,1)<average
desmat = zeros(length(FA_subjs),2);
for mm=1:size(desmat)
    if behav_mean(mm,1)~=0
        if behav_mean(mm,1)>all_behav_mean
            desmat(mm,1)= 1; %get 1 in desmat if the value is above the average
        elseif behav_mean(mm,1)<all_behav_mean
            desmat(mm,1) = 0; %get 0 in desmat if the value is above the average
            if behav_mean(mm,1)<all_behav_mean
                desmat(mm,2)= 1; %get 1 in desmat if the value is above the average
            elseif behav_mean(ii,1)>all_behav_mean
                desmat(mm,2) = 0; %get 0 in desmat if the value is above the average
            end
        end
    end
end

%% IN-HOUSE CLUSTER-ANALYSIS

%starting up some of the functions by LB for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl);

%% PROPER SECTION FOR T-TESTS

%%% RUNNING T-STAT ANALYSIS (IN-HOUSE SCRIPT) USING all_FA_skeletonised as input %%%
%1mm FA mask in TBSS
fname = '/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/tbss_bis/stats/all_FA_skeletonised.nii.gz';
skel = nii.load(fname); %Load all_FA_skeletonised file (4D file)
SS = size(skel(:,:,:,1));
%YOU NEED: desmat, a design matrix previously calculated (IT MUST BE COHERENT WITH THE STATISTICS THAT YOU WANT TO RUN)
ind = find(skel(:,:,:,1) ~= 0); %getting indices of non-0 voxels (indeces that only belong to brain areas, not empty space)
[i1,i2,i3] = ind2sub(size(skel(:,:,:,1)),ind); %reshaping the vector into a 3D matrix shaped as the original data (skel) = get one vector for each dimension

%% Loop for t-tests

P = zeros(length(i1),1);
TVAL = zeros(length(i1),1);
for nn = 1:length(ind) %Loop across voxels
    a = squeeze(skel(i1(nn),i2(nn),i3(nn),(desmat(:,1)==1))); %take 1st group of subjects
    b = squeeze(skel(i1(nn),i2(nn),i3(nn),(desmat(:,2)==1))); %take 2nd group of subjects
    %think if using ranksum or ttest (doesn't really change much..)
    %[p,~,stats] = ranksum(a,b);
    [~,p,~,stats] = ttest2(a,b); %compute t-test between the two groups
    P(nn,1) = p; %store p-values in P
    TVAL(nn,1) = stats.tstat; %store t-values in TVAL
    disp([num2str(nn) '/' num2str(length(ind))])
end
%reshaping them into 3D matrices
TVAL2 = zeros(SS(1),SS(2),SS(3));
P2 = zeros(SS(1),SS(2),SS(3));
for ii =1:length(i1)
    TVAL2(i1(ii),i2(ii),i3(ii)) = TVAL(ii);
    P2(i1(ii),i2(ii),i3(ii)) = P(ii);
end
% %saving image (nifti)
% fname = '/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/tbss_2mm/stats/mean_FA_skeleton_mask.nii.gz';
% aalnii = load_nii(fname);
% aalnii.img = TVAL2;
% save_nii(aalnii,'jes_inhousttest_tval.nii.gz');

%% PERMUTATION FUNCTION FOR THE CLUSTER OF COMPUTERS (PARALLEL COMPUTING - AARHUS UNVIERSITY)

cogmeas='PR_WM_SP'; %name of the cognitive function considered - TO BE CHANGED
jobs_vector = [0 0 0 1]; %specify which jobs you want to run: 1) and 2) together first, then 3) and then 4).

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/DTI_Portis');
%function inputs
S = [];
S.V_p = []; %data in 3D (p-val)
S.V_t = []; %data in 3D (t-val)
S.p_thresh = 0.05; %defining threshold for p-values
dum = num2str(S.p_thresh); %to get name with p_thresh inside..
S.thresh_mc = 0.001; %%defining threshold for Monte-Carlo simulations
S.perm_numb = 1000; %number of permutations (for now less than 100 or multiples of 100. If e.g. 954, the function will anyway compute 1000 permutations. If 923, it will compute 900 (so it rounds it..). This is not ideal nor elegant, but for now that is good enough..)
S.sizemass = 'size'; %cluster 'size' (=n. voxels forming the cluster) or 'mass' (=sum of the t-value of the voxels forming the cluster)
S.outdir = ['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/tbss_bis/stats/' cogmeas '_' dum(3:end)]; %output directory where to store results (should be a different one for each analysis you run)
%1mm FA mask
S.maskFA = '/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/tbss_bis/stats/mean_FA_skeleton_mask.nii.gz'; %mask FA created by FSL TBSS
%2mm FA mask
% S.maskFA = '/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/tbss_2mm/stats/mean_FA_skeleton_mask.nii.gz'; %mask FA created by FSL TBSS
S.analysis_name = ['PVAL_' dum(3:end)]; %name of the output file with the results stored in
S.clsopt.queue = 1; S.clsopt.slot = 1; S.clsopt.scheduler = 'cluster'; %options for cluster (parallel computing)
S.fun_run = jobs_vector; %specifying which steps you want to run 1)original cluster, 2)permutations, 3)putting everything together, 4)plotting
%actual function        
DTI_cluster_perm_2groups_LBPD(S);

%% %%% %%% %%% %%

%% PROBABILISTIC TRACTOGRAPHY - FSL
%BedpostX (using markov chain Monte Carlo) (slow, about 15 hours)

path_out='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/'; %set the path to the folder where you will save your converted data
path1=dir([path_out '/00*']); %query the content of the folder "raw": it contains one folder for each subject
% A) Create and input directory for BedpostX (DTI001_BedpostX in this case)
% for each subject
for ii=1:length(path1) %loop across all subjects
    mkdir([path1(1).folder '/' path1(ii).name], [path1(ii).name '_BedpostX']);  %create a new folder for each subject (ii) and name it 'subjectID_Bedpostx'
end

%%  B) Add the following files to the folder
% 1) "data.nii.gz" which is the renamed version of "eddy_unwarped_images.nii.gz"
% 2) "nodif_brain_mask.nii.gz" mask (you already have it)
% 3) "bvals" obtained from initial nifti convertion (here we renamed "20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bval" as "bvals"
% 4) "bvecs" obtained from initial nifti convertion (here we renamed "20180123_102505cmrrmbep2ddiffmultidirAPSubjectNo0001s014a001.bvec" as "bvecs"

path_out='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/'; %set the path to the folder where you will save your converted data
path1=dir([path_out '/00*']); %query the content of the folder "raw": it contains one folder for each subject
for gg=2:length(path1)
    path_bed=dir([path1(1).folder '/' path1(gg).name '/' path1(gg).name '_BedpostX' ]); %path to BedpostX folder
    path_datanii=[path1(1).folder '/' path1(gg).name '/eddy_unwarped_images.nii.gz']; %path(gg) to eddy unwarped images
    path_nodif=[path1(1).folder '/' path1(gg).name '/hifi_nodif_brain_mask.nii.gz']; %path to nodif_brain_mask
    path_sub=dir([path1(1).folder '/' path1(gg).name '/2*']); %path to each subject's main folder and relative content
    path_bvals=[path_sub(1).folder '/' path_sub(1).name];%path to each subject's .bval file
    path_bvecs= [path_sub(1).folder '/' path_sub(2).name];%path to each subject's .bvec file
    copyfile(path_datanii,[path_bed(1).folder '/data.nii.gz']); %copy eddy_unwarped_images and rename it as "data.nii.gz"
    delete(path_datanii); %remove eddy_unwarped_images from subject's main folder
    %remove oiginal eddy
    copyfile(path_nodif,[path_bed(1).folder '/nodif_brain_mask.nii.gz']); %copy hifi_nodif_brain_mask and rename it as nodif_brain_mask
    copyfile(path_bvals,[path_bed(1).folder '/bvals']); %copy .bval file and rename it as "bvals"
    copyfile(path_bvecs,[path_bed(1).folder '/bvecs']); %copy .bvec file and rename it as "bvecs"
    disp(gg)
end

%%  C) calling terminal command with "bedpostx" operating in the desired directory (bedpostx) - slow

path_out='/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/';
path1=dir([path_out '00*']); %query the content of the folder "raw": it contains one folder for each subject
for gg = 1:length(path1)
    path_bed_run = [path1(1).folder '/' path1(gg).name '/' path1(gg).name '_BedpostX' ]; %path to BedpostX folder
    %cd(path_bed); %change directory for each subject
    cmd = ['submit_to_cluster -q long.q -n 1 -p MINDLAB2017_MEG-LearningBach "/usr/local/fsl/bin/bedpostx ' path_bed_run '"']; %submit BedpostX to cluster
    system(cmd)
    disp(gg)
end
% cmd = 'bedpostx /projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/DTI001/DTI001_BedpostX';
%bedpostx command
%directory that you created (/projects/MINDLAB2017_MEG-LearningBach/scratch/DTI_Portis/DTI001/DTI001_BedpostX)
% system(cmd)

%% CREATING INDEPENDENT FILES FOR EACH AAL REGION

%%% THIS SHOULD BE RUN ONLY ONCE AND THEN THE MASKS CAN BE USED FOR DIFFERENT PROJECTS %%%

% creating indepdent image files with AAL ROIs
% command to split one ROI (80) from AAL
% cmd = 'fslmaths aal_1mm.nii.gz -thr 80 -uthr 80 -bin 80.nii.gz';
% system(cmd)
for ii = 1:90
    cmd = ['fslmaths /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/DTI001/Tractography/prova_nazi/AAL_templates/aal_2mm_try.nii.gz -thr ' num2str(ii) ' -uthr ' num2str(ii) ' -bin /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/DTI001/Tractography/prova_nazi/AAL_templates/AAL_2mm_90ROIs/' num2str(ii) '.nii.gz'];
    system(cmd)
end
%then, we manually moved such files into the following directory:
% /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/AAL_2mm_90ROIs

%% REGISTRATION AND PROBTRACKX

%registration after bedpostx
path = dir('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/00*'); %path to subjects
for ii = 3:length(path) %over subjects
    %going to subject's directory
    cd([path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX'])
    copyfile('nodif_brain_mask.nii.gz','nodif_brain.nii.gz') %copying mask file (just to rename it..)
    %actual command line for registration with flirt (FSL)
    cmd = 'flirt -in nodif_brain.nii.gz -ref /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_T1_2mm_brain.nii.gz -omat diff2stand.mat';
    system(cmd)
    %moving the output file into the bedpostX xfms subfolder
    copyfile('diff2stand.mat',[path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/diff2stand.mat'])
    load([path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/diff2stand.mat'],'-ascii') %loading the matrix computed by flirt
    %then, since we need to get the inverse matrix, it seems that we can simply compute it, as follows:
    stand2diff = inv(diff2stand);
    save([path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/stand2diff.mat'],'stand2diff','-ascii') %saving the inverse matrix
    disp(ii)
end

%% ACTUAL TRACTOGRAPHY (probtrackx - FSL)

mask_label=1; %flag for mask type: 1=3559 parcels; 2=2mm AAL
streamlines=1; %set streamlines number (1000 streamlines for AAL 2mm)

%1000 streamlines since it seems to give (almost) the same results as 5000 streamlines (the most usual option..)
path = dir('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/00*'); %path to subjects
for ii = 2 %59:length(path) %over subjects
    if mask_label ==1
        mask = ['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/masks.txt']; %3559 parcels
        outdir = [path(ii).folder '/' path(ii).name '/Tractography/parce3559_' num2str(streamlines) 'stream']; %output directory
    elseif mask_label==2
        mask = ['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/masks.txt']; %2-mm AAL
        outdir = [path(ii).folder '/' path(ii).name '/Tractography/AAL90_' num2str(streamlines) 'stream']; %output directory
    end
    xfm_mat = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/stand2diff.mat']; %registration matrix
    invxfm_mat = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/xfms/diff2stand.mat']; %registration matrix (inverse)
    mergedk = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/merged']; %merged files from bedpostX
    nodifmask = [path(ii).folder '/' path(ii).name '/' path(ii).name '_BedpostX.bedpostX/nodif_brain_mask']; %brain mask from bedpostX
    %actual line for probtrackx - FSL
    cmd = ['submit_to_cluster -q long.q -n 2 -p MINDLAB2017_MEG-LearningBach "/usr/local/fsl/bin/probtrackx2 --network -x ' mask ' -l --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P ' num2str(streamlines) ' --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --xfm=' xfm_mat ' --invxfm=' invxfm_mat ' --forcedir --opd -s ' mergedk ' -m ' nodifmask ' --dir=' outdir '"'];
    system(cmd)
end

%% WORKING ON AAL CONNECTIVITY MATRIX (1000 STREAMLINES)

%starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl); % add to path LBPD and other OSL functions

%% Defining different groups of subjects for each cognitive measure

%GENERAL SETTINGS
%loading xlsx file with all behavioral data you may be interested in for MINLABD2017-MEG_LearningBach
shiet = 2; %select the proper excel sheet (1 = background; 2 = WAIS-IV; 3 = MET; 4 = BDI; 5 = GOLDSMITH; 6 = MEG behavioral task)
[~,~,raw] = xlsread('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/tbss_bis/BehavioralMeasuresLearningBach.xlsx',shiet);
%getting index of participants whose AAL connectivity matrix was actually calculated
FA_subjs = dir('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/tbss_bis/origdata/0*');

%% WM-PR-SP settings

% once again, not the most elegant code, but it works..
% %building a design matrix for a two-sample unpaired t-test
% behav_ind = 19; %selecting the behavioral index (check the raw variable to get which index you need; e.g. WM = 24)
% desmat = zeros(length(FA_subjs),2);
% for ii = 1:size(desmat)
%     if ~isnan(raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind}) %if nan just ignore and leave 0 (that should not take that participant into account in any of the groups
%         desmat(ii,1) = raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind}; %24 = WM
%         if raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind} == 0 %getting 1s for the second column if the value for the first column is 0
%             desmat(ii,2) = 1;
%         end
%     end
% end

% Defining different groups of subjects for average across cognitive measures
behav_ind = [17 22 27]; %select the behavioral index of the measures (check the raw variable to get which index you need; e.g. PR,WM,SP=[17 22 27])
behav_mean = zeros(length(FA_subjs),1); %Initialize vector where you will store the mean values across the behavioural measures of interest
%Compute mean across conditions, for every subject
for ii = 1:size(behav_mean)
    if ~isnan(raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind(:,1)}) %if nan just ignore and leave 0 (that should not take that participant into account in any of the groups
        behav_mean(ii,1)=(raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind(:,1)}+raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind(:,2)}+raw{str2double(FA_subjs(ii).name(1:4)) + 1,behav_ind(:,3)})/length(behav_ind); %Compute the average across the three cognitive measures
    end 
end
all_behav_mean=mean(behav_mean(behav_mean~=0)); %Get mean value of cognitive measures across subjects
%Fill desmat based on the average values: 1st col:
%1=behav_mean(mm,1)>average, 0=behav_mean(mm,1)<average
desmat = zeros(length(FA_subjs),2);
for mm=1:size(desmat)
    if behav_mean(mm,1)~=0
        if behav_mean(mm,1)>all_behav_mean
            desmat(mm,1)= 1; %get 1 in desmat if the value is above the average
        elseif behav_mean(mm,1)<all_behav_mean
            desmat(mm,1) = 0; %get 0 in desmat if the value is above the average
            if behav_mean(mm,1)<all_behav_mean
                desmat(mm,2)= 1; %get 1 in desmat if the value is above the average
            elseif behav_mean(ii,1)>all_behav_mean
                desmat(mm,2) = 0; %get 0 in desmat if the value is above the average
            end
        end
    end
end

%% COMPUTING STRUCTURAL CONNECTIVITY MATRIX FOR EACH SUBJECT (the plotting option here is not very much meaningful anymore)

% 1) compute matrices for each subject
% 2) choose if plotting them separately or averaging over subject
single_subj_plot = 0; % 1 for single subject plot; 0 for 2 groups of subjects (e.g. high and low WM)
compute_mat = 0; % 1 for computing the matrices; 0 for loading the previously computed matrices

path = dir('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/00*'); %path to subjects
order = [1:2:90 90:-2:2]; %to reshape the matrix symmetricaly
if single_subj_plot == 1
    for ii = 1:length(path) %over subjects
        if compute_mat == 1
            %reading output matrix from "probtrackx2" as a table
            M = readtable([path(ii).folder '/' path(ii).name '/Tractography/AAL90_1000stream/fdt_network_matrix']);
            M = table2array(M); %converting table to double
            M2 = M(1:90,1:90); %extracting only first 90 rows and columns (it seems that we get by default one additional NaN column)
            %corrections
            % 1) averaging of ROI1-ROI2 and ROI2-ROI1
            % 2) dividing by sum of sizes of ROIs
            cnt = 1; %intialize counter
            M4 = zeros(90); %preallocate space for matrix
            for pp = 1:90 %over ROI1
                %loading image ROI1
                ROI1 = load_nii(['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/AAL_2mm_90ROIs/' num2str(pp) '.nii.gz']);
                sizeROI1 = length(find(ROI1.img==1)); %calculating size of ROI1 (number of voxels)
                cnt = cnt + 1;
                for oo = cnt:90 %over ROI2
                    dummyelv = (M2(pp,oo) + M2(oo,pp))/2; %average of ROI1-ROI2 and ROI2-ROI1 (outputted from tractography)
                    %loading image ROI2
                    ROI2 = load_nii(['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/AAL_2mm_90ROIs/' num2str(oo) '.nii.gz']);
                    sizeROI2 = length(find(ROI2.img==1)); %calculating size of ROI2 (number of voxels)
                    M4(pp,oo) = dummyelv/(sizeROI1 + sizeROI2); %dividing previous average by sum of sizes of the two ROIs that are connected (storing it in the upper triangle)
                    M4(oo,pp) = M4(pp,oo); %storing it also in the lower triangle
                end
                disp(['Subject = ' path(ii).name ' - ROI1 = ' num2str(pp)])
            end
            %since we provided the AAL ROIs in LRLRLR order, we reshape them to get the symmetric order (LLLRRR) for nicer and more usual visualization purposes
            M3 = M4(order,order);
            %saving matrix
            save(['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Graph_Theory_AAL90Matrix/Subj_' path(ii).name '.mat'],'M3')
        else
            load(['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Graph_Theory_AAL90Matrix/Subj_' path(ii).name '.mat']); %loading matrix
        end
        %plotting connectivity matrix
        figure
        imagesc(M3)
        colorbar
%         caxis([0 12e+5])
        set(gcf,'color','w')
        %changing colormap of figures (red-bue with white for 0 values)
        x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
        colormap(bluewhitered_PD(0,x))
    end
else
    g_H = zeros(90,90,length(find(desmat(:,1)==1)));
    g_L = zeros(90,90,length(find(desmat(:,2)==1)));
    cntH = 0;
    cntL = 0;
    for ii = 1:length(path)
        if sum(desmat(ii,:)) ~= 0 %if subject belongs to one of the two groups
            load(['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Graph_Theory_AAL90Matrix/Subj_' path(ii).name '.mat']); %loading matrix
            if desmat(ii,1) == 1
                cntH = cntH + 1;
                g_H(:,:,cntH) = M3;
            else
                cntL = cntL + 1;
                g_L(:,:,cntL) = M3;
            end
        end
        disp(ii)
    end
end
%actual plotting
MH = mean(g_H,3);
ML = mean(g_L,3);
%matrix 1
figure
imagesc(MH)
colorbar
caxis([0 12e+5])
set(gcf,'color','w')
%changing colormap of figures (red-bue with white for 0 values)
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))
%matrix 2
figure
imagesc(ML)
colorbar
caxis([0 12e+5])
set(gcf,'color','w')
%changing colormap of figures (red-bue with white for 0 values)
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))
%difference matrix
DM = MH - ML;
figure
imagesc(DM)
colorbar
caxis([-10e+4 10e+4])
set(gcf,'color','w')
%changing colormap of figures (red-bue with white for 0 values)
x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
colormap(bluewhitered_PD(0,x))
load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat');
labsymm = lab(order,:);

%% %%% %%% %%% %%

%% GRAPH THEORY ON STRUCTURAL CONNECTIVITY MATRICES

clear
label_cogmeas = 4; %1=PR, 2=WM, 3=SP, 4=PR+WM+SP, 5=BDI
BDI_nodepr = 4; %value to be used for thresholding healthy subjects (i.e. 2 means subjects with BDI no more than 2 are healthy..)
label_test = 2; %1=Degree, 2=modularity intra and extra connectivity; 0 = no
brain_plot = 0; %1 = brain plot for label_test (degree and connector/provincial hubs); 0 = no
violins = 1; %1 = violins plot for label_test (degree and connector/provincial hubs); 0 = no
% label_ttest = 0; %1=t-test, 0=no
label_additiona_graph_measures = 0; %additional measures: 1 = char path length; 2 = glob eff; 3 = nod eff; 4 = loc eff; 5 = modularity
label_density = 0; %0 = no; otherwise percentage of values we want to remove for thresholding (for density) (e.g. 10 = 10%)
label_modularity_brainplot = 0; %1 to plot modularity in the brain
intra_con = 1; %1 for intra subnetwork connections, 0 inter connectivity
label_schem_mod = 0; %1 for schemaball plotting of modularity
mod_perm_label = 0; %1 for MCS (1000 permutations) on modularity (establishing whether the brain is more "modulable" than a random graph with its elements; 0 for not
% inter_con = 0; %1 for inter subnetwork connections, 0 no
% perc_intra = 100; %percentage of intra subnetwork connections to be plotted
limitt_mod = [0.1 5]; %scaling values for modularity plotting
label_matrix_plot = 0; %1 for plotting the connectivity matrix (as a matrix and as a schemaball)

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/DTI_Portis'); %adding path to function
close all
COL2 = [0.8 0 0; 0 0 0.8; 1 0 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0]; %vector with colors
%%% 1) DEFINING DESMAT (DIFFERENT GROUPS OF SUBJECTS)
%Defining different groups of subjects for each cognitive measure
%loading xlsx file with all behavioral data you may be interested in for MINLABD2017-MEG_LearningBach
if label_cogmeas ~= 5
    shiet = 2; %select the proper excel sheet (1 = background; 2 = WAIS-IV; 3 = MET; 4 = BDI; 5 = GOLDSMITH; 6 = MEG behavioral task)
else
    shiet = 4; %select the proper excel sheet (1 = background; 2 = WAIS-IV; 3 = MET; 4 = BDI; 5 = GOLDSMITH; 6 = MEG behavioral task)
end
[~,~,raw] = xlsread('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/tbss_bis/BehavioralMeasuresLearningBach.xlsx',shiet);
%getting index of participants whose AAL connectivity matrix was actually calculated
FA_subjs = dir('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Graph_Theory_AAL90Matrix/Subj*'); %list of subjects
%WM-PR-SP settings
%building a design matrix for a two-sample unpaired t-test
if label_cogmeas == 1 %PR
    behav_ind = 19; %selecting the behavioral index (check the raw variable to get which index you need; e.g. WM = 24)
elseif label_cogmeas == 2 %WM
    behav_ind = 24; %selecting the behavioral index (check the raw variable to get which index you need; e.g. WM = 24)
elseif label_cogmeas == 3 %SP
    behav_ind = 29; %selecting the behavioral index (check the raw variable to get which index you need; e.g. WM = 24)
elseif label_cogmeas == 5 %BDI
    behav_ind = 2; %selecting the behavioral index (check the raw variable to get which index you need; e.g. WM = 24)
end
desmat = zeros(length(FA_subjs),2);
if label_cogmeas < 4 %IQ independent measures
    for ii = 1:size(desmat)
        if ~isnan(raw{str2double(FA_subjs(ii).name(6:9)) + 1,behav_ind}) %if nan just ignore and leave 0 (that should not take that participant into account in any of the groups)
            desmat(ii,1) = raw{str2double(FA_subjs(ii).name(6:9)) + 1,behav_ind}; %24 = WM
            if raw{str2double(FA_subjs(ii).name(6:9)) + 1,behav_ind} == 0 %getting 1s for the second column if the value for the first column is 0
                desmat(ii,2) = 1;
            end
        end
    end
end
if label_cogmeas == 4 %IQ all measures together
    %Defining different groups of subjects based on the average values of cognitive measures
    behav_ind = [17 22 27]; %select the behavioral index of the measures (check the raw variable to get which index you need; e.g. PR,WM,SP=[17 22 27])
    behav_mean = zeros(length(FA_subjs),1); %Initialize vector where you will store the mean values across the behavioural measures of interest
    %Compute mean across conditions, for every subject
    for ii = 1:size(behav_mean)
        if ~isnan(raw{str2double(FA_subjs(ii).name(6:9)) + 1,behav_ind(:,1)}) %if nan just ignore and leave 0 (that should not take that participant into account in any of the groups
            behav_mean(ii,1)=(raw{str2double(FA_subjs(ii).name(6:9)) + 1,behav_ind(:,1)}+raw{str2double(FA_subjs(ii).name(6:9)) + 1,behav_ind(:,2)}+raw{str2double(FA_subjs(ii).name(6:9)) + 1,behav_ind(:,3)})/length(behav_ind); %Compute the average across the three cognitive measures
        end
    end
    all_behav_mean=mean(behav_mean(behav_mean~=0)); %Get mean value of cognitive measures across subjects
    %Fill desmat based on the average values: 1st col:
    %1=behav_mean(mm,1)>average, 0=behav_mean(mm,1)<average
    %     desmat = zeros(length(FA_subjs),2);
    for mm=1:size(desmat)
        if behav_mean(mm,1)~=0
            if behav_mean(mm,1)>all_behav_mean
                desmat(mm,1)= 1; %get 1 in desmat if the value is above the average
            elseif behav_mean(mm,1)<all_behav_mean
                desmat(mm,1) = 0; %get 0 in desmat if the value is above the average
                if behav_mean(mm,1)<all_behav_mean
                    desmat(mm,2)= 1; %get 1 in desmat if the value is above the average
                elseif behav_mean(ii,1)>all_behav_mean
                    desmat(mm,2) = 0; %get 0 in desmat if the value is above the average
                end
            end
        end
    end
end
if label_cogmeas == 5 %BDI
    %Defining different groups according to BDI (col1=LOW, col2=HIGH)
    for ii = 1:size(desmat)
        if ~isnan(raw{str2double(FA_subjs(ii).name(6:9)) + 1,behav_ind})
            if (raw{str2double(FA_subjs(ii).name(6:9)) + 1,behav_ind}<=BDI_nodepr)
                desmat(ii,1) = 1; %1 in the 1st column for BDI<=2 (NO tendency to depression)
                desmat(ii,2) = 0;
            elseif (raw{str2double(FA_subjs(ii).name(6:9)) + 1,behav_ind}>=7)   %getting 1s for the second column if the value for the first column is higher than 7
                desmat(ii,2) = 1; % 1 in the 2nd column for BDI>7 (tendency to depression)
                desmat(ii,1) = 0;
            end
        end
    end
end
%%% 2) GRAPH THEORY
%%%%%% GRAPH THEORY APPLIED TO DTI %%%%%%
%Graph Theory - Degree and degree-based hubs, modularity  
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/DTI_Portis/GraphTheory/BCT/2019_03_03_BCT'); %add BCT functions to path (FFM)
louv_l = 0; %1 for Louvain; 0 for modularity
list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Graph_Theory_AAL90Matrix/Subj*'); %list of subjects
%initialising matrices for several measures of graph theory
degree = zeros(length(list),90);
hubs=zeros(length(list),90);
provincial=zeros(length(list),90);
connector=zeros(length(list),90);
optstruct= zeros(90,length(list));
maxmod=zeros(length(list),1);
permnum= 1000;%num permutations
nodesum=zeros(length(list),90);
MM = zeros(90,90,length(list));
%actual computation
for ii = 1:length(list) %over subjects
    %loading matrix with strenght of connections (AAL, same ROIs used for distance, of course..)
    load([list(ii).folder '/' list(ii).name]); %Load matrix (M3) containing the strength of the connections across regions of interest
    MM(:,:,ii) = M3;
    degree(ii,:) = sum(M3); %compute degree of connectivity (weighted)
    %according to previous DTI-graph theory papers
    meandegreen=mean(degree(ii,:))/mean(degree(ii,:)); %mean degree, normalized = 1 (done just for the sake of formality)
    degreestd=std(degree(ii,:)./mean(degree(ii,:))); %std of degree, normalized
    degreen=degree(ii,:)./mean(degree(ii,:)); %normalized degree
    hubs(ii,:)=double(degreen>(meandegreen+degreestd));%compute degree-based hubs
    if louv_l == 1
        optstructP= zeros(90,permnum);
        maxmodP=zeros(permnum,1);
        for pp=1:permnum
            %                 [optstructP(:,pp), maxmodP(pp,1)] = modularity_und(M3);
            [optstructP(:,pp), maxmodP(pp,1)] = community_louvain(M3); %Compute modularity with Louvain algorithm
            %         disp(pp)
        end
        [maxval, maxind]=max(maxmodP); %getting maximum value showing the best modularity over the 1000 permutations
        maxmod(ii,1) =maxval; %storing such value for each subject
        optstruct(:,ii)=optstructP(:,maxind); %storing modularity corresponding to that value
    else
        optstruct(:,ii) = modularity_und(M3); %Newmann modularity (deterministic)
    end
    %strength of subnetworks as outputted by Louvain modularity
    for hh=1:90 %loop over the nodes (initial AAL areas)
        %       nodesum(:,hh)=degree(ii,hh)+degree(:,hh); %sum every connection with all the other connections, for each area (areas=columns)
        subnethh=optstruct(hh,ii); %Subnetwork to which the node hh belongs
        subidx=find(optstruct(:,ii)==subnethh); %Index of the nodes within the subnetwork
        nodesum(ii,hh)=(sum(M3(hh,subidx)))./sum(M3(hh,:)); %sum of the connections between node (hh) and all the
        %other nodes of its community, divided by all connections hh has with all the other nodes
    end
    %provincial and connector hubs based on modularity
    meanode=mean(nodesum(ii,:)); %Compute mean of the hubs (provincial and connectors), for each subject
    stdnode=std(nodesum(ii,:)); %Compute std of the hubs (provincial and connectors)
    provincial(ii,:)=double(nodesum(ii,:)>(meanode+stdnode));%compute provincial hubs (=values above mean value + std)
    connector(ii,:)=double(nodesum(ii,:)<(meanode-stdnode));%compute connector hubs (=values below mean value - std)
    %showing progressive subject ID
    disp(ii)
end
% save(['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Graph_Theory_AAL90Matrix/' cog_meas '_GT.mat'], 'degree', 'nodesum', 'hubs', 'provincial', 'connector') %saving the GT measures for each cognitive measure
if label_test == 1 %computing MCS for degree and connector/provincial hubs
    clear coldum
    coldum(1) = 'r'; coldum(2) = 'b'; coldum(3) = 'r';
    % PERMUTATION TEST (DEGREE)
    addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/DTI_Portis'); %adding path to function
    S = [];
    S.thr = [];
    S.data{1,1} = degree(desmat(:,1)==1,:); %degrees of subjects with high PR
    S.data{1,2} = degree(desmat(:,2)==1,:); %degrees of subjects with low PR
    S.thr = std(median(S.data{1,1},1) - median(S.data{1,2},1)); %keeping only extreme data (leave empty [] for using all data)
%     S.thr = 20; %keeping only extreme data (leave empty [] for using all data)
    S.permtype = 1; %choose the permutation type: 1 for permuting the degrees of each couple of regions (high vs low PR) independently; 2 for permuting all degrees together
    S.permnum = 10000;
    disp('degree results')
    [p_val_pos, p_val_neg, posidx, negidx, diff1] = DTI_GT_MCS(S);
    disp(['p_val_pos: ' num2str(p_val_pos) ' - p_val_neg: ' num2str(p_val_neg)])
    %loading AAL labels
    if ~exist('lab','var')
        load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat');
        order = [1:2:90 90:-2:2]; %to reshape the matrix symmetricaly
        lab2 = lab(order,:);
    end
    disp('pos nodes')
    ROIspos = lab2(posidx',:)
    valuespos = diff1(posidx)';
    disp('neg nodes')
    ROIsneg = lab2(negidx',:)
    valuesneg = diff1(negidx)';
    save('DTI_degree.mat','ROIspos','ROIsneg','p_val_pos','p_val_neg','valuespos','valuesneg')
    if violins == 1
        %violins
        po1 = (median(S.data{1,1},1))'; %median over subjects (group1)
        po2 = (median(S.data{1,2},1))'; %median over subjects (group2)
        datt{1,1} = po1; datt{1,2} = po2;
        cl = COL2(1:2, :); %same..
        figure
        h = rm_raincloud2(datt',cl)
        grid minor
        set(gcf,'color','w')
        set(gcf,'Position',[200,200,400,550])
        xlim([0 1100])
        saveas(gcf,['DTI_Degree.svg'])
        %difference (violin plot)
        poo = (po1 - po2); %+ (mean(po1)+mean(po2))/2;
        datt2{1,1} = poo;
        figure
        h = rm_raincloud2(datt2',COL2(3, :))
        xlim([-100 100])
        grid minor
        set(gcf,'color','w')
        set(gcf,'Position',[200,200,400,550])
        saveas(gcf,['DTI_Degree_diff.svg'])
    end
    %PLOTTING DEGREE IN THE BRAIN (IN THE CENTROID LOCATIONS)
    if brain_plot == 1
%         posidx = negidx; %to get BDI depressed people.. if BDI label is set 
        lo = (squeeze(median(S.data{1,1},1)) - squeeze(median(S.data{1,2},1)))'; %storing difference between median of the degree of the two groups
        S.data{1,3} = lo(posidx');
        S.data{1,1} = (median(S.data{1,1},1))'; %median over subjects (group1)
        S.data{1,2} = (median(S.data{1,2},1))'; %median over subjects (group2)
        rmin12 = min(cat(1,S.data{1},S.data{2})); %for later scaling..
        rmax12 = max(cat(1,S.data{1},S.data{2}));
        posidx = find(posidx==1)';
        for cc = 1:3 %group1 - group2 - difference of their medians
            limitt = [2,30]; %size of centroids for the 2 groups
            openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
            load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI_RC_AAL_2mm.mat')
            %scaling on the basis of connections strength
            dataa = S.data{cc};
            if cc == 3 %scaling difference  between groups
                rmin = min(dataa);
                rmax = max(dataa);
                %creating nifti image for WorkBench
                vector = zeros(90,1);
                vector(posidx) = dataa;
                create_AALnifti(vector,['DTI_degree_contr_' num2str(cc) '.nii'],1);
            else %otherwise scaling together group1 and group2
                rmin = rmin12; rmax = rmax12;
                %creating nifti image for WorkBench
                vector = dataa;
                create_AALnifti(vector,['DTI_degree_contr_' num2str(cc) '.nii'],1);
            end
            %scaling connection strengths with extreme values (limitt) requested by
            size_con = (dataa-rmin)./(rmax-rmin).*(limitt(2)-limitt(1)) + limitt(1);
            for gg = 1:length(dataa) %over AAL ROIs
                hold on
                if cc == 3
                    plot3(MNI_RC(posidx(gg),1), MNI_RC(posidx(gg),2), MNI_RC(posidx(gg),3), '.', 'Color', coldum(cc), 'MarkerSize', size_con(gg)); %centroid of ROI ii
                else
                    plot3(MNI_RC(gg,1), MNI_RC(gg,2), MNI_RC(gg,3), '.', 'Color', coldum(cc), 'MarkerSize', size_con(gg)); %centroid of ROI ii
                end
            end
            rotate3d on; axis off; axis vis3d; axis equal
            set(gcf,'color','w')
            title(['degree - contrast ' num2str(cc)])
            view([0 90])
%             export_fig(['/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021/Dots_top_DTI_degree_contr_' num2str(cc) '.png'])
%             saveas(gcf,['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Images_Paper_IQ/Dots_top_DTI_degree_contr_' num2str(cc) '.svg'])
            view([-90 0])
%             export_fig(['/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021/Dots_left_DTI_degree_contr_' num2str(cc) '.png'])
%             saveas(gcf,['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Images_Paper_IQ/Dots_left_DTI_degree_contr_' num2str(cc) '.svg'])
        end
    end
else
    %PERMUTATION TEST (NODESUM - INTRASUBNETWORK CONNECTIVITY DIVIDED BY ALL CONNECTIVITY FOR EACH NODE AT A TIME)
    coldum(1) = 'r'; coldum(2) = 'b'; coldum(3) = 'r';
    S = [];
    S.data{1,1} = nodesum(desmat(:,1)==1,:); %degrees of subjects with high PR
    S.data{1,2} = nodesum(desmat(:,2)==1,:); %degrees of subjects with low PR
    S.thr = []; %keeping only extreme data (leave empty [] for using all data)
    S.thr = std(median(S.data{1,1},1) - median(S.data{1,2},1)); %keeping only extreme data (leave empty [] for using all data)
    % S.thr = std(mean(S.data{1,1},1) - mean(S.data{1,2},1)); %keeping only extreme data (leave empty [] for using all data)
    S.permtype = 1; %choose the permutation type: 1 for permuting the degrees of each couple of regions (high vs low PR) independently; 2 for permuting all degrees together
    S.permnum = 10000;
    disp('connector/provincial results')
    [p_val_pos, p_val_neg, posidx2, negidx2, diff2] = DTI_GT_MCS(S);
    disp(['p_val_pos: ' num2str(p_val_pos) ' - p_val_neg: ' num2str(p_val_neg)])
    %loading AAL labels
    if ~exist('lab','var')
        load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/AAL_labels.mat');
        order = [1:2:90 90:-2:2]; %to reshape the matrix symmetricaly
        lab2 = lab(order,:);
    end
    disp('pos nodes')
    ROIspos = lab2(posidx2',:)
    valuespos = diff2(posidx2)';
    disp('neg nodes')
    ROIsneg = lab2(negidx2',:)
    valuesneg = diff2(negidx2)';
    save('DTI_connector.mat','ROIspos','ROIsneg','p_val_pos','p_val_neg','valuespos','valuesneg')    
    if violins == 1
        %violins
        po1 = (median(S.data{1,1},1))'; %median over subjects (group1)
        po2 = (median(S.data{1,2},1))'; %median over subjects (group2)
        datt{1,1} = po1; datt{1,2} = po2;
        cl = COL2(1:2, :); %same..
        figure
        h = rm_raincloud2(datt',cl)
        xlim([0.2 1.2])
        grid minor
        set(gcf,'color','w')
        set(gcf,'Position',[200,200,400,550])
        %         saveas(gcf,['DTI_Connector.svg'])
        %difference (violin plot)
        poo = (po1 - po2); %+ (mean(po1)+mean(po2))/2;
        datt2{1,1} = poo;
        figure
        h = rm_raincloud2(datt2',COL2(3, :))
        xlim([-0.25 0.25])
        grid minor
        set(gcf,'color','w')
        set(gcf,'Position',[200,200,400,550])
        %         saveas(gcf,['DTI_Connector_diff.svg'])
    end    
    %PLOTTING CONNECTOR HUBS IN THE BRAIN (IN THE CENTROID LOCATIONS)
    if brain_plot == 1
        if p_val_pos < 0.05 || p_val_neg < 0.05
            if length(find(negidx2==1)) < length(find(posidx2==1)) %taking the side  of the distribution where the difference was larger
                negidx2 = posidx2;
            end
            lo = (squeeze(median(S.data{1,1},1)) - squeeze(median(S.data{1,2},1)))'; %storing difference between median of the degree of the two groups
            S.data{1,3} = lo(negidx2');
            S.data{1,1} = (median(S.data{1,1},1))'; %median over subjects (group1)
            S.data{1,2} = (median(S.data{1,2},1))'; %median over subjects (group2)
            rmin12 = min(cat(1,S.data{1},S.data{2})); %for later scaling..
            rmax12 = max(cat(1,S.data{1},S.data{2}));
            negidx2 = find(negidx2==1)';
            for cc = 1:3 %group1 - group2 - difference of their medians
                limitt = [2,30]; %size of centroids for the 2 groups
                openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
                load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI_RC_AAL_2mm.mat')
                %scaling on the basis of connections strength
                dataa = abs(S.data{cc});
                if cc == 3 %scaling difference between groups
                    %             dataa = abs(dataa); %since we are interested in connector hubs (which are negative)
                    rmin = min(dataa);
                    rmax = max(dataa);
                    %scaling connection strengths with extreme values (limitt) requested by
                    size_con = (dataa-rmin)./(rmax-rmin).*(limitt(2)-limitt(1)) + limitt(1);
                else %otherwise scaling together group1 and group2
                    rmin = rmin12; rmax = rmax12;
                    %scaling connection strengths with extreme values (limitt) requested by
                    size_con = (dataa-rmin)./(rmax-rmin).*(limitt(2)-limitt(1)) + limitt(1);
                    size_con2 = limitt(1) + (limitt(2) - size_con);
                end
                if cc == 3
                    %creating nifti image for WorkBench
                    vector = zeros(90,1);
                    vector(negidx2) = size_con;
                    create_AALnifti(vector,['DTI_connectorhubs_contr_' num2str(cc) '.nii'],1);
                else
                    %creating nifti image for WorkBench
                    vector = size_con2;
                    create_AALnifti(vector,['DTI_connectorhubs_contr_' num2str(cc) '.nii'],1);
                end
                for gg = 1:length(dataa) %over AAL ROIs
                    hold on
                    if cc == 3
                        plot3(MNI_RC(negidx2(gg),1), MNI_RC(negidx2(gg),2), MNI_RC(negidx2(gg),3), '.', 'Color', coldum(cc), 'MarkerSize', size_con(gg)); %centroid of ROI ii
                    else
                        plot3(MNI_RC(gg,1), MNI_RC(gg,2), MNI_RC(gg,3), '.', 'Color', coldum(cc), 'MarkerSize', size_con2(gg)); %centroid of ROI ii
                    end
                end
                rotate3d on; axis off; axis vis3d; axis equal
                set(gcf,'color','w')
                title(['connector hubs - contrast ' num2str(cc)])
                view([0 90])
%                 saveas(gcf,['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Images_Paper_IQ/Dots_top_DTI_connector_contr_' num2str(cc) '.svg'])
                view([-90 0])
%                 saveas(gcf,['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Images_Paper_IQ/Dots_left_DTI_connector_contr_' num2str(cc) '.svg'])
            end
        else
            disp('no significant connector hubs..')
        end
    end
end
% Density
if label_density ~= 0
    list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Graph_Theory_AAL90Matrix/Subj*'); %list of subjects
    Mones=ones(90);
    Mcat=[]; %initialize matrix which will contain the concatenated matrices of every subject
    perc=label_density; %percentage of values we want to remove for thresholding
    for ii = 1:length(list) %over subjects
        %loading matrix with strenght of connections (AAL, same ROIs used for distance, of course..)
        load([list(ii).folder '/' list(ii).name]); %Load matrix (M3) containing the strength of the connections across regions of interest
        Mtriu=M3(find(triu(Mones,1)~=0));
        Mcat=cat(1,Mcat, Mtriu);
    end
    thresh1=sort(Mcat); %sort Mcat
    indthresh1=round(length(thresh1)/100*perc); %find indeces corresponding to 1% of Mcat
    thresh1val=thresh1(indthresh1); %values corresponding to 1% of Mcat
    density = zeros(length(list),1);
    for ii = 1:length(list) %over subjects
        %loading matrix with strenght of connections (AAL, same ROIs used for distance, of course..)
        load([list(ii).folder '/' list(ii).name]); %Load matrix (M3) containing the strength of the connections across regions of interest
        M4=zeros(90);
        M4(M3>thresh1val)=M3(M3>thresh1val);
        density(ii,1) = density_und(M4); %compute degree of connectivity
    end
%     figure
%     scatter(1:68,density) %plot
    %ylim([-100 100])
    highPR = density(desmat(:,1)==1,1); %density of subjects with high PR
    lowPR = density(desmat(:,2)==1,1); %density of subjects with low PR
    figure
    scatter(1:length(highPR), highPR, 'o', 'r');
    hold on
    scatter(1:length(lowPR), lowPR, 'o', 'b');
    %t-test
    [~,p,~,stats] = ttest2(highPR,lowPR);
    title(['density - pval ' num2str(p) ' - tval ' num2str(stats.tstat)])
end
%Distance (distM), characteristich path length (CPL), global (GE), nodal (NE) and local efficiency (LE), modularity (maybe)
if label_additiona_graph_measures == 1 %additional measures: 1 = char path length; 2 = glob eff; 3 = nod eff; 4 = loc eff; 5 = modularity
    list = dir('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Graph_Theory_AAL90Matrix/Subj*'); %list of subjects
    Mdiag = diag(1:90);
    idx = find(Mdiag~=0);
    CPL = zeros(length(list),1);
    GE = zeros(length(list),1);
    NE = zeros(length(list),90);
    LE = zeros(length(list),1);
    MOD = zeros(length(list),1);
    for ii = 1:length(list) %over subjects
        %loading matrix with strenght of connections (AAL, same ROIs used for distance, of course..)
        load([list(ii).folder '/' list(ii).name]); %Load matrix (M3) containing the strength of the connections across regions of interest
        mM3=M3/(max(max(M3))); %Divide every value of M3 by the greatest value within M3 -> get weighted values between 0 and 1
        iM3=1-mM3; %create the inverse of mM3 (strongest connections = shortest paths)
        [distM, edges]=distance_wei(iM3); %Compute distance matrix containing the shortest path length between each pair of nodes
        distM(idx) = NaN; %assigning NaNs to diagonal
        CPL(ii,1) = nanmean(nanmean(distM)); %characteristic path length
        GE(ii,1) = nanmean(nanmean(1./distM)); %global efficiency
        NE(ii,:) = nanmean(1./distM);%nodal efficiency
        LE(ii,1) = mean(NE(ii,:)); %local efficiency
        [~,modd] = modularity_und(M3); %Newmann modularity (deterministic)
        MOD(ii,1) = modd;
        disp(ii)
    end
    %plotting
    highPR = CPL(desmat(:,1)==1,1); %characteristic path length of subjects with high PR
    lowPR = CPL(desmat(:,2)==1,1); %characteristic path length of subjects with low PR
    figure
    scatter(1:length(highPR), highPR, 'o', 'r');
    hold on
    scatter(1:length(lowPR), lowPR, 'o', 'b');
    [~,p,~,stats] = ttest2(highPR,lowPR);
    title(['CPL - pval ' num2str(p) ' - tval ' num2str(stats.tstat)])
    % Global efficiency
    highPR = GE(desmat(:,1)==1,1); %Global efficiency of subjects with high PR
    lowPR = GE(desmat(:,2)==1,1); %Global efficiency of subjects with low PR
    figure
    scatter(1:length(highPR), highPR, 'o', 'r');
    hold on
    scatter(1:length(lowPR), lowPR, 'o', 'b');
    [~,p,~,stats] = ttest2(highPR,lowPR);
    title(['GE - pval ' num2str(p) ' - tval ' num2str(stats.tstat)])
    %Nodal efficiency
    highPR = NE(desmat(:,1)==1,:); %nodal efficiency of subjects with high PR
    lowPR = NE(desmat(:,2)==1,:); %nodal efficiency of subjects with low PR
    diff2=mean(highPR,1)-mean(lowPR,1); %Difference between highPR and lowPR
    figure
    scatter(1:90, diff2);
    [~,p,~,stats] = ttest2(highPR,lowPR);
    disp('NE')
    lab2(find(p<0.05),:)
    %     title(['NE - pval ' num2str(p)])
    title('NE')
    % Local efficiency
    highPR = LE(desmat(:,1)==1,1); %Local efficiency of subjects with high PR
    lowPR = LE(desmat(:,2)==1,1); %Local efficiency of subjects with low PR
    figure
    scatter(1:length(highPR), highPR, 'o', 'r');
    hold on
    scatter(1:length(lowPR), lowPR, 'o', 'b');
    [~,p,~,stats] = ttest2(highPR,lowPR);
    title(['LE - pval ' num2str(p) ' - tval ' num2str(stats.tstat)])
    %Modularity
    highPR = MOD(desmat(:,1)==1,1); %characteristic path length of subjects with high PR
    lowPR = MOD(desmat(:,2)==1,1); %characteristic path length of subjects with low PR
    figure
    scatter(1:length(highPR), highPR, 'o', 'r');
    hold on
    scatter(1:length(lowPR), lowPR, 'o', 'b');
    [~,p,~,stats] = ttest2(highPR,lowPR);
    title(['mod - pval ' num2str(p) ' - tval ' num2str(stats.tstat)])
end
if label_modularity_brainplot == 1
    %loading MNI coordinates of AAL 2-mm centroids
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/MNI_RC.mat')
    %getting mean over groups of subjects and then computing modularity
    MM1 = mean(MM(:,:,desmat(:,1)==1),3); %group 1
    [com1,modl1] = modularity_und(MM1); %Newmann modularity (deterministic)
    COMM = zeros(length(com1),2);
    COMM(:,1) = com1; %storing ROIs communities ID for group 1
    MM2 = mean(MM(:,:,desmat(:,2)==1),3); %group 2
    [com1,modl2] = modularity_und(MM2); %Newmann modularity (deterministic)
    COMM(:,2) = com1; %storing ROIs communities ID for group 2
    %permutations to establish whether the brain is more "modulable" than a random graph with its same elements
    if mod_perm_label == 1
        %group1
        [~,m] = size(MM1);
        permut = 1000;
        MODP1 = zeros(permut,1);
        disp('computing permutations ofr MCS on modularity - group 1')
        for pp = 1:permut
            f = find(triu(ones(m),1) == 1); %find indexes of elements of the upper triangle of the square matrix
            idx_dummy = randperm(length(f)); %create a permuted array from 1 to length(f)
            idx_2 = f(idx_dummy); %shuffle the indexes in f according to the order in idx_dummy
            r_dummy = zeros(m*m,1); %initialise vector
            r_dummy(f) = MM1(idx_2); %taking element in matrix data with index idx_2 (shuffled indexes of the upper triangle)
            %rebuild 2D matrix with actual data
            r3 = reshape(r_dummy,[m,m]);
            %ricreate the simmetric matrix for calculating more easily the sums
            r3 = (r3 + r3') - eye(size(r3,1)).*diag(r3);
            %calculating segregation of the randomly rearranged matrix
            [~,mod_glob_r3] = modularity_und(r3); %only global value of segregation
            MODP1(pp) = mod_glob_r3;
%             disp(pp)
        end
        pval1 = length(find(modl1<MODP1))/permut;
        disp(['modularity p-value group 1 = ' num2str(pval1)])
        %group1
        [~,m] = size(MM2);
        permut = 1000;
        MODP2 = zeros(permut,1);
        disp('computing permutations ofr MCS on modularity - group 2')
        for pp = 1:permut
            f = find(triu(ones(m),1) == 1); %find indexes of elements of the upper triangle of the square matrix
            idx_dummy = randperm(length(f)); %create a permuted array from 1 to length(f)
            idx_2 = f(idx_dummy); %shuffle the indexes in f according to the order in idx_dummy
            r_dummy = zeros(m*m,1); %initialise vector
            r_dummy(f) = MM2(idx_2); %taking element in matrix data with index idx_2 (shuffled indexes of the upper triangle)
            %rebuild 2D matrix with actual data
            r3 = reshape(r_dummy,[m,m]);
            %ricreate the simmetric matrix for calculating more easily the sums
            r3 = (r3 + r3') - eye(size(r3,1)).*diag(r3);
            %calculating segregation of the randomly rearranged matrix
            [~,mod_glob_r3] = modularity_und(r3); %only global value of segregation
            MODP2(pp) = mod_glob_r3;
%             disp(pp)
        end
        pval2 = length(find(modl2<MODP2))/permut;
        disp(['modularity p-value group 2 = ' num2str(pval2)])
    end
    %elaborated way to match the color of most similar communities between groups
    if max(COMM(:,2)) > max(COMM(:,1)) %working finding communities of group with less communities in group with more communities
        labcom1 = 1;
        labcom2 = 2;
    else
        labcom1 = 2;
        labcom2 = 1;
    end
    matchh = zeros(max(COMM(:,labcom1)),2);
    for jj = 1:max(COMM(:,labcom1))
        dum = COMM(find(COMM(:,labcom1)==jj),labcom2); %getting elements of group2 with indices of group1 == community jj
        a = unique(dum); %getting unique values
        PP = zeros(1,length(a));
        for pp = 1:length(a) %over communities in dum
            PP(pp) = length(find(dum==a(pp))); %finding occurrencies of unique values (communities)
        end
        [~,ipp] = max(PP); %best match between communities in group1 and group2
        matchh(jj,1) = jj; %storing community for group1
        matchh(jj,2) = a(ipp); %with best corresponding community of group2
    end
    a4 = 1:max(COMM(:,labcom1));
     if labcom1 == 2 %matching communities color here
        %color codes using letters
%         COL{1} = ['r','m','k','g','y','c']; %assigning colors to communities of group1
%         COL{2} = COL{1}(matchh(:,2)); %colors to community of group2 from match between the two groups
%         COL{2}  = cat(2,COL{2},COL{1}(find(sum(double(a4==matchh(:,2)),1)==0)));
        %color codes using numbers (triplets)
        COL{1} = [0 0.6 1; 1 0.4 0; 1 0 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0]; %vector with colors
        COL{2} = COL{1}(matchh(:,2),:); %colors to community of group2 from match between the two groups
        COL{2}  = cat(1,COL{2},COL{1}(find(sum(double(a4==matchh(:,2)),1)==0),:));
    else
        %color codes using letters
%         COL{2} = ['r','m','k','g','y','c']; %assigning colors to communities of group2
%         COL{1} = COL{2}(matchh(:,2)); %colors to community of group1 from match between the two groups
%         COL{1}  = cat(2,COL{1},COL{2}(find(sum(double(a4==matchh(:,2)),1)==0)));
        %color codes using numbers (triplets)
        COL{2} = [0 0.6 1; 1 0.4 0; 1 0 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0]; %vector with colors
        COL{1} = COL{2}(matchh(:,2),:); %colors to community of group2 from match between the two groups
        COL{1}  = cat(1,COL{1},COL{2}(find(sum(double(a4==matchh(:,2)),1)==0),:));
    end
    for ii = 1:2 %over experimental groups
        S = [];
        S.intra_con = intra_con; %1 for intra-connectivity plots
        if intra_con == 1
            inter_con = 0;
        else
            inter_con = 1;
        end
        S.inter_con = inter_con; %1 for inter-connectivity plots
        S.perc_intra = 50; %percentage of the connections that you want to plot in each community (intra)
        S.perc_inter = 3; %percentage of the connections that you want to plot between each couple of communities (inter)
        S.MNI_RC = MNI_RC;
        S.conn_matrix = mean(MM(:,:,desmat(:,ii)==1),3); %connectivity matrix
        S.communities = COMM(:,ii);
        S.colors_mod = COL{ii};
        S.col_dot = 'b';
        S.limits = limitt_mod;
        S.title = ['group ' num2str(ii)];
        %actual plotting function
        GT_modul_plot_LBPD(S)

        %printing optimal community structure in excel files (high IQ)
        PD = cell(91,max(COMM(:,1)));
        for ll = 1:max(COMM(:,1))
            PD{1,ll} = ['Module ' num2str(ll)];
            PD(2:size(lab2((COMM(:,1)==ll),:),1)+1,ll) = cellstr(lab2((COMM(:,1)==ll),:));
        end
        PDnh = cell2table(PD); %remove the possible empty cell
        writetable(PDnh,['DTI_High_IQ.xlsx'],'Sheet',ll) %printing excel file
        %printing optimal community structure in excel files (average IQ)
        PD = cell(91,max(COMM(:,2)));
        for ll = 1:max(COMM(:,2))
            PD{1,ll} = ['Module ' num2str(ll)];
            PD(2:size(lab2((COMM(:,2)==ll),:),1)+1,ll) = cellstr(lab2((COMM(:,2)==ll),:));
        end
        PDna = cell2table(PD); %remove the possible empty cell
        writetable(PDna,['DTI_Average_IQ.xlsx'],'Sheet',ll) %printing excel file
        %saving images
%         if intra_con == 1
%             view([0 90])
%             saveas(gcf,['IntraMod_DTI_group_' num2str(ii) '_top.svg'])
%             view([-90 0])
%             saveas(gcf,['IntraMod_DTI_group_' num2str(ii) '_left.svg'])
%         else
%             view([0 90])
%             saveas(gcf,['InterMod_DTI_group_' num2str(ii) '_top.svg'])
%             view([-90 0])
%             saveas(gcf,['InterMod_DTI_group_' num2str(ii) '_left.svg'])
%         end
        if label_schem_mod == 1
            %testing schemaball with modularity
            figure
            order = [1:2:90 90:-2:2];
            load('/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat'); %loading labels
            lab2 = lab(order,:);
            r = S.conn_matrix./max(max(S.conn_matrix)); %scaling matrix
            perc_intra = 50; %TO BE DEFINED
            perc_inter = 3;
            limitt = [0.01 3]; %TO BE DEFINED
            %actual function
            schemaball_modularity_LBPD(r, lab2, S.communities, COL{ii}, perc_intra, perc_inter, limitt)
            title(['group ' num2str(ii)])
            set(gcf,'color','w')
%             export_fig(['DTI_schemball_group_' num2str(ii) '.eps'])
        end
    end
    %CENTROIDS OF SUBNETWORKS (FROM MODULARITY)
    %%% PROBABLY RELEVANT, maybe in supplementary (at the end we did not report it.. not that relevant..)
%     centroids_l = 0;
%     if centroids_l == 1
%         MS = [30,10]; %size of centroids for the 2 groups
%         openfig('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/BrainTemplate_GT.fig')
%         load('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI_RC_AAL_2mm.mat')
%         for gg = 1:2 %over experimental groups
%             for cc = 1:max(COMM(:,gg)) %over communities of group gg
%                 ind = find(COMM(:,gg)==cc); %getting indices of community cc of group gg
%                 centroid = mean(MNI_RC(ind,:),1); %MNI 3D coordinates of centroid of community cc of group gg
%                 hold on
%                 plot3(centroid(1), centroid(2), centroid(3), '.', 'Color', COL{gg}(cc,:), 'MarkerSize', MS(gg)); %centroid of ROI ii
%             end
%         end
%         rotate3d on; axis off; axis vis3d; axis equal
%         set(gcf,'color','w')
%         title(['centroids of communities'])
%         view([0 90])
% %         saveas(gcf,['/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021/Modularity/Centroids_DTI_AllSubjs_top.svg'])
%         view([-90 0])
% %         saveas(gcf,['/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021/Modularity/Centroids_DTI_AllSubjs_left.svg'])
%     end
end
%plotting average FC matrix
if label_matrix_plot == 1    
    Mdiag = diag(1:90);
    idx = find(Mdiag~=0);
    MM2 = mean(MM,3);    
    MM2(idx) = NaN; %NaNs to diagonal
    figure
    imagesc(MM2)
    colorbar
    % clims = [0 1];
%     caxis(clims)
%     caxis([0 160])
    set(gcf,'color','w')
    %changing colormap of figures (red-bue with white for 0 values)
    x = []; x.bottom = [0 0 0.5]; x.botmiddle = [0 0.5 1]; x.middle = [1 1 1]; x.topmiddle = [1 0 0]; x.top = [0.6 0 0]; %red - blue
    colormap(bluewhitered_PD(0,x))
    export_fig('/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021/Matrix_DTI_AllSubjs.eps')
    %schemaball plotting    
    S = [];
    S.MATF = MM2; %connectivity matrix ROIs x ROIs (can be both binary or non-binary) (double)
    S.symmetric_l = 1; %1 if MATF is LLLRRR; 0 if MATF is LRLRLR
    S.outpath = '/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021'; %path where saving images (characters)
    S.p = []; %parcellation object (created by using OSL)
    S.thr_cortex = []; %  strongest % of connectivity to be plotted (it is usually a quite high number, e.g. 0.993). If empty [] default is 0.993.
    S.name_gif = []; %name of gif/video (character). Leave empty [] if you do not want to plot the connectivity in the brain (e.g. you want only the schemaball).
    S.frame_vect = []; %vector with frames you want. If empty [] default is [15,90], another good option could be [1,15,59,90,130]
    S.fr_spec = []; %cell array with characters specifing the name of the requested frames. E.g. [15,90] corresponds to fr_spec = {'Posterior_L_R';'Frontal_R_L'} [1,15,59,90] corresponds to fr_spec = {'Orig_Angle';'Posterior_L_R';'Hem_R';'Frontal_R_L';'Hem_L'}
    S.schball_l = 1; %1 for having the schemaball; 0 for not having it.
    S.extr = []; %minimum and maximum values to be used for normalizing the matrix to be submitted to schemaball function.
        %Leave empty [] for scaling the matrix on its max(abs value).
    S.lab_parc_ph = '/projects/MINDLAB2017_MEG-LearningBach/scripts/scripts_osl_learningbach/LEiDA/AAL_labels.mat'; %path to labels. It must be a .mat file with a character array named 'lab' (characters).
    S.colbacksch = []; %background color for shcemaball. If empty [], default is 'w', white.
    S.name_sch = [];
    S.schemplth = 100; %percentage of the strongest connections to be plotted in schemaball (5% to be submitted as '5').
    S.name_sch = ['DTI_allSubjs']; %name for saved schemaball image. Leave empty [] for not saving it.
    %actual function
    FC_Brain_Schemaball_plotting_LBPD_D(S)
    %plotting connectivity in the brain (using the modularity plotting with all ROIs belonging to 1 community)
    %loading MNI coordinates of AAL 2-mm centroids
    load('/scratch7/MINDLAB2017_MEG-LearningBach/Leonardo/MNI_RC.mat')
    S = [];
    S.intra_con = 1; %1 for intra-connectivity plots
    S.inter_con = 0; %1 for inter-connectivity plots
    S.perc_intra = 10; %percentage of the connections that you want to plot in each community (intra)
    S.perc_inter = 0; %percentage of the connections that you want to plot between each couple of communities (inter)
    S.MNI_RC = MNI_RC;
    S.conn_matrix = MM2; %connectivity matrix
    S.communities = ones(90,1); %assigning all ROIs to just 1 community
    S.colors_mod = [0.7 0 0]; %; 1 0.4 0; 1 0 0; 0.6 0.2 0.6; 0 0 0.6; 0 0 0];
    S.col_dot = [0 0 0.5];
    S.limits = [0.1 5];
    S.title = ['group ' num2str(ii)];
    %actual plotting function
    GT_modul_plot_LBPD(S)
    view([0 90])
%     saveas(gcf,['/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021/Full_DTI_AllSubjs_top.svg'])
    view([-90 0])
%     saveas(gcf,['/aux/MINDLAB2020_MEG-AuditoryPatternRecognition/IQ_GT_Portis/April_06_2021/Full_DTI_AllSubjs_left.svg'])
end

%%
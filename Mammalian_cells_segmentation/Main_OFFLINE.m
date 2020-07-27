clear all;
close all;
clc;
exit_flag=1;

% This MAIN code is needed to set up the segmentation
% algorithm. The user can define the  .mat file from
% which extracting cropping variables. This .mat file is originally
% generated directly from the microscope when running experiments and
% contains the crop variables that are saved as 'CropRect<i>' where
% i is the position of interest, and 'CropRectBCKG' to define the 
% region where background fluorescence should be calculated. 
% Crop variables can also be manually generated and saved in a file called 'cropfiles.mat'  
% with the name of 'crop' and 'CropRectBCKG'. If the
% user wants to load crops from 'cropfiles.mat' s/he must first rename the
% variable 'crop' into 'CropRect<i>', as this is the format accepted by the
% CROPfunction.
% The user can define 'channel numbers'; in our example, they are
% associated with names of raw images from the microscope; for example, the
% raw phase contrast image of Position 3 is called
% 'Position003...--C00.jpg', where C00 identifies the channel associated to
% the phase contrast. Here a list of channels in our experiments is shown.

% Rex1-GFPd2 mESCs (cell cluster segmentation):
% C00 = Phase Contrast
% C01 = GFP
% C02 = H2B
% C03 = Blue Dye

% dual-reporter mESCs (single cell tracking):
% C00 = Phase Contrast
% C01 = mCherry
% C02 = GFP
% C03 = H2B
% C04 = Blue Dye

% Please read also the "Mammalian cell segmentation" document

while (exit_flag)
    prompt = 'Are you running a tracking experiment? (Answer Y for yes and N for no): \n ';
    temp=input(prompt, 's');
    
    if strcmp(temp,'Y') || strcmp(temp,'y')
        tracking_flag=1;
  
        exit_flag=0;
        
    elseif strcmp(temp,'N') || strcmp(temp,'n')
        tracking_flag=0;
        exit_flag=0;
       
    else
        fprintf('Typing Error!\n');
    end
    fprintf('\n');
end

prompt = 'Insert the number of channels for the experiment (2, 3 or 4): \n ';
n_Channels = input(prompt);fprintf('\n');

exit_flag = 1;
if n_Channels>2
    while (exit_flag)
        prompt = 'Do you have a Blue Dye Channel? (Answer Y for yes and N for no): \n ';
        Blue_dye_flag_temp=input(prompt, 's');
        if strcmp(Blue_dye_flag_temp,'Y') || strcmp(Blue_dye_flag_temp,'y')
            Blue_dye_flag=1;
            prompt = 'Which is the Channel associated to the blue dye? (Answer with anumber, i.e.4) \n ';
            Parameters.fluoeval.blue_dye_channel = input(prompt);
            exit_flag=0;
            
        elseif strcmp(Blue_dye_flag_temp,'N') || strcmp(Blue_dye_flag_temp,'n')
            Blue_dye_flag=0;
            exit_flag=0;
            
        else
            fprintf('Typing Error!\n');
        end
    end
    
    fprintf('\n');
    exit_flag=1;
    
    while (exit_flag)
        prompt = 'Do you have a Channel for the Nucleo tag? (Answer Y for yes and N for no): \n ';
        IRFP_flag_temp=input(prompt, 's');
        if strcmp(IRFP_flag_temp,'Y') || strcmp(IRFP_flag_temp,'y')
            prompt = 'Which is the Channel for the Nucleo tag? (Answer with anumber, i.e.3) \n ';
            Parameters.fluoeval.nucleo_tag_channel = input(prompt);
            IRFP_flag=1;
            exit_flag=0;
            
        elseif strcmp(IRFP_flag_temp,'N') || strcmp(IRFP_flag_temp,'n')
            IRFP_flag=0;
            exit_flag=0;
            
        else
            fprintf('Typing Error!\n');
        end
        fprintf('\n');
    end
    
else
    Blue_dye_flag=0;
    IRFP_flag=0;
end
exit_flag = 1;

while (exit_flag)
    prompt = 'Do you need a crop for your frames? (Answer Y for yes and N for no): \n ';
    crop_temp=input(prompt, 's');
    
    if strcmp(crop_temp,'Y') || strcmp(crop_temp,'y')
        crop_flag=1;
        prompt = 'Do you want to draw the crop? (Answer Y for yes and N for no): \n IF NO, THE CROP WILL BE UPLOADED FROM A MAT FILE \n';
        crop_ext_temp=input(prompt, 's');
        
        if strcmp(crop_ext_temp,'Y') || strcmp(crop_ext_temp,'y')
            crop_drawn = 1;
            
        elseif strcmp(crop_ext_temp,'N') || strcmp(crop_ext_temp,'n')
            prompt = 'Insert the name of the .mat file from which extracting the crop variable \n';
            Parameters.name_matfile = input(prompt,'s');
            crop_drawn = 0;
        end
        
        exit_flag=0;
        
    elseif strcmp(crop_temp,'N') || strcmp(crop_temp,'n')
        crop_flag = 0;
        crop_ext = 0;
        exit_flag=0;
        
    else
        fprintf('Typing Error!\n');
    end
    fprintf('\n');
end

prompt = 'Which number is the channel you want to use to evaluate the mask? (insert the number): \n ';
mask_channel=input(prompt);fprintf('\n');

prompt = 'Which number is the channel you want to use to evaluate the fluorescence? (insert the number): \n ';
fluorescence_channel=input(prompt);fprintf('\n');

prompt = 'Insert the number of timeframes to analyse for the experiment: \n ';
n_timeframes = input(prompt);fprintf('\n');

prompt = 'Insert the position of interest for the experiment (can be the controlled one, or simply a position you want to analyse offline): \n ';
position_of_interest = input(prompt);

prompt = 'Insert the position of the background (if the background is computed as the ANTIMASK insert 0): \n ';
position_of_background = input(prompt);fprintf('\n');

prompt = 'Insert the file name for the VideoMaker: \n ';
video_name = input(prompt,'s');fprintf('\n');

prompt = 'Insert the directory of this Matlab folder plus a slash ( \\ or /) \n according to the path method of your operative system \n (i.e. C:Users\\ ... \\Name_matlab_folder\\ ): \n ';
dir_mat_folder = input(prompt,'s');fprintf('\n');

prompt = 'Insert the name of the folder containing the images to be analysed plus a slash ( \\ or /) \n according to the path method of your operative system \n (i.e. Mark_and_Find_001\\ for Windows, Mark_and_Find_001/ for MacOS): \n ';
img_folder = input(prompt,'s');fprintf('\n');

Parameters.n_timeframes = n_timeframes;% This is the number of hours associated to the experiment
Parameters.position_of_interest = position_of_interest;
Parameters.crop_flag = crop_flag;
Parameters.crop_drawn = crop_drawn;
Parameters.mat_figures_and_data_folder = 'mat_figures_and_data_folder';
Parameters.tracking_flag = tracking_flag;

Parameters.fluoeval.n_Channels = n_Channels;
Parameters.fluoeval.position_of_background = position_of_background;
Parameters.fluoeval.Blue_dye_flag = Blue_dye_flag;
Parameters.fluoeval.IRFP_flag = IRFP_flag;
Parameters.fluoeval.mask_channel = mask_channel;
Parameters.fluoeval.fluorescence_channel = fluorescence_channel;
Parameters.fluoeval.cell_of_interest = 5; %we require to identify a label for a 
                                            %single cell, in order
                                            % to compute a fluorescence
                                            % vector associated. This lable
                                            % is generally not available
                                            % before a first run of the
                                            % code. The authors then
                                            % suggest to have a first run
                                            % without the computation of
                                            % fluorescence, just to detect
                                            % labels associated to single
                                            % cells.

Parameters.Video.video_name = video_name;

%PATH VARIABLES
Parameters.Paths.directory_mat_folder = dir_mat_folder;
Parameters.Paths.img_folder = img_folder;

% INSERT THE NAME OF THE EXPERIMENT FOR THE MAT FILE
save('Parameters_for_name_of_the_experiment.mat', 'Parameters' );
% save('Parameters_for_Cluster_segmentation.mat', 'Parameters');
% save('Parameters_for_test_SingleCells.mat', 'Parameters');

%% OFFLINE FUNCTION CALL
global CELLS
CELLS = cell(3,1);
CELLS{3,1} = 1; 

[outputArg1,Parameters] = OFFLINE_function(Parameters);
save('results_for_name_of_the_experiment.mat','outputArg1', 'Parameters' );
% save('Results_Parameters_for_Cluster_segmentation.mat','outputArg1', 'Parameters', 'CELLS' );

%% PLOT

plot(0:60:(Parameters.n_timeframes-1)*60,outputArg1.Fluorescence_Channel_A_A(:,2),'color','k','linewidth',1.5)
xlim([0 Parameters.n_timeframes]*60);
ylim([min(outputArg1.Fluorescence_Channel_A_A(:,2)) max(outputArg1.Fluorescence_Channel_A_A(:,2))]);
ylabel('Fluorescence [A.u.]','fontsize',16);
xlabel('time [min]','fontsize',16);
set(gca,'Fontsize',16);
hold on
% plot([5*60 5*60],[min(outputArg1.Fluorescence_Channel_A_A(:,2)) max(outputArg1.Fluorescence_Channel_A_A(:,2))],'--','linewidth',1.5,'Color',[1 0 0])
% plot([15*60 15*60],[min(outputArg1.Fluorescence_Channel_A_A(:,2)) max(outputArg1.Fluorescence_Channel_A_A(:,2))],'--','linewidth',1.5,'Color',[1 0 0])

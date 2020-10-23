function [Output,Parameters] = OFFLINE_function(Parameters)

close all
clc
warning off;
global CELLS

%% Definition of parameters
n_CHANNELS = Parameters.fluoeval.n_Channels;
n_timeframes = Parameters.n_timeframes;
selected_BCKGfield = Parameters.fluoeval.position_of_background; %if 0, background as antimask
tracking_flag = Parameters.tracking_flag;

%% Fluorescence vectors
Output.Absolute_Fluorescence_AVG = [];
Output.Fluorescence_Channel_A_A = [];
Output.Nuclei_Normalized_Fluorescence_Channel_A_A = [];
Output.Average_Background_fluorescence = [];
Output.MASK = {};
%% CROP FUNCTION call
Parameters = CROPfunction(Parameters);

%% MAIN LOOP with mask and fluorescence computation
for i=1:n_timeframes
    % Function defining addresses associated to images
    [cropped_image,fileNames, fileNamesBCKGRND] = fileNames_addressFunction(i, Parameters);
       
    for j= 1:n_CHANNELS
        Parameters.fluoeval.CHANNEL{j} = fileNames{j};
        if (selected_BCKGfield~=0)
            Parameters.fluoeval.BCKGRND_CHANNEL{j} = fileNamesBCKGRND{j};
        end
    end
    
    if (tracking_flag==1)
    % single cells mask n track
    [XYpos,Radii,orientation, num_el, mask, singleMASK]=Mask_n_Track(Parameters.fluoeval.CHANNEL{Parameters.fluoeval.mask_channel},...
        Parameters.fluoeval.CHANNEL{Parameters.fluoeval.nucleo_tag_channel},...
        Parameters.fluoeval.crop, i);
    % fluorescence computation
        [Output_fluo_eval] = fluorescence_evaluation(Parameters.fluoeval, singleMASK(Parameters.fluoeval.cell_of_interest).pixel);

    else
        % cell clusters mask
        [mask, Output.area_mask(i)] = mask_function(Parameters.fluoeval.CHANNEL{Parameters.fluoeval.mask_channel},Parameters.fluoeval.crop);
        [Output_fluo_eval] = fluorescence_evaluation(Parameters.fluoeval, mask);
        Output.MASK{i} = mask;
    end
    % saving fluorescence values
    Output.Absolute_Fluorescence_AVG = [Output.Absolute_Fluorescence_AVG ; Output_fluo_eval.Absolute_Fluorescence_AVG];
    Output.Fluorescence_Channel_A_A = [ Output.Fluorescence_Channel_A_A ; Output_fluo_eval.Fluorescence_Channel_A_A ];
    
    if (Parameters.fluoeval.IRFP_flag~=0)
        Output.Nuclei_Normalized_Fluorescence_Channel_A_A = [ Output.Nuclei_Normalized_Fluorescence_Channel_A_A ; Output_fluo_eval.Nuclei_Normalized_Fluorescence_Channel_A_A ];
    end
    Output.Average_Background_fluorescence = [ Output.Average_Background_fluorescence ; Output_fluo_eval.Average_Background_fluorescence ];
    
    
end
% VIDEO GENERATOR
Generate_Video(Parameters,Output);
% saving current workspace
curr_filename = strcat('field_',num2str(Parameters.position_of_interest));
save(curr_filename);
end

%% FUNCTION FOR THE CROP
function [Parameters] = CROPfunction(Parameters)
%Input: stack of parameters that are set up in the main file.
%Output: stack of parameters updated with the crop variable.
source = Parameters.Paths.directory_mat_folder ;
imageDir_exp = Parameters.Paths.img_folder ;
position_of_interest =Parameters.position_of_interest;
selected_BCKGfield = Parameters.fluoeval.position_of_background;

if (Parameters.crop_flag==1) 
    if (Parameters.crop_drawn==1) 
        % names of images must be changed according to the experimental
        % setup and users standards
        name_position_ctrld = strcat('Position', ...
            pad(num2str(position_of_interest), 3, 'left', '0'), ...
            '--t', pad(num2str(0), 2, 'left', '0'),...
            '--C', pad(num2str(0), 2, 'left', '0'), '.tif');
        %-----------------------------------------------------------
        image = strcat(source,imageDir_exp,name_position_ctrld);
        [~,crop] = imcrop(imadjust(imread(image)));
        save('cropfiles.mat','crop'); 
        Parameters.fluoeval.crop = crop;
        if (Parameters.fluoeval.position_of_background~=0)
            % names of images must be changed according to the experimental
            % setup and users standards
            name_BCKGRND_position = strcat('Position', ...
                pad(num2str(selected_BCKGfield), 3, 'left', '0'), ...
                '--t', pad(num2str(0), 2, 'left', '0'),...
                '--C', pad(num2str(0), 2, 'left', '0'), '.tif');
            %-------------------------------------------------------------
            imageBCKG = strcat(source,imageDir_exp,name_BCKGRND_position);
            [~,CropRectBCKG]=imcrop(imadjust(imread(imageBCKG)));
            save('cropfiles.mat','CropRectBCKG','-append');
            Parameters.fluoeval.CropRectBCKG = CropRectBCKG;
        end
    else
        crop_var = strcat('CropRect',num2str(position_of_interest));
        temp_struct = load(Parameters.name_matfile,crop_var);
        crop = temp_struct.(crop_var);
        Parameters.fluoeval.crop = crop;
        
        if (selected_BCKGfield~=0)
            temp_struct = load(Parameters.name_matfile,'CropRectBCKG');
            CropRectBCKG = temp_struct.('CropRectBCKG');
            Parameters.fluoeval.CropRectBCKG = CropRectBCKG;
        end
    end
else
    Parameters.fluoeval.crop = 0;
end

end

%% FUNCTION TO READ CROPPED IMAGES AND GET IMAGES PATHS
% Images names must be changed according to the user's standard
function [cropped_image,fileNames, fileNamesBCKGRND] = fileNames_addressFunction(i, Parameters)
% Input: a struct of parameters, and the current
% time 'i'. Output: a n*m matrix (cropped_image), where n*m depends on the crop size, the complete
% address to the position of interest (fileNames) and the complete address
% to the background position (fileNamesBCKGRND)
fileNames = {};
uncropped_image = {};
cropped_image = {};
fileNamesBCKGRND = {};

for type_i = 1:Parameters.fluoeval.n_Channels
    % names of images must be changed according to the experimental setup
    ctrl_position = strcat('Position', ...
        pad(num2str(Parameters.position_of_interest), 3, 'left', '0'), ...
        '--t', pad(num2str(i), 2, 'left', '0'),...
        '--C', pad(num2str(type_i - 1), 2, 'left', '0'), '.tif');
    %---------------------------------------------------------------
    ctrl_position = fullfile(Parameters.Paths.img_folder , ctrl_position);
    fileNames{type_i} = fullfile(Parameters.Paths.directory_mat_folder , ctrl_position);
    uncropped_image{type_i} =imread(strcat(Parameters.Paths.directory_mat_folder ,ctrl_position));
    
    if (Parameters.crop_flag~=0)
        cropped_image{type_i}=imcrop(uncropped_image{type_i},Parameters.fluoeval.crop);
    else
        cropped_image{type_i}=uncropped_image{type_i};
    end
    
    if (Parameters.fluoeval.position_of_background~=0)
        % names of images must be changed according to the experimental setup
        BCKGRND_position = strcat('Position', ...
            pad(num2str(Parameters.fluoeval.position_of_background), 3, 'left', '0'), ...
            '--t', pad(num2str(i), 2, 'left', '0'),...
            '--C', pad(num2str(type_i - 1), 2, 'left', '0'), '.tif');
        %---------------------------------------------------------------
        BCKGRND_position = fullfile(Parameters.Paths.img_folder , BCKGRND_position);
        fileNamesBCKGRND{type_i} = fullfile(Parameters.Paths.directory_mat_folder, BCKGRND_position);
        
    end
    
end
end


%% MASK FUNCTION -USED FOR CELL CLUSTERS EXPERIMENTS

function [mask,areamask] = mask_function(Image,crop)
%Input: image of interest (Image), depending on the channel used to compute the
%mask, and the crop size n*m (crop).
% Outputs: n*m mask matrix  (mask) and associated area ( areamask).

Image_temp1=imread(Image);
if (crop~=0)
    Image_temp2=imcrop(Image_temp1,crop);
else
    Image_temp2 = Image_temp1;
end
Img = imadjust(Image_temp2);

%CUSTOM PARAMETER
filter_size = 3;
morphology_size =12;
% filter_size = 5;
% morphology_size =25;
%---------------

% filtering functions:
BW1 = uint16(Img);
% I = medfilt2(BW1,[filter_size filter_size]);
I = BW1;
I2 = adapthisteq(I);
BW2 = imbinarize(I2);
BW3 = imdilate(BW2, strel('disk',morphology_size));
BW4 = imfill(BW3,'holes');
BWfinal = imerode(BW4,strel('disk',morphology_size));

% BW4 = BW2;
% BW = imerode(BW4,strel('disk',filter_size));
% BWfinal = imdilate(BW, strel('disk',morphology_size));

if (numel(find(BWfinal))/numel(BWfinal)>=.95)
    BWfinal = ones(size(Image));
end

mask = BWfinal;

areamask=sum(sum(mask));

end
%% VIDEO GENERATOR
function [] = Generate_Video(Parameters,Output)
global CELLS

n_timeframes = Parameters.n_timeframes;
profile = 'MPEG-4';
video_name = Parameters.Video.video_name ;
time_lapse = VideoWriter(strcat(video_name,'_field_',num2str(Parameters.position_of_interest)),profile);
time_lapse.FrameRate = 2;
time_lapse.Quality = 100;
open(time_lapse);%open file for writing video data
fig1 = figure;
set(fig1, 'Position', [5, 5, 1500, 800])

n_CHANNELS = Parameters.fluoeval.n_Channels;
n_timeframes = Parameters.n_timeframes;
selected_BCKGfield = Parameters.fluoeval.position_of_background; %if 0, background as antimask
tracking_flag = Parameters.tracking_flag;

for i=1:n_timeframes

    [cropped_image,fileNames, fileNamesBCKGRND] = fileNames_addressFunction(i, Parameters);
    
    %% FLUO EVAL
    
    for j= 1:n_CHANNELS
        Parameters.fluoeval.CHANNEL{j} = fileNames{j};
        if (selected_BCKGfield~=0)
            Parameters.fluoeval.BCKGRND_CHANNEL{j} = fileNamesBCKGRND{j};
        end
    end
    if (Parameters.crop_flag)
        Ph_image = imcrop(imread(Parameters.fluoeval.CHANNEL{1}),Parameters.fluoeval.crop);
        Im_Green = imcrop(imread(Parameters.fluoeval.CHANNEL{3}),Parameters.fluoeval.crop);
        BackGround = imcrop(imread(Parameters.fluoeval.BCKGRND_CHANNEL{3}),Parameters.fluoeval.CropRectBCKG);
    else
        Ph_image = imread(Parameters.fluoeval.CHANNEL{1});
        Im_Green = (imread(Parameters.fluoeval.CHANNEL{3}));
        BackGround = (imread(Parameters.fluoeval.BCKGRND_CHANNEL{3}));
    end

    
    %% Video
    if tracking_flag==1
        CENTERS =[];
        for a=1:length(CELLS{1,1})
            CENTERS = [CENTERS;CELLS{1,1}(a).center(i,1) CELLS{1,1}(a).center(i,2)];
        end
        label = [];
        label = [label; CELLS{1,1}.label];
        
        figure(fig1)
        set(fig1, 'Position', [5, 5, 1500, 800])
        
        %% Plots of the cropped images
        
        title(strcat('Frame n. ',num2str(i-1)));
        
        subplot(2,3,3);imshow(imadjust(imread(Parameters.fluoeval.CHANNEL{Parameters.fluoeval.nucleo_tag_channel})),[]);title(strcat('Nucleo Tag Frame n. ',num2str(i-1)));
        subplot(2,3,2);imshow(imadjust(imread(Parameters.fluoeval.CHANNEL{Parameters.fluoeval.fluorescence_channel})));title(strcat('mCherry Fluo Frame n. ',num2str(i-1)));
        subplot(2,3,1);
        RGB = insertText(imadjust(Ph_image),CENTERS,label,'FontSize',30);
        imshow(RGB),title(strcat('Phase Contrast Frame n. ',num2str(i-1)));
        title(strcat('Phase Contrast Frame n. ',num2str(i-1)));
        
        subplot(2,1,2);
        
        plot((0:i-1)*60,Output.Fluorescence_Channel_A_A (1:i,Parameters.fluoeval.fluorescence_channel),'color','k','linewidth',1.5);
        hold on;
        xlim([0 (n_timeframes-1)*60]);
        title(strcat('Frame n. ',num2str(i-1)));
        ylabel('Fluorescence [A.u.]','fontsize',16);
        ylim([min(Output.Fluorescence_Channel_A_A (:,Parameters.fluoeval.fluorescence_channel)) max(Output.Fluorescence_Channel_A_A (:,Parameters.fluoeval.fluorescence_channel))]);
        xlabel('time [min]');
        legend( 'Cell < 5 > fluorescence' , 'fontsize',12);
        hold off;
    else
        title(strcat('Frame n. ',num2str(i-1)));
        
        subplot(2,4,1);imshow( cropped_image{1},[]);title(strcat('Phase Contrast Frame n. ',num2str(i-1)));
        subplot(2,4,2);imshow(imadjust(cropped_image{Parameters.fluoeval.fluorescence_channel }));title(strcat('Green Fluo Frame n.',num2str(i-1)));
        subplot(2,4,4); imshow((Output.MASK{i})); title(strcat('Mask Frame n.', num2str(i-1)));
        subplot(2,4,3); imshow(imadjust(cropped_image{Parameters.fluoeval.blue_dye_channel})); title(strcat('Dye Fluo Frame n.', num2str(i-1)));%.*uint16(mask)
        
        subplot(2,1,2);
        
        plot((0:i-1)*60,Output.Fluorescence_Channel_A_A (1:i,Parameters.fluoeval.fluorescence_channel),'color','k','linewidth',1.5);
        hold on;
        xlim([0 (n_timeframes-1)*60]);
%         ylim([20 80])
        ylim([min(Output.Fluorescence_Channel_A_A (:,Parameters.fluoeval.fluorescence_channel)) max(Output.Fluorescence_Channel_A_A (:,Parameters.fluoeval.fluorescence_channel))]);
        title(strcat('Frame n. ',num2str(i-1)));
        ylabel('Fluorescence [A.u.]','fontsize',16);
        xlabel('time [min]');
        hold off;
    end
    
    
    pause(.2);
    frame = getframe(gcf);
    writeVideo(time_lapse,frame);
    pause(.1);
    writeVideo(time_lapse,getframe(fig1));
    
end
close(time_lapse);
end



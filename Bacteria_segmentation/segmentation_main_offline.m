%% This script calculates the fluorescence of single cells and the average fluorescence of population  
%This algorithm works using phase contrast and fluorescence images.
%The phase contrast image is used to calculate a binary mask of single
%cells, and the fluorescence image is used to calculate the fluorescence of
%the mask.
%The user will obtain a video with fluorescence quantification and the final mask.
%A mat file is created to save the variables and quantifications.

%Channel list for bacteria experiment:
% C00 = Phase Contrast
% C01 = GFP

%% Clear previous data
close all;
clear all;
clc;

%% Initial conditions
calculate_fluorescence = 1; % If the user wants to compute the fluorescence 
plot_results =1; % If the user wants to plot results
make_movie = 1; % If the user wants to generate a movie 

%% This is the path where the raw images are located
path='/Volumes/MicaliLab/Constantinos/Instrument data/Microscope DMi8/Experiment data/20240307/Tiff/';

%% Load images 
%Position001 and *C00 are the strings associated to the name of images the
% algorithm will analyse and .tif is the format of image acquired 
%from the microscope. The extension should be changed if needed.

% %%%Loading all the phase contrast data
PhContrast = dir(strcat(path,'20240307_Prot-Aux+Aux-Aux_Cm4uM_Pos1_1-281_p.tiff'));
PhContrast = table2struct(sortrows(struct2table(PhContrast),{'name'},{'ascend'}));

%Loading the green field data 
GreenField = dir(strcat(path,'20240307_Prot-Aux+Aux-Aux_Cm4uM_Pos1_1-281_g.tiff'));
GreenField = table2struct(sortrows(struct2table(GreenField),{'name'},{'ascend'}));

%% Cropping the zone of the image we want to segment
%This crop is done manually by the user
disp('draw the crop rectangle to select cells in the frame...');
[~,rect] = imcrop(imadjust(imread(strcat(path,PhContrast(1).name))));

%% Cropping the zone of the image we want to use to measure background fluorescence
%This crop is done manually by the user
disp('draw the crop rectangle to measure the dye in the frame...');
[~,rectInput] = imcrop(imadjust(imread(strcat(path,PhContrast(1).name))));
close all; %close the opened images 

%% Setting up the number of images we want to segment
init_frame = 1;                %setting the number of the initial fram
end_frame = length(PhContrast); 

%%Our control experiments have two phases called preinitialisation and
%%initialisation. These phases are used to calculate the minimum and maximum
%%fluorescence that will be used to normalize fluorescence values between 0 and 1.

%%Max_frame_end and min_frame_end define the range of images used to set
%%the normalization of fluorescence amplitude 
max_frame_end=120;   %end frame of the initialization (Max value of fluo)
min_frame_end=60;    %end frame of the preinitialization (Min value of fluo)

%%Number of images included in the initialization and preinitialization
%%phases
numbPreInitpics=60;  
numbinitpics=60;

%Initialization of variables
cell_number = zeros(1,end_frame);
elapsed_time = zeros(1,end_frame);
cell_masks = struct();
Sampling_time=5;

%% Segmentation 
%%This routine is carried out for every phase contrast image to calculate the binary mask. 
%%This function is an independent code 'segmentation_GF.m'
for i = init_frame:end_frame
disp('analysing picture: ');                                  
disp(num2str(i));
image = imcrop(imread(strcat(path,PhContrast(i).name)),rect); %Crop the image 
[elapsed_time(i),cell_masks(i).mask,cell_number(i),~] = segmentation_GF(image,2,0); %Obtain the masks and cell numbers (Segmentation code)
clc;% clear the screen
end

%% Threshold Evaluation
%%A threshold value was set to include only cells whose fluorescence was within 
%%a maximum range and to discard fluorescence from outlier cells. The threshold 
%%was calculated as the fluorescence value of the cells exhibiting the highest 
%%fluorescence in a set of images (60 in this example) acquired before the control experiment.

if calculate_fluorescence %If the user wants to compute the fluorecence
    
    %Variables initialization
    cell_fluorescenceInit = struct();
    Obj_Fluo  = [];
    background_fluorescence = zeros(1,end_frame);
    
    %%Loop to calculate the fluorescence to set a fluorescence threshold 
    for i = init_frame:end_frame %for each time step
        fprintf('Evaluating fluorecence of image at time %.2f of %.2f (Threshold) \n',(i)*Sampling_time,(end_frame-1)*Sampling_time);
        image = imcrop(imread(strcat(path,GreenField(i).name)),rect); %Cropping the image for the cells 
        ImageInput = imcrop(imread(strcat(path,GreenField(i).name)),rectInput); %Cropping the image  for the background 
        [cell_fluorescenceInit(i).val,background_fluorescence(i), cell_fluorescenceInit(i).average, cell_fluorescenceInit(i).std_dev,L] = fluorescence_eval_GF_Init(image,ImageInput,cell_masks(i).mask); %Computing the cell fluorecence with respect to the background
        %Threshold variables
        a=cell2mat({L(1:length(L)).MeanIntensity}); 
        b = max(a);
        Obj_Fluo  = [Obj_Fluo;b];
        clc;
    end
    
end

av_fluoInit=cell2mat({cell_fluorescenceInit(1:length(cell_fluorescenceInit)).average});

%% [numbPreInitpics+1:numbPreInitpics+numbinitpics] define the set of images used to calculated the threshold. 
th = max(Obj_Fluo(numbPreInitpics+1:numbPreInitpics+numbinitpics));  %Evalutating the threshold (cells brighter than this value are not taken into account in the evaluation of the fluo)

%% Fluorescence evaluation
if calculate_fluorescence          % If the user wants to compute the fluorescence
    %Initialization of the variables
    cell_fluorescence = struct();
    cell_fluorescence1 = struct();
    background_fluorescence = zeros(1,end_frame);
    fluorescence = zeros(1,end_frame);
       
    for i = 1:end_frame
        fprintf('Evaluating fluorecence of image at time %.2f of %.2f (Cells Fluo) \n',(i)*Sampling_time,(end_frame-1)*Sampling_time);
        image = imcrop(imread(strcat(path,GreenField(i).name)),rect); %Crop the image of the cells
        imageInput = imcrop(imread(strcat(path,GreenField(i).name)),rectInput); %Crop the image for the background
        [cell_fluorescence(i).val,background_fluorescence(i), cell_fluorescence(i).average,cell_fluorescence(i).std_dev,L,dObjs,fluorescence(i),cell_fluorescence(i).median] = fluorescence_eval_GF(image,imageInput,cell_masks(i).mask,th); %Compute the fluo of all cells with respect to the background
        clc;
    end
end
  av_fluo=cell2mat({cell_fluorescence(1:length(cell_fluorescence)).average}); %Average fluorescence of all cells
  stdav_fluo=cell2mat({cell_fluorescence(1:length(cell_fluorescence)).std_dev}); %Standard deviation of the fluorescence
  med_fluo=cell2mat({cell_fluorescence(1:length(cell_fluorescence)).median}); %Median value of the fluorescence
  
%% Making a video to visualise the fluorescence output in the microfluidic device
if make_movie&&calculate_fluorescence   
    
    %Starting video recording
    exp_name = datestr(PhContrast(1).datenum,'ddmmyyyy');
    time_lapse = VideoWriter(strcat('movie_exp',exp_name,'bacteria segmentation.avi'));
    time_lapse.FrameRate = 10;
    time_lapse.Quality = 100;
    open(time_lapse);
    fig1 = figure;
    
    max_fluo=mean(av_fluo(max_frame_end-6:max_frame_end));
    min_fluo=mean(av_fluo(min_frame_end-6:min_frame_end));
    av_fluo_norm=(av_fluo-min_fluo)./(max_fluo-min_fluo);
   
    for i = init_frame:end_frame -1
            PCimage = imcrop(imread(strcat(path,PhContrast(i).name)),rect);
            GFimage = imcrop(imread(strcat(path,GreenField(i).name)),rect);
            GFimage1 = imcrop(imread(strcat(path,GreenField(i).name)),rect);
            image = imcrop(imread(strcat(path,GreenField(i).name)),rect);
            PCimage = interp2(single(PCimage),1,'cubic');
            PCimage = uint16(PCimage);
            PCimage = reshapeHist(PCimage);
            GFimage2show = reshapeHist(GFimage1);
            GFimage = reshapeHist(GFimage).*uint16(cell_masks(i).mask);
            edges = zeros(size(PCimage));
            M = cell_masks(i).mask;
            M = interp2(single(M),1,'nearest');
            M = logical(M);
            edges(bwperim(M)) = 2^16-1;
            PC_rgb = cat(3, PCimage, PCimage, PCimage);
            GF_rgb = cat(3, zeros(size(GFimage)),GFimage2show,zeros(size(GFimage)));
            edges_rgb = cat(3, zeros(size(edges)),zeros(size(edges)),edges);
                         
            figure(fig1)
            set(fig1, 'Position',  [5, 5, 1500, 800]);
            str = strcat('exp - ','TEST',', time: ',num2str((i)*5),' min.');
            subplot(3,3,1);imshow(PCimage,[]);title(strcat('Phase Contrast Frame n. ',num2str(i-1)));
            subplot(3,3,2);imshow(GFimage);title(strcat('Green Fluo Frame n. ',num2str(i-1)))
            M = cell_masks(i).mask;
            M = logical(M);
            subplot(3,3,3);imshow(M);title(strcat('Mask Frame n. ',num2str(i-1)))
            
            subplot(3,1,2);
            title(strcat('Frame n. ',num2str(i-1)));
            hold on;
            plot((init_frame:i)*Sampling_time,av_fluo(init_frame:i),'color','k','linewidth',1.5);
            xlim([init_frame end_frame]*Sampling_time);
            ylim([min(av_fluo) max(av_fluo)]);
            ylabel('Fluorescence [A.u.]','fontsize',15);
            xlabel('Time [minutes]');
           
            
            subplot(3,1,3)
            hold on;
            plot((init_frame:i)*Sampling_time,cell_number(init_frame:i),'k','linewidth',2);
            xlim([init_frame end_frame]*Sampling_time);
            ylabel('Cell count','fontsize',15);
            xlabel('Time [minutes]');
               
    hold off;
    pause(.2);
    frame = getframe(gcf);
    writeVideo(time_lapse,frame);        
    end
    close(time_lapse);
end

close all;
close(time_lapse);
% filename = strcat('exp_',exp_name,'bacteria_segmentation.mat');
% save(filename);

function [] = Generate_Video(Parameters,Output)

n_timeframes = Parameters.n_timeframes;
profile = 'MPEG-4';
video_name = Parameters.Video.video_name ;
time_lapse = VideoWriter(strcat(video_name,'_field_',num2str(Parameters.position_of_interest)),profile);
time_lapse.FrameRate = 2;
time_lapse.Quality = 100;
open(time_lapse);%open file for writing video data
fig1 = figure;
set(fig1, 'Position', [5, 5, 1500, 800])

for i=1:n_timeframes

    [~,fileNames, fileNamesBCKGRND] = fileNames_addressFunction(i, Parameters);
    
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
    
    pause(.2);
    frame = getframe(gcf);
    writeVideo(time_lapse,frame);
    pause(.1);
    writeVideo(time_lapse,getframe(fig1));
    
end
close(time_lapse);
end

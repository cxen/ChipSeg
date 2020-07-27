function [Output] = fluorescence_evaluation(Input_data,BW_MASK)

%% reading images, also background ones if required, or building background
for i=1:Input_data.n_Channels
   Frame_from_Channel{i} = imread(Input_data.CHANNEL{i});
   if  (Input_data.position_of_background~=0)
       BCKGRND_Frame_from_Channel{i} = imread(Input_data.BCKGRND_CHANNEL{i});
   end
end

%% cropping images, if required
if (Input_data.crop~=0)
    for i=1:Input_data.n_Channels
        Frame_from_Channel{i} = imcrop(Frame_from_Channel{i},Input_data.crop);
        if  (Input_data.position_of_background~=0)
            BCKGRND_Frame_from_Channel{i} = imcrop(BCKGRND_Frame_from_Channel{i},Input_data.CropRectBCKG);
        else
            BCKGRND_from_antimask{i}=Frame_from_Channel{i}.*uint16(~BW_MASK);
        end
    end    
end

%% application of mask 
for i=1:Input_data.n_Channels
   Masked_Frame_from_Channel{i} = Frame_from_Channel{i}.*uint16(BW_MASK);
   
end

%% BACKGROUNDS A SEPARATE FIELD: AVERAGES ON CROPPED IMGS
for i=1:Input_data.n_Channels
    if (Input_data.position_of_background~=0)
        AVG_measure_from_BCKGRND_Frame_Channel(i) = mean(mean(BCKGRND_Frame_from_Channel{i}));

    else
        AVG_measure_from_antimask_BCKGRND_Channel(i) = mean(BCKGRND_from_antimask{i}(BCKGRND_from_antimask{i}~=0));

    end
end

%% FLUO AVERAGES ON MASKED IMGS 
for i=1:Input_data.n_Channels
  AVG_measure_from_MaskedFrame_Channel(i) = mean(Masked_Frame_from_Channel{i}(Masked_Frame_from_Channel{i}~=0));
end

%%  FLUO EVALUATION as AVG - AVG

for i=1:Input_data.n_Channels    
    if (Input_data.position_of_background~=0)
        Fluorescence_Channel_A_A(i)  = AVG_measure_from_MaskedFrame_Channel(i) -  AVG_measure_from_BCKGRND_Frame_Channel(i); 
    else
        Fluorescence_Channel_A_A(i)  = AVG_measure_from_MaskedFrame_Channel(i) -  AVG_measure_from_antimask_BCKGRND_Channel(i);
    end
end

%% OUTPUT FLUO AS NORMALIZED FLUO
if (Input_data.IRFP_flag~=0)
    for i=1:Input_data.n_Channels
        Nuclei_Normalized_Fluorescence_Channel_A_A(i)  = Fluorescence_Channel_A_A(i)/Fluorescence_Channel_A_A(Input_data.nucleo_tag_channel);
    end   
    Output.Nuclei_Normalized_Fluorescence_Channel_A_A = Nuclei_Normalized_Fluorescence_Channel_A_A;

end

%% Outputs
Output.Absolute_Fluorescence_AVG = AVG_measure_from_MaskedFrame_Channel; %absolute fluo as avg of pixels
Output.Fluorescence_Channel_A_A = Fluorescence_Channel_A_A; %avg of masked fluo frame - avg background frame
if(Input_data.position_of_background~=0)
    Output.Average_Background_fluorescence = AVG_measure_from_BCKGRND_Frame_Channel;
else
    Output.Average_Background_fluorescence = AVG_measure_from_antimask_BCKGRND_Channel;
end
end




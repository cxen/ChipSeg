%%This function is called to calculate the fluorescence to set the threshold.
%%Input: image,ImageInput,mask
%%Output: cell_fluorescence,background_fluorescence, fluorescence average, fluorescence std_dev, and
%%L that contains the fluorescence of every element in mask input 

function [cell_fluorescence,background_fluorescence, average, std_dev,L] = fluorescence_eval_GF_Init(image,ImageInput,mask)

background_fluorescence =  mean(mean(ImageInput(ImageInput~=0)));

CC = bwconncomp(mask,4);
L = regionprops(CC,image,'pixelidxlist','MeanIntensity');
cell_fluorescence = zeros(1,length(L));
    for index=1:length(L)
       FilledImage = ~ones(size(mask)); 
       FilledImage(L(index).PixelIdxList) = 1;
       cell_fluorescence(index) =  sum(sum(image.*uint16(FilledImage)))/(sum(sum(FilledImage)));
    end
cell_fluorescence = cell_fluorescence - background_fluorescence;
cell_fluorescence(cell_fluorescence<0) = 0;
average = mean(cell_fluorescence);
std_dev = std(cell_fluorescence);
return
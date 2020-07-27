% This function is called to calculate the fluorescence for each image.
% Input: image,ImageInput,mask
% Output: cell_fluorescence,background_fluorescence, fluorescence average, fluorescence std_dev,
% L contains the fluorescence of every element in mask input 

function [cell_fluorescence,background_fluorescence, average,std_dev,L,deletedObjs,fluorescence,med] = fluorescence_eval_GF(image,ImageInput,mask,th)

fluorescence = sum(sum(image.*uint16(~mask)))/(sum(sum(~mask)));
background_fluorescence =  mean(mean(ImageInput(ImageInput~=0)));

CC = bwconncomp(mask,4);
L = regionprops(CC,image,'pixelidxlist','MeanIntensity','Centroid','Area');
idx = find([L.MeanIntensity]<=th);
idy = find([L.MeanIntensity]>th);
LLL =L(idy);

LL=L(idx);
%Use L instead of LL to debug the fluorescence levels over the threshold 
deletedObjs = length(L)-length(LL);
cell_fluorescence = zeros(1,length(LL));
cell_fluorescence1 = zeros(1,length(LL));

    for index=1:length(LL)
       FilledImage = ~ones(size(mask)); 
       FilledImage(LL(index).PixelIdxList) = 1;
       cell_fluorescence(index) =  sum(sum(image.*uint16(FilledImage)))/(sum(sum(FilledImage)));
    end

cell_fluorescence = cell_fluorescence - background_fluorescence;

cell_fluorescence(cell_fluorescence<0) = 0;
average = mean(cell_fluorescence);
med = median(cell_fluorescence);
std_dev = std(cell_fluorescence);
return
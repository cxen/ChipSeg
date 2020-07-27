%% Adjust the image to its max and min brightness 
% The histogram of img2 is stretched to the min and max
% intensity value of img1. %This function uses imadjust. 
function img2 = reshapeHist(img1)
y = double(reshape(img1,[1,numel(img1)]));
low_in = min(y)/(2^16-1);
high_in = max(y)/(2^16-1);
img2 = imadjust(img1,[low_in high_in],[]);
end
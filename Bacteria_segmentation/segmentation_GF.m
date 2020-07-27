%% segmentation_GF calculates the fluorescence at every image.
% Input: image (Phase contrast image), K, DEBUG (1-activated or 0-inactivated (default value))  
% Output: elapse_time, BW_final (binary mask), cell_number and
% L contains the fluorescence of every element 

function [elapsed_time,BWfinal,cell_number,L] = segmentation_GF(image,K,DEBUG)
tic;
%% Custom parameters: K, filter_size, morphology_size and l_threshold_size 
%The user will define these parameters manually. Their value should be 
%tuned and it depends on the image quality.
filter_size = 3;
morphology_size = 20;
l_threshold_size = 15;

%% interpolation
old_image = image;
image = interp2(single(image),K,'cubic');
image = uint16(image);


%% macroBlobs performs global thresholding of Otsu algorithm.
%The output will contain the area in which cells are located
%I is the output image (filtered) and BW global is a binary mask 
[BWglobal,I] = macroBlobs(image,filter_size,morphology_size);

%% localThresholding performs global thresholding of Otsu algorithm.
%BWlocal is a binay mask which contains the individual objects (single cells) 
%This mask may contain small objects which are not cells and that
%come from fluorescence levels in the background.
BWlocal = localThresholding(I,BWglobal,filter_size,l_threshold_size);

%% segmentationRefinement separates bigger connected objects at convex points 
%BWrefined is the resulting mask after the connected elements (BWbig) are separated and
%small objects are removed
[BWrefined,BWbig] = segmentationRefinement(BWlocal,K);

%% sampling the selection mask
%BWsampled is the resized mask after the refinement process 
BWsampled = BWrefined(1:2^K:end,1:2^K:end);

%% smallObjectsRemoval will delete the small objects around the cells  
%BWfinal is the final mask
BWfinal = smallObjectsRemoval(BWsampled);

CC = bwconncomp(BWfinal,4);
cell_number = CC.NumObjects;
L = regionprops(CC,'centroid','convexhull','PixelIdxList');

elapsed_time = toc; %The tic and toc functions work together to measure elapsed time. 

%% DEBUG is active when we put 1 in the function segmentation_GF(image,K,DEBUG) in the main code
% fig 1 to fig 8 shows the debugging of the segmentation process 

if DEBUG
    fig1 = figure;
    figure(fig1)
    title('Binarization after global thresholding');
    hold on;
    imshowpair(image,BWglobal,'montage')
    
    fig2 = figure;
    figure(fig2)
    title('Binarization after local thresholding');
    hold on;
    imshowpair(image,BWlocal,'montage')
    
    fig3 = figure;
    figure(fig3)
    title('Global vs local thresholding');
    hold on;
    imshowpair(BWglobal,BWlocal);
    
    fig4 = figure;
    figure(fig4)
    title('Binarization refined');
    hold on;
    imshowpair(image,BWrefined,'montage')
    
    fig5 = figure;
    figure(fig5)
    title('Big objects being segmented');
    hold on;
    imshow(BWbig);    
   
    fig6 = figure;
    figure(fig6)
    title('Refinement vs local thresholding');
    hold on;
    imshowpair(BWlocal,BWrefined);
   
    fig6 = figure;
    figure(fig6)
    title('Segmentation result');
    hold on;
    imshowpair(old_image,old_image.*uint16(BWfinal),'montage');
    
    fig8=figure;
    figure(fig8);
    title('Segmentation result');
    hold on;
    imshow(old_image,[],'InitialMagnification','fit');

    for i=1:length(L)
    
    figure(fig8)
    hold on;
     plot(L(i).Centroid(1),L(i).Centroid(2),'y*')
     plot(L(i).ConvexHull(:,1),L(i).ConvexHull(:,2),'r-','linewidth',0.5)
    end
    elapsed_time = toc;
end

return

function [BWout,I] = macroBlobs(image,filter_size,morphology_size)
% %----median filter
I = wiener2(image,[filter_size^2 filter_size^2]);
I = medfilt2(I,[filter_size filter_size]);
% %----contrast enhancement
I2 = adapthisteq(I);
% %----Otsu binarization
level = graythresh(I2);
BW = imbinarize(I2,level);
% %----cell areas selection
BW1 = imdilate(BW, strel('disk',morphology_size));
BW3 = imfill(BW1,'holes');
% BWout contains the macro blobs where cells are located
BWout = imerode(BW3,strel('disk',morphology_size)); %(optional)
if (numel(find(BWout))/numel(BWout)>=.95)
    BWout = ones(size(image));
end
return

function BWout = localThresholding(I,BW,filter_size,l_threshold_size)
CC = bwconncomp(BW);
L = regionprops(CC,'boundingbox','area');
% BW7 is composed of all the local areas analysed locally (cells in white)
BW7 = ~(ones(size(I)));
for i=1:length(L)
    rect = ceil(L(i).BoundingBox);
    temp = imcrop(I,rect);%cropping the filtered image (for improved results)
    tBW = imcrop(BW,rect);
    tempBW = ~niblack(temp, [l_threshold_size l_threshold_size]);
    tempBW(~tBW) = 0;
    if (rect(3)~=size(tBW,2))
        rect(3)=size(tBW,2);
    end
    if (rect(4)~=size(tBW,1))
        rect(4)=size(tBW,1);
    end
          BW7(rect(2):rect(2)+rect(4)-1,rect(1):rect(1)+rect(3)-1) = tempBW;
end
BW8 = bwareafilt(BW7,[10 Inf],4);
% Mask BW8 gaussian filtering -> background pixels have a higher intensity value
% than those belonging to cells
BW9 = reshape(BW8,[1 numel(BW8)]);
pic = I.*uint16(BW8);
pic = reshape(pic,[1 numel(pic)]);
idx = pic~=0;
[y,x] = hist(double(pic(idx)),100);
f = fit(x.',y.','gauss1');
cut_off = f.b1+f.c1;
idx2 = pic>cut_off;
BW9(idx2) = 0;
BW9 = reshape(BW9,size(BW8));
% smoothing BW9's edges
BW9 = imfill(BW9,'holes');
BWout = imopen(BW9,strel('disk',3));
return

function [BWout,BW10] = segmentationRefinement(BW,K)
CC2 = bwconncomp(BW,4);
L2 = regionprops(CC2,'area','centroid');
cell_area = zeros(1,length(L2));
    for i=1:length(L2)
        cell_area(i) = L2(i).Area;
    end
% % poisson fitting
lambda = poissfit(cell_area);
BW10 = bwareafilt(BW,[2.*lambda Inf],4);
% BW10 contains all the connetted elements that need to be further
% segmented
BWout = (BW & ~BW10) | BW10&~bwperim(BW10);
max_iter = 5;
iter = 0;
BW12 = bwareafilt(BWout,[2.*lambda Inf],4);
    while (sum(sum(BW12))>0&&iter < max_iter)
        iter = iter+1;
        BWout = (BWout & ~BW12) | convex_points(BW12,K);
        BW12 = bwareafilt(BWout,[2.*lambda Inf],4);
    end
return

function BWout = smallObjectsRemoval(BW)
CC = bwconncomp(BW,4);
L = regionprops(CC,'area');
cell_area = zeros(1,length(L));
    for i=1:length(L)
        cell_area(i) = L(i).Area;
    end
lambda = poissfit(cell_area); 
BWout = bwareafilt(BW,[.5*lambda Inf],4);
return

%%This function calculates convex points of the final mask
function bw = convex_points(BW,K)
clc;
bw = ~ones(size(BW));
L2 = regionprops(bwconncomp(BW),'filledimage','conveximage','pixelidxlist','boundingbox');
for index=1:length(L2)
FilledImage = ~ones(size(BW)); ConvexImage = ~ones(size(BW));
FilledImage(L2(index).PixelIdxList) = 1;
rect = floor(L2(index).BoundingBox);
ConvexImage(rect(2):rect(2)+rect(4)-1,rect(1):rect(1)+rect(3)-1) = L2(index).ConvexImage;
%isolating the areas between the convexhull and the cells;
ProcessImage = ConvexImage&~FilledImage;
ProcessImage = bwareafilt(ProcessImage,[4 Inf],8);
cc2 = bwconncomp(ProcessImage,8);
l2 = regionprops(cc2,'all');
if length(l2)>1
    regions_area = zeros(1,length(l2));
    for i=1:length(l2)
        regions_area(i) = l2(i).Area;
    end
    [~,idx] = max(regions_area);

%isolating the edge points of all the selected areas
    perims = struc([]);
    for i=1:length(l2)
        temp = ~ones(size(ProcessImage));
        temp(l2(i).PixelIdxList) = 1;
        [perims(i).y,perims(i).x] = find(bwperim(temp));
        perims(i).perim = bwperim(temp);
    end

    closest_pts.x1 =0;closest_pts.y1 =0;closest_pts.x2 =0;closest_pts.y2 =0;closest_pts.dist =Inf;
    %1st set of points
    A = [perims(idx).x,perims(idx).y];
    indices = 1:length(l2);
    indices(indices==idx) = [];
    for j=1:length(indices)
        %2nd set of points
        B = [perims(indices(j)).x,perims(indices(j)).y];
        D2 = pdist2(A,B);
        dist = min(min(D2));
            if dist<closest_pts.dist
               closest_pts.dist = dist;
               [r,c] = find(D2==dist);
               closest_pts.x1 = A(r(1),1);
               closest_pts.x2 = B(c(1),1);
               closest_pts.y1 = A(r(1),2);
               closest_pts.y2 = B(c(1),2);
            end
    end
    Line = Drawline(ProcessImage, closest_pts.y1,closest_pts.x1,closest_pts.y2, closest_pts.x2,K);
    bw = bw | (FilledImage&~Line);
else
    bw = bw | FilledImage;
end
end
return

%%This function draws a line to separate connected elements in BWbig
function img = Drawline(img, X0, Y0, X1, Y1,K)
img = ~ones(size(img));
min_boundary = 1;
max_boundary_x = size(img,2);
x0 = (X0-K:X0+K);x1 = (X1-K:X1+K);
x0(x0>max_boundary_x) = max_boundary_x;
x0(x0<min_boundary) = min_boundary;
x1(x1>max_boundary_x) = max_boundary_x;
x1(x1<min_boundary) = min_boundary;
for i = 1:min(length(x0),length(x1))
    for n = 0:(1/round(sqrt((x1(i)-x0(i))^2 + (Y1-Y0)^2))):1 
        xn = ceil(x0(i) +(x1(i) - x0(i))*n); 
        yn = ceil(Y0 +(Y1 - Y0)*n); 
        img(xn,yn) = 1; 
    end
end
return

function output = niblack(image, varargin)
% Initialization
numvarargs = length(varargin);          % only want 4 optional inputs at most
if numvarargs > 4
    error('myfuns:somefun2Alt:TooManyInputs', ...
     'Possible parameters are: (image, [m n], k, offset, padding)');
end
 
optargs = {[3 3] -0.2 0 'replicate'};   % set defaults
 
optargs(1:numvarargs) = varargin;        % use memorable variable names
[window, k, offset, padding] = optargs{:};

if ~ismatrix(image)
    error('The input image must be a two-dimensional array.');
end

% Convert to double
image = double(image);

% Mean value
mean = averagefilter(image, window, padding);

% Standard deviation
meanSquare = averagefilter(image.^2, window, padding);
deviation = (meanSquare - mean.^2).^0.5;

% Initialize the output
output = zeros(size(image));

% Niblack
output(image > mean + k * deviation - offset) = 1;
return

function image=averagefilter(image, varargin)
% Parameter checking.
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end
 
optargs = {[3 3] 0};            % set defaults for optional inputs
optargs(1:numvarargs) = varargin;
[window, padding] = optargs{:}; % use memorable variable names
m = window(1);
n = window(2);

if ~mod(m,2) 
    m = m-1; 
end       % check for even window sizes
if ~mod(n,2) 
    n = n-1; 
end

% Initialization.
[rows,columns] = size(image);   % size of the image

% Pad the image.
imageP  = padarray(image, [(m+1)/2 (n+1)/2], padding, 'pre');
imagePP = padarray(imageP, [(m-1)/2 (n-1)/2], padding, 'post');

% Always use double because uint8 would be too small.
imageD = double(imagePP);

% Matrix 't' is the sum of numbers on the left and above the current cell.
t = cumsum(cumsum(imageD),2);

% Calculate the mean values from the look up table 't'.
imageI = t(1+m:rows+m, 1+n:columns+n) + t(1:rows, 1:columns)...
    - t(1+m:rows+m, 1:columns) - t(1:rows, 1+n:columns+n);

imageI = imageI/(m*n);

% Return matrix in the original type class.
% image = cast(imageI, class(image));
image = cast(imageI, 'like', image);
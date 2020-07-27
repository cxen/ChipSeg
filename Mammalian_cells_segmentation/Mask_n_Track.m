
function [XYpos,Radii,orientation, num_el, totalMask, singleMASK]=Mask_n_Track(Image_for_mask, Im_H2B, crop, curr_time)
%read and crop
Image_for_mask=imread(Image_for_mask);
Im_H2B=imread(Im_H2B);

if (crop~=0)
    Image_for_mask=imcrop(Image_for_mask,crop);
    Im_H2B=imcrop(Im_H2B,crop);
end
% CUSTOM PARAMETERS
radius = 3;
nhood = ones(4,2);
lower_bound_avg_area_percent = 0.05;
Min_area = 80;
%-----------

% image filtering
Image = imadjust(Image_for_mask);
Ih2b = adapthisteq(Im_H2B);

bw_1 = imbinarize(Image);
bw_h2b = imbinarize(Ih2b, graythresh(Ih2b));
bw2_h2b = imclose(bw_h2b, nhood);


se = strel('disk',radius);
bw2 = imopen(bw_1, nhood);
bw2 = imerode(bw2, se);

bw3 = imfill(bw2,'holes');
bw3_h2b = imfill(bw2_h2b,'holes');

s_temp = regionprops(bw3,{'Area'});
AVGarea = 0;
for i=1:length(s_temp)
   AVGarea =  AVGarea + s_temp(i).Area;
end
AVGarea = AVGarea/length(s_temp);

bw4 = bwareaopen(bw3,floor(lower_bound_avg_area_percent*AVGarea));
bw4_h2b = bwareaopen(bw3_h2b,floor(Min_area));

[BWh,n_nuclei] =bwlabel(bw4_h2b,8);
[BW,num_el] =bwlabel(bw4,8);

s = regionprops(BWh,{'Area','Centroid','MajorAxisLength','MinorAxisLength','Orientation'});

XYpos = [];
Radii = [];
orientation = [];

for q = 1:length(s)
    XYpos = [XYpos; [s(q).Centroid(1) s(q).Centroid(2)]];   % centres positions
    Radii = [Radii; [s(q).MajorAxisLength s(q).MinorAxisLength]];    % radii lenghts
    orientation = [orientation; s(q).Orientation];   % orientations
end
% Mask function, needed for tracking
[totalMask,singleMASK] = MaskGenerator( BWh, XYpos, Radii, orientation);

% tracking function
tracking_fun(singleMASK, Radii, XYpos, orientation, curr_time, n_nuclei)

end

%%  Ellipse masks GENERATOR
% function that generates elliptic masks needed for tracking single cells
function [totalMask,singleMASK]=MaskGenerator( mask_rough, XYpos, Radii, orientation)
   
    ImSize=size(mask_rough);
    Image = imadjust(mask_rough);
    falseMASK=false(ImSize);  
    totalMask = falseMASK;
    
    singleMASK.pixel=[];
    
    t = linspace(0,2*pi,50);
    
    for indexOBJ = 1:length(Radii)
               
        a = Radii(indexOBJ,1)/2;%major axis
        b = Radii(indexOBJ,2)/2;%minor axis
        Xc = XYpos(indexOBJ,1);%x coordinate centre
        Yc = XYpos(indexOBJ,2);%y coordinate centre
        phi = deg2rad(-orientation(indexOBJ));
                
        plineX = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi); 
        plineY = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi); 
        
        M = poly2mask(plineX, plineY, ImSize(1), ImSize(2));  
     
        croppedIMG = imcrop(Image, [Xc-a Yc-a 2*a 2*a]);  

        BW = imfill(imdilate(imbinarize(croppedIMG),strel('disk',2)),'holes');   
        L = bwlabel(BW, 8); 
        S = regionprops(L, 'Area'); 
        [~, IND] = max(cat(1, S.Area)); 

        if ~isempty(IND)    
            
            startX = round(Yc-a);
            startY = round(Xc-a);
        
            countX = size(BW,1);   
            countY = size(BW,2);    
           
            if startX<=0
                startX = 1; 
            end
            
            if startY<=0
                startY = 1; 
            end
            
            M(startX:startX+countX-1,startY:startY+countY-1) = M(startX:startX+countX-1,startY:startY+countY-1) & BW;   
            M = imdilate(M,strel('disk',3));
            
            totalMask = totalMask | M;  
            singleMASK(indexOBJ).pixel = falseMASK | M; 
            
        else   
            singleMASK(indexOBJ).pixel = falseMASK; 
        end
        
    end
    
end

function [] = tracking_fun(singleMASK, Radii, XYpos, orientation, curr_time, num_el)

global CELLS

MAX_NUCLEI_DIMENSION = 50;

Cells = CELLS{1,1};
RemovedCELLS = CELLS{2,1};
idN = CELLS{3,1};

        if (num_el~=0)
            
            if curr_time==1              
                for k=1:num_el
                    POS = 1;  
                    Cells(k).label = idN;
                    idN = idN+1;
                    [Cells,idN] = SaveStructEllipse(Cells,Radii,XYpos,orientation,singleMASK(k).pixel,idN,curr_time,1,POS,k,k);
                 
                    for q=1:length(Cells)
                        if norm([Cells(q).center(end,1) Cells(q).center(end,2)] - [Cells(k).center(end,1) Cells(k).center(end,2)]) < MAX_NUCLEI_DIMENSION
                            Cells(k).nucleusLABEL = 1;
                            Cells(k).originalCELL = Cells(q).label;
                        end
                    end
                    
                end
                
            else              
                oldXYpos = extractOBJS(Cells);
                [m] = TraceOldNew(oldXYpos,XYpos);               
                [Cells,RemovedCELLS,idN] = updateOBJS( m, Cells, RemovedCELLS, XYpos, Radii, orientation, idN,curr_time,1,singleMASK);  

            end
            
        else
            
            RemovedCELLS = [RemovedCELLS Cells];
            
        end
           
    
    CELLS{1,1} = Cells;
    CELLS{2,1} = RemovedCELLS;
    CELLS{3,1} = idN;
end
    
    %%                             SAVE STRUCT
% -------------------------------------------------------------------------
% Function that saves all information in a structure
% -------------------------------------------------------------------------

function [Temp,idN]=SaveStructEllipse(Temp,radius,center,orientation,Mask,idN,frame,field,POS,k,c)  
    Temp(c).rad(POS,1)=radius(k,1);
    Temp(c).rad(POS,2)=radius(k,2);
    Temp(c).center(POS,1)=center(k,1);
    Temp(c).center(POS,2)=center(k,2);
    Temp(c).orientation(POS)=orientation(k);
    Temp(c).frame(POS)=frame;
    Temp(c).field=field;
    Temp(c).nucleusLABEL=0;  
    Temp(c).originalCELL=Temp(c).label;
    Temp(c).mask{POS}=Mask;

    
end

function [M] = TraceOldNew(oldP, newP)

N1 = size(oldP,1);

N2 = size(newP,1);

COST = zeros(N1,N2);

M = zeros(N1,N2);

for i = 1:N1    
    
    for j = 1:N2   
        
        COST(i,j) = norm(oldP(i,:) - newP(j,:));   % cost function: distance between two cells
        
    end
    
end
for i= 1:N1
    [c] = find(COST(i,:)==min(COST(i,:)));
    M(i,c) = 1;
end


M = uint16(M);  

end

%%                            UPDATE CELL OBJETS
% -------------------------------------------------------------------------

function [cells,removedCELLS,idN]=updateOBJS(m,cells,removedCELLS,XYpos,Radius,orientation,idN,frame,field,NUCLEIsingleMASK)

MAX_CIRCLE_DISPLACEMENT=50;
MAX_NUCLEI_DIMENSION=100;
M = size(m,1);
updatedINDEXES = [];
toREADD=[];

for k = 1:M
    INDEX = find(m(k,:) > 0); 
    
    if length(INDEX) == 1
        
        Z = INDEX(1);
        POS = size(cells(k).center,1)+1;   
        [cells,idN]=SaveStructEllipse(cells,Radius,XYpos,orientation,NUCLEIsingleMASK(Z).pixel,idN,frame,field,POS,Z,k);
        
        updatedINDEXES(end+1) = k;  
        
    elseif length(INDEX) > 1   
        
        NORMS = [];
        
        for h = 1:length(INDEX)
            NORMS(h) = norm([cells(k).center(end,1) cells(k).center(end,2)] - XYpos(INDEX(h),:));  
        end
        
        [~,IDS] = sort(NORMS); 
        sortedINDEXES = INDEX(IDS); 
        
        POS = size(cells(k).center,1)+1;   
        
        [cells,idN]=SaveStructEllipse(cells,Radius,XYpos,orientation,NUCLEIsingleMASK(sortedINDEXES(1)).pixel,idN,frame,field,POS,sortedINDEXES(1),k);
        
        updatedINDEXES(end+1) = k;
        
        for D = 2:length(sortedINDEXES) 
            
            S = length(cells);
            POS=1;
            
            for l=1:length(removedCELLS)
                
                DIST(l)= norm([removedCELLS(l).center(end,1) removedCELLS(l).center(end,2)] - XYpos(sortedINDEXES(D),:));
                
                if DIST(l)<MAX_CIRCLE_DISPLACEMENT
                    
                    cells(S+1)=removedCELLS(l);
                    toREADD=[toREADD l];
                    POS = size(cells(S+1).center,1)+1;
                    
                end
                
            end
            
            if S+1>length(cells)
                cells(S+1).label=idN;
                idN=idN+1;
            end
            
            [cells,idN]=SaveStructEllipse(cells,Radius,XYpos,orientation,NUCLEIsingleMASK(sortedINDEXES(D)).pixel,idN,frame,field,POS,sortedINDEXES(D),S+1);
            
            if norm([cells(k).center(end,1) cells(k).center(end,2)] - [cells(S+1).center(end,1) cells(S+1).center(end,2)]) < MAX_NUCLEI_DIMENSION
                cells(S+1).nucleusLABEL=1;
                cells(S+1).originalCELL=cells(k).label;
            end
            
            updatedINDEXES(end+1) = S+1;
            
        end
        
    end
    
end

removedCELLS(toREADD)=[];

toREMOVE = [];

for h = 1:length(cells)
    
    if isempty(find(updatedINDEXES==h,1))   
        
        toREMOVE = [toREMOVE h];   
        
    end
    
end

for h = 1:length(toREMOVE)
    
    removedCELLS = [removedCELLS cells(toREMOVE(h))]; 
    
end

cells(toREMOVE) = [];

end

%% Extracting previous object centroids
function OBJS = extractOBJS(trackedOBJS)

for q = 1:length(trackedOBJS)
    
    OBJS(q,1) = trackedOBJS(q).center(end,1);
    
    OBJS(q,2) = trackedOBJS(q).center(end,2);
    
end

end

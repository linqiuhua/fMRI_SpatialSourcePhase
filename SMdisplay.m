%*****************************************************************
% Description: show the saptial maps. 
%     The function is adapted from the GIFT toolbox
%     (http://mialab.mrn.org/software/gift/index.html) and some functions
%     and data involved are provided in the folder named "SMshow"
% Usage:
%   SMdisplay(ic,mask_ind,Postion,kk)
% Input: 
    % ic£ºthe spatial maps with dimension 1 * length(mask_ind)
	% mask_ind: the index for voxels to show
    % Postion: the slices to show e.g.,[89,70,85]
    % kk: the selection of colormap, the details are introduced in SMshow\icatb_getColormap2.m
% Output: 
    % figure for saptial maps
% Date: March 2017
% Author: Yue Qiu
%*****************************************************************
    
function SMdisplay(ic,mask_ind,Postion,kk)
cd SMshow
icaDIM = [181,217,181]; structDIM = [181,217,181];
S = zeros(53,63,46); S(mask_ind) = ic;
load niisetV.mat V 
sm = reshape(S,[53 63 46 1]);
icatb_write_nifti_data('2DSMshow.nii', V, sm);
load parameters_returnResizedImage_ch2betr.mat
compFiles = '2DSMshow.nii';
[structuralImage, icasig, coords, HInfo, text_left_right] = icatb_returnResizedImage(structFile, compFiles, ...
    anatomicalPlane, slices_in_mm, dataType, complexInfo, 1); %, mask_file
% icasig = SS(:,:,:,ii);
icasig = reshape(icasig, size(icasig, 1), prod(structDIM)); 
% apply display parameters
imageValues = 1;
threshValue = 0;
convertToZ = 0;

returnValue = imageValues ; 
[icasig] = icatb_applyDispParameters(icasig, convertToZ, returnValue, threshValue, ...
    structDIM, HInfo);
icasig = reshape(icasig, [size(icasig,1), HInfo.DIM(1), HInfo.DIM(2), HInfo.DIM(3)]);
tempIC = reshape(icasig, 1, icaDIM(1), icaDIM(2), icaDIM(3));

mm = max(max(max(max(tempIC(:,:,:,:)))));
nn = min(min(min(min(tempIC(:,:,:,:)))));
maxz=max(abs(mm),abs(nn));

%  [y,i,j,z]=ind2sub(size(tempIC),find(tempIC==mm));
%  Postion = [i,j,z]; 
[im, maxICAIM, minICAIM, minInterval, maxInterval] = icatb_overlayImages(tempIC, structuralImage, icaDIM, ...
    structDIM, imageValues,maxz);
% get unit color
unitColor = (icatb_range(im))/256;

% reshape imDIM
imDIM = structDIM;
im = squeeze(reshape(im, imDIM(1), imDIM(2), imDIM(3)));
% color maps
cm = icatb_getColormap2(kk,imageValues,1);
%cm(1:64,:)=[255*ones(64,1),55*ones(64,1),2*ones(64,1)]/255; %mask

% get max value for scaling
maxICA = maxz;%maxICAIM
minICA = -maxz;%minICAIM
if minICA >= 0
    icaCLIM = [minICA maxICA];
    % Special case when the maxICA is equal to zero
elseif maxICA == 10^-8 & minICA < 0
    icaCLIM = [minICA 0];
else
    tempMax = max([abs(maxICA) abs(minICA)]);
    icaCLIM = [-tempMax tempMax];
end

% reshape image so that it doesn't look "squashed"
im = reshape(im, structDIM(1), structDIM(2), structDIM(3));


% redraw callback

%-------------------------------------------------------------------
%CODE TO Draw figure
%-------------------------------------------------------------------

%get global defaults
icatb_defaults;
global COLORLIST;
% Screen Color Defaults
global BG_COLOR;
global BG2_COLOR;
global BUTTON_COLOR;
global BUTTON_FONT_COLOR;
global FONT_COLOR;
global AXES_COLOR;

% Fonts
global UI_FONTNAME;
global UI_FONTUNITS;
global UI_FS;
global LEGENDCOLOR;
global ASPECT_RATIO_FOR_SQUARE_FIGURE;
global DETRENDNUMBER;

%--load data and grab variables
% minICATC = data.minICATC;
%GraphicsHandle = data.GraphicsHandle;
realPos =  Postion;
icaPos = Postion;
% voxel origin
% voxelOrigin = data.voxelOrigin;
% unitColor = data.unitColor;
% im = data.im;
% fmriHInfo.DIM = structDIM; 
% fmriFiles = data.fmriFiles;


% numOfSub = size(fmriFiles,2);
% numOfSess = size(fmriFiles(1).ses,2);

% [realWorldPos, voxelCoord] = getCoord(fmriHInfo.V(1), real_world_coords, realPos);
% [realWorldPos, voxelCoord] = getCoord(fmriHInfo.V(1).mat, real_world_coords, pixelPos);
realWorldPos = [0 0 0]; voxelCoord =Postion;

icasig = squeeze(icasig);
% icaCLIM = data.icaCLIM;
voxelIntensity = icasig(realPos(1), realPos(2), realPos(3));


%display current position
%disp(['Current Voxel Position ', num2str(voxelCoord), ' [x y z]. Value is ', num2str(voxelIntensity)]);

% Print real world coordinates
%disp(['Real world coordinate is ', num2str(realWorldPos(1)), ' ', num2str(realWorldPos(2)), ' ', num2str(realWorldPos(3))]);

%--save image
xdim = size(im,1);
ydim = size(im,2);
zdim = size(im,3);

%--get max and min of image
maxIM = max(max(max(im)));
minIM = min(min(min(im)));

signal = im(icaPos(1),icaPos(2),icaPos(3));
diffInterval = icatb_range([maxInterval abs(minInterval)]);
diffInterval = diffInterval / 2;
if(diffInterval - signal <=0)
    voxelColor=minInterval+unitColor;
else
    voxelColor=maxInterval-unitColor;
end

% im(icaPos(1), icaPos(2),icaPos(3))=voxelColor;

sliceXY=reshape(im(:,:,icaPos(3)),size(im,1),size(im,2));
sliceXZ=reshape(im(:,icaPos(2),:),size(im,1),size(im,3));
sliceYZ=reshape(im(icaPos(1),:,:),size(im,2),size(im,3));


%scaledA from icasig;
% maxAbsValue = max([max(max(max(icasig))) abs(min(min(min(icasig))))]);
% A=A;
% % Scaling with positive number
% %scaledVal = A * (abs(icasig(icaPos(1),icaPos(2),icaPos(3)))/maxAbsValue);
% scaledVal = A * (icasig(icaPos(1),icaPos(2),icaPos(3))/maxAbsValue);


%figure out max of y-axis and min y-axis for plotting ica timecourse
% maxA=max(A);
% minA=min(A);

newColor = cm;
% figure handle
figLabel = '  ';
GraphicsHandle = icatb_getGraphics(figLabel, 'graphics', 'orthoviewer');
set(GraphicsHandle, 'units', 'normalized' ,'Toolbar', 'figure', 'menubar', 'figure');
optionsH = uimenu('parent', GraphicsHandle, 'label', 'Options');
optionsSubH = uimenu(optionsH, 'label', 'Set Voxel Position', 'callback', {@setVoxelPosCallback, GraphicsHandle});

%----------------------------------------------------
%Draw and scale image in Quad 1

%use dimensions of y and x axis to keep aspect ratio
% ie makes axis square within allocated space
% VOX = data.structHInfo.VOX;
axQuad1Pos = [.35 .40 .36 .20]; handles = 2;
dim = size(sliceXZ);
ImageAxis = axes( 'units', 'normalized', 'position', axQuad1Pos, 'color', [1 1 1]);
axis(ImageAxis, 'off');
% using image function as it is general function than imagesc
origData.axisHandle1 = image(rot90(sliceXZ), 'parent', ImageAxis, 'CDataMapping', 'scaled');
set(ImageAxis, 'clim', [minInterval 2*maxInterval]);
%origData.axisHandle1 = imagesc(rot90(sliceXZ), [minInterval 2*maxInterval]); % replace with image function
imageAxisPos = get(ImageAxis, 'position');
yAxisRatio = imageAxisPos(4)/imageAxisPos(3);
xAxisRatio = imageAxisPos(3)/imageAxisPos(4);
if(yAxisRatio>1)
    yAxisRatio = 1;
else
    xAxisRatio = 1;
end
R= ASPECT_RATIO_FOR_SQUARE_FIGURE;
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*yAxisRatio imageAxisPos(4)*xAxisRatio];
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*R(1) imageAxisPos(4)*R(2)];
set(ImageAxis, 'position', imageAxisPos);


%use actual dimensions of image to keep aspect ratio
imageAxisPos = get(ImageAxis,'position');
imageXRatio = dim(2)/dim(1);
imageYRatio = dim(1)/dim(2);
if(imageXRatio>1)
    imageXRatio = 1;
else
    imageYRatio = 1;
end
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*imageYRatio imageAxisPos(4)*imageXRatio];
set(ImageAxis, 'position', imageAxisPos);


%use voxel dimensions of image to keep aspect ratio
voxYAxisRatio = 1;
voxXAxisRatio = 1;
R=ASPECT_RATIO_FOR_SQUARE_FIGURE;
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*voxYAxisRatio imageAxisPos(4)*voxXAxisRatio];
set(ImageAxis, 'position', imageAxisPos);

%keeping the same aspect ratio as calculated above enlarge image to fit
% in to its allocated space
enlargeScalar = min([axQuad1Pos(3) axQuad1Pos(4)])/max([imageAxisPos(3) imageAxisPos(4)]);
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*enlargeScalar imageAxisPos(4)*enlargeScalar];
set(ImageAxis,'position',imageAxisPos);
colormap(newColor);
quad(1).xCord =[ imageAxisPos(1), imageAxisPos(1)+imageAxisPos(3),imageAxisPos(1)+imageAxisPos(3),imageAxisPos(1)];
quad(1).yCord =[ imageAxisPos(2)+imageAxisPos(4), imageAxisPos(2)+imageAxisPos(4),imageAxisPos(2),imageAxisPos(2)];
set(ImageAxis,'XTickLabel',[]);
set(ImageAxis,'YTickLabel',[]);

%--------------------------------------------------------
%Draw and scale image in Quad 2


%use dimensions of y and x axis to keep aspect ratio
% ie makes axis square within allocated space
axQuad2Pos = [.1 .40 .36 .20];
dim = size(sliceYZ);
ImageAxis = axes( 'units', 'normalized', 'position', axQuad2Pos, 'color', [1 1 1]);
axis(ImageAxis, 'off');
%origData.axisHandle1=imagesc(rot90(sliceYZ),[minInterval 2*maxInterval]);
% using image as it is more general than the imagesc
origData.axisHandle1 = image(rot90(sliceYZ), 'parent', ImageAxis, 'CDataMapping', 'scaled');
set(ImageAxis, 'clim', [minInterval 2*maxInterval]); % set the axis clim property
imageAxisPos = get(ImageAxis, 'position');
yAxisRatio = imageAxisPos(4)/imageAxisPos(3);
xAxisRatio = imageAxisPos(3)/imageAxisPos(4);
if(yAxisRatio>1)
    yAxisRatio = 1;
else
    xAxisRatio = 1;
end
R= ASPECT_RATIO_FOR_SQUARE_FIGURE;
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*yAxisRatio imageAxisPos(4)*xAxisRatio];
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*R(1) imageAxisPos(4)*R(2)];
set(ImageAxis,'position',imageAxisPos);
axQuad2Pos = imageAxisPos;

%use actual dimensions of image to keep aspect ratio
imageAxisPos = get(ImageAxis, 'position');
imageXRatio = dim(2)/dim(1);
imageYRatio = dim(1)/dim(2);
if(imageXRatio>1)
    imageXRatio = 1;
else
    imageYRatio = 1;
end
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*imageYRatio imageAxisPos(4)*imageXRatio];
set(ImageAxis,'position',imageAxisPos);

%use voxel dimensions of image to keep aspect ratio
voxYAxisRatio = 1;
voxXAxisRatio = 1;
if(voxYAxisRatio > 1)
    voxYAxisRatio =1;
else
    voxXAxisRatio =1;
end
R=ASPECT_RATIO_FOR_SQUARE_FIGURE;
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*voxYAxisRatio imageAxisPos(4)*voxXAxisRatio];
set(ImageAxis,'position',imageAxisPos);

%keeping the same aspect ratio as calculated above enlarge image to fit
% in to its allocated space
enlargeScalar = min([axQuad2Pos(3) axQuad2Pos(4)])/max([imageAxisPos(3) imageAxisPos(4)]);
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*enlargeScalar imageAxisPos(4)*enlargeScalar];
set(ImageAxis,'position',imageAxisPos);
colormap(newColor);
quad(2).xCord =[ imageAxisPos(1), imageAxisPos(1)+imageAxisPos(3),imageAxisPos(1)+imageAxisPos(3),imageAxisPos(1)];
quad(2).yCord =[ imageAxisPos(2)+imageAxisPos(4), imageAxisPos(2)+imageAxisPos(4),imageAxisPos(2),imageAxisPos(2)];
set(ImageAxis,'XTickLabel',[]);
set(ImageAxis,'YTickLabel',[]);


%Quad 3

%use dimensions of y and x axis to keep aspect ratio
% ie makes axis square within allocated space
axQuad3Pos = [.6 .4 .36 .20];
dim = size(sliceXY);
ImageAxis = axes( 'units','normalized','position',axQuad3Pos,'color',[1 1 1]);
axis(ImageAxis, 'off');
% origData.axisHandle1=imagesc(rot90(sliceXY),[minInterval 2*maxInterval]);
% use image as it is more general than imagesc
origData.axisHandle1 = image(rot90(sliceXY), 'parent', ImageAxis, 'CDataMapping', 'scaled');
set(ImageAxis, 'clim', [minInterval 2*maxInterval]); % set the axis clim property
imageAxisPos = get(ImageAxis,'position');
yAxisRatio = imageAxisPos(4)/imageAxisPos(3);
xAxisRatio = imageAxisPos(3)/imageAxisPos(4);
if(yAxisRatio>1)
    yAxisRatio = 1;
else
    xAxisRatio = 1;
end
R= ASPECT_RATIO_FOR_SQUARE_FIGURE;
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*yAxisRatio imageAxisPos(4)*xAxisRatio];
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*R(1) imageAxisPos(4)*R(2)];
set(ImageAxis,'position',imageAxisPos);
axQuadPos = imageAxisPos;

%use actual dimensions of image to keep aspect ratio
imageAxisPos = get(ImageAxis,'position');
imageXRatio = dim(2)/dim(1);
imageYRatio = dim(1)/dim(2);
if(imageXRatio>1)
    imageXRatio = 1;
else
    imageYRatio = 1;
end
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*imageYRatio imageAxisPos(4)*imageXRatio];
set(ImageAxis,'position',imageAxisPos);


%use voxel dimensions of image to keep aspect ratio
voxYAxisRatio = 1;
voxXAxisRatio = 1;
if(voxYAxisRatio > 1)
    voxYAxisRatio =1;
else
    voxXAxisRatio =1;
end
R=ASPECT_RATIO_FOR_SQUARE_FIGURE;
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*voxYAxisRatio imageAxisPos(4)*voxXAxisRatio];
set(ImageAxis,'position',imageAxisPos);

%keeping the same aspect ratio as calculated above enlarge image to fit
% in to its allocated space
enlargeScalar = min([axQuad1Pos(3) axQuad1Pos(4)])/max([imageAxisPos(3) imageAxisPos(4)]);
imageAxisPos = [imageAxisPos(1) imageAxisPos(2) imageAxisPos(3)*enlargeScalar imageAxisPos(4)*enlargeScalar];
set(ImageAxis,'position',imageAxisPos);
colormap(newColor);
quad(3).xCord =[ imageAxisPos(1), imageAxisPos(1)+imageAxisPos(3),imageAxisPos(1)+imageAxisPos(3),imageAxisPos(1)];
quad(3).yCord =[ imageAxisPos(2)+imageAxisPos(4), imageAxisPos(2)+imageAxisPos(4),imageAxisPos(2),imageAxisPos(2)];
set(ImageAxis,'XTickLabel',[]);
set(ImageAxis,'YTickLabel',[]);

%colorbar
colorbarPos = [.6 .4 .36 .20];
drawnow;
cbAxis = axes( 'units','normalized', 'position', colorbarPos, 'color', [1 1 1]);
ColorBarHandle = colorbar; %(cbAxis);
colorPos = get(ColorBarHandle, 'position');
axis off;

childH = get(ColorBarHandle,'Children');
%set(childH, 'YData', [minInterval 2*maxInterval]);
imInd = strmatch('image', lower(get(childH, 'Type')), 'exact');
set(childH(imInd), 'YData', [minInterval 2*maxInterval]);
set(ColorBarHandle, 'YLim', [minInterval maxInterval]);
set(ColorBarHandle, 'YColor', FONT_COLOR);
set(ColorBarHandle, 'YTick', [minInterval maxInterval]);
icaCLIM = round(icaCLIM*100)/100;
set(ColorBarHandle, 'YTickLabel', [icaCLIM(1) icaCLIM(2)]);
%colorPos=get(ColorBarHandle,'position');
set(ColorBarHandle, 'position', [colorPos(1)-.04,colorPos(2)+.03,colorPos(3),colorPos(4)]);
cd ..

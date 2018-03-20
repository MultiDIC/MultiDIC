function [ROImask] = selectROI(ImPaths, varargin)
%% function for selecting the ROI on speckle images in step2
% options:
%[] = selectROI(IMset);
%[] = selectROI(IMset,nROI);
%
% function for selecting the region of interest (ROI) for stereo DIC (2 views)
% requirements: GIBBON toolbox
% The idea is to draw a polygon around the ROI, which is the overlapping
% area that is visible in both views. The ROI is selected on the reference
% image, but in order to select it properly, it is advisable to look at
% both the reference and the "deformed" views.
%
% INPUT:
% * IMset - a 2nX1 cell array containing paths 2n grayscale images. The first n
% images are from camera A (the "reference" camera), and the last n images
% are from camera B (the "deformed" camera). The first image in the set is
% considered as the reference image, on which the reference grid of points
% is defined, and all the correlated points and consequent displacements
% and strains, are relative to this image.
% * nROI (optional, default=1) - number of ROIs (closed polygons). must be
% integer>=1; If nROI is empty [] or not given then it is set to the default(=1).

%%
nVars = length(varargin);

switch nVars
    case 0
        nROI=1;
    case 1
        nROI=varargin{1};
        if isempty(nROI), nROI=1; end
    otherwise
        error('Wrong number of input arguments');
end

% number of frames
nCur=numel(ImPaths);

% load images
ImSet=cell(nCur,1);
for ii=1:nCur
   ImSet{ii}=imread(ImPaths{ii}); 
end


nImages=nCur/2; %Number of animation steps
if  rem(nImages,1)~=0
    error('Number of images in the set should be even');
end

polyOK='No';
while ~strcmp(polyOK,'Yes')
    hf=figure('Name','Select ROI');
    hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';

    suptitle(['Draw ' num2str(nROI) ' ROI polygons on the reference images. double-click on the polygon to finish']);
    
%     subplot(2,2,2);
%     imshow(ImSet{nImages}); hold on; axis on;
%     title('Reference camera last image');
%     axis off
%     
%     subplot(2,2,4);
%     imshow(ImSet{2*nImages}); hold on; axis on;
%     title('Deformed camera last image');
%     axis off
%         
%     subplot(1,2,1);
    imshow(ImSet{1}); hold on; axis on;
%     title('Reference image');
    axis off
    
    for iROI=1:nROI
        hp(iROI) = impoly;
        wait(hp(iROI));
    end
    polyOK = questdlg('polygon mask OK? (if you select "No", this GUI will start over)', 'polygon mask OK?', 'Yes', 'No', 'Yes');
end

for iROI=1:nROI
    nodes = getPosition(hp(iROI));
    ROImask(:,:,iROI) = poly2mask(nodes(:,1),nodes(:,2),size(ImSet{1},1),size(ImSet{1},2));
end

ROImask=logical(sum(ROImask,3));



end
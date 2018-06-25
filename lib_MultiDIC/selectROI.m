function [ROImask] = selectROI(IM, varargin)
%% function for selecting the ROI on speckle images in step2
% options:
%[] = selectROI(IM);
%[] = selectROI(IM,nROI);
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


polyOK='No';
while ~strcmp(polyOK,'Yes')
    hf=figure('Name','Select ROI');
    hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';

    gtitle(['Draw ' num2str(nROI) ' ROI polygons on the reference images. double-click on the polygon to finish'],20);

    imshow(IM); hold on; axis on;
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
    ROImask(:,:,iROI) = poly2mask(nodes(:,1),nodes(:,2),size(IM,1),size(IM,2));
end

ROImask=logical(sum(ROImask,3));



end

%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>
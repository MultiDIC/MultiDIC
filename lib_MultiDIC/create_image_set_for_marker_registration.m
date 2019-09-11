function [ImPaths,ImSet]=create_image_set_for_marker_registration(folderPaths)
%% function for creating the image set of stereo pair in STEP2.
% this function takes the folder names where the reference and current
% image sets are located, and returns cell arrays with the image paths
% images in the folder must be named with the number after an underscore
% images numbers have to be increasing monotonically but don't have to be
% sequential. For example, you can have im_01,im_04,im_05.
% In addition, image numbers for both cameras must match.
%
% INPUTS:
% * folderPaths: 1-by-2 cell array, containing in each
%   cell the path to the speckle images. The camera indices are retrieved
%   from the folder names
% * processedImagePath: create a processed image folder, change the image path and save the processed image there
%
% OUTPUTS:
% * ImPaths: 2*Nimages cell array of paths to the speckle images
% * ImSet: 2*Nimages cell array of the images to transfer to Ncorr

%%
imageSetRef1=imageSet(folderPaths{1,1}); % image set for reference camera
imageSetRef2=imageSet(folderPaths{1,2}); % image set for reference camera

imageSetDef1=imageSet(folderPaths{2,1}); % image set for reference camera
imageSetDef2=imageSet(folderPaths{2,2}); % image set for reference camera

% image paths
imagePathsRef1=imageSetRef1.ImageLocation{1};
imagePathsRef2=imageSetRef2.ImageLocation{1};
imagePathsDef1=imageSetDef1.ImageLocation{1};
imagePathsDef2=imageSetDef2.ImageLocation{1};

ImPaths={imagePathsRef1; imagePathsDef1; imagePathsRef2; imagePathsDef2} ;
ImSet=cell(4,1);

%%  read original image, turn to gray and undistort if necessary, write if necessary

% for reference camera
for ii=1:4
    IM=imread(ImPaths{ii}); % read original image
    
    % if image is RGB, create a processed image folder, change the image path and save the processed image there
    if size(IM,3)==3
%         folderName=[processedImagePath  '\' icam1Str];
%         warning('off','MATLAB:MKDIR:DirectoryExists');
%         mkdir(folderName);
%         imName=[folderName '\' name ext];
%         
        % turn to gray if RGB
        IM=rgb2gray(IM);
        
%         % change path to processed image
%         ImPaths{ii}=imName; 
%         % write to path
%         imwrite(IM,imName,'Quality',100); 
    end
    
    ImSet{ii}=IM; % save image to set
end



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

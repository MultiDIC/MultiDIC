function [ImPaths,ImSet]=createDICimageSet(folderPaths,processedImagePath)
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
imageSet1=imageSet(folderPaths{1}); % image set for reference camera
nImages1=imageSet1.Count; % number of images for reference camera
imageSet2=imageSet(folderPaths{2}); % image set for second camera
nImages2=imageSet2.Count; % number of images for second camera
if nImages1==nImages2 % check that both sets have same number of images
    nImages=nImages1; % total number of images in the set (including ref)
else
    error('number of images for both cameras must be equal');
end

% extract camera numbers
icam1Cell=strsplit(folderPaths{1},'\');
icam1Str=icam1Cell{end};
icam1=str2double(icam1Str);
icam2Cell=strsplit(folderPaths{2},'\');
icam2Str=icam2Cell{end};
icam2=str2double(icam2Str);

% image paths
imagePaths1=imageSet1.ImageLocation;
imagePaths2=imageSet2.ImageLocation;
ImPaths=cell(2*nImages,1);
ImSet=cell(2*nImages,1);

%% find image numbers
imNums=zeros(1,nImages);
for ii=1:nImages
    [~,name1,~]=fileparts(imagePaths1{ii});
    [~,name2,~]=fileparts(imagePaths2{ii});
    imNumCell1=strsplit(name1,'_');
    imNumCell2=strsplit(name2,'_');
    imNum1=str2double(imNumCell1{2});
    imNum2=str2double(imNumCell2{2});
    if imNum1==imNum2
        imNums(ii)=imNum1;
    else
        error('Image numbers from both cameras must match');
    end
end

%%  read original image, turn to gray and undistort if necessary, write if necessary

% for reference camera
for ii=1:nImages
%     [~,name,ext]=fileparts(imagePaths1{ii});
    ImPaths{ii}=imagePaths1{ii};
    IM=imread(imagePaths1{ii}); % read original image
    
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

% for second camera
for ii=1:nImages
%     [~,name,ext]=fileparts(imagePaths2{ii});
    ImPaths{nImages+ii}=imagePaths2{ii};
    IM=imread(imagePaths2{ii}); % read original image
    
    % if image is RGB or if undistortion required, create a processed image folder, change the image path and save the processed image there
    if size(IM,3)==3
%         folderName=[processedImagePath  '\' icam2Str];
%         warning('off','MATLAB:MKDIR:DirectoryExists');
%         mkdir(folderName);
%         imName=[folderName '\' name ext];
%         
        % turn to gray if RGB
        IM=rgb2gray(IM);
        
%         % change path to processed image
%         ImPaths{nImages+ii}=imName;
%         % write to path
%         imwrite(IM,imName,'Quality',100);
    end
    
    ImSet{nImages+ii}=IM; % save image to set
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

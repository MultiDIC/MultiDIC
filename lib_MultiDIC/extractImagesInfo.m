function [CBimagesInfo]=extractImagesInfo(imagesFolderPath)
%% function for extracting the image info from the path selected by the user in STEP0
%
% INPUTS:
% * imagesFolderPath: 1-by-NumberOfCameras cell array, containing in each
%   cell the path to the folder contatining the checkerboard images. The
%   camera index is retrieved from the folder name (_Ncam)
%
% OUTPUTS:
% * CBimagesInfo: structure with fields:
%       * icam: camera index (2 last digits of folder name)
%       * I: a matrix of all CB images (Height -by- Width -by- 3 -by- Nimages)
%       * imageFileNames: 1 -by- NumberOfCameras cell array with the image locations
%       * Nimages: number of images in folder
%       * imageType: file extension
%       * imageSize: in pixels

%%
pathStringCell=strsplit(imagesFolderPath,filesep);
icamString=pathStringCell{end};
icamCell=strsplit(icamString,'_');
icam=str2num(icamCell{end});
images = imageSet(imagesFolderPath);
imageFileNames = images.ImageLocation;
[~,~,imageType]=fileparts(imageFileNames{1});
Nimages=images.Count;

I1=imread(images.ImageLocation{1});

if size(size(I1),2)==2
    imageSize3=[size(I1) 1];
else
    imageSize3=size(I1);
end
imageSize=imageSize3(1,1:2);

I=zeros([imageSize3 Nimages],'uint8');
for ip=1:Nimages
    I(:,:,:,ip)=imread(images.ImageLocation{ip});
end

CBimagesInfo.icam=icam;
CBimagesInfo.I=I;
CBimagesInfo.imageFileNames=imageFileNames;
CBimagesInfo.Nimages=Nimages;
CBimagesInfo.imageType=imageType;
CBimagesInfo.imageSize=imageSize;

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
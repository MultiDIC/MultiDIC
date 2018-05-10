function []=createCheckerBoardImage(Nrows,Ncols,squareSize,resolution,path)
%% function for writing a checker board image to be printed for camera clibration in step0. 
% The function creates the image and saves it in the desired path.

% INPUT:
% * Nrows: number of rows (uneven)
% * Ncols: number of columns (even)
% * squareSize: size of each square (in meters)
% * resolution: image resolution (pixels per meter)
% * path: path where the image will be saves (including file name)

%%
squareSizePixel=squareSize*resolution;
CB = checkerboardBW(squareSizePixel,Nrows,Ncols);
figure; imshow(CB);
imwrite(CB,[path '\CB_' num2str(Nrows) '_' num2str(Nrows) '_' num2str(squareSize*1000) '.png'],'png','ResolutionUnit','meter','XResolution',resolution); %XResolution=pixels per meter

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
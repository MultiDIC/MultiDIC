function []=plotAllCameraImages(CBimagesInfo)
%% function for asking the user for saving options for the camera parmaters in STEP0.
% Including overwriting options
%
% INPUTS:
% * CBimagesInfo: structure created in STEP0 and contains the image paths
%
% OUTPUTS:
% * a figure of all the checkerboard images

%%

I=CBimagesInfo.I;
icam=CBimagesInfo.icam;
imageSize=CBimagesInfo.imageSize;

figure('name',['Camera ' num2str(icam) ', all checkerboard images'],'units','normalized','outerposition',[.1 .1 .8 .8]);

Nimages=size(I,4);
Nsubplots=numSubplots(Nimages);
for ip=1:Nimages
    subtightplot(Nsubplots(1),Nsubplots(2),ip);
    imshow(I(:,:,:,ip)); hold on
    text(.1*imageSize(2),.1*imageSize(1),num2str(ip),'color','b','backgroundcolor','w');
%     title(['Image ' num2str(ip)]);
    drawnow
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
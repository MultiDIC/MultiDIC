function [hf] = anim8_DIC_images(ImPaths,ImSet)
%% function for plotting animated figure of all speckle images, for reviweing before selecting ROI, in step 2
% anim8_DIC_images
% [] = anim8_DIC_images(ImPaths);
%
% function for plotting the images chosen for stereo DIC (2 views)
% requirements: GIBBON toolbox
%
% INPUT: 
% * ImPaths - a 2nX1 cell array containing 2n paths to grayscale images. The first n
% images are from camera A (the "reference" camera), and the last n images
% are from camera B (the "deformed" camera). The first image in the set is
% considered as the reference image, on which the reference grid of points
% is defined, and all the correlated points and consequent displacements
% and strains, are relative to this image.

%%
nCur=numel(ImPaths); % number of frames

nSteps=nCur/2; %Number of animation steps
if  rem(nSteps,1)~=0
    error('Number of images in the set should be even');
end

nImCurSet=repmat(1:nSteps,1,2);

%%
hf=cFigure;
hf.Name='Preview Images';
hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';

animStruct=struct; % create animation structure

% show 1st image from view 1 (Reference Camera)
subplot(1,2,1)
hp1=imagesc(repmat(ImSet{1},1,1,3)); hold on
pbaspect([size(ImSet{1},2) size(ImSet{1},1) 1]);
hs1=title(['Cur ' num2str(1) ' (Cam Ref, frame ' num2str(nImCurSet(1)) ')']);
axis off

% show 1st image from view 2 ("Deformed" Camera)
subplot(1,2,2)
hp2=imagesc(repmat(ImSet{1+nCur/2},1,1,3)); hold on
pbaspect([size(ImSet{1+nCur/2},2) size(ImSet{1+nCur/2},1) 1]);
hs2=title(['Cur ' num2str(1) ' (Cam Def, frame ' num2str(nImCurSet(1)) ')']);
axis off

%Create the "time" (frame) vector
animStruct.Time=linspace(0,1,nSteps);

%Populate the animStruct
for ii=1:1:nSteps      
    TitleNow2=['Cur ' num2str(ii) ' (Cam Ref, frame ' num2str(nImCurSet(ii)) ')'];
    TitleNow3=['Cur ' num2str(ii+nCur/2) ' (Cam Def, frame ' num2str(nImCurSet(ii+nCur/2)) ')'];
    
    %Set entries in animation structure
    animStruct.Handles{ii}=[hp1,hp2,hs1,hs2]; %Handles of objects to animate
    animStruct.Props{ii}={'CData','CData','String','String'}; %Properties of objects to animate
    animStruct.Set{ii}={repmat(ImSet{ii},1,1,3),repmat(ImSet{ii+nCur/2},1,1,3),TitleNow2,TitleNow3}; %Property values for to set in order to animate
   
end

gtitle('Image set for DIC. Press the play button to preview all images, press the stop button to stop animation, and then press any key in the command window to continue',14);

% Start |anim8| gui
anim8(hf,animStruct);


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
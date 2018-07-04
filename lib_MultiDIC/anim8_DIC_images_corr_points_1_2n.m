function [] = anim8_DIC_images_corr_points_1_2n(ImSet,DIC_2Dpair_results,varargin)
%% function for plotting 2D-DIC results imported from Ncorr in step 2
% called inside plotNcorrPairResults
% plotting the images chosen for stereo DIC (2 views) with the
% correlated points results plotted on top, colored as their correlation
% coefficient.
% on the left side the reference image, and on the right side the current
% images from both views
% requirements: GIBBON toolbox
%
% calling options:
% [] = anim8_DIC_images(IMset,DIC_2Dpair_results);
% [] = anim8_DIC_images(IMset,DIC_2Dpair_results,CorCoeffCutOff,CorCoeffDispMax);
%
% INPUT:
% * IMset - a 2nX1 cell array containing 2n grayscale images. The first n
% images are from camera A (the "reference" camera), and the last n images
% are from camera B (the "deformed" camera). The first image in the set is
% considered as the reference image, on which the reference grid of points
% is defined, and all the correlated points and consequent displacements
% and strains, are relative to this image.
% * DIC_2Dpair_results - containig the correlated points, correlation
% coefficients, faces..
% * optional: CorCoeffCutOff - - maximal correlation coefficient to plot
% points
% * optional: CorCoeffDispMax - maximal correlation coefficient in colorbar

%%
nCur=numel(ImSet); % number of frames

nSteps=nCur/2; %Number of animation steps
if  rem(nSteps,1)~=0
    error('Number of images in the set should be even');
end

%%
Points=DIC_2Dpair_results.Points;
CorCoeffVec=DIC_2Dpair_results.CorCoeffVec;
nCamRef=DIC_2Dpair_results.nCamRef;
nCamDef=DIC_2Dpair_results.nCamDef;
nImages=DIC_2Dpair_results.nImages;

nVars = length(varargin);
switch nVars
    case 2
        CorCoeffCutOff=varargin{1};
        if isnan(CorCoeffCutOff)
            CorCoeffCutOff=max(max([CorCoeffVec{:}]));
        end
        CorCoeffDispMax=varargin{2};
        if isnan(CorCoeffDispMax)
            CorCoeffDispMax=max(max([CorCoeffVec{:}]));
        end
    case 1
        CorCoeffCutOff=varargin{1};
        if isnan(CorCoeffCutOff)
            CorCoeffCutOff=max(max([CorCoeffVec{:}]));
        end
        CorCoeffDispMax=max(max([CorCoeffVec{:}]));
    case 0
        CorCoeffCutOff=max(max([CorCoeffVec{:}]));
        CorCoeffDispMax=max(max([CorCoeffVec{:}]));
    otherwise
        error('Wrong number of input arguments');
end

%% cut out point with large correlation coefficient
for ii=1:2*nImages
    CorCoeffVec{ii}(CorCoeffVec{ii}>CorCoeffCutOff)=NaN;
end

%%
hf=cFigure;
hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';

ii=1;
% Ref
subplot(1,2,1)
hp1=imagesc(repmat(ImSet{ii},1,1,3)); hold on;
P=Points{ii}(~isnan(CorCoeffVec{ii}),:);
hp2=scatter(P(:,1),P(:,2),6,CorCoeffVec{ii}(~isnan(CorCoeffVec{ii})),'+');
colormap jet
pbaspect([size(ImSet{ii},2) size(ImSet{ii},1) 1])
title(['Ref (Cam ' num2str(nCamRef) ' frame ' num2str(1) ')']);
hc1=colorbar; caxis([0 CorCoeffDispMax])
title(hc1, 'Corr-Coeff')
axis off

% Cur
subplot(1,2,2)
hp3=imagesc(repmat(ImSet{ii},1,1,3)); hold on
hp4=scatter(P(:,1),P(:,2),6,CorCoeffVec{ii}(~isnan(CorCoeffVec{ii})),'+');
colormap jet
pbaspect([size(ImSet{ii},2) size(ImSet{ii},1) 1])
hs2=title(['Cur ' num2str(ii) ' (Cam ' num2str(nCamRef) ' frame ' num2str(1) ')']);
hc2=colorbar;
caxis([0 CorCoeffDispMax])
title(hc2, 'Corr-Coeff')
axis off

drawnow

%Create the time vector
animStruct.Time=linspace(0,1,2*nImages);

for ii=1:2*nImages
    xNow1=Points{1}(~isnan(CorCoeffVec{ii}),1);
    yNow1=Points{1}(~isnan(CorCoeffVec{ii}),2);
    xNow2=Points{ii}(~isnan(CorCoeffVec{ii}),1);
    yNow2=Points{ii}(~isnan(CorCoeffVec{ii}),2);
    
    cNow=CorCoeffVec{ii}(~isnan(CorCoeffVec{ii}));
    
    if ii<=nImages
        nCamCur=nCamRef;
        nImCur=ii;
    else
        nCamCur=nCamDef;
        nImCur=ii-nImages;
    end
    TitleNow=['Cur ' num2str(ii) ' (Cam ' num2str(nCamCur) ' frame ' num2str(nImCur) ')'];
    
    %Set entries in animation structure
    animStruct.Handles{ii}=[hp3,hp2,hp2,hp2,hp4,hp4,hp4,hs2]; %Handles of objects to animate
    animStruct.Props{ii}={'CData','XData','YData','CData','XData','YData','CData','String'}; %Properties of objects to animate
    animStruct.Set{ii}={repmat(ImSet{ii},1,1,3),xNow1,yNow1,cNow,xNow2,yNow2,cNow,TitleNow}; %Property values for to set in order to animate
    
end

anim8(hf,animStruct);

addColorbarLimitsButton(hf);
addColormapButton(hf);

end

%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% 
% Modified by Rana Odabas 2018
% 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>
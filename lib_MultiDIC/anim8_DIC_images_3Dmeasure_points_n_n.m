function [] = anim8_DIC_images_3Dmeasure_points_n_n(ImSet,DIC2DpairResults,DIC3DpairResults,pointMeasureStr,varargin)
%% function for plotting 2D-DIC results imported from Ncorr in step 2
% called inside plotNcorrPairResults
% plotting the images chosen for stereo DIC (2 views) with the
% correlated points results plotted on top, colored as their correlation
% coefficient.
% on the left side the images from the reference camera (reference image and current images), and on the right side the 
% images from the deformed camera
% requirements: GIBBON toolbox
%
% calling options:
% [] = anim8_DIC_images_3Dmeasure_points_n_n(IMset,DIC2DpairResults,DIC3DpairResults,pointMeasureStr);
% [] = anim8_DIC_images_3Dmeasure_points_n_n(IMset,DIC2DpairResults,DIC3DpairResults,pointMeasureStr,optStruct);
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
Points=DIC2DpairResults.Points;
CorCoeffVec=DIC2DpairResults.CorCoeffVec;
nCamRef=DIC2DpairResults.nCamRef;
nCamDef=DIC2DpairResults.nCamDef;
nImages=DIC2DpairResults.nImages;

switch nargin
    case 4 % in case no results were entered
        optStruct=struct;
    case 5
        optStruct=varargin{1};
    otherwise
        error('wrong number of input arguments');
end


%% cut out point with large correlation coefficient
if ~isfield(optStruct,'CorCoeffCutOff')
    CorCoeffCutOff=max(cell2mat(CorCoeffVec));
else
    CorCoeffCutOff=optStruct.CorCoeffCutOff;
end

for ii=1:2*nImages
    CorCoeffVec{ii}(CorCoeffVec{ii}>CorCoeffCutOff)=NaN;   
end

%%
switch pointMeasureStr
    case {'DispMgn'}
        PC=DIC3DpairResults.Disp.DispMgn;      
        colorBarLogic=1;
        cMap='jet';
    case {'DispX'}
        for ii=1:size(DIC3DpairResults.Disp.DispVec,1)
            PC{ii}=DIC3DpairResults.Disp.DispVec{ii}(:,1);
        end
        cMap='jet';
    case {'DispY'}
        for ii=1:size(DIC3DpairResults.Disp.DispVec,1)
            PC{ii}=DIC3DpairResults.Disp.DispVec{ii}(:,2);
        end
        cMap='jet';
    case {'DispZ'}
        for ii=1:size(DIC3DpairResults.Disp.DispVec,1)
            PC{ii}=DIC3DpairResults.Disp.DispVec{ii}(:,3);
        end
        cMap='jet';
    otherwise
        error('unexpected face measure string. plots not created');
        
end

%%

if ~isfield(optStruct,'PClimits')    
    PCmin=0;
    PCmax=0;
    for ii=1:nImages
        PCmin=min([min(PC{ii}(~isnan(CorCoeffVec{ii}))) PCmin]);
        PCmax=max([max(PC{ii}(~isnan(CorCoeffVec{ii}))) PCmax]);
    end
    PClimits=[PCmin PCmax];
else
    PClimits=optStruct.PClimits;
end


%%
hf=cFigure;
hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';

% Ref
ii=1;
subplot(1,2,1)
hp1=imagesc(repmat(ImSet{ii},1,1,3)); hold on;
P=Points{ii}(~isnan(CorCoeffVec{ii}),:);
hp2=scatter(P(:,1),P(:,2),6,PC{ii}(~isnan(CorCoeffVec{ii})),'+');
pbaspect([size(ImSet{ii},2) size(ImSet{ii},1) 1])
hs1=title(['Ref (Cam ' num2str(nCamRef) ' frame ' num2str(1) ')']);
hc1=colorbar; 
caxis(PClimits)
title(hc1, pointMeasureStr)
hc1.FontSize=16;
axis off

% Cur
ii=nImages+1;
subplot(1,2,2)
hp3=imagesc(repmat(ImSet{ii},1,1,3)); hold on
P=Points{ii}(~isnan(CorCoeffVec{ii}),:);
hp4=scatter(P(:,1),P(:,2),6,PC{ii-nImages}(~isnan(CorCoeffVec{ii})),'+');
colormap jet
pbaspect([size(ImSet{ii},2) size(ImSet{ii},1) 1])
hs2=title(['Cur ' num2str(ii) ' (Cam ' num2str(nCamDef) ' frame ' num2str(1) ')']);
hc2=colorbar; 
caxis(PClimits)
title(hc2, pointMeasureStr)
hc2.FontSize=16;
axis off

colormap(cMap);
drawnow

%Create the time vector
animStruct.Time=linspace(0,1,nImages);

for ii=1:nImages  
    xNow1=Points{ii}(~isnan(CorCoeffVec{ii}),1);
    yNow1=Points{ii}(~isnan(CorCoeffVec{ii}),2);
    xNow2=Points{ii+nImages}(~isnan(CorCoeffVec{ii+nImages}),1);
    yNow2=Points{ii+nImages}(~isnan(CorCoeffVec{ii+nImages}),2);
    
    cNow1=PC{ii}(~isnan(CorCoeffVec{ii}));
    cNow2=PC{ii}(~isnan(CorCoeffVec{ii+nImages}));
    
    TitleNow1=['Cur ' num2str(ii) ' (Cam ' num2str(nCamRef) ' frame ' num2str(ii) ')'];
    TitleNow2=['Cur ' num2str(ii) ' (Cam ' num2str(nCamDef) ' frame ' num2str(ii) ')'];
    
   %Set entries in animation structure
    animStruct.Handles{ii}=[hp1,hp3,hp2,hp2,hp2,hp4,hp4,hp4,hs1,hs2]; %Handles of objects to animate
    animStruct.Props{ii}={'CData','CData','XData','YData','CData','XData','YData','CData','String','String'}; %Properties of objects to animate
    animStruct.Set{ii}={repmat(ImSet{ii},1,1,3),repmat(ImSet{ii+nImages},1,1,3),xNow1,yNow1,cNow1,xNow2,yNow2,cNow2,TitleNow1,TitleNow2}; %Property values for to set in order to animate
   
end

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

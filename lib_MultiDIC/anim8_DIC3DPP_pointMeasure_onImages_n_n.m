function [] = anim8_DIC3DPP_pointMeasure_onImages_n_n(DIC3DPPresults,pairIndex,pointMeasureStr,varargin)
%% function for plotting 3D-DIC post-processing results from step 4 projected on the 2D images
% plotting the images with the 3D point measure results plotted on top
% on the left side the images from the reference camera (reference image and current images), and on the right side the
% images from the deformed camera
%
%
% calling options:
% [] = anim8_DIC3D_faceMeasures_onImages_n_n(DIC3DPPresults,pairIndex,pointMeasureString);
% [] = anim8_DIC3D_faceMeasures_onImages_n_n(DIC3DPPresults,pairIndex,pointMeasureString,optStruct);
%
% INPUT:

%%
Points=DIC3DPPresults.DIC2Dinfo{pairIndex}.Points;
nCamRef=DIC3DPPresults.DIC2Dinfo{pairIndex}.nCamRef;
nCamDef=DIC3DPPresults.DIC2Dinfo{pairIndex}.nCamDef;
nImages=DIC3DPPresults.DIC2Dinfo{pairIndex}.nImages;
currentPointLogic=DIC3DPPresults.PointPairInds==pairIndex;
for ii=1:nImages
    CorCoeffVec{ii}=DIC3DPPresults.corrComb{ii}(currentPointLogic,:);
end

ImPaths=DIC3DPPresults.DIC2Dinfo{pairIndex}.ImPaths;

% if ImPaths are not valid (for example if using on another computer, ask
% user to provide a new folder where all the images are located.
try
   ImSet{1}=imread(ImPaths{1});
catch
    newPath=uigetdir([],'Path to images is invalid. Please provide the correct path to the images (the folder containing all camera folders with the processed gray images)');
    for ii=1:length(ImPaths)
        [a,b,c]=fileparts(ImPaths{ii});
        as = strsplit(a,{'\','/'});
        ImPaths{ii}=[newPath '\' as{end} '\' b c];
    end
end

for ii=1:2*nImages
    ImSet{ii}=imread(ImPaths{ii});
    
    if size(ImSet{ii},3)==3
        ImSet{ii}=rgb2gray(ImSet{ii});     
    end
end

switch nargin
    case 3 % in case no results were entered
        optStruct=struct;
    case 4
        optStruct=varargin{1};
    otherwise
        error('wrong number of input arguments');
end

%% cut out point with large correlation coefficient
if ~isfield(optStruct,'CorCoeffCutOff')
    CorCoeffCutOff=max(max(cell2mat(CorCoeffVec)));
else
    CorCoeffCutOff=optStruct.CorCoeffCutOff;
end

if ~isfield(optStruct,'logicRBM')
    logicRBM=0;
else
    logicRBM=optStruct.logicRBM;
end

for ii=1:nImages
    CorCoeffVec{ii}(CorCoeffVec{ii}>CorCoeffCutOff)=NaN;   
end

%%
switch pointMeasureStr
    case {'DispMgn'}
        if logicRBM
            for ii=1:nImages
                PC{ii}=DIC3DPPresults.Disp.DispMgn_ARBM{ii}(currentPointLogic,:);
            end
        else
            for ii=1:nImages
                PC{ii}=DIC3DPPresults.Disp.DispMgn{ii}(currentPointLogic,:);
            end
        end
        cMap='parula';
    case {'DispX'}
        for ii=1:nImages
            PC{ii}=DIC3DPPresults.Disp.DispVec{ii}(currentPointLogic,1);
        end
        cMap='coldwarm';
    case {'DispY'}
        for ii=1:nImages
            PC{ii}=DIC3DPPresults.Disp.DispVec{ii}(currentPointLogic,2);
        end
        cMap='coldwarm';
    case {'DispZ'}
        for ii=1:nImages
            PC{ii}=DIC3DPPresults.Disp.DispVec{ii}(currentPointLogic,3);
        end
        cMap='coldwarm';
    otherwise
        error('unexpected face measure string. plots not created');
        
end

%%

if ~isfield(optStruct,'PClimits') 
    switch pointMeasureStr
        case {'DispMgn'}
            PCmin=0;
            PCmax=0;
            for ii=1:nImages
                PCmax=max([max(PC{ii}(~isnan(CorCoeffVec{ii}))) PCmax]);
            end
            PClimits=[PCmin PCmax];
        case {'DispX','DispY','DispZ'}
            PCmax=0;
            for ii=1:nImages
                PCmax=max([max(abs(PC{ii}(~isnan(CorCoeffVec{ii})))) PCmax]);
            end
            PClimits=[-PCmax PCmax];
            
    end   
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
P=Points{ii}(~isnan(CorCoeffVec{ii-nImages}),:);
hp4=scatter(P(:,1),P(:,2),6,PC{ii-nImages}(~isnan(CorCoeffVec{ii-nImages})),'+');
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
    xNow2=Points{ii+nImages}(~isnan(CorCoeffVec{ii}),1);
    yNow2=Points{ii+nImages}(~isnan(CorCoeffVec{ii}),2);
    
    cNow1=PC{ii}(~isnan(CorCoeffVec{ii}));
    cNow2=PC{ii}(~isnan(CorCoeffVec{ii}));
    
    TitleNow1=['Cur ' num2str(ii) ' (Cam ' num2str(nCamRef) ' frame ' num2str(ii) ')'];
    TitleNow2=['Cur ' num2str(ii) ' (Cam ' num2str(nCamDef) ' frame ' num2str(ii) ')'];
    
   %Set entries in animation structure
    animStruct.Handles{ii}=[hp1,hp3,hp2,hp2,hp2,hp4,hp4,hp4,hs1,hs2]; %Handles of objects to animate
    animStruct.Props{ii}={'CData','CData','XData','YData','CData','XData','YData','CData','String','String'}; %Properties of objects to animate
    animStruct.Set{ii}={repmat(ImSet{ii},1,1,3),repmat(ImSet{ii+nImages},1,1,3),xNow1,yNow1,cNow1,xNow2,yNow2,cNow2,TitleNow1,TitleNow2}; %Property values for to set in order to animate
   
end

anim8(hf,animStruct);
addFigureButtons;

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

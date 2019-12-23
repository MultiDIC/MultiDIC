function [] = anim8_DIC3DPP_faceMeasure_onImages_oneCamera(DIC3DPPresults,cameraIndex,faceMeasureString,varargin)
%% function for plotting 3D-DIC results on images from one camera (1 or 2 pairs)
% plotting the images with the 3D face measure results plotted on top
%
% calling options:
% [] = anim8_DIC3DPP_faceMeasure_onImages_oneCamera(DIC3DPPresults,pairIndex,faceMeasureString);
% [] = anim8_DIC3DPP_faceMeasure_onImages_oneCamera(DIC3DPPresults,pairIndex,faceMeasureString,optStruct);
%
% INPUT:
%%
%% find pairs containing cameraIndex
logicPairIndex=DIC3DPPresults.pairIndices==cameraIndex;
numPairs=nnz(logicPairIndex);
pairIndices=find(any(logicPairIndex,2));
for ii=1:numPairs
    pairCamera(ii,1)=find(logicPairIndex(pairIndices(ii),:));
end

nImages=DIC3DPPresults.DIC2Dinfo{pairIndices(1)}.nImages;

%%
Points1=DIC3DPPresults.DIC2Dinfo{pairIndices(1)}.Points((pairCamera(1)-1)*nImages+(1:nImages));
Points2=DIC3DPPresults.DIC2Dinfo{pairIndices(2)}.Points((pairCamera(2)-1)*nImages+(1:nImages));
currentFacesLogic1=DIC3DPPresults.FacePairInds==pairIndices(1);
currentPointIndex1=find(DIC3DPPresults.PointPairInds==pairIndices(1));
currentFacesLogic2=DIC3DPPresults.FacePairInds==pairIndices(2);
currentPointIndex2=find(DIC3DPPresults.PointPairInds==pairIndices(2));

firstCurrentPointIndex1=currentPointIndex1(1);
firstCurrentPointIndex2=currentPointIndex2(1);

F1=DIC3DPPresults.Faces(currentFacesLogic1,:);
F1=F1-firstCurrentPointIndex1+1;
F2=DIC3DPPresults.Faces(currentFacesLogic2,:);
F2=F2-firstCurrentPointIndex2+1;


FaceCorr1=cell(nImages,1);
FaceCorr2=cell(nImages,1);
for ii=1:nImages
    FaceCorr1{ii}=DIC3DPPresults.FaceCorrComb{ii}(currentFacesLogic1,:);
    FaceCorr2{ii}=DIC3DPPresults.FaceCorrComb{ii}(currentFacesLogic2,:);
end
ImSet=cell(nImages,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImPaths=DIC3DPPresults.DIC2Dinfo{pairIndices(1)}.ImPaths;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:nImages
    ImSet{ii}=imread(ImPaths{(pairCamera(1)-1)*nImages+ii});
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
    CorCoeffCutOff1=max(max(cell2mat(FaceCorr1)));
    CorCoeffCutOff2=max(max(cell2mat(FaceCorr2)));
else
    CorCoeffCutOff1=optStruct.CorCoeffCutOff;
    CorCoeffCutOff2=optStruct.CorCoeffCutOff;
end

for ii=1:nImages
    FaceCorr1{ii}(FaceCorr1{ii}>CorCoeffCutOff1)=NaN;
    FaceCorr2{ii}(FaceCorr2{ii}>CorCoeffCutOff2)=NaN;
end

%%
FC1=cell(nImages,1);
FC2=cell(nImages,1);
switch faceMeasureString
    case {'J','Lamda1','Lamda2'}
        for ii=1:nImages
            FC1{ii}=DIC3DPPresults.Deform.(faceMeasureString){ii}(currentFacesLogic1,:);
            FC2{ii}=DIC3DPPresults.Deform.(faceMeasureString){ii}(currentFacesLogic2,:);
        end
        cMap=coldwarm;
        if ~isfield(optStruct,'FClimits')
            FCmax=0;
            for ii=1:nImages
                FCmax=max([max(abs(FC1{ii}(~isnan(FaceCorr1{ii}))-1)) max(abs(FC2{ii}(~isnan(FaceCorr2{ii}))-1)) FCmax]);
            end
            FClimits=[1-FCmax 1+FCmax];
        else
            FClimits=optStruct.FClimits;
        end
    case {'Emgn','emgn','Eeq','eeq','EShearMax','eShearMax'}
        for ii=1:nImages
            FC1{ii}=DIC3DPPresults.Deform.(faceMeasureString){ii}(currentFacesLogic1,:);
            FC2{ii}=DIC3DPPresults.Deform.(faceMeasureString){ii}(currentFacesLogic2,:);
        end
        cMap='parula';
        if ~isfield(optStruct,'FClimits')
            FCmax=0;
            for ii=1:nImages
                FCmax=max([max(FC1{ii}(~isnan(FaceCorr1{ii}))) max(FC2{ii}(~isnan(FaceCorr2{ii}))) FCmax]);
            end
            FClimits=[0 FCmax];
        else
            FClimits=optStruct.FClimits;
        end
    case {'Epc1','Epc2','epc1','epc2'}
        for ii=1:nImages
            FC1{ii}=DIC3DPPresults.Deform.(faceMeasureString){ii}(currentFacesLogic1,:);
            FC2{ii}=DIC3DPPresults.Deform.(faceMeasureString){ii}(currentFacesLogic2,:);
        end
        cMap=coldwarm;
        if ~isfield(optStruct,'FClimits')
            FCmax=0;
            for ii=1:nImages
                FCmax=max([max(abs(FC1{ii}(~isnan(FaceCorr1{ii})))) max(abs(FC2{ii}(~isnan(FaceCorr2{ii})))) FCmax]);
            end
            FClimits=[-FCmax FCmax];
        else
            FClimits=optStruct.FClimits;
        end
    otherwise
        error('unexpected face measure string. plots not created');
end

for ii=1:nImages
    FC1{ii}(isnan(FaceCorr1{ii}))=NaN;
    FC2{ii}(isnan(FaceCorr2{ii}))=NaN;
end

%% cut out point with extreme face color values
if ~isfield(optStruct,'dataLimits')
    dataLimits=[min(min([cell2mat(FC1); cell2mat(FC2)])) max(max([cell2mat(FC1); cell2mat(FC2)]))];
else
    dataLimits=optStruct.dataLimits;
end

for ii=1:nImages
    FC1{ii}(FC1{ii}<dataLimits(1) | FC1{ii}>dataLimits(2))=NaN;
    FC2{ii}(FC2{ii}<dataLimits(1) | FC2{ii}>dataLimits(2))=NaN;
end

%% plot
hf=cFigure;
hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';

% Ref
ii=1;
hp1=imagesc(repmat(ImSet{ii},1,1,3)); hold on;
hp2=gpatch(F1,Points1{ii},FC1{ii},'none',0.5); hold on;
hp3=gpatch(F2,Points2{ii},FC2{ii},'none',0.5);
pbaspect([size(ImSet{ii},2) size(ImSet{ii},1) 1])
hs1=title(['Cam ' num2str(cameraIndex) ' frame ' num2str(1)]);
colormap(cMap);
hc1=colorbar;
caxis(FClimits)
title(hc1, faceMeasureString);
hc1.FontSize=16;
axis off
drawnow

%Create the time vector
animStruct.Time=linspace(0,1,nImages);

for ii=1:nImages
    Pnow1=Points1{ii};
    Pnow2=Points2{ii};
    
    cNow1=FC1{ii};
    cNow2=FC2{ii};
    
    TitleNow=['Cam ' num2str(cameraIndex) ' frame ' num2str(ii)];
    
    %Set entries in animation structure
    animStruct.Handles{ii}=[hp1,hp2,hp2,hp3,hp3,hs1]; %Handles of objects to animate
    animStruct.Props{ii}={'CData','Vertices','CData','Vertices','CData','String'}; %Properties of objects to animate
    animStruct.Set{ii}={repmat(ImSet{ii},1,1,3),Pnow1,cNow1,Pnow2,cNow2,TitleNow}; %Property values for to set in order to animate
    
end

anim8(hf,animStruct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% cut out point with large correlation coefficient
% if ~isfield(optStruct,'CorCoeffCutOff')
%     CorCoeffCutOff=max([cell2mat(FaceCorrL); cell2mat(FaceCorrR)]);
% else
%     CorCoeffCutOff=optStruct.CorCoeffCutOff;
% end
% 
% for ii=1:nImages
%     FaceCorrL{ii}(FaceCorrL{ii}>CorCoeffCutOff)=NaN;
%     FaceCorrR{ii}(FaceCorrR{ii}>CorCoeffCutOff)=NaN;
% end
% 
% %%
% switch faceMeasureString
%     case {'J','Lamda1','Lamda2'}
%         FCL=DIC3DpairResultsL.Deform.(faceMeasureString);
%         FCR=DIC3DpairResultsR.Deform.(faceMeasureString);
%         cMap=coldwarm;
%         if ~isfield(optStruct,'FClimits')
%             FCmax=0;
%             for ii=1:nImages
%                 FCmax=max([max(abs(FCL{ii}(~isnan(FaceCorrL{ii}))-1)) max(abs(FCR{ii}(~isnan(FaceCorrR{ii}))-1)) FCmax]);
%             end
%             FClimits=[1-FCmax 1+FCmax];
%         else
%             FClimits=optStruct.FClimits;
%         end
%     case {'Emgn','emgn'}
%         FCL=DIC3DpairResultsL.Deform.(faceMeasureString);
%         FCR=DIC3DpairResultsR.Deform.(faceMeasureString);        cMap='parula';
%         if ~isfield(optStruct,'FClimits')
%             FCmax=0;
%             for ii=1:nImages
%                 FCmax=max([max(FCL{ii}(~isnan(FaceCorrL{ii}))) max(FCR{ii}(~isnan(FaceCorrR{ii}))) FCmax]);
%             end
%             FClimits=[0 FCmax];
%         else
%             FClimits=optStruct.FClimits;
%         end
%     case {'Epc1','Epc2','epc1','epc2'}
%         FCL=DIC3DpairResultsL.Deform.(faceMeasureString);
%         FCR=DIC3DpairResultsR.Deform.(faceMeasureString);        cMap=coldwarm;
%         if ~isfield(optStruct,'FClimits')
%             FCmax=0;
%             for ii=1:nImages
%                 FCmax=max([max(abs(FCL{ii}(~isnan(FaceCorrL{ii})))) max(abs(FCR{ii}(~isnan(FaceCorrR{ii})))) FCmax]);
%             end
%             FClimits=[-FCmax FCmax];
%         else
%             FClimits=optStruct.FClimits;
%         end
%     otherwise
%         error('unexpected face measure string. plots not created');      
% end
% 
% for ii=1:nImages
%     FCL{ii}(isnan(FaceCorrL{ii}))=NaN;
%     FCR{ii}(isnan(FaceCorrR{ii}))=NaN;
% end
% %%
% hf=figure; hold all;
% hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';
% 
% ii=1;
% hp1=imagesc(repmat(ImSet{ii},1,1,3)); hold on; axis ij
% hp2=gpatch(FL,PointsL{1},FCL{ii},'none',0.5);
% hp3=gpatch(FR,PointsR{1},FCR{ii},'none',0.5);
% pbaspect([size(ImSet{ii},2) size(ImSet{ii},1) 1])
% hs1=title(['Cam ' num2str(nCam) ' frame ' num2str(1)]);
% colormap(cMap);
% hc1=colorbar;
% caxis(FClimits)
% title(hc1, faceMeasureString);
% hc1.FontSize=16;
% axis off
% drawnow
% 
% %Create the time vector
% animStruct.Time=linspace(0,1,nImages);
% 
% for ii=1:nImages
%     PnowL=PointsL{ii};
%     PnowR=PointsR{ii};
%     
%     cNowL=FCL{ii};
%     cNowR=FCR{ii};
%     
%     TitleNow=['Cam ' num2str(nCam) ' frame ' num2str(ii)];
%     
%     %Set entries in animation structure
%     animStruct.Handles{ii}=[hp1,hp2,hp2,hp3,hp3,hs1]; %Handles of objects to animate
%     animStruct.Props{ii}={'CData','CData','Vertices','CData','Vertices','String'}; %Properties of objects to animate
%     animStruct.Set{ii}={repmat(ImSet{ii},1,1,3),cNowL,PnowL,cNowR,PnowR,TitleNow}; %Property values for to set in order to animate
%     
% end
% 
% anim8(hf,animStruct);

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
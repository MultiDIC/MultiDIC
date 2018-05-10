function [] = anim8_DIC_image_3Dmeasure_faces_2n(ImSet,DIC2DpairResults,DIC3DpairResults,faceMeasureString,varargin)
%% function for plotting 2D-DIC results imported from Ncorr in step 2
% called inside plotNcorrPairResults
% plotting the images chosen for stereo DIC (2 views) with the
% correlated points results plotted on top, colored as their correlation
% coefficient.
% on the left side the images from the reference camera (reference image and current images), and on the right side the
% images from the deformed camera
% requirements: GIBBON toolbox
%
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
DIC2DpairResultsL=DIC2DpairResults{1};
DIC2DpairResultsR=DIC2DpairResults{2};
DIC3DpairResultsL=DIC3DpairResults{1};
DIC3DpairResultsR=DIC3DpairResults{2};

nImages=DIC2DpairResultsL.nImages;
nCam=DIC2DpairResultsL.nCamDef;

PointsL=DIC2DpairResultsL.Points(nImages+1:end);
PointsR=DIC2DpairResultsR.Points(1:nImages);
FL=DIC2DpairResultsL.Faces;
FR=DIC2DpairResultsR.Faces;
FaceCorrL=DIC3DpairResultsL.FaceCorrComb;
FaceCorrR=DIC3DpairResultsR.FaceCorrComb;

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
    CorCoeffCutOff=max([cell2mat(FaceCorrL); cell2mat(FaceCorrR)]);
else
    CorCoeffCutOff=optStruct.CorCoeffCutOff;
end

for ii=1:nImages
    FaceCorrL{ii}(FaceCorrL{ii}>CorCoeffCutOff)=NaN;
    FaceCorrR{ii}(FaceCorrR{ii}>CorCoeffCutOff)=NaN;
end

%%
switch faceMeasureString
    case {'J','Lamda1','Lamda2'}
        FCL=DIC3DpairResultsL.Deform.(faceMeasureString);
        FCR=DIC3DpairResultsR.Deform.(faceMeasureString);
        cMap=coldwarm;
        if ~isfield(optStruct,'FClimits')
            FCmax=0;
            for ii=1:nImages
                FCmax=max([max(abs(FCL{ii}(~isnan(FaceCorrL{ii}))-1)) max(abs(FCR{ii}(~isnan(FaceCorrR{ii}))-1)) FCmax]);
            end
            FClimits=[1-FCmax 1+FCmax];
        else
            FClimits=optStruct.FClimits;
        end
    case {'Emgn','emgn'}
        FCL=DIC3DpairResultsL.Deform.(faceMeasureString);
        FCR=DIC3DpairResultsR.Deform.(faceMeasureString);        cMap='parula';
        if ~isfield(optStruct,'FClimits')
            FCmax=0;
            for ii=1:nImages
                FCmax=max([max(FCL{ii}(~isnan(FaceCorrL{ii}))) max(FCR{ii}(~isnan(FaceCorrR{ii}))) FCmax]);
            end
            FClimits=[0 FCmax];
        else
            FClimits=optStruct.FClimits;
        end
    case {'Epc1','Epc2','epc1','epc2'}
        FCL=DIC3DpairResultsL.Deform.(faceMeasureString);
        FCR=DIC3DpairResultsR.Deform.(faceMeasureString);        cMap=coldwarm;
        if ~isfield(optStruct,'FClimits')
            FCmax=0;
            for ii=1:nImages
                FCmax=max([max(abs(FCL{ii}(~isnan(FaceCorrL{ii})))) max(abs(FCR{ii}(~isnan(FaceCorrR{ii})))) FCmax]);
            end
            FClimits=[-FCmax FCmax];
        else
            FClimits=optStruct.FClimits;
        end
    otherwise
        error('unexpected face measure string. plots not created');      
end

for ii=1:nImages
    FCL{ii}(isnan(FaceCorrL{ii}))=NaN;
    FCR{ii}(isnan(FaceCorrR{ii}))=NaN;
end
%%
hf=figure; hold all;
hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';

ii=1;
hp1=imagesc(repmat(ImSet{ii},1,1,3)); hold on; axis ij
hp2=gpatch(FL,PointsL{1},FCL{ii},'none',0.5);
hp3=gpatch(FR,PointsR{1},FCR{ii},'none',0.5);
pbaspect([size(ImSet{ii},2) size(ImSet{ii},1) 1])
hs1=title(['Cam ' num2str(nCam) ' frame ' num2str(1)]);
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
    PnowL=PointsL{ii};
    PnowR=PointsR{ii};
    
    cNowL=FCL{ii};
    cNowR=FCR{ii};
    
    TitleNow=['Cam ' num2str(nCam) ' frame ' num2str(ii)];
    
    %Set entries in animation structure
    animStruct.Handles{ii}=[hp1,hp2,hp2,hp3,hp3,hs1]; %Handles of objects to animate
    animStruct.Props{ii}={'CData','CData','Vertices','CData','Vertices','String'}; %Properties of objects to animate
    animStruct.Set{ii}={repmat(ImSet{ii},1,1,3),cNowL,PnowL,cNowR,PnowR,TitleNow}; %Property values for to set in order to animate
    
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
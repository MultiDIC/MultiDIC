function [] = anim8_DIC_image_3Dmeasure_points_2n(ImSet,DIC2DpairResults,DIC3DpairResults,faceMeasureString,varargin)
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
CorCoeffVecL=DIC2DpairResultsL.CorCoeffVec(nImages+1:end);
CorCoeffVecR=DIC2DpairResultsR.CorCoeffVec(1:nImages);

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
    CorCoeffCutOff=max([cell2mat(CorCoeffVecL); cell2mat(CorCoeffVecR)]);
else
    CorCoeffCutOff=optStruct.CorCoeffCutOff;
end

for ii=1:nImages
    CorCoeffVecL{ii}(CorCoeffVecL{ii}>CorCoeffCutOff)=NaN;   
    CorCoeffVecR{ii}(CorCoeffVecR{ii}>CorCoeffCutOff)=NaN;  
end


%%
switch faceMeasureString
    case {'DispMgn'}
        PCL=DIC3DpairResultsL.Disp.DispMgn;  
        PCR=DIC3DpairResultsR.Disp.DispMgn;  
        cMap='jet';
     case {'DispX'}
         for ii=1:nImages
             PCL{ii,1}=DIC3DpairResultsL.Disp.DispVec{ii}(:,1);
             PCR{ii,1}=DIC3DpairResultsR.Disp.DispVec{ii}(:,1);
         end
         cMap='jet';
    case {'DispY'}
        for ii=1:nImages
            PCL{ii,1}=DIC3DpairResultsL.Disp.DispVec{ii}(:,2);
            PCR{ii,1}=DIC3DpairResultsR.Disp.DispVec{ii}(:,2);
        end
        cMap='jet';
    case {'DispZ'}
        for ii=1:nImages
            PCL{ii,1}=DIC3DpairResultsL.Disp.DispVec{ii}(:,3);
            PCR{ii,1}=DIC3DpairResultsR.Disp.DispVec{ii}(:,3);
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
        PCmin=min([min(PCL{ii}(~isnan(CorCoeffVecL{ii}))) min(PCL{ii}(~isnan(CorCoeffVecL{ii}))) PCmin]);
        PCmax=max([max(PCR{ii}(~isnan(CorCoeffVecR{ii}))) max(PCR{ii}(~isnan(CorCoeffVecR{ii}))) PCmax]);
    end
    PClimits=[PCmin PCmax];
else
    PClimits=optStruct.PClimits;
end

%%
hf=figure; hold all;
hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';

ii=1;
hp1=imagesc(repmat(ImSet{ii},1,1,3)); hold on; axis ij
PL=PointsL{ii}(~isnan(CorCoeffVecL{ii}),:);
PR=PointsR{ii}(~isnan(CorCoeffVecR{ii}),:);
hp2=scatter(PL(:,1),PL(:,2),6,PCL{ii}(~isnan(CorCoeffVecL{ii})),'+');
hp3=scatter(PR(:,1),PR(:,2),6,PCR{ii}(~isnan(CorCoeffVecR{ii})),'+');
pbaspect([size(ImSet{ii},2) size(ImSet{ii},1) 1])
hs1=title(['Cam ' num2str(nCam) ' frame ' num2str(1)]);
colormap(cMap);
hc1=colorbar;
caxis(PClimits)
title(hc1, faceMeasureString);
hc1.FontSize=16;
axis off
drawnow

%Create the time vector
animStruct.Time=linspace(0,1,nImages);

for ii=1:nImages
    xNowL=PointsL{ii}(~isnan(CorCoeffVecL{ii}),1);
    yNowL=PointsL{ii}(~isnan(CorCoeffVecL{ii}),2);
    xNowR=PointsR{ii}(~isnan(CorCoeffVecR{ii}),1);
    yNowR=PointsR{ii}(~isnan(CorCoeffVecR{ii}),2);
    
    cNowL=PCL{ii};
    cNowR=PCR{ii};
    
    TitleNow=['Cam ' num2str(nCam) ' frame ' num2str(ii)];
    
    %Set entries in animation structure
    animStruct.Handles{ii}=[hp1,hp2,hp2,hp2,hp3,hp3,hp3,hs1]; %Handles of objects to animate
    animStruct.Props{ii}={'CData','XData','YData','CData','XData','YData','CData','String'}; %Properties of objects to animate
    animStruct.Set{ii}={repmat(ImSet{ii},1,1,3),xNowL,yNowL,cNowL,xNowR,yNowR,cNowR,TitleNow}; %Property values for to set in order to animate
    
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
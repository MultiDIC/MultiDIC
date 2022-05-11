function []=anim8_DIC3DPP_faceMeasure(DIC3DPPresults,faceMeasureString,RBMlogic,varargin)
%% function for plotting 3D-DIC results of face measures in STEP4.
% plotting 3D surfaces from camera pairs, animation changing
% with time, and the faces colored according to faceMeasureString
% this function is called in plotMultiDICPairResults
%
% Options:
% anim8_DIC3DPP_faceMeasure(D(DIC3DPPresults,faceMeasureString,RBMlogic)
% anim8_DIC3DPP_faceMeasure((DIC3DPPresults,faceMeasureString,RBMlogic,optStruct)
%
% Inputs:
% * DIC3DAllPairsResults
% * faceMeasureString: can be any of the following:
%   'J','Emgn','emgn','Epc1','Epc2','epc1','epc2','dispMgn','dispX','dispY','dispZ','FaceColors','FaceIsoInd','pairInd','Lamda1','Lamda2'
% * optStruct: optional structure for plotting options which may include any of the following fields:
%   - smoothLogic: logical variable for smoothing (true)/not smoothing (false) the face measure
%   - FaceAlpha: transparacy of the faces (scalar between 0 and 1, where zero is transparent and 1 is opaque)
%   - colorBarLimits: a 2x1 scalar vector for the colobar limits. if not set, it's automatic
%   - dataLimits: a 2x1 scalar vector for the data limits of the face measure. if a face measure is outside these limits, it is set to NaN. if not set no face is set to NaN
%   - colorMap
%   - zDirection: 1 for z up and -1 for z down
%   - lineColor: line color for the mesh. can be for example 'b','k','none',etc...
%   - supTitleString=faceMeasureString;

%% Assign plot options
Narg=numel(varargin);
switch Narg
    case 1
        optStruct=varargin{1};
    case 0
        optStruct=struct;
    otherwise
        ('wrong number of input arguments');
end

% complete the struct fields
if ~isfield(optStruct,'smoothLogic')
    optStruct.smoothLogic=0;
end
if ~isfield(optStruct,'dataLimits')
    optStruct.dataLimits=[-inf inf];
end
if ~isfield(optStruct,'zDirection') % 1 or -1
    optStruct.zDirection=1;
end

if ~isfield(optStruct,'maxCorrCoeff')
    optStruct.maxCorrCoeff=[];
end

%%
nFrames=numel(DIC3DPPresults.Points3D);

%% Assign the right face measure into FC
FC=cell(nFrames,1);
switch faceMeasureString
    
    case {'FaceColors'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.(faceMeasureString);
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[0 255];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='gray';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Color texture from image';
    case 'FacePairInds'
        for it=1:nFrames
            FC{it}=DIC3DPPresults.FacePairInds;
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[1 max(DIC3DPPresults.FacePairInds)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='gjet';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=0.8;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='k';
        end
        faceMeasureTitle='Camera-pair index';
    case {'J'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            Jmax=max(abs(FCmat(:)-1));
            if Jmax>1
                optStruct.colorBarLimits=[0 2];
            else
                optStruct.colorBarLimits=[1-prctile(abs(FCmat(:)-1),100) 1+prctile(abs(FCmat(:)-1),100)];
            end
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap=0.8*coldwarm;
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Surface Area Change (Dilatation J)';
    case {'Lamda1'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            Lmax=max(abs(FCmat(:)-1));
            if Lmax>1
                optStruct.colorBarLimits=[0 2];
            else
                optStruct.colorBarLimits=[1-prctile(abs(FCmat(:)-1),100) 1+prctile(abs(FCmat(:)-1),100)];
            end
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap=0.8*coldwarm;
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='1st principal stretch';
    case {'Lamda2'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            Lmax=max(abs(FCmat(:)-1));
            if Lmax>1
                optStruct.colorBarLimits=[0 2];
            else
                optStruct.colorBarLimits=[1-prctile(abs(FCmat(:)-1),100) 1+prctile(abs(FCmat(:)-1),100)];
            end
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap=0.8*coldwarm;
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='2nd principal stretch';
    case {'Emgn'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[0 prctile(FCmat(:),100)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='parula';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Lagrangian strain norm';
    case {'emgn'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[0 prctile(FCmat(:),100)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='parula';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Eulerian strain norm';
    case {'Eeq'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[0 prctile(FCmat(:),100)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='parula';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Equivalent Lagrangian strain';
    case {'eeq'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[0 prctile(FCmat(:),100)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='parula';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Equivalent Eulerian strain';
    case {'EShearMax'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[0 prctile(FCmat(:),100)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='parula';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Max Lagrangian shear strain';
    case {'eShearMax'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[0 prctile(FCmat(:),100)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='parula';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Max Eulerian shear strain';
    case {'Epc1'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[-prctile(abs(FCmat(:)),100) prctile(abs(FCmat(:)),100)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap=0.8*coldwarm;
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='1st principal Lagrangian strain';
    case {'Epc2'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[-prctile(abs(FCmat(:)),100) prctile(abs(FCmat(:)),100)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap=0.8*coldwarm;
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='2nd principal Lagrangian strain';
    case {'epc1'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[-prctile(abs(FCmat(:)),100) prctile(abs(FCmat(:)),100)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap=0.8*coldwarm;
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='1st principal Eulerian strain';
    case {'epc2'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.Deform.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[-prctile(abs(FCmat(:)),100) prctile(abs(FCmat(:)),100)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap=0.8*coldwarm;
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='2nd principal Eulerian strain';
    case {'DispMgn'}
        for it=1:nFrames
            Fnow=DIC3DPPresults.Faces;
            if RBMlogic
                dispNow=DIC3DPPresults.Disp.DispMgn_ARBM{it}; % point measure
            else
                dispNow=DIC3DPPresults.Disp.DispMgn{it}; % point measure
            end
            FC{it}=mean(dispNow(Fnow),2); % turn into face measure
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[0 prctile(FCmat(:),100)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='parula';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Displacement magnitude';
    case {'DispX'}
        for it=1:nFrames
            Fnow=DIC3DPPresults.Faces;
            if RBMlogic
                dispNow=DIC3DPPresults.Disp.DispVec_ARBM{it}(:,1); % point measure
            else
                dispNow=DIC3DPPresults.Disp.DispVec{it}(:,1); % point measure
            end
            FC{it}=mean(dispNow(Fnow),2); % turn into face measure
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[min(FCmat(:)) max(FCmat(:))];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='parula';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='X Displacement';
    case {'DispY'}
        for it=1:nFrames
            Fnow=DIC3DPPresults.Faces;
            if RBMlogic
                dispNow=DIC3DPPresults.Disp.DispVec_ARBM{it}(:,2); % point measure
            else
                dispNow=DIC3DPPresults.Disp.DispVec{it}(:,2); % point measure
            end
            FC{it}=mean(dispNow(Fnow),2); % turn into face measure
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[min(FCmat(:)) max(FCmat(:))];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='parula';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Y Displacement';
    case {'DispZ'}
        for it=1:nFrames
            Fnow=DIC3DPPresults.Faces;
            if RBMlogic
                dispNow=DIC3DPPresults.Disp.DispVec_ARBM{it}(:,3); % point measure
            else
                dispNow=DIC3DPPresults.Disp.DispVec{it}(:,3); % point measure
            end
            FC{it}=mean(dispNow(Fnow),2); % turn into face measure
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[min(FCmat(:)) max(FCmat(:))];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='parula';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Z Displacement';
    case {'FaceIsoInd'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[0 1];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='parula';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Triangular face isotropy (regularity) index ';
    case {'FaceCorrComb'}
        for it=1:nFrames
            FC{it}=DIC3DPPresults.(faceMeasureString){it};
            if ~isempty(optStruct.maxCorrCoeff)
                corrNow=DIC3DPPresults.FaceCorrComb{it};
                FC{it}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            end
        end
        FCmat = cell2mat(FC);
        if ~isfield(optStruct,'colorBarLimits')
            optStruct.colorBarLimits=[0 prctile(FCmat(:),100)];
        end
        colorBarLogic=1;
        if ~isfield(optStruct,'colorMap')
            optStruct.colorMap='parula';
        end
        if ~isfield(optStruct,'FaceAlpha')
            optStruct.FaceAlpha=1;
        end
        if ~isfield(optStruct,'lineColor') % 'none' or 'k'
            optStruct.lineColor='none';
        end
        faceMeasureTitle='Combined correlation coefficient';
        
    otherwise
        error('unexpected face measure string. plots not created');
        
end

%% Plot
[xl,yl,zl]=axesLimits(DIC3DPPresults.Points3D);

animStruct=struct;

hf=cFigure;
hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';

axisGeom;
ax=gca;
ax.CameraUpVector=[0 0 optStruct.zDirection];

colormap(optStruct.colorMap);
if colorBarLogic
    if strcmp(faceMeasureString,'FacePairInds')
        icolorbar(optStruct.colorBarLimits);
    else
        colorbar;
        if optStruct.colorBarLimits(2)>optStruct.colorBarLimits(2)
            caxis(optStruct.colorBarLimits);
        else
            optStruct.colorBarLimits=[-eps eps];
        end
    end
end

suptitle(faceMeasureTitle);
% axis off
% camlight headlight

it=1;
Fnow=DIC3DPPresults.Faces;
if RBMlogic
    Pnow=DIC3DPPresults.Points3D_ARBM{it};
else
    Pnow=DIC3DPPresults.Points3D{it};
end
CFnow=FC{it};
if optStruct.smoothLogic
    [CFnow]=patchSmoothFaceMeasure(Fnow,Pnow,CFnow);
end
CFnow(CFnow<optStruct.dataLimits(1))=NaN;
CFnow(CFnow>optStruct.dataLimits(2))=NaN;
hp=gpatch(Fnow,Pnow,CFnow,optStruct.lineColor,optStruct.FaceAlpha); hold on


h_ax=gca;
h_ax.XLim = xl; h_ax.YLim = yl; h_ax.ZLim = zl;

animStruct.Time=1:nFrames;
animStruct.Handles=cell(1,nFrames);
animStruct.Props=cell(1,nFrames);
animStruct.Set=cell(1,nFrames);


for it=1:nFrames
    animStruct.Handles{it}=[];
    animStruct.Props{it}=cell(1,2);
    animStruct.Set{it}=cell(1,2);
    
    Fnow=DIC3DPPresults.Faces;
    if RBMlogic
        Pnow=DIC3DPPresults.Points3D_ARBM{it};
    else
        Pnow=DIC3DPPresults.Points3D{it};
    end
    CFnow=FC{it};
    if optStruct.smoothLogic
        [CFnow]=patchSmoothFaceMeasure(Fnow,Pnow,CFnow);
    end
    CFnow(CFnow<optStruct.dataLimits(1))=NaN;
    CFnow(CFnow>optStruct.dataLimits(2))=NaN;
    animStruct.Handles{it}=[animStruct.Handles{it} hp hp]; %Handles of objects to animate
    animStruct.Props{it}{1}='CData';
    animStruct.Props{it}{2}='Vertices'; %Properties of objects to animate
    animStruct.Set{it}{1}=CFnow;
    animStruct.Set{it}{2}=Pnow; %Property values for to set in order to animate
    
    h_ax.XLim = xl; h_ax.YLim = yl; h_ax.ZLim = zl;
    
end
anim8(hf,animStruct);

addColorbarLimitsButton(hf);
addColormapButton(hf);
addEdgeColorButton(hf);
addFaceAlphaButton(hf);
addLightButton(hf);
addAmbientStrengthButton(hf);
addDiffuseStrengthButton(hf);
addSpecularStrengthButton(hf);
addFaceLightingButton(hf);


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
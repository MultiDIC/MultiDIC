function []=anim8_DIC3DPP_pointMeasure(DIC3DPPresults,pointMeasureString,RBMlogic,varargin)
%% function for plotting 3D-DIC results of point measures in STEP4.
% plotting 3D points from camera pairs, animation changing
% with time, and the points colored according to pointMeasureString
% this function is called in plotMultiDICPairResults
%
% Options:
% anim8_DIC3DPP_pointMeasure(DIC3DPPresults,pointMeasureString,RBMlogic)
% anim8_DIC3DPP_pointMeasure(DIC3DPPresults,pointMeasureString,RBMlogic,optStruct)
%
% Inputs:
% * DIC3DPPresults
% * pointMeasureString: can be any of the following:
%   'corrComb','DispMgn','dispMgn','dispX','dispY','dispZ','pairInd'
% * optStruct: optional structure for plotting options which may include any of the following fields:
%   - smoothLogic: logical variable for smoothing (true)/not smoothing (false) the face measure
%   - alphaVal: transparacy of the faces (scalar between 0 and 1, where zero is transparent and 1 is opaque)
%   - colorBarLimits: a 2x1 scalar vector for the colobar limits. if not set, it's automatic
%   - dataLimits: a 2x1 scalar vector for the data limits of the face measure. if a face measure is outside these limits, it is set to NaN. if not set no face is set to NaN
%   - colorMap
%   - zDirection: 1 for z up and -1 for z down
%   - lineColor: line color for the mesh. can be for example 'b','k','none',etc...
%   - supTitleString=faceMeasureString;

%%
Narg=numel(varargin);
switch Narg
    case 1
        optStruct=varargin{1};
    case 0
        optStruct=struct;
    otherwise
        ('wrong number of input arguments');
end
%%
nFrames=numel(DIC3DPPresults.Points3D);

%% Assign the right point measure into PC
PC=cell(nFrames);
switch pointMeasureString
    case {'corrComb'}
        for it=1:nFrames
            PC{it}=DIC3DPPresults.corrComb{it};
        end
        PCmat = cell2mat(PC);
        PCmax=max(PCmat(:));
        PClimits=[0 PCmax];
        colorBarLogic=1;
        cMap='parula';
        pointMeasureTitle='Combined Correlation Coefficient';
    case {'DispMgn'}
        for it=1:nFrames
            if RBMlogic
                PC{it}=DIC3DPPresults.Disp.DispMgn_ARBM{it};
            else
                PC{it}=DIC3DPPresults.Disp.DispMgn{it};
            end
        end
        PCmat = cell2mat(PC);
        PClimits=[0 max(PCmat(:))];
        colorBarLogic=1;
        cMap='parula';
        pointMeasureTitle='Displacement Magnitude';
    case {'DispX'}
        for it=1:nFrames
            if RBMlogic
                PC{it}=DIC3DPPresults.Disp.DispVec_ARBM{it}(:,1);
            else
                PC{it}=DIC3DPPresults.Disp.DispVec{it}(:,1);
            end
        end
        PCmat = cell2mat(PC);
        PClimits=[min(PCmat(:)) max(PCmat(:))];
        colorBarLogic=1;
        cMap='parula';
        pointMeasureTitle='X Displacement';
    case {'DispY'}
        for it=1:nFrames
            if RBMlogic
                PC{it}=DIC3DPPresults.Disp.DispVec_ARBM{it}(:,2);
            else
                PC{it}=DIC3DPPresults.Disp.DispVec{it}(:,2);
            end
        end
        PCmat = cell2mat(PC);
        PClimits=[min(PCmat(:)) max(PCmat(:))];
        colorBarLogic=1;
        cMap='parula';
        pointMeasureTitle='Y Displacement';
    case {'DispZ'}
        for it=1:nFrames
            if RBMlogic
                PC{it}=DIC3DPPresults.Disp.DispVec_ARBM{it}(:,3);
            else
                PC{it}=DIC3DPPresults.Disp.DispVec{it}(:,3);
            end
        end
        PCmat = cell2mat(PC);
        PClimits=[min(PCmat(:)) max(PCmat(:))];
        colorBarLogic=1;
        cMap='parula';
        pointMeasureTitle='Z Displacement';
    case {'pairInd'}
        PCtemp=DIC3DPPresults.PointPairInds;
        for it=1:nFrames
            PC{it}=PCtemp;
        end
        PClimits=[1 max(PCtemp)];
        colorBarLogic=1;
        cMap='gjet';
        pointMeasureTitle='Camera-pair index';
    otherwise
        error('unexpected face measure string. plots not created');
        
end

%% Assign plot options


% complete the struct fields
if ~isfield(optStruct,'smoothLogic')
    optStruct.smoothLogic=0;
end
if ~isfield(optStruct,'colorBarLimits')
    optStruct.colorBarLimits=PClimits;
end
if ~isfield(optStruct,'dataLimits')
    optStruct.dataLimits=PClimits;
end
if ~isfield(optStruct,'colorMap')
    optStruct.colorMap=cMap;
end
if ~isfield(optStruct,'zDirection') % 1 or -1
    optStruct.zDirection=1;
end
if ~isfield(optStruct,'supTitleString')
    optStruct.supTitleString=pointMeasureTitle;
end
if ~isfield(optStruct,'maxCorrCoeff')
    optStruct.maxCorrCoeff=[];
end

%% Plot
[xl,yl,zl]=axesLimits(DIC3DPPresults.Points3D);

animStruct=struct;

hf=cFigure;
hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';

axisGeom;

suptitle(optStruct.supTitleString);
% axis off
% camlight headlight

it=1;
if RBMlogic
    Pnow=DIC3DPPresults.Points3D_ARBM{it};
else
    Pnow=DIC3DPPresults.Points3D{it};
end

if ~isempty(optStruct.maxCorrCoeff)
    corrNow=DIC3DPPresults.corrComb{it};
    Pnow(corrNow>optStruct.maxCorrCoeff,:)=NaN;
end
Cnow=PC{it};

%%%%%%% add point smoothing
%     if optStruct.smoothLogic
%         [Cnow]=triSmoothFaceMeasure(Cnow,Fnow,Pnow,[],[]);
%     end

Cnow(Cnow<optStruct.dataLimits(1))=NaN;
Cnow(Cnow>optStruct.dataLimits(2))=NaN;

hp=scatter3(Pnow(:,1),Pnow(:,2),Pnow(:,3),[],Cnow,'+'); hold on

axisGeom;

colormap(cMap);
if colorBarLogic
    if strcmp(pointMeasureString,'pairInd')
        icolorbar(optStruct.colorBarLimits);
    else
        colorbar;
        caxis(optStruct.colorBarLimits);
    end
end

h_ax=gca;
h_ax.XLim = xl; h_ax.YLim = yl; h_ax.ZLim = zl;
h_ax.CameraUpVector=[0 0 optStruct.zDirection];

animStruct.Time=1:nFrames;
animStruct.Handles=cell(1,nFrames);
animStruct.Props=cell(1,nFrames);
animStruct.Set=cell(1,nFrames);

for it=1:nFrames
    animStruct.Handles{it}=[];
    animStruct.Props{it}=cell(1,4);
    animStruct.Set{it}=cell(1,4);
    
    if RBMlogic
        Pnow=DIC3DPPresults.Points3D_ARBM{it};
    else
        Pnow=DIC3DPPresults.Points3D{it};
    end
    
    if ~isempty(optStruct.maxCorrCoeff)
        corrNow=DIC3DPPresults.corrComb{it};
        Pnow(corrNow>optStruct.maxCorrCoeff,:)=NaN;
    end
    
    Cnow=PC{it};
    Cnow(Cnow<optStruct.dataLimits(1))=NaN;
    Cnow(Cnow>optStruct.dataLimits(2))=NaN;
    
    animStruct.Handles{it}=[animStruct.Handles{it} hp hp hp hp]; %Handles of objects to animate
    
    animStruct.Props{it}{1}='XData';
    animStruct.Props{it}{2}='YData';
    animStruct.Props{it}{3}='ZData';
    animStruct.Props{it}{4}='CData';
    
    animStruct.Set{it}{1}=Pnow(:,1);
    animStruct.Set{it}{2}=Pnow(:,2);
    animStruct.Set{it}{3}=Pnow(:,3);
    animStruct.Set{it}{4}=Cnow;
    %         %Property values for to set in order to animate
    
    h_ax.XLim = xl; h_ax.YLim = yl; h_ax.ZLim = zl;
    
end
anim8(hf,animStruct);

addColorbarLimitsButton(hf);
addColormapButton(hf); % only applicable for plots without designated colors

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
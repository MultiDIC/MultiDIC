function []=anim8_DIC_3D_pairs_pointMeasure(DIC3DAllPairsResults,pointMeasureString,varargin)
%% function for plotting 3D-DIC results of point measures in STEP3.
% plotting 3D points from camera pairs, animation changing
% with time, and the points colored according to pointMeasureString
% this function is called in plotMultiDICPairResults
%
% Options:
% anim8_DIC_3D_pairs_pointMeasure(DIC_3Dallpairs_results,pointMeasureString)
% anim8_DIC_3D_pairs_pointMeasure(DIC_3Dallpairs_results,pointMeasureString,optStruct)
% 
% Inputs:
% * DIC3DAllPairsResults
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
%   - TitleString=faceMeasureString;

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
nPairs=numel(DIC3DAllPairsResults);
nFrames=numel(DIC3DAllPairsResults{1}.Points3D);

%% Assign the right face measure into FC
PC=cell(nPairs,nFrames);
switch pointMeasureString
    case {'corrComb'}
        for ip=1:nPairs
            for it=1:nFrames
                PC{ip,it}=DIC3DAllPairsResults{ip}.corrComb{it};
            end
        end
        PCmat = cell2mat(PC);
        PCmax=max(PCmat(:));
        PClimits=[0 PCmax];
        colorBarLogic=1;
        cMap='parula';
        legendLogic=0;
    case {'DispMgn'}
        for ip=1:nPairs
            for it=1:nFrames
                PC{ip,it}=DIC3DAllPairsResults{ip}.Disp.DispMgnTransformed{it};
            end
        end
        PCmat = cell2mat(PC);
        PClimits=[0 max(PCmat(:))];
        colorBarLogic=1;
        cMap='parula';
        legendLogic=0;
    case {'DispX'}
        for ip=1:nPairs
            for it=1:nFrames
                PC{ip,it}=DIC3DAllPairsResults{ip}.Disp.DispVecTransformed{it}(:,1);
            end
        end
        PCmat = cell2mat(PC);
        PClimits=[min(PCmat(:)) max(PCmat(:))];
        colorBarLogic=1;
        cMap='parula';
        legendLogic=0;
    case {'DispY'}
        for ip=1:nPairs
            for it=1:nFrames
                PC{ip,it}=DIC3DAllPairsResults{ip}.Disp.DispVecTransformed{it}(:,2);
            end
        end
        PCmat = cell2mat(PC);
        PClimits=[min(PCmat(:)) max(PCmat(:))];
        colorBarLogic=1;
        cMap='parula';
        legendLogic=0;
    case {'DispZ'}
        for ip=1:nPairs
            for it=1:nFrames
                PC{ip,it}=DIC3DAllPairsResults{ip}.Disp.DispVecTransformed{it}(:,3);
            end
        end
        PCmat = cell2mat(PC);
        PClimits=[min(PCmat(:)) max(PCmat(:))];
        colorBarLogic=1;
        cMap='parula';
        legendLogic=0;
    case {'pairInd'}
        for ip=1:nPairs
            for it=1:nFrames
                PC{ip,it}=ip*ones(size(DIC3DAllPairsResults{ip}.Points3D{it},1),1);
            end
        end
        PClimits=[1 nPairs];
        colorBarLogic=0;
        cMap='gjet';
        legendLogic=1;
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
if ~isfield(optStruct,'TitleString')
    optStruct.TitleString=pointMeasureString;
end
if ~isfield(optStruct,'maxCorrCoeff')
    optStruct.maxCorrCoeff=[];
end

%% Plot
DIC_3Djoined_results=joinPairs(DIC3DAllPairsResults);
[xl,yl,zl]=axesLimits(DIC_3Djoined_results.Points3D);

animStruct=struct;

hf=cFigure;
hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';

axisGeom;

gtitle(optStruct.TitleString,20);
% axis off
% camlight headlight

it=1;
for ip=1:nPairs
    Pnow=DIC3DAllPairsResults{ip}.Points3Dtransformed{it};
    if ~isempty(optStruct.maxCorrCoeff)
        corrNow=DIC3DAllPairsResults{ip}.corrComb{it};
        Pnow(corrNow>optStruct.maxCorrCoeff,:)=NaN;
    end
    Cnow=PC{ip,it};
    
    %%%%%%% add point smoothing
    %     if optStruct.smoothLogic
    %         [Cnow]=triSmoothFaceMeasure(Cnow,Fnow,Pnow,[],[]);
    %     end
    
    Cnow(Cnow<optStruct.dataLimits(1))=NaN;
    Cnow(Cnow>optStruct.dataLimits(2))=NaN;
    
    hp(ip)=scatter3(Pnow(:,1),Pnow(:,2),Pnow(:,3),[],Cnow,'+'); hold on
    
    axisGeom;
    
end

colormap(cMap);
if colorBarLogic
    cbh=colorbar;
    caxis(optStruct.colorBarLimits);
end

if legendLogic
    pairStrings=arrayfun(@num2str,1:nPairs,'unif',0);
    legendStrings=strcat('pair',pairStrings);
    legend(legendStrings);
end

h_ax=gca;
h_ax.XLim = xl; h_ax.YLim = yl; h_ax.ZLim = zl;
h_ax.CameraUpVector=[0 0 optStruct.zDirection];

animStruct.Time=1:nFrames;
animStruct.Handles=cell(1,nFrames);
animStruct.Props=cell(1,nFrames);
animStruct.Set=cell(1,nFrames);

ic=1;
for it=1:nFrames
    animStruct.Handles{ic}=[];
    animStruct.Props{ic}=cell(1,4*nPairs);
    animStruct.Set{ic}=cell(1,4*nPairs);
    
    for ip=1:nPairs
        Pnow=DIC3DAllPairsResults{ip}.Points3Dtransformed{it};
        if ~isempty(optStruct.maxCorrCoeff)
            corrNow=DIC3DAllPairsResults{ip}.corrComb{it};
            Pnow(corrNow>optStruct.maxCorrCoeff,:)=NaN;
        end
        
        Cnow=PC{ip,it};
        
        Cnow(Cnow<optStruct.dataLimits(1))=NaN;
        Cnow(Cnow>optStruct.dataLimits(2))=NaN;
        
        animStruct.Handles{ic}=[animStruct.Handles{ic} hp(ip) hp(ip) hp(ip) hp(ip)]; %Handles of objects to animate
        
        animStruct.Props{ic}{4*ip-3}='XData';
        animStruct.Props{ic}{4*ip-2}='YData';
        animStruct.Props{ic}{4*ip-1}='ZData';
        animStruct.Props{ic}{4*ip-0}='CData';
        
        animStruct.Set{ic}{4*ip-3}=Pnow(:,1);
        animStruct.Set{ic}{4*ip-2}=Pnow(:,2);
        animStruct.Set{ic}{4*ip-1}=Pnow(:,3);
        animStruct.Set{ic}{4*ip-0}=Cnow;
        %         %Property values for to set in order to animate
    end
    
    h_ax.XLim = xl; h_ax.YLim = yl; h_ax.ZLim = zl;
    
    ic=ic+1;
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
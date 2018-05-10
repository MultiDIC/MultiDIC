function []=anim8_DIC_3D_pairs_faceMeasure_direction(DIC3DAllPairsResults,faceMeasureString,varargin)
%% function for plotting 3D-DIC results of face measures as color + direction as arrows in STEP3.
% plotting 3D surfaces from camera pairs, animation changing
% with time, and the faces colored according to faceMeasureString 
% this function is called in plotMultiDICPairResults
%
% Options:
% anim8_DIC_3D_pairs_faceMeasure(DIC_3Dallpairs_results,faceMeasureString)
% anim8_DIC_3D_pairs_faceMeasure(DIC_3Dallpairs_results,faceMeasureString,optStruct)
% 
% Inputs:
% * DIC3DAllPairsResults
% * faceMeasureString: can be any of the following:
%   'Epc1','Epc2','epc1','epc2','Lamda1','Lamda2', or a combination of two  in a cell array to plot side by side
% * optStruct: optional structure for plotting options which may include any of the following fields:
%   - smoothLogic: logical variable for smoothing (true)/not smoothing (false) the face measure 
%   - alphaVal: transparacy of the faces (scalar between 0 and 1, where zero is transparent and 1 is opaque) 
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
if ~isfield(optStruct,'alphaVal')
    optStruct.alphaVal=1;
end
if ~isfield(optStruct,'dataLimits')
    optStruct.dataLimits=[-inf inf];
end
if ~isfield(optStruct,'zDirection') % 1 or -1
    optStruct.zDirection=1;
end
if ~isfield(optStruct,'lineColor') % 'none' or 'k'
    optStruct.lineColor='none';
end
if ~isfield(optStruct,'maxCorrCoeff')
    optStruct.maxCorrCoeff=[];
end

%%
nFrames=numel(DIC3DAllPairsResults{1}.Points3D);

quiverScaleFactor=20;

DIC_3Djoined_results=joinPairs(DIC3DAllPairsResults);
[xl,yl,zl]=axesLimits(DIC_3Djoined_results.Points3Dtransformed);

%% Assign the right face measure into FC

if iscell(faceMeasureString)
    faceMeasureCell=faceMeasureString;
    nStrains=numel(faceMeasureString);
    switch nStrains
        case 1
            FC=cell(nFrames,nStrains);
            D=cell(nFrames,nStrains);
            Ds=cell(nFrames,nStrains);
            Vc=cell(nFrames,nStrains);
        case 2
            FC=cell(nFrames,nStrains);
            D=cell(nFrames,nStrains);
            Ds=cell(nFrames,nStrains);
            Vc=cell(nFrames,nStrains);
        otherwise
            error('wrong number of cells in faceMeasureString cell array (must be 1 or 2)');
    end
elseif ischar(faceMeasureString)
    nStrains=1;
    faceMeasureCell=cell(1);
    faceMeasureCell{1}=faceMeasureString;
    FC=cell(nFrames,1);
    D=cell(nFrames,1);
    Ds=cell(nFrames,1);
    Vc=cell(nFrames,1);
else
    error('wrong face measure (second input variable)');
end

if ~isfield(optStruct,'supTitleString') % 'none' or 'k'
    optStruct.supTitleString=faceMeasureCell;
end

directionStringCell=cell(nStrains,1);
for is=1:nStrains
    switch faceMeasureCell{is}
        case 'Epc1'
            directionStringCell{is}='Epc1vecCur';
        case 'Epc2'
            directionStringCell{is}='Epc2vecCur';
        case 'epc1'
            directionStringCell{is}='epc1vec';
        case 'epc2'
            directionStringCell{is}='epc2vec';
        case 'Lamda1'
            directionStringCell{is}='Epc1vecCur';
        case 'Lamda2'
            directionStringCell{is}='Epc2vecCur';
        otherwise
            error('unexpected face measure string. plots not created');
    end
    for it=1:nFrames
        FC{it,is}=DIC_3Djoined_results.Deform.(faceMeasureCell{is}){it}; % face color (strain)
        D{it,is}=DIC_3Djoined_results.Deform.(directionStringCell{is}){it}; % direction (unit vector)
        Vc{it,is}=DIC_3Djoined_results.FaceCentroidsTransformed{it};
        if ~isempty(optStruct.maxCorrCoeff)
            corrNow=DIC_3Djoined_results.FaceCorrComb{it};
            FC{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            D{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
            Vc{it,is}(corrNow>optStruct.maxCorrCoeff,:)=NaN;
        end
        switch faceMeasureCell{is}
            case {'Epc1','Epc2','epc1','epc2'}
                Ds{it,is}=quiverScaleFactor*FC{it,is}.*D{it,is}; % direction with magnitude (scaled vector)
            case {'Lamda1','Lamda2'}
                Ds{it,is}=quiverScaleFactor*(FC{it,is}-1).*D{it,is}; % direction with magnitude (scaled vector)
        end

    end
end
FCmat = cell2mat(FC);
if ~isfield(optStruct,'colorBarLimits')
    switch faceMeasureCell{is}
        case {'Epc1','Epc2','epc1','epc2'}
            optStruct.colorBarLimits=[-prctile(abs(FCmat(:)),100) prctile(abs(FCmat(:)),100)];
        case {'Lamda1','Lamda2'}
            optStruct.colorBarLimits=[1-prctile(abs(FCmat(:)-1),100) 1+prctile(abs(FCmat(:)-1),100)];
    end
end
colorBarLogic=1;
if ~isfield(optStruct,'colorMap')
    optStruct.colorMap=.8*coldwarm;
end


%% Plot

animStruct=struct;

hf=cFigure;
hf.Units='normalized'; hf.OuterPosition=[.05 .05 .9 .9]; hf.Units='pixels';

for is=1:nStrains
    subplot(1,nStrains,is);
    
    axisGeom;     
    ax=gca;     
    ax.CameraUpVector=[0 0 optStruct.zDirection];
    colormap(optStruct.colorMap);
    if colorBarLogic
        colorbar;
        caxis(optStruct.colorBarLimits);
    end
    
    title(optStruct.supTitleString{is});
    axis off
    % camlight headlight
    
    it=1;
    Fnow=DIC_3Djoined_results.Faces;
    Pnow=DIC_3Djoined_results.Points3Dtransformed{it};
    Vnow=Vc{it,is};
    FCnow=FC{it,is};
    Dnow=Ds{it,is};
    if optStruct.smoothLogic
        [FCnow]=triSmoothFaceMeasure(FCnow,Fnow,Pnow,[],[]);
    end
    FCnow(FCnow<optStruct.dataLimits(1))=NaN;
    FCnow(FCnow>optStruct.dataLimits(2))=NaN;

    hp(is)=gpatch(Fnow,Pnow,FCnow,optStruct.lineColor,optStruct.alphaVal); hold on
    hq(is)=quiver3(Vnow(:,1),Vnow(:,2),Vnow(:,3),Dnow(:,1),Dnow(:,2),Dnow(:,3),0,'Color',.2*[1 1 1],'ShowArrowHead','off','AutoScale','on'); hold on;
%         
    h_ax=gca;
    h_ax.XLim = xl; h_ax.YLim = yl; h_ax.ZLim = zl;
    
end

%% fill in the animstruct

animStruct.Time=1:nFrames;
animStruct.Handles=cell(1,nFrames);
animStruct.Props=cell(1,nFrames);
animStruct.Set=cell(1,nFrames);


for it=1:nFrames
    animStruct.Handles{it}=[];
    animStruct.Props{it}=cell(1,8*nStrains);
    animStruct.Set{it}=cell(1,8*nStrains);
    Pnow=DIC_3Djoined_results.Points3Dtransformed{it};
    
    for is=1:nStrains
        FCnow=FC{it,is};
        if optStruct.smoothLogic
            [FCnow]=triSmoothFaceMeasure(FCnow,Fnow,Pnow,[],[]);
        end
        FCnow(FCnow<optStruct.dataLimits(1))=NaN;
        FCnow(FCnow>optStruct.dataLimits(2))=NaN;
        
        Fnow=DIC_3Djoined_results.Faces;
        Pnow=DIC_3Djoined_results.Points3D{it};
        Vnow=Vc{it,is};
        Dnow=Ds{it,is};
    
        animStruct.Handles{it}=[animStruct.Handles{it} hp(is) hp(is) hq(is) hq(is) hq(is) hq(is) hq(is) hq(is)]; %Handles of objects to animate (add one every pair)
        
        animStruct.Props{it}{1+8*(is-1)}='CData';
        animStruct.Props{it}{2+8*(is-1)}='Vertices'; %Properties of objects to animate
        animStruct.Props{it}{3+8*(is-1)}='XData'; %Properties of objects to animate
        animStruct.Props{it}{4+8*(is-1)}='YData'; %Properties of objects to animate
        animStruct.Props{it}{5+8*(is-1)}='ZData'; %Properties of objects to animate
        animStruct.Props{it}{6+8*(is-1)}='UData'; %Properties of objects to animate
        animStruct.Props{it}{7+8*(is-1)}='VData'; %Properties of objects to animate
        animStruct.Props{it}{8+8*(is-1)}='WData'; %Properties of objects to animate
        
        animStruct.Set{it}{1+8*(is-1)}=FCnow;
        animStruct.Set{it}{2+8*(is-1)}=Pnow; %Property values for to set in order to animate
        animStruct.Set{it}{3+8*(is-1)}=Vnow(:,1); %Property values for to set in order to animate
        animStruct.Set{it}{4+8*(is-1)}=Vnow(:,2); %Property values for to set in order to animate
        animStruct.Set{it}{5+8*(is-1)}=Vnow(:,3); %Property values for to set in order to animate
        animStruct.Set{it}{6+8*(is-1)}=Dnow(:,1); %Property values for to set in order to animate
        animStruct.Set{it}{7+8*(is-1)}=Dnow(:,2); %Property values for to set in order to animate
        animStruct.Set{it}{8+8*(is-1)}=Dnow(:,3); %Property values for to set in order to animate
        
        h_ax.XLim = xl; h_ax.YLim = yl; h_ax.ZLim = zl;
        
    end
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
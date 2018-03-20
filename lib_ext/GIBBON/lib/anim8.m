function hf=anim8(varargin)

%% Parse input

switch nargin
    case 2
        hf=varargin{1};
        animStruct=varargin{2};
end

%% Visualisation settings

% fontColor='w';
fontSize=15;
% cMap=gjet(250);

scrollBarWidth=50;

%% Defining slider

animTime=animStruct.Time(:);

minT=1;
maxT=numel(animTime);
w=maxT-minT;
sliceIndexI=1; %Initial index
tickSizeMajor_I=ceil(w/20);

%% Initialize display

figure(hf); drawnow;

%Initialize slider
jSlider = javax.swing.JSlider(minT,maxT);
javacomponent(jSlider,[0,0,round(hf.Position(3)),scrollBarWidth]);
set(jSlider, 'MajorTickSpacing',tickSizeMajor_I, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@updateViewFunc,{hf,jSlider}},'Orientation',jSlider.HORIZONTAL);

%% Set resize function

set(hf,'ResizeFcn',{@setScrollSizeFunc,{hf,scrollBarWidth,jSlider}});

%% Initialize figure callbacks

set(hf,'KeyPressFcn', {@figKeyPressFunc,{hf}});%'WindowButtonDownFcn', {@figMouseDown,{hf}},'WindowButtonUpFcn', {@mouseup,hf});

%% Initialise buttons

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
iconPath=fullfile(toolboxPath,'icons');

% hb = findall(hf,'Type','uitoolbar');
hb = uitoolbar(hf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Help button

%get icon
D=importdata(fullfile(iconPath,'help.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Help','CData',S,'Tag','help_button','ClickedCallback',@helpFunc,'Separator','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Play button

%get icon 1
D=importdata(fullfile(iconPath,'play.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
iconPlay=S;

%get icon 2
D=importdata(fullfile(iconPath,'stop.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
iconStop=S;

% Create a uitoggletool in the toolbar
hPlay=uitoggletool(hb,'TooltipString','Play','CData',iconPlay,'Tag','play_button','ClickedCallback',{@playFunc,{hf}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time button

%get icon
D=importdata(fullfile(iconPath,'time.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Time','CData',S,'Tag','time_button','ClickedCallback',{@timeFunc,{hf}}); %,'Separator','on'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cycle button

%get icon
D=importdata(fullfile(iconPath,'cycle.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uitoggletool in the toolbar
hCycle=uitoggletool(hb,'TooltipString','Cycle forward-backward','CData',S,'Tag','cycle_button');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save animation button

%get icon
D=importdata(fullfile(iconPath,'saveAnimation.jpg'));
S=double(D);
S=S-min(S(:));
S=S./max(S(:));
S(S==1)=NaN;
if size(S,3)==1
    S=repmat(S,[1 1 3]);
end
% Create a uipushtool in the toolbar
uipushtool(hb,'TooltipString','Save .gif animation','CData',S,'Tag','saveAnimation_button','ClickedCallback',{@saveAnimationFunc,{hf}}); %,'Separator','on'

%% Text fields

% title(sprintf('%6.16e ',T),'FontSize',hf.UserData.anim8.fontSize);
% title(sprintf('%f',T),'FontSize',hf.UserData.anim8.fontSize);

hTextTime = uicontrol(hf,'Style','text',...
    'String',[' Time: ',sprintf('%f',animTime(1))],...
    'Position',[0 scrollBarWidth round(hf.Position(3)) round(scrollBarWidth/2)],'BackgroundColor',[1 1 1],'HorizontalAlignment','Left'); %hf.Color

%% Set figure UserData

hf.UserData.anim8.sliderHandles={jSlider};
% hf.UserData.anim8.ButtonHandles.Sample=hSample;
hf.UserData.anim8.animStruct=animStruct;
hf.UserData.anim8.fontSize=fontSize;
% animHandles=animStruct.Handles; %Handles of objects to animate
% animProps=animStruct.Props; %Properties of objects to animate
% animSets=animStruct.Set; %Property values for to set in order to animate
% animTime=animStruct.Time;
hf.UserData.anim8.playDir=1;
hf.UserData.anim8.ButtonHandles.Play=hPlay;
hf.UserData.anim8.ButtonHandles.hCycle=hCycle;
hf.UserData.anim8.ButtonHandles.hTextTime=hTextTime;
hf.UserData.anim8.pauseTime=3/numel(animStruct.Time);%mean(diff(animStruct.Time));
hf.UserData.anim8.shiftMag=1;%ceil(numel(animStruct.Time)/20);
hf.UserData.anim8.icons.play=iconPlay;
hf.UserData.anim8.icons.stop=iconStop;

% Store current settings
hf.UserData.anim8.WindowButtonDownFcn=hf.WindowButtonDownFcn;
hf.UserData.anim8.WindowButtonUpFcn=hf.WindowButtonUpFcn;
hf.UserData.anim8.KeyPressFcn=hf.KeyPressFcn;
hf.UserData.anim8.WindowScrollWheelFcn=hf.WindowScrollWheelFcn;
hf.UserData.anim8.BusyAction=hf.BusyAction;

%Export figure widget settings
hf.UserData.efw.defaultPath=fullfile(cd,'efw');
hf.UserData.efw.imName=['figure',num2str(get(hf,'Number'))];
hf.UserData.efw.imExt='png';
hf.UserData.efw.imRes='50';
hf.UserData.efw.exportFigOpt='-nocrop';
hf.UserData.efw.exportGifOpt='1';

%% Initialize slider locations
set(jSlider,'Value',sliceIndexI);

%%
drawnow;
end

%% Scroll bar resizing

function setScrollSizeFunc(~,~,inputCell)
hf=inputCell{1};
scrollBarWidth=inputCell{2};

jSlider=inputCell{3};
javacomponent(jSlider,[0,0,round(hf.Position(3)), scrollBarWidth]);

set(hf.UserData.anim8.ButtonHandles.hTextTime,'Position',[0 scrollBarWidth round(hf.Position(3)) round(scrollBarWidth/2)]);
end

%% Help

function helpFunc(~,~)

msgText={'* Play button or press the space bar -> Start animation',...
    '* Stop button or press the space bar -> Stop animation',...
    '* Time button -> Change time stepping settings',...
    '* Cycle button -> Toggle forward / forward-backward mode',...
    '* Save button -> Export snap shorts (e.g. jpg) and animated .gif files',...
    '* Press the v key to activate the view control widget (vcw), note that vcw changes key press functions until vcw is deactivated',...
    };
helpdlg(msgText,'Help information');

end

%% play

function playFunc(~,~,inputCell)
hf=inputCell{1};
shiftMag=hf.UserData.anim8.shiftMag;

set(hf.UserData.anim8.ButtonHandles.Play,'CData',hf.UserData.anim8.icons.stop,'TooltipString','Stop');

while strcmp(get(hf.UserData.anim8.ButtonHandles.Play,'State'),'on')
    tic
    if strcmp(get(hf.UserData.anim8.ButtonHandles.hCycle,'State'),'on')
        
        jSlider=hf.UserData.anim8.sliderHandles{1};
        
        sliderValue=get(jSlider,'Value');
        sliderValueNew=sliderValue+(shiftMag*hf.UserData.anim8.playDir);
        
        sliderMax=get(jSlider,'Maximum');
        sliderMin=get(jSlider,'Minimum');
        
        if sliderValueNew<sliderMin
            hf.UserData.anim8.playDir=hf.UserData.anim8.playDir*-1;
            sliderValueNew=sliderValue+(shiftMag*hf.UserData.anim8.playDir);
        elseif sliderValueNew>sliderMax
            hf.UserData.anim8.playDir=hf.UserData.anim8.playDir*-1;
            sliderValueNew=sliderValue+(shiftMag*hf.UserData.anim8.playDir);
        end
        
        set(jSlider,'Value',sliderValueNew);
        
    else
        shiftSlider(hf.UserData.anim8.sliderHandles{1},shiftMag);
    end
    
    
    drawnow;
    
    t=toc;
    if (hf.UserData.anim8.pauseTime-t)>0
        pause(hf.UserData.anim8.pauseTime-t);
    end
end

set(hf.UserData.anim8.ButtonHandles.Play,'CData',hf.UserData.anim8.icons.play,'TooltipString','Play');

end

%% Time
function timeFunc(~,~,inputCell)
hf=inputCell{1};

prompt = {'Enter pause time:','Enter step size:'};
dlg_title = 'Time stepping settings';
defaultOptions = {num2str(hf.UserData.anim8.pauseTime),num2str(hf.UserData.anim8.shiftMag)};
s=40+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);
Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);
if ~isempty(Q)
    if ~isempty(Q{1})
        hf.UserData.anim8.pauseTime=str2double(Q{1});
    end
    if ~isempty(Q{2})
        hf.UserData.anim8.shiftMag=str2double(Q{2});
    end
end

end

%% Save Animation

function saveAnimationFunc(~,~,inputCell)
hf=inputCell{1};
jSlider=hf.UserData.anim8.sliderHandles{1};

%%

defStruct=hf.UserData.efw;
prompt = {'Save path (leave empty to browse to desired folder instead):',...
    'Image name:','Image extension (i.e. png, jpg, bmp, or tif):',...
    'Image resolution (e.g. 120):',...
    'Extra export_fig options (comma seperated, no spaces e.g. -nocrop,-transparent,-painters):',...
    'Export gif option'};
dlg_title = 'Export Gif Widget (see: help efw and help export_fig)';
defaultOptions = {defStruct.defaultPath,defStruct.imName,defStruct.imExt,defStruct.imRes,defStruct.exportFigOpt,hf.UserData.efw.exportGifOpt};

s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);

Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);

if ~isempty(Q)
    if isempty(Q{1})
        Q{1}=uigetdir(defStruct.defaultPath,'Select save path');
        if Q{1}==0
            return;
        end
    end
    
    if isempty(Q{2})
        error('Empty input. Please enter a file name');
    end
    
    if ~exist(Q{1},'dir') %create output folder if it does not exist already
        mkdir(Q{1});
    end
    
    if isempty(Q{3})
        warning('Empty input. No image format specified, using jpg files.');
        Q{3}='jpg';
    end
    
    fileName=fullfile(Q{1},Q{2});
    exportGifCell{1,1}=fileName;
    
    stringSet=Q{3}; %The image extension
    stringNoSpaces=regexprep(stringSet,'[^\w'']',''); %Remove potential extra spaces
    
    if ~strcmp(stringNoSpaces(1),'-') %If first character is not '-'
        stringNoSpaces=['-',stringNoSpaces]; %Add '-' to start, e.g. 'jpg' becomes '-jpg'
    end
    
    %Check format validaty and keep if valid
    if any(strcmp(stringNoSpaces,{'-png','-jpg','-tiff','-bmp'}))
        exportGifCell{1,end+1}=stringNoSpaces; %Add to input list
    else
        error('Wrong image format requested');
    end
    
    figRes=['-r',Q{4}];
    exportGifCell{1,end+1}=figRes;
    
    if ~isempty(Q{5})
        stringSet=Q{5}; %The set of potentially multiple options
        stringSetSep = strsplit(stringSet,',');
        for q=1:1:numel(stringSetSep)
            stringNoSpaces=regexprep(stringSetSep{q},'[^\w'']',''); %Remove potential extra spaces
            if ~strcmp(stringNoSpaces(1),'-') %If first character is not '-'
                stringNoSpaces=['-',stringNoSpaces]; %Add '-' to start, e.g. 'jpg' becomes '-jpg'
            end
            exportGifCell{1,end+1}=stringNoSpaces; %Add to input list
        end
    end
    
    if ~isempty(Q{6})
        exportGifOpt=Q{6};
    end
    
    fileNameGif=exportGifCell{1,1};
    exportGifCellSub=exportGifCell;
    
    c=1;
    stepRange=1:hf.UserData.anim8.shiftMag:numel(hf.UserData.anim8.animStruct.Time);    
    for q=stepRange
        set(jSlider,'Value',q);
        fileNameNow=[fileNameGif,'_',num2str(q)];
        exportGifCellSub{1,1}=fileNameNow;
        figure(hf);
        export_fig(exportGifCellSub{:});
        gifStruct.FileNames{c}=[fileNameNow,'.',exportGifCell{1,2}(2:end)];
        c=c+1;
    end
    
    if strcmp(exportGifOpt,'1')
        %Add reverse path
        if strcmp(get(hf.UserData.anim8.ButtonHandles.hCycle,'State'),'on')
            numFiles=numel(gifStruct.FileNames);
            if numFiles>2
                for q=(numFiles-1):-1:2
                    gifStruct.FileNames{end+1}=gifStruct.FileNames{q};
                end
            end
        end
        
        gifStruct.DelayTime=hf.UserData.anim8.pauseTime;
        gifStruct.FileNameGif=fileNameGif;
        
        exportGif(gifStruct);
    end
    
    %Override defaults
    defStruct.defaultPath=Q{1};
    defStruct.imName=Q{2};
    defStruct.imExt=Q{3};
    defStruct.imRes=Q{4};
    defStruct.exportFigOpt=Q{5};
    defStruct.efw.exportGifOpt=Q{6};
    
    hf.UserData.efw=defStruct;
    
end

end

%% Figure key press

function figKeyPressFunc(~,eventData,inputCell)

hf=inputCell{1}; %Figure handle

% step = 1;
% if ismember('shift', eventData.Modifier)
%     step = -step; %Make negative while shift is down
% end
%
% if ismember('control', eventData.Modifier)
%     step = step * 4; %Increase speed
% end

% Key input options
switch eventData.Key
    case {'leftarrow','downarrow'}
        shiftSlider(hf.UserData.anim8.sliderHandles{1},-1*hf.UserData.anim8.shiftMag);
    case {'rightarrow','uparrow'}
        shiftSlider(hf.UserData.anim8.sliderHandles{1},1*hf.UserData.anim8.shiftMag);
    case 'v' % Activate vcw
        set(hf.UserData.cFigure.Handles.vcw,'State','On');
    case'space'
        if strcmp(get(hf.UserData.anim8.ButtonHandles.Play,'State'),'off')
            set(hf.UserData.anim8.ButtonHandles.Play,'State','on');
            playFunc([],[],{hf});
        else
            set(hf.UserData.anim8.ButtonHandles.Play,'State','off');
        end
end

end

%% updateViewFunc

function updateViewFunc(~,~,inputCell)

hf=inputCell{1};
jSlider=inputCell{2};
animStruct=hf.UserData.anim8.animStruct;

sliderValue=get(jSlider,'Value');

T=animStruct.Time(sliderValue);
H=animStruct.Handles{sliderValue};%e.g. [hp,hp]; %Handles of objects to animate
P=animStruct.Props{sliderValue};% e.g. {'FaceColor','Vertices'}; %Properties of objects to animate
S=animStruct.Set{sliderValue};%e.g.{c(q,:),V}; %Property values for to set in order to animate

set(hf.UserData.anim8.ButtonHandles.hTextTime,'String',[' Time: ',sprintf('%f',T)],'FontSize',hf.UserData.anim8.fontSize);

for q=1:1:numel(H)
    h=H(q); % Current graphics handle
    p=P{q}; % Current graphics property
    s=S{q}; % Current property setting
    h.(p)=s;% Setting the property
end

end

%% Shift slider

function shiftSlider(jSlider,shiftMag)

sliderValue=get(jSlider,'Value');
sliderValueNew=sliderValue+shiftMag;

sliderMax=get(jSlider,'Maximum');
sliderMin=get(jSlider,'Minimum');

if sliderValueNew<sliderMin
    set(jSlider,'Value',sliderMax);
elseif sliderValueNew>sliderMax
    set(jSlider,'Value',1);
else
    set(jSlider,'Value',sliderValueNew);
end

end

 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2017  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

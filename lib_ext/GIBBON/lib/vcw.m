function [varargout]=vcw(varargin)

% function vcw(hf,buttonOpt)
% ------------------------------------------------------------------------
% vcw, View Control Widget, Allows users to manipulate a view in 3D using
% key and button presses
%
%   vcw(hf,buttonOpt)
%
% Allows the user to rotate, pan and zoom a figure using key presses and
% mouse gestures. Additionally, press q to quit the widget, r to reset the
% axes and escape to close the figure. This function is non-blocking, but
% fixes axes aspect ratios.
%
% IN:
%   hf - Handle of the figure to be manipulated (default: gcf).
%   buttonOpt - 4x1 cell array indicating the function to associate with
%             each mouse button (left to right) and the scroll action.
%             Functions can be any of:
%                'rot' - Rotate about x and y axes of viewer's coordinate
%                        frame
%                'rotz' - Rotate about z axis of viewer's coordinate frame
%                'zoom' - Zoom (change canera view angle)
%                'zoomz' - Move along z axis of viewer's coordinate frame
%                'pan' - Pan
%                '' - Don't use that button
%             Default: {'pan','rot','zoomz','zoomz'};).
%
% This code was inspired by the fcw function by Oliver Woodford (which was
% based on Torsten Vogel's view3d function, which was in turn inspired by
% rotate3d from The MathWorks, Inc.).
%
% See modification log below to see how vcw function differs from the fcw
% function. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2015/04/15 %Copied from fcw and renamed to vcw due to planned revision
% 2015/04/15 %Added: 1) handing of colorbars (bug in fcw when view(2) is
% used combined with panning which induced zooming and panning), 2) overobj
% axes selection so that the current axes is determined based on mouse
% pointer location for most functions, 3) A toggle button for activation
% and deactivation in the figure toolbar, 4) ability to start vcw before
% objects are plotted, 5) "proper" closure of the vcw widget, in fcw the q
% button did not exit the keyDown functions such as panning etc. Now the
% quit action deactivates the widget, 6) Uppon activation of the vcw widget
% the plotting and default view manipulation tools and buttons are disabled
% (to avoid interference with vcw), 7) Added "linked" mode by using ALT
% button to alter views for all axes in figure uppon keypress, 8) Altered
% keypress functions and behaviour with SHIFT, also added i to display help
% information for the vcw function.
% 2015/04/20 Added to GIBBON toolbox
% 2015/04/22 Added JavaFrame handling of ALT related mnemonics
% 2015/04/28 Fixed behaviour for repated vcw; commands (only generate a
% single vcw button even if vcw is called multiple times). 
% 2015/04/28 Fixed behaviour for figures without axes. I.e. vcw will only
% start if an axis is present. 
% 2016/01/13 Added that clipping is turned off
%
% TO DO: 1) Improved handling of colorbars. Currently requires colorbar
% locations to be set to 'manual' for vcw. However this causes the figure
% to rescale/adjust after deactivation/activation of vcw. It would be best
% if the colorbar locations settings could remain constant. 2) Proper
% restoring of all figure properties. Currently defaults are set manually
% which could remove user defined figure features.
% 2) Turn back the clipping type after exciting vcw
%------------------------------------------------------------------------

%% Parse input arguments
switch nargin
    case 0
        hf = gcf;
        buttonOpt = {'pan','rot','zoomz','zoomz'};
    case 1
        hf=varargin{1};
        buttonOpt = {'pan','rot','zoomz','zoomz'};
    case 2
        hf=varargin{1};
        buttonOpt=varargin{2};
    otherwise
        error('Wrong number of input arguments');
end

if ~ishandle(hf)
    buttonOpt=hf;
    hf = gcf;
end

%% Check axis limits

h = findobj(hf, 'Type', 'axes', '-depth', 1)'; %All axis handles
if ~isempty(h)
    for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
        axis(h);
        
        xLim=get(h,'xlim');
        yLim=get(h,'ylim');
        zLim=get(h,'zlim');
        
        wx=abs(diff(xlim));
        wy=abs(diff(ylim));
        wz=abs(diff(zlim));
        
        w_max=max([wx wy wz]);
        min_w=0.1;
        if w_max<min_w
            w_max=min_w;
        end
        
        w_min=w_max/10;
        if w_min<min_w
            w_min=min_w;
        end
        w_add=[-w_min w_min]/2;       
        
        if wx<w_min
            set(h,'xlim',xLim+w_add);
        end
        
        if wy<w_min
            set(h,'ylim',yLim+w_add);
        end
        
        if wz<w_min
            set(h,'zlim',zLim+w_add);
        end              
    end
end


%% Initialise button and button/keypress wait
hb = findall(hf,'Type','uitoolbar');

%Check for presence of a vcw button
hp = findobj(hb,'Tag','tBar');
if isempty(hp) %If vcw button is not present create one and wait for key/button press
    
    % Build icon
    s=[ NaN,NaN,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;...
        NaN,NaN,NaN,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;...
        NaN,NaN,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,NaN,NaN;...
        NaN,1  ,1  ,1  ,NaN,1  ,NaN,NaN,NaN,NaN,1  ,1  ,1  ,NaN,NaN,NaN;...
        NaN,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,NaN,NaN;...
        1,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,NaN;...
        1,1  ,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,1  ,NaN,NaN,1  ,1  ,NaN;...
        1,1  ,NaN,NaN,NaN,1  ,1  ,NaN,NaN,1  ,1  ,NaN,NaN,NaN,1  ,1  ;...
        1,1  ,NaN,NaN,NaN,NaN,1  ,NaN,NaN,1  ,NaN,NaN,NaN,NaN,1  ,1  ;...
        NaN,1  ,1  ,NaN,NaN,NaN,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,1  ,1  ;...
        NaN,1  ,1  ,NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,NaN,NaN,NaN,1  ,1  ;...
        NaN,NaN,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,NaN;...
        NaN,NaN,NaN,1  ,1  ,1  ,NaN,NaN,NaN,NaN,1  ,NaN,1  ,1  ,1  ,NaN;...
        NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,NaN,NaN;...
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,NaN,NaN,NaN;...
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,NaN,NaN];
    
    S=zeros(16,16,3);
    S(:,:,1)=0.7.*s;
    S(:,:,2)=0.25.*s;
    S(:,:,3)=0.05.*s;
    
    % Create a uipushtool in the toolbar
    hp=uitoggletool(hb,'TooltipString','Activate View Control Widget (or enter v)','CData',S,'Tag','tBar','Separator','on');
    set(hp,'OnCallback',{@start_vcw_toggle,{hf,buttonOpt,hp}});
    set(hp,'OffCallback',{@quit_vcw_toggle,{hf,buttonOpt,hp}});

    %% Wait for start using key-press
    
    set(hf,'KeyPressFcn', {@keyPress_wait,buttonOpt,hp,hf},'BusyAction','cancel');
    
end

%% Outputs
switch nargout
    case 1
        varargout{1}=hp;
end

end

%%
function start_vcw_toggle(hObject,callbackdata,indputCell)
start_vcw(indputCell{1},indputCell{2},indputCell{3});
end

function keyPress_wait(src,eventData,buttonOpt,hp,hf)

cax = overobj2('axes');

if isempty(cax)
%     cax=gca; %this gets current axis or if none exists creates one
cax = get(hf, 'CurrentAxes');
else
    axes(cax)
end

if isempty(cax)
    return
end

%Key actions
switch eventData.Key
    case {'v'} %Start vcw        
        start_vcw(hf,buttonOpt,hp);
end
end
%%
function start_vcw(hf,buttonOpt,hp)

% Store current settings
hf.UserData.WindowButtonDownFcn=hf.WindowButtonDownFcn;
hf.UserData.WindowButtonUpFcn=hf.WindowButtonUpFcn;
hf.UserData.KeyPressFcn=hf.KeyPressFcn;
hf.UserData.WindowScrollWheelFcn=hf.WindowScrollWheelFcn;
hf.UserData.BusyAction=hf.BusyAction;

checkAxisLimits(hf);

set(hp,'State','On');
set(hp,'TooltipString','Dectivate View Control Widget (or enter v)');

% Clear any visualization modes we might be in
pan(hf, 'off');
zoom(hf, 'off');
rotate3d(hf, 'off');

%Quick fix for colorbars
H=findobj(hf,'Type','colorbar'); %Handle set
figUserDataStruct=get(hf,'UserData');
if isempty(figUserDataStruct)
    figUserDataStruct.colorbarLocSet=get(H,'Location');
    set(hf,'UserData',figUserDataStruct);
end
colorbarLocSet(hf,'manual');

% Disable Plottools Buttons and Exploration Buttons
initialState.toolbar = findobj(allchild(hf),'flat','Type','uitoolbar');
if ~isempty(initialState.toolbar)
    initialState.ptButtons = [uigettool(initialState.toolbar,'Plottools.PlottoolsOff'), ...
        uigettool(initialState.toolbar,'Plottools.PlottoolsOn'),...
        uigettool(initialState.toolbar,'Exploration.Rotate'), ...
        uigettool(initialState.toolbar,'Exploration.Pan'),...
        uigettool(initialState.toolbar,'Exploration.ZoomIn'),...
        uigettool(initialState.toolbar,'Exploration.ZoomOut'),...
        ];
    initialState.ptState = get (initialState.ptButtons,'Enable');
    set (initialState.ptButtons,'Enable','off');
end

% For each set of axes
cax = get(hf,'CurrentAxes');
% cax=gca; %this gets current axis or if none exists creates one
if isempty(cax)
    set(hp,'State','Off');
    set(hp,'TooltipString','Activate View Control Widget (or enter v)');
    return
end
h = findobj(hf, 'Type', 'axes', '-depth', 1)'; %All axis handles

if ~isempty(h)
    for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
        % Set everything to manual
        set(h, 'CameraViewAngleMode', 'manual', 'CameraTargetMode', 'manual', 'CameraPositionMode', 'manual');
        % Store the camera viewpoint
        axes(h); axis vis3d;
        caxUserDataStruct.defaultView=camview(h);
        set(h, 'UserData',caxUserDataStruct);
        %Turn clipping off
        set(h,'Clipping','off');   
    end
    axes(cax);
else
    set(hp,'State','Off');
    set(hp,'TooltipString','Activate View Control Widget (or enter v)');
    return
end

% Initialize the callbacks
set(hf, 'WindowButtonDownFcn', {@mousedown, {str2func(['vcw_' buttonOpt{1}]), str2func(['vcw_' buttonOpt{2}]), str2func(['vcw_' buttonOpt{3}])},hf}, ...
    'WindowButtonUpFcn', {@mouseup,hf}, ...
    'KeyPressFcn', {@keypress,buttonOpt,hp,hf}, ...
    'WindowScrollWheelFcn', {@scroll, str2func(['vcw_' buttonOpt{4}])}, ...
    'BusyAction', 'cancel');

end

%%
function keypress(src, eventData,buttonOpt,hp,hf)

cax = overobj2('axes');
if isempty(cax)
    cax = get(hf, 'CurrentAxes');
else
    axes(cax)
end
if isempty(cax)
    return;
end

% checkAxisLimits(hf);

step = 1;
if ismember('shift', eventData.Modifier)
    step = -step; %Make negative while shift is down
end

if ismember('control', eventData.Modifier)
    step = step * 4; %Increase speed
end

mnemOff=1;
if ismember('alt', eventData.Modifier)
    linkedOn=1;
    % Try to turn off menu Mnemonics
    try
        warning off; %Stop jframe warning
        jFrame = get(handle(hf),'JavaFrame');
        jMenuBar=jFrame.fHG2Client.getMenuBar;
        for q=0:1:jMenuBar.getComponentCount-1
            jComp=jMenuBar.getComponent(q);
            jComp.setMnemonic(' ');
        end
        mnemOff=1;
        warning on;
    catch
        %Remove toolbar when ALT is pressed (QUICK FIX)
        mnemOff=0;
        set(hf,'MenuBar','none');
        t = uitoolbar;
        set(t,'Tag','emptyBar_vcw');
    end
else
    linkedOn=0;
end

% Key input options
switch eventData.Key
    case 'leftarrow' % Pan left
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_pan([], [step 0], h, hf);
            end
        else
            vcw_pan([], [step 0], cax, hf);
        end
    case 'rightarrow' % Pan right
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_pan([], [-step 0], h, hf);
            end
        else
            vcw_pan([], [-step 0], cax, hf);
        end
    case 'downarrow' % Pan down
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_pan([], [0 step], h, hf);
            end
        else
            vcw_pan([], [0 step], cax, hf);
        end
    case 'uparrow' % Pan up
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_pan([], [0 -step], h, hf);
            end
        else
            vcw_pan([], [0 -step], cax, hf);
        end
    case 'x' % Rotate around x
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_rot([], [0 step], h, hf);
            end
        else
            vcw_rot([], [0 step], cax, hf);
        end
    case 'y' % Rotate around y
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_rot([], [step 0], h, hf);
            end
        else
            vcw_rot([], [step 0], cax, hf);
        end
    case 'z' % Rotate around z
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_rotz([], [0 step], h, hf);
            end
        else
            vcw_rotz([], [0 step], cax, hf);
        end
    case 'm' % Magnify/zoom (positive or negative)
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                vcw_zoom([], [0 -step],h, hf);
            end
        else
            vcw_zoom([], [0 -step], cax, hf);
        end
    case 't' % top view
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(0,90);
            end
        else
            view(0,90);
        end
    case 'b' % bottom view
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(0,-90);
            end
        else
            view(0,-90);
        end
    case 'f' % front view
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(-90,0);
            end
        else
            view(-90,0);
        end
    case 'h' % back view
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(90,0);
            end
        else
            view(90,0);
        end
    case 'l' % left view
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(0,0);
            end
        else
            view(0,0);
        end
    case 'r' % right view
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                axes(h);view(-180,0);
            end
        else
            view(-180,0);
        end
    case 's' % Store current views as default/initial views (return to it using 'd')
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                caxUserDataStruct=get(h,'UserData');
                caxUserDataStruct.defaultView=camview(h);
                set(h,'UserData',caxUserDataStruct);
            end
        else
            caxUserDataStruct=get(cax,'UserData');
            caxUserDataStruct.defaultView=camview(cax);
            set(cax,'UserData',caxUserDataStruct);
        end
    case 'd' % Reset all the axes to default
        if linkedOn==1
            for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
                caxUserDataStruct=get(h,'UserData');
                camview(h,caxUserDataStruct.defaultView);
            end
        else
            caxUserDataStruct=get(cax,'UserData');
            camview(cax,caxUserDataStruct.defaultView);
        end
    case 'i' % Show help
        msgText={'Key input options:',...
            '------------------------------------------------------------------------',...
            'Arrow keys = panning',...
            'x = Rotate around viewer x-axis',...
            'y = Rotate around viewer y-axis',...
            'z = Rotate around viewer z-axis',...
            'm = Zoom/magnify',...
            'Hold down SHIFT to change direction of change for key input based rotation/pan/zoom',...
            'Hold down CTRL to use x4 speed of change  for key input based rotation/pan/zoom',...
            'Hold down ALT to link manipulations for all figure axes for key input based rotation/pan/zoom',...
            'f,h,t,b,l,r = Set front, hind, top, bottom, left, or right view respectively',...
            's = Store current view states as default (return to default using d)',...
            'd = Restore/reset to default view',...
            'v = activate/deactivate vcw mode',...
            'i = Display help information',...
            'The key inputs and mouse scroll work in the axis defined by mouse pointer location (overobj)',...
            'The mouse inputs (other than scroll) work in current axis (e.g. gca) irrespective of point location',...
            '------------------------------------------------------------------------',...
            };
        helpButton = questdlg(msgText,'Help for vcw','OK','OK');
    case 'v' % Quit vcw mode
        quit_vcw(hf,buttonOpt,hp);
    otherwise
%         quit_vcw(hf,buttonOpt,hp);
end

if mnemOff==0
    set(hf,'MenuBar','figure');
    h1 = findobj(hf,'Tag','emptyBar_vcw');
    if ~isempty(h1)
        delete(h1);
    end
end
end

%%
function mousedown(src, eventData, funcs, hf)

% Get the button pressed
% cax = overobj2('axes');

cax = get(hf, 'CurrentAxes');
if isempty(cax)
    return;
end

% checkAxisLimits(hf);

switch get(hf, 'SelectionType')
    case 'extend' % Middle button
        method = funcs{2};
    case 'alt' % Right hand button
        method = funcs{3};
    case 'open' % Double click
        caxUserDataStruct=get(cax,'UserData');
        camview(cax,caxUserDataStruct.defaultView);
        return;
    otherwise
        method = funcs{1};
end

% Set the cursor
switch func2str(method)
    case {'vcw_zoom', 'vcw_zoomz'}
        shape=[ 2   2   2   2   2   2   2   2   2   2 NaN NaN NaN NaN NaN NaN  ;
            2   1   1   1   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN  ;
            2   1   2   2   2   2   2   2   2   2 NaN NaN NaN NaN NaN NaN  ;
            2   1   2   1   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN  ;
            2   1   2   1   1   1   1   1   2 NaN NaN NaN NaN NaN NaN NaN  ;
            2   1   2   1   1   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN  ;
            2   1   2   1   1   1   1   1   2 NaN NaN NaN   2   2   2   2  ;
            2   1   2   1   1   2   1   1   1   2 NaN   2   1   2   1   2  ;
            2   1   2   1   2 NaN   2   1   1   1   2   1   1   2   1   2  ;
            2   2   2   2 NaN NaN NaN   2   1   1   1   1   1   2   1   2  ;
            NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   1   2   1   2  ;
            NaN NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   2   1   2  ;
            NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   1   2   1   2  ;
            NaN NaN NaN NaN NaN NaN   2   2   2   2   2   2   2   2   1   2  ;
            NaN NaN NaN NaN NaN NaN   2   1   1   1   1   1   1   1   1   2  ;
            NaN NaN NaN NaN NaN NaN   2   2   2   2   2   2   2   2   2   2  ];
    case 'vcw_pan'
        shape=[ NaN NaN NaN NaN NaN NaN NaN   2   2 NaN NaN NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN   2   1   1   1   1   2 NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN   1   1   1   1   1   1 NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
            NaN NaN   2   1 NaN NaN   2   1   1   2 NaN NaN   1   2 NaN NaN ;
            NaN   2   1   1   2   2   2   1   1   2   2   2   1   1   2 NaN ;
            2   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2 ;
            2   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2 ;
            NaN   2   1   1   2   2   2   1   1   2   2   2   1   1   2 NaN ;
            NaN NaN   2   1 NaN NaN   2   1   1   2 NaN NaN   1   2 NaN NaN ;
            NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN   1   1   1   1   1   1 NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN   2   1   1   1   1   2 NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN NaN   2   1   1   2 NaN NaN NaN NaN NaN NaN ;
            NaN NaN NaN NaN NaN NaN NaN   2   2 NaN NaN NaN NaN NaN NaN NaN ];
    case {'vcw_rotz', 'vcw_rot'}
        % Rotate
        shape=[ NaN NaN NaN   2   2   2   2   2 NaN   2   2 NaN NaN NaN NaN NaN ;
            NaN NaN NaN   1   1   1   1   1   2   1   1   2 NaN NaN NaN NaN ;
            NaN NaN NaN   2   1   1   1   1   2   1   1   1   2 NaN NaN NaN ;
            NaN NaN   2   1   1   1   1   1   2   2   1   1   1   2 NaN NaN ;
            NaN   2   1   1   1   2   1   1   2 NaN NaN   2   1   1   2 NaN ;
            NaN   2   1   1   2 NaN   2   1   2 NaN NaN   2   1   1   2 NaN ;
            2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
            2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
            2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
            2   1   1   2 NaN NaN NaN NaN NaN NaN NaN NaN   2   1   1   2 ;
            NaN   2   1   1   2 NaN NaN   2   1   2 NaN   2   1   1   2 NaN ;
            NaN   2   1   1   2 NaN NaN   2   1   1   2   1   1   1   2 NaN ;
            NaN NaN   2   1   1   1   2   2   1   1   1   1   1   2 NaN NaN ;
            NaN NaN NaN   2   1   1   1   2   1   1   1   1   2 NaN NaN NaN ;
            NaN NaN NaN NaN   2   1   1   2   1   1   1   1   1 NaN NaN NaN ;
            NaN NaN NaN NaN NaN   2   2 NaN   2   2   2   2   2 NaN NaN NaN ];
    otherwise
        return
end

% Record where the pointer is
global VCW_POS
VCW_POS = get(0, 'PointerLocation');
% Set the cursor and callback
set(hf, 'Pointer', 'custom', 'pointershapecdata', shape, 'WindowButtonMotionFcn', {method, cax,hf});

end

%%
function mouseup(src, eventData,hf)
% Clear the cursor and callback
set(hf, 'WindowButtonMotionFcn', '', 'Pointer', 'arrow');
end

%%
function scroll(src, eventData, func, hf)
% Get the axes handle
cax = overobj2('axes');
if isempty(cax)
    cax = get(hf, 'CurrentAxes');
else
    axes(cax)
end
if isempty(cax)
    return;
end

% Call the scroll function
func([], [0 -10*eventData.VerticalScrollCount], cax);
end

%%

function d = check_vals(s, d)
% Check the inputs to the manipulation methods are valid
global VCW_POS
if ~isempty(s)
    % Return the mouse pointers displacement
    new_pt = get(0, 'PointerLocation');
    d = VCW_POS - new_pt;
    VCW_POS = new_pt;
end
end

%%

% function setLightPos(hf,cax)
% 
% hLights=findobj(cax,'Type','light');
% 
% cameraPositionsNew=cax.CameraPosition;
% 
% cameraPositionsOld=hf.UserData.CameraPosition;
% 
% a=vecnormalize(cameraPositionsOld);
% b=vecnormalize(cameraPositionsNew);
% [R]=vecAngle2Rot(acos(dot(a,b)),vecnormalize(cross(a,b)));
% 
% % [F,V,~]=quiver3Dpatch(0,0,0,a(1),a(2),a(3),[],[2 2]);
% % patch('Faces',F,'Vertices',V,'FaceColor','r');
% % 
% % c=b*R
% % 
% % [F,V,~]=quiver3Dpatch(0,0,0,c(1),c(2),c(3),[],[2 2]);
% % patch('Faces',F,'Vertices',V,'FaceColor','none','EdgeColor','g');
% 
% for q=1:1:numel(hLights)    
%     hLights(q).Position=hLights(q).Position*R';   
%     drawnow;
% end
% 
% end

% Figure manipulation functions
function vcw_rot(s, d, cax, hf)
d = check_vals(s, d);
% hf.UserData.CameraPosition=cax.CameraPosition;
try
    % Rotate XYt            
    camorbit(cax, d(1), d(2), 'camera', [0 0 1]);              
catch
    % Error, so release mouse down
    mouseup([],[],hf);
end
end

function vcw_rotz(s, d, cax, hf)    
%     hf.UserData.CameraPosition=cax.CameraPosition;    
d = check_vals(s, d);
try
    % Rotate Z
    camroll(cax, d(2));        
catch
    % Error, so release mouse down
    mouseup([],[],hf);
end
end

function vcw_zoom(s, d, cax, hf)
d = check_vals(s, d);
% Zoom
d = (1 - 0.01 * sign(d(2))) ^ abs(d(2));
try
    camzoom(cax, d);
catch
    % Error, so release mouse down
    mouseup([],[],hf);
end
end

function vcw_zoomz(s, d, cax, hf)
d = check_vals(s, d);
% Zoom by moving towards the camera
d = (1 - 0.01 * sign(d(2))) ^ abs(d(2)) - 1;
try
    camdolly(cax, 0, 0, d, 'fixtarget', 'camera');
catch
    % Error, so release mouse down
    mouseup([],[],hf);
end
end

function vcw_pan(s, d, cax, hf)
d = check_vals(s, d);
try
    % Pan
    camdolly(cax, d(1), d(2), 0, 'movetarget', 'pixels');
catch
    % Error, so release mouse down
    mouseup([],[],hf);
end
end

function colorbarLocSet(hf,locOpt)
H=findobj(hf,'Type','colorbar'); %Handle set
for q=1:1:numel(H)
    if isa(locOpt,'cell')
        set(H(q),'Location',locOpt{q});
    else
        set(H(q),'Location',locOpt);
    end
end
end


function checkAxisLimits(hf)

h = findobj(hf, 'Type', 'axes', '-depth', 1)'; %All axis handles
if ~isempty(h)
    for h = findobj(hf, 'Type', 'axes', '-depth', 1)'
        axis(h);
        
        xLim=get(h,'xlim');
        yLim=get(h,'ylim');
        zLim=get(h,'zlim');
        
        wx=abs(diff(xlim));
        wy=abs(diff(ylim));
        wz=abs(diff(zlim));
        
        w_max=max([wx wy wz]);
        min_w=1e-3;
        if w_max<min_w
            w_max=min_w;
        end
        
        w_min=w_max/10;
        if w_min<min_w
            w_min=min_w;
        end
        w_add=[-w_min w_min]/2;       
        
        if wx<w_min
            set(h,'xlim',xLim+w_add);
        end
        
        if wy<w_min
            set(h,'ylim',yLim+w_add);
        end
        
        if wz<w_min
            set(h,'zlim',zLim+w_add);
        end              
    end
    drawnow; 
end

end

%%
function quit_vcw_toggle(~,~,indputCell)
quit_vcw(indputCell{1},indputCell{2},indputCell{3})
end

%%

function quit_vcw(hf,buttonOpt,hp)

% Restore colorbar state
figUserDataStruct=get(hf,'UserData');
if isfield(figUserDataStruct,'colorbarLocSet')
    colorbarLocSet(hf,figUserDataStruct.colorbarLocSet);
end
% Enable Plottools Buttons and Exploration Buttons
initialState.toolbar = findobj(allchild(hf),'flat','Type','uitoolbar');
if ~isempty(initialState.toolbar)
    initialState.ptButtons = [uigettool(initialState.toolbar,'Plottools.PlottoolsOff'), ...
        uigettool(initialState.toolbar,'Plottools.PlottoolsOn'),...
        uigettool(initialState.toolbar,'Exploration.Rotate'), ...
        uigettool(initialState.toolbar,'Exploration.Pan'),...
        uigettool(initialState.toolbar,'Exploration.ZoomIn'),...
        uigettool(initialState.toolbar,'Exploration.ZoomOut'),...
        ];
    initialState.ptState = get (initialState.ptButtons,'Enable');
    set (initialState.ptButtons,'Enable','on');
end

%         hp=findobj(hf,'Tag','tBar');
set(hp,'State','Off');
set(hp,'TooltipString','Activate View Control Widget (or enter v)');

% Restore figure settings except for key press (if empty) to allow for
% reactivation with v key
hf.WindowButtonDownFcn=hf.UserData.WindowButtonDownFcn;
hf.WindowButtonUpFcn=hf.UserData.WindowButtonUpFcn;
if isempty(hf.UserData.KeyPressFcn)
    hf.KeyPressFcn={@keyPress_wait,buttonOpt,hp,hf};%[];%hf.UserData.KeyPressFcn;
else
    hf.KeyPressFcn=hf.UserData.KeyPressFcn;
end
hf.WindowScrollWheelFcn=hf.UserData.WindowScrollWheelFcn;
hf.BusyAction=hf.UserData.BusyAction;
hf.MenuBar='figure';

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

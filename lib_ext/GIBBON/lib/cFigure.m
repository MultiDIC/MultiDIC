function [varargout]=cFigure(varargin)

% function [h]=cFigure(figStruct)
% ------------------------------------------------------------------------
% Creates a custom figure using the input structure figStruct. The cFigure
% function provides easy control of background color, the color
% definitions, the figure window size (e.g. near maximal), and enables
% figure property cloning. It also allows users to create hidden figures
% which can be made visible for instance using the mfv command.
%
% The content of figStruct may follow all properties of a normal figure
% i.e. such that figStruct=figure. Which could lead to (amonst other
% properties):
%
% figStruct=figure
%
%   Figure (1) with properties:
%
%       Number: 10
%         Name: ''
%        Color: [0.9400 0.9400 0.9400]
%     Position: [680 558 560 420]
%        Units: 'pixels'
%
% figStruct used to be a handle in which case its use in this function
% involves the set command e.g.: set(h,'outerPosition',[a b c d]);
% For newer MATLAB versions however the cFigure function uses a different
% but equivalent syntax i.e.:
% h.outerPosition=[a b c d];
%
% Some additional fields can be added that are not normally part of the
% figure property set: ColorDef and ScreenOffset.
% ColorDef sets the color definition which is either 'white' or 'black'.
% This allows the user to select a dark background and appropriately set
% the colorscheme for it e.g. for a black background:
%         figStruct.ColorDef='black';
%         figStruct.Color='k';
% Where the Color property sets the figure background color while the
% ColorDef property sets the colorscheme used (of axes etc.).
% By default cFigure creates figures that are the full screensize but
% reduced 10% away from the edges. The spacing between the figure window
% and the screen edges is set by the figStruct.ScreenOffset property. The
% units are pixels.
%
% See also: figure, set, get, colordef, mfv, scf
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2014/11/25 %Created
% 2015/04/15 %Added vcw functionality
% 2018/02/02 %Fixed bug in relation to groot units (e.g. figure size is
% wrong if units are not pixels). 
%------------------------------------------------------------------------

%% Parse input and set defaults

%Force groot units to be pixels
graphicalRoot=groot;
grootUnits=graphicalRoot.Units;
if ~strcmp(grootUnits,'pixels')
    graphicalRoot.Units='pixels';
end

switch nargin
    case 0
        %Create defaults
        figStruct.Visible='on';
        figStruct.ColorDef='white';
        figStruct.Color='w';
%         figStruct.Colormap=gjet(250);
        screenSizeGroot = get(groot,'ScreenSize');
        figStruct.ScreenOffset=round(max(screenSizeGroot)*0.1); %i.e. figures are spaced around 10% of the sreensize from the edges
        vcwOpt={'pan','rot','zoomz','zoomz'};
        efwOpt=1;
    case 1
        %Use custom
        figStruct=varargin{1};
        
        %Use defaults where nothing is provided
        if ~isfield(figStruct,'Visible')
            figStruct.Visible='on';
        end
        
        if ~isfield(figStruct,'ColorDef')
            figStruct.ColorDef='white';
        end
        
        if ~isfield(figStruct,'Color')
            figStruct.Color='w';
        end
        
        if ~isfield(figStruct,'ScreenOffset')
            screenSizeGroot = get(groot,'ScreenSize');
            figStruct.ScreenOffset=round(max(screenSizeGroot)*0.1); %i.e. figures are spaced around 10% of the sreensize from the edges
        end
        
        if ~isfield(figStruct,'vcw')
            vcwOpt={'pan','rot','zoomz','zoomz'};
        else
            vcwOpt=figStruct.vcw;
            figStruct=rmfield(figStruct,'vcw'); %Remove field from structure array
        end
        
        if ~isfield(figStruct,'efw')
            efwOpt=1;
        else
            efwOpt=figStruct.efw;
            figStruct=rmfield(figStruct,'efw'); %Remove field from structure array
        end
        
end

% if ~isfield(figStruct,'Clipping');
%     figStruct.Clipping='off';
% end

%%

isOld=verLessThan('matlab', '8.4.0.150421 (R2014b)');

%% Create a hidden figure

hf = figure('Visible', 'off'); %create an invisible figure

%% Setcolor definition and associated defaults

hf=colordef(hf,figStruct.ColorDef); %Update figure handle
figStruct=rmfield(figStruct,'ColorDef'); %Remove field from structure array

%% Set figure size

if isfield(figStruct,'ScreenOffset')
    screenSizeGroot = get(groot,'ScreenSize');
    screenSizeGroot=screenSizeGroot(3:4); % width, height
    figSizeEdgeOffset=figStruct.ScreenOffset; % Figure offset from border
    figSize=screenSizeGroot-figSizeEdgeOffset; % width, height
    
    if isOld
        set(hf,'units','pixels');
        set(hf,'outerPosition',[(screenSizeGroot(1)-figSize(1))/2 (screenSizeGroot(2)-figSize(2))/2 figSize(1) figSize(2)]); % left bottom width height
    else
        hf.Units='pixels';
        hf.Position=[(screenSizeGroot(1)-figSize(1))/2 (screenSizeGroot(2)-figSize(2))/2 figSize(1) figSize(2)]; % left bottom width height
    end
    
    figStruct=rmfield(figStruct,'ScreenOffset'); %Remove field from structure array
end

%% Parse remaining figure properties

% Note: This is where figure becomes visible if figStruct.Visible='on'

fieldSet = fieldnames(figStruct); % Cell containing all structure field names
for q=1:1:numel(fieldSet)
    fieldNameCurrent=fieldSet{q};
    try
        if isOld
            set(hf,fieldNameCurrent,figStruct.(fieldNameCurrent));
        else
            hf.(fieldNameCurrent)=figStruct.(fieldNameCurrent);
        end
    catch errorMsg
        rethrow(errorMsg); %likely false option
    end
end

%% Check for activation of vcw

if isa(vcwOpt,'cell') %Allow enabling of vcw mode    
    hp=vcw(hf,vcwOpt);
    hf.UserData.cFigure.Handles.vcw=hp;
end

%% Check for activation of efw
if efwOpt
    efw; 
end

%%
if nargout>0
    varargout{1}=hf;
end

%%
% Reset groot units if a change was needed
if ~strcmp(grootUnits,'pixels')
    graphicalRoot.Units=grootUnits;
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
% Copyright (C) 2018  Kevin Mattheus Moerman
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

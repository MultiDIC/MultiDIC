function [varargout]=icolorbar(varargin)

% function [h]=icolorbar(cLim,h)
%-------------------------------------------------------------------------
% 
%
% Change log: 
% 2018/05/10 Fixed bug in relation to older matlab versions which do not
% support h.Colormap=c for axis specific color maps. Older versions use
% colormap(h,c) instead. 
% 2019/08/14 Added colorbar handle output
%-------------------------------------------------------------------------

%% parse input

switch nargin
    case 0 
        h=[]; 
        cLim=caxis;
    case 1
        cLim=varargin{1};
        h=[];
    case 2
        cLim=varargin{1};
        h=varargin{2};
end

if isempty(cLim)
    cLim=caxis;
end

if isempty(h)
    %Get current figure or open a new one
    if isempty(findobj('type','figure'))
        hf=cFigure; %Open a new figure
    else
        hf=gcf; %Get the current figure
    end    
   figure(hf);       
   if isempty(hf.CurrentAxes) %isempty(get(hf,'CurrentAxes'))
       h=axes;
   else
       h=gca;
   end
else    
    switch class(h)
        case {'matlab.ui.Figure'} %A figure handle
            hf=h;
            if isempty(hf.CurrentAxes) %isempty(get(hf,'CurrentAxes'))
                h=axes;
            else
                h=gca;
            end         
        case {'matlab.graphics.axis.Axes'} %An axes handle            
            hf=h.Parent; %Get figure handle
        otherwise
            error('Handle does not represent an axis or figure handle');
    end
end

%% Get colormap
if isprop(h,'Colormap') %Available in MATLAB 2018  
    c=h.Colormap;
    matlabNew=1;
else %Older MATLAB versions
    c=colormap(h);
    matlabNew=0;
end

%% Set limits and tick labels
caxis([cLim(1)-0.5 cLim(2)+0.5]);
hc=colorbar; 
hc.Ticks=cLim(1):1:cLim(2);
cn=resampleColormap(c,numel(hc.Ticks));
if matlabNew==1 %axis colormap
    h.Colormap=cn;
else %Figure colormap
    colormap(h,cn);
end

%% Collect output
switch nargout
    case 1
        varargout{1}=h;
    case 2
        varargout{1}=h;
        varargout{2}=hc;
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
% Copyright (C) 2019  Kevin Mattheus Moerman
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

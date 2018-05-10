function varargout = roipolyFit(varargin)
% This is a modification of the original MATLAB roipoly function, with these slight
% modifications:
% * when using imshow, the size is constrained to fit the figure window
% size
% used in STEP1 to select the polygon mask on calibration images
%
%ROIPOLY Select polygonal region of interest.
%   Use ROIPOLY to select a polygonal region of interest within an
%   image. ROIPOLY returns a binary image that you can use as a mask for
%   masked filtering.
%   
%   BW = ROIPOLY creates an interactive polygon tool, associated with the
%   image displayed in the current figure, called the target image. You
%   place the tool interactively, using the mouse to specify the vertices
%   of the polygon. Double-click to add a final vertex to the polygon and
%   close the polygon. Right-click to close the polygon without adding a
%   vertex. You can adjust the position of the polygon and individual
%   vertices in the polygon by clicking and dragging.
%
%   To add new vertices, position the pointer along an edge of the polygon
%   and press the "A" key. The pointer changes shape. Left-click to add a
%   vertex at the specified position.
%
%   After positioning and sizing the polygon, create the mask image, BW, by
%   either double-clicking over the polygon or choosing Create Mask from
%   the tool's context menu. BW is a binary image the same size as I with
%   0's outside the region of interest and 1's inside.
%
%   To delete the polygon, press Backspace, Escape or Delete, or choose the
%   Cancel option from the context menu. If the polygon is deleted, all
%   return values are set to empty.
%
%   BW = ROIPOLY(I) displays the image I and creates an interactive polygon
%   tool associated with that image.
%
%   BW = ROIPOLY(I,C,R) returns the region of interest specified by the
%   polygon described by vectors C and R. C and R, which must be the same
%   size, specify the column and row coordinates of the vertices.
%
%   BW = ROIPOLY(x,y,I,xi,yi) uses the vectors x and y to establish a
%   nondefault spatial coordinate system. xi and yi are equal-length
%   vectors that specify the coordinates of vertices as locations in this
%   coordinate system.
%
%   [BW,xi,yi] = ROIPOLY(...) returns the polygon coordinates in xi and yi.
%   Note that ROIPOLY always produces a closed polygon. If the points
%   specified describe a closed polygon (i.e., if the last pair of
%   coordinates is identical to the first pair), the length of xi and yi is
%   equal to the number of points specified. If the points specified do not
%   describe a closed polygon, ROIPOLY adds a final point having the same
%   coordinates as the first point. In this case the length of xi and yi is
%   one greater than the number of points specified.
%
%   [x,y,BW,xi,yi] = ROIPOLY(...) returns the XData and YData in x and y;
%   the mask image in BW; and the polygon coordinates in xi and yi.
%
%   If ROIPOLY is called with no output arguments, the resulting image is
%   displayed in a new figure.
%
%   Class Support
%   -------------
%   The input image I can be uint8, uint16, int16, single or double.  The
%   output image BW is logical. All other inputs and outputs are double.
%
%   Remarks
%   -------
%   For any of the ROIPOLY syntaxes, you can replace the input image I with
%   two arguments, M and N, that specify the row and column dimensions of an
%   arbitrary image. If you specify M and N with an interactive form of
%   ROIPOLY, an M-by-N black image is displayed, and you use the mouse to
%   specify a polygon with this image.
%
%   Example
%   -------
%       I = imread('eight.tif');
%       c = [222 272 300 270 221 194];
%       r = [21 21 75 121 121 75];
%       BW = roipoly(I,c,r);
%       figure, imshow(I), figure, imshow(BW)
%
%   See also IMPOLY, POLY2MASK, REGIONFILL, ROIFILT2, ROICOLOR.

%   Copyright 1993-2014 The MathWorks, Inc.

[xdata,ydata,num_rows,num_cols,xi,yi,placement_cancelled] = parse_inputs(varargin{:});

% return empty if user cancels operation
if placement_cancelled
    varargout = repmat({[]},nargout,1);
    return;
end

if length(xi)~=length(yi)
    error(message('images:roipoly:xiyiMustBeSameLength')); 
end

% Make sure polygon is closed.
if (~isempty(xi))
    if ( xi(1) ~= xi(end) || yi(1) ~= yi(end) )
        xi = [xi;xi(1)]; 
        yi = [yi;yi(1)];
    end
end
% Transform xi,yi into pixel coordinates.
roix = axes2pix(num_cols, xdata, xi);
roiy = axes2pix(num_rows, ydata, yi);

d = poly2mask(roix, roiy, num_rows, num_cols);

switch nargout
case 0
    figure
    imshow(d,'XData',xdata,'YData',ydata,'InitialMagnification', 'fit');
    
case 1
    varargout{1} = d;
    
case 2
    varargout{1} = d;
    varargout{2} = xi;
    
case 3
    varargout{1} = d;
    varargout{2} = xi;
    varargout{3} = yi;
    
case 4
    varargout{1} = xdata;
    varargout{2} = ydata;
    varargout{3} = d;
    varargout{4} = xi;
    
case 5
    varargout{1} = xdata;
    varargout{2} = ydata;
    varargout{3} = d;
    varargout{4} = xi;
    varargout{5} = yi;
    
otherwise
    error(message('images:roipoly:tooManyOutputArgs'));
    
end

end % roipoly

%%%
%%% parse_inputs
%%%

%--------------------------------------------------------
function [x,y,nrows,ncols,xi,yi,placement_cancelled] = parse_inputs(varargin)

% placement_cancelled only applies to interactive syntaxes. Assume placement_cancelled is false for initialization.
placement_cancelled = false;

cmenu_text = getString(message('images:roiContextMenuUIString:createMaskContextMenuLabel'));

switch nargin

case 0, 
    % ROIPOLY
    
    % verify we have a target image
    hFig = get(0,'CurrentFigure');
    hAx  = get(hFig,'CurrentAxes');
    hIm = findobj(hAx, 'Type', 'image');
    if isempty(hIm)
        error(message('images:roipoly:noImage'))
    end
    
    %  Get information from the current figure
    [x,y,a,hasimage] = getimage;
    if ~hasimage,
        error(message('images:roipoly:needImageInFigure'));
    end
    [xi,yi,placement_cancelled] = createWaitModePolygon(gca,cmenu_text);
    nrows = size(a,1);
    ncols = size(a,2);
    
case 1
    % ROIPOLY(A)
    a = varargin{1};
    nrows = size(a,1);
    ncols = size(a,2);
    x = [1 ncols];
    y = [1 nrows];
    imshow(a, 'InitialMagnification', 'fit');
    [xi,yi,placement_cancelled] = createWaitModePolygon(gca,cmenu_text);
    
case 2
    % ROIPOLY(M,N)
    nrows = varargin{1};
    ncols = varargin{2};
    a = repmat(uint8(0), nrows, ncols);
    x = [1 ncols];
    y = [1 nrows];
    imshow(a, 'InitialMagnification', 'fit');
    [xi,yi,placement_cancelled] = createWaitModePolygon(gca,cmenu_text);
    
case 3,
    % SYNTAX: roipoly(A,xi,yi)
    a = varargin{1};
    nrows = size(a,1);
    ncols = size(a,2);
    xi = varargin{2}(:);
    yi = varargin{3}(:);
    x = [1 ncols]; y = [1 nrows];

case 4,
    % SYNTAX: roipoly(m,n,xi,yi)
    nrows = varargin{1}; 
    ncols = varargin{2};
    xi = varargin{3}(:);
    yi = varargin{4}(:);
    x = [1 ncols]; y = [1 nrows];
    
case 5,
    % SYNTAX: roipoly(x,y,A,xi,yi)
    x = varargin{1}; 
    y = varargin{2}; 
    a = varargin{3};
    xi = varargin{4}(:); 
    yi = varargin{5}(:);
    nrows = size(a,1);
    ncols = size(a,2);
    x = [x(1) x(end)];
    y = [y(1) y(end)];
    
case 6,
    % SYNTAX: roipoly(x,y,m,n,xi,yi)
    x = varargin{1}; 
    y = varargin{2}; 
    nrows = varargin{3};
    ncols = varargin{4};
    xi = varargin{5}(:); 
    yi = varargin{6}(:);
    x = [x(1) x(end)];
    y = [y(1) y(end)];
    
otherwise,
    error(message('images:roipoly:invalidInputArgs'));

end

xi = cast_to_double(xi);
yi = cast_to_double(yi);
x = cast_to_double(x);
y = cast_to_double(y);
nrows= cast_to_double(nrows);
ncols = cast_to_double(ncols);

end % parse_inputs

%%%
% cast_to_double
%%%

%-----------------------------
function a = cast_to_double(a)
  if ~isa(a,'double')
    a = double(a);
  end
end    
    
%% 
function [xi, yi, placement_cancelled] = createWaitModePolygon(h_ax,finished_cmenu_text)
%createWaitModePolygon creates a blocking instance of impoly.
% [XI,YI,PLACEMENT_CANCELLED] =
% createWaitModePolygon(H_AX,FINISHED_CMENU_TEXT) creates a blocking
% instance of impoly parented to H_AX. The context menu description of the
% end of placement gesture receives the label FINISHED_CMENU_TEXT. Once the
% polygon is placed, it returns [XI,YI] vertices of where the polygon was
% oriented and PLACEMENT_CANCELLED which signals whether placement was
% cancelled.

%   Copyright 2007-2009 The MathWorks, Inc.

% initialize to "identity" mask
[xi,yi] = deal([]);

% interactively place tool
h_poly = iptui.roiPolygon(h_ax,finished_cmenu_text);
placement_cancelled = isempty(h_poly);
if placement_cancelled
    return;
end

pos = wait(h_poly);

abortedWaitMode = isempty(pos);
if abortedWaitMode
    placement_cancelled = true;
else
    xi = pos(:,1);
    yi = pos(:,2);
end
% We are done. Delete the polygon to clean up. Use isvalid to account for
% Cancel context menu item which will have already deleted the impoly
% instance.
if isvalid(h_poly)
    h_poly.delete();
end

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
    

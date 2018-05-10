function I = checkerboardBW(varargin)
%% This is a modification of Matlab's checkerboard function, to produce black&white tiles instead of with gray.
% to be used in Step0, and in the function createCheckerBoardImage
%
%CHECKERBOARD Create checkerboard image.
%   I = CHECKERBOARD creates a checkerboard image composed of squares that
%   have 10 pixels per side. The light squares on the left half of the
%   checkerboard are white. The light squares on the right half of the
%   checkerboard are gray.
%
%   I = CHECKERBOARD(N) creates a checkerboard where each square has N
%   pixels per side.
%
%   I = CHECKERBOARD(N,P,Q) creates a rectangular checkerboard. There are P
%   rows (uneven) and Q columns (even).
%
%   Input-output specs
%   ------------------ 
%   N,P,Q:    scalar positive integers 
%
%   I:        real double 2D matrix
%
%%
[n, p, q] = ParseInputs(varargin{:});

black = zeros(n);
white = ones(n);
tile = [black white; white black];

I = repmat(tile,floor(p/2),q/2);
lastRow=repmat([black  white],1,q/2);
I=[I; lastRow];


% % make right half plane have light gray tiles
% ncols = size(I,2);
% midcol = ncols/2 + 1; 
% I(:,midcol:ncols) = I(:,midcol:ncols) - .3;
% I(I<0) = 0;

%-------------------------------
% Function  ParseInputs
%
function [n, p, q] = ParseInputs(varargin)

% defaults
n = 10;
p = 4;
q = p;

narginchk(0,3);

varNames={'N', 'P', 'Q'};
for x = 1:1:length(varargin)
    validateattributes(varargin{x}, {'numeric'},...
                  {'integer' 'real' 'positive' 'scalar'}, ...
                  mfilename,varNames{x},x);
end

switch nargin
  case 0
    % I = CHECKERBOARD
    return;
    
  case 1
    % I = CHECKERBOARD(N)
    n = varargin{1};

  case 2
    % I = CHECKERBOARD(N,P)
    n = varargin{1};
    p = varargin{2};
    q = p;

  case 3
    % I = CHECKERBOARD(N,P,Q)
    n = varargin{1};
    p = varargin{2};    
    q = varargin{3};    
end

n = double(n);
p = double(p);
q = double(q);

%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>
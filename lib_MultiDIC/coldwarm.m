function [cMap]=coldwarm(varargin)
% function [cMap]=coldwarm
% function [cMap]=coldwarm(n)
% ------------------------------------------------------------------------
% Creates the colormap data for n levels for the cold and warm colormap. Low
% values define a cold blue color while high values define a warm/hot
% color. The 0 value is white. Therefore this map is suitable for data centered around 0
%------------------------------------------------------------------------

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

% cFire=[0 0 0; 1 0 0; 1 1 0; 1 1 1]; %Simple version
cFire=flipud([0.526,0.0720,0.00200;0.620,0.100,0.00300;...
       0.693,0.131,0.0110;0.770,0.172,0.0190;0.837,0.217,0.0320;0.894,0.261,0.0430;...
       0.954,0.315,0.0520;0.993,0.379,0.0630;0.998,0.484,0.0780;1,0.574,0.110;...
       1,0.644,0.152;1,0.711,0.211;1,0.773,0.274;1,0.837,0.322;1,0.914,0.369;...
       1,0.953,0.443;1,0.982,0.539;0.999,1,0.646;0.985,1,0.772;1,1,1]);
cMap=[rot90(cFire,2); cFire(2:end,:)];

[cMap]=resampleColormap(cMap,n);
 
%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>
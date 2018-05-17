function [varargout]=maxnumel(varargin)

% function [n,numBytesType]=maxnumel(a,numFreeBytes)

%% Parse input

switch nargin
    case 1
        a=varargin{1};
        numFreeBytes=[];
    case 2
        a=varargin{1};
        numFreeBytes=varargin{2};
end
%%

if isempty(numFreeBytes)
    %Max variable size available
    [numFreeBytes]=freeMemory;
end

%Size of input variable
sizA=size(a);
whos_struct=whos('a');
numBytesType=whos_struct.bytes;

%Max number of input type available
n=floor(numFreeBytes./numBytesType);

varargout{1}=n;
switch nargout    
    case 2
        varargout{2}=numBytesType;
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

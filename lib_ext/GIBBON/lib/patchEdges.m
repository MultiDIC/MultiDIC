function varargout=patchEdges(varargin)

% function [E,indUni1,indUni2]=patchEdges(F,uniOpt)
% -----------------------------------------------------------------------
% E=patchEdges(F,uniOpt)
% Uses the input faces array F to compute an edge array E. If uniOpt==1
% then the output array contain unique edges (irrespective of node order
% such that e.g. [4 1] and [1 4] are seen as the same edge). If uniOpt~=1
% then the double edges are maintained. If only one input is provided it is
% assumed to represent F and the default, uniOpt=0, is used. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/03/17
% 2016/01/19 Added additional output in relation to unique operation
%------------------------------------------------------------------------

%% PARSE INPUT
F=varargin{1};
switch nargin
    case 1
        uniOpt=0;
    case 2
        uniOpt=varargin{2};
    otherwise
        error('Wrong number of input arguments');
end
%% DERIVE NON-UNIQUE EDGES MATRIX
E1=F';
E2=F(:,[2:end 1])';
E=[E1(:) E2(:)];

%% REMOVE DOUBLE ENTRIES IF DESIRED

if uniOpt==1
    Es=sort(E,2); %Sorted so [1 4] and [4 1] are seen as the same edge
    [~,indUni1,indUni2]=unique(Es,'rows'); %Get indices for unique edges
    E=E(indUni1,:);    
end

%% Collect output

varargout{1}=E; 

if nargout>1 && uniOpt==1
    varargout{2}=indUni1;
    varargout{3}=indUni2;
elseif nargout>1 && uniOpt~=1
    error('Multiple outputs only available if uniOpt=1');
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

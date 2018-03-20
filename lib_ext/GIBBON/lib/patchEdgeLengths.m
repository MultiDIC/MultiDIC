function [D]=patchEdgeLengths(F,V)

% function [D]=patchEdgeLengths(F,V)
% -----------------------------------------------------------------------
% Computes the edge lengths (D) for the patch data specified by the faces
% (F) and vertices (V) arrays. If size(F,2)>2 it is assumed that F indeed
% represents faces. If however size(F,2)==2 it is instead assumed that F is
% an array representing edges. As such it skips the computation of the
% edges array. The edges array used is non-unique by default. See the
% |patchEdges| function for more details if the lengths of a unique set of
% edges is desired. 
%
%
% See also: |patchEdges|
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/03/17
%------------------------------------------------------------------------

%%

%Derive edge array
if size(F,2)>2 %The input is assumed to represent faces hence an edge array is derived
    E=patchEdges(F);
else %It is assumed that the input array represents an edges array
    E=F; 
end

%Derive edge vertex arrays
V_E1=V(E(:,1),:);
V_E2=V(E(:,2),:);

%Derive difference vectors
VD=(V_E1-V_E2);

%Compute the edge lengths
D=sqrt(sum(VD.^2,2));

 
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

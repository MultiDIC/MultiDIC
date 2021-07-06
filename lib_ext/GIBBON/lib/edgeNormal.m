function [varargout]=edgeNormal(F,V)

% function [NE,VE]=edgeNormal(F,V)
%-------------------------------------------------------------------------
% 
% To do: Expand to add corner normals
%-------------------------------------------------------------------------

%%

%Create edge-face indices
CE=(1:size(F,1))'; %Initialize face index set
CE=[CE(:,ones(size(F,2),1))]'; %Replicate
CE=CE(:); %Force as column

%Get edge description
E=patchEdges(F,0); 

%Get edge vectors
V_edge=vecnormalize(V(E(:,2),:)-V(E(:,1),:)); %The normalized edge vectors
N_face=patchNormal(F,V); %Face normals
N_face=N_face(CE,:); %replicate using edge-face indices
NE=vecnormalize(cross(V_edge,N_face,2)); %Calculate edge normals based on cross product

%%
varargout{1}=NE;
if nargout>1
    %Compute central edge coordinates if requested
    VE=patchCentre(E,V);
    varargout{2}=VE;
end

%%
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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

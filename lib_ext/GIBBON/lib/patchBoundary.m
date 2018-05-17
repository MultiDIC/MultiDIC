function [varargout]=patchBoundary(F,V)

%function [Eb,E,indBoundary]=patchBoundary(F,V)

%%
%Get non-unique edges
E1=F';
E2=F(:,[2:end 1])';
E=[E1(:) E2(:)];

%Get boundary indices
[indBoundary]=tesBoundary(E,V);

%Boundary edges
Eb=E(indBoundary,:);

%Output
varargout{1}=Eb; %Boundary edges
varargout{2}=E; %All edges
varargout{3}=indBoundary; %Indices for boundary edges
 
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

function [varargout]=pointSetPrincipalDir(X)


%% 

%Centre on own mean
MU=mean(X,1); %Point set mean
X=X-MU(ones(size(X,1),1),:); %Centre points around mean

%Compute singular value decomposition to get principal directions
[U,S,V]=svd(X,0); 

%% Collect output
switch nargout
    case 1
        varargout{1}=V; 
    case 2
        varargout{1}=V;
        varargout{2}=S; 
    case 3
        varargout{1}=V;
        varargout{2}=S;
        varargout{3}=U; 
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

function [D]=distND(V1,V2)

D=zeros(size(V1,1),size(V2,1));

for q=1:1:size(V1,2) %For all dimensions
    A=V1(:,q); %Coordinates of first set
    B=V2(:,q); %Coordinates of second set
    
    AB=A(:,ones(1,size(B,1)));
    BA=B(:,ones(1,size(A,1)))';
    
    D=D+(AB-BA).^2;
end
D=sqrt(D);
 
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

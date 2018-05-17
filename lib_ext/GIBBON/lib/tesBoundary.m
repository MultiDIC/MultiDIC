function [indBoundary]=tesBoundary(F,V)

if numel(V)==1
    numPoints=V;
else
    numPoints=size(V,1);
end

Fbs=sort(F,2);

sizVirt=numPoints*ones(1,size(F,2));

ind_F=subMat2ind(sizVirt,Fbs);

[~,indUni1,~]=unique(Fbs,'rows'); %Get indices for unique faces
F_uni=F(indUni1,:);

ind_F_uni=ind_F(indUni1,:);

ind=1:1:size(F,1);
ind=ind(~ismember(ind,indUni1));
ind_Fb_cut=ind_F(ind,:);
L_uni=~ismember(ind_F_uni,ind_Fb_cut);

indBoundary=indUni1(L_uni,:);
 
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

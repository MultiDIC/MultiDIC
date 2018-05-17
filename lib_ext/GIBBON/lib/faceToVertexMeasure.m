function [CV]=faceToVertexMeasure(F,V,CF)

[IND_F]=tesIND(F,V,0);

L=IND_F>0;

CV=ones(size(V,1),size(CF,2));
for q=1:1:size(CF,2)
    cf=CF(:,q);        
    cv=nan(size(IND_F));     
    cv(L)=cf(IND_F(L));
    cv=nanmean(cv,2);
    CV(:,q)=cv; 
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

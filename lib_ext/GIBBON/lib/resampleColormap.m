function [cmap_i]=resampleColormap(cmap,n)

% function [cmap_i]=resampleColormap(cmap,n)
% ------------------------------------------------------------------------
%
% This function resamples the input colormap cmap using n steps. Resampling
% is based on linear interpolation. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

ind=(1:1:size(cmap,1))';
ind_i=linspace(1,size(cmap,1),n)';
cmap_i=zeros(n,size(cmap,2));

%Interpolate color data
for q=1:1:size(cmap,2)
    cmap_i(:,q)=interp1(ind,cmap(:,q),ind_i,'linear');
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

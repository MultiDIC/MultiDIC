function [V_intersect,L_intersect,T] = triangleRayIntersection (V_ori,R,V,F,optStruct)

%
% Ray/triangle intersection using the algorithm proposed by M�ller and
% Trumbore (1997)
% 
% Based on implementation by: Jarek Tuszynski (jaroslaw.w.tuszynski@saic.com)
%
%
% Algorithm:
%  Function solves
%        |t|
%    M * |u| = (o-v0)
%        |v|
%  for [t; u; v] where M = [-d, v1-v0, v2-v0]. u,v are barycentric coordinates
%  and t - the distance from the ray origin in |d| units
%  ray/triangle intersect if u>=0, v>=0 and u+v<=1
%
% Note:
%  The algorithm is able to solve several types of problems:
%  * many faces / single ray  intersection
%  * one  face  / many   rays intersection
%  * one  face  / one    ray  intersection
%  * many faces / many   rays intersection
%  In order to allow that to happen all imput arrays are expected in Nx3
%  format, where N is number of vertices or rays. In most cases number of
%  vertices is different than number of rays, so one of the imputs will
%  have to be cloned to have the right size. Use "repmat(A,size(B,1),1)".
%
% Input (all arrays in in Nx3 format, where N is number of vertices or rays):
%  * orig : ray's origin
%  * dir  : ray's direction
%  * vert0, vert1, vert2: vertices of the triangle
%  * options: aditional customization options
%    * options.triangle - 'one sided' or 'two sided' (default) - how to treat
%        triangles. In 'one sided' version only intersections in single
%        direction are counted and intersections with back facing
%           tringles are ignored
%    * options.ray - 'ray' (default) or 'segment' - how to treat ray as an
%        infinite line (ray) or as line segment defined by a vector
%    * option.border - controls border handling. If 'normal'(default)
%        border points are included, but can be easily lost due to
%        rounding errors. If option.border='inclusive' borders points are
%        included, with a margin of option.eps. If option.border='exclusive'
%        borders points are excluded, with margin of option.eps.
%    * options.epsilon (default = 1e-5)
%
% Output:
%   * Intersect - boolean array of length N
%   * t   - distance from the ray origin to the intersection point in |dir|
%   * u,v - barycentric coordinates of the intersection point units
%
% Based on:
%  *"Fast, minimum storage ray-triangle intersection". Tomas Muller and
%    Ben Trumbore. Journal of Graphics Tools, 2(1):21--28, 1997.
%    http://www.graphics.cornell.edu/pubs/1997/MT97.pdf
%  * http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/raytri/
%  * http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/raytri/raytri.c
%
% Kevin Matthaus Moerman (kevinmoerman@hotmail.com)
% 2013/12/18: Changed for patch type input and Cartesian coordinate output
% TO DO: Full validation, help adjustment according to changes, generation
% of DEMO file

%%

%Get triangle vertices
V0=V(F(:,1),:);
V1=V(F(:,2),:);
V2=V(F(:,3),:);

%% verify that inputs are in correct format

if (size(V_ori,1)==3 && size(V_ori,2)~=3);
    V_ori =V_ori';
end

if (size(R,1)==3 && size(R,2)~=3);
    R =R';
end

if (size(V0,1)==3 && size(V0,2)~=3);
    V0=V0';
end

if (size(V1,1)==3 && size(V1,2)~=3);
    V1=V1';
end

if (size(V2,1)==3 && size(V2,2)~=3);
    V2=V2';
end

if (any(size(V_ori)~=size(V0)) || ...
        any(size(V_ori)~=size(V1)) || ...
        any(size(V_ori)~=size(V2)) || ...
        any(size(V_ori)~=size(R  )) )
    error('All input vectors have to be of the same size.');
end

if (size(V_ori,2)~=3)
    error('All input vectors have to be in Nx3 format.');
end

%% Read options from optStruct

%Create defaults
if (nargin<5)    
    optStruct.eps      = 1e-5;
    optStruct.triangle = 'two sided';
    optStruct.ray      = 'ray';
    optStruct.border   = 'normal';
end

%% initialize default output
L_intersect = false(size(V_ori,1),1);
T=nan(size(V_ori,1),1); 
U_bar=T;
V_bar=T;
V_intersect=nan(size(V_ori));

%% Find faces parallel to the ray
E1 = V1-V0;          % find vectors for two edges sharing vert0
E2 = V2-V0;
Tvec  = V_ori -V0;          % distance from vert0 to ray origin
Pvec  = cross(R, E2,2);  % begin calculating determinant - also used to calculate U parameter

detMat   = sum(E1.*Pvec,2);   % determinant of the matrix M = dot(E1,Pvec)

logicParalel = (abs(detMat)<optStruct.eps);    % if determinant is near zero then ray lies in the plane of the triangle

if all(logicParalel) % if all parallel than no intersections
    return;
end

switch optStruct.border
    case 'normal'
        zeroLim=0.0;
    case 'inclusive'
        zeroLim=optStruct.eps;
    case 'exclusive'
        zeroLim=-optStruct.eps;
end

%% Different behavior depending on one or two sided triangles
if strcmpi(optStruct.triangle,'two sided')          % treats triangles as two sided    
    
    detMat(logicParalel) = 1;                       % change to avoid division by zero
    
    U_bar  = sum(Tvec.*Pvec,2)./detMat;             % calculate U parameter used to test bounds
    L_ok = (~logicParalel & U_bar>=-zeroLim & U_bar<=1.0+zeroLim);% mask which allows performing next 2 operations only when needed
    
    if ~any(L_ok) % if all ray/plane intersections are outside the triangle than no intersections
        return;
    end
    
    Qvec = cross(Tvec(L_ok,:), E1(L_ok,:),2); % prepare to test V parameter
    V_bar(L_ok,:) = sum(R(L_ok,:).*Qvec,2) ./ detMat(L_ok,:);  % calculate V parameter used to test bounds    
    L_intersect = (V_bar>=-zeroLim & U_bar+V_bar<=1.0+zeroLim & L_ok);
    
%     if (nargout==1 && strcmpi(ray,'ray'));
%         disp('RAY1');
%         return;
%     end
    
    T(L_ok,:) = sum(E2(L_ok,:).*Qvec,2)./detMat(L_ok,:);
    
    if ~(strcmpi(optStruct.ray,'ray'))
%         disp('RAY2');
        L_intersect = (L_intersect & T>=-zeroLim & T<=1.0+zeroLim);
    end
    
else % treats triangles as one sided
    
    U_bar = sum(Tvec.*Pvec,2);                   % calculate U parameter used to test bounds
    L_ok = (detMat>optStruct.eps & U_bar>=0.0 & U_bar<=detMat);        % mask which allows performing next 2 operations only when needed
    
    if ~any(L_ok) % if all ray/plane intersections are outside the triangle than no intersections
        return;
    end
    
    Qvec = cross(Tvec(L_ok,:), E1(L_ok,:),2); % prepare to test V parameter
    V_bar(L_ok,:) = sum(R(L_ok,:).*Qvec,2);        % calculate V parameter used to test bounds
    L_intersect = (detMat>optStruct.eps & U_bar>=-zeroLim & V_bar>=-zeroLim & U_bar+V_bar<=detMat*(1+zeroLim));
    
%     if (nargout==1 && strcmpi(ray,'ray'));
%         return;
%     end
    
    T(L_ok,:)  = sum(E2(L_ok,:).*Qvec,2);
    inv_det = zeros(size(detMat));
    inv_det(L_ok,:) = 1./detMat(L_ok,:);
    T = T.*inv_det;  % calculate t - distance from origin to the intersection in |d| units
    U_bar = U_bar.*inv_det;
    V_bar = V_bar.*inv_det;
    
    if ~(strcmpi(optStruct.ray,'ray'))
        L_intersect = (L_intersect & T>=-zeroLim & T<=1.0+zeroLim); % intersection between origin and destination
    end
end

V_intersect=((1-U_bar(:,ones(1,3))-V_bar(:,ones(1,3))).*V0) + (U_bar(:,ones(1,3)).*V1) + (V_bar(:,ones(1,3)).*V2);
V_intersect(~L_intersect,:)=NaN;
T(~L_intersect,:)=NaN;
 
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

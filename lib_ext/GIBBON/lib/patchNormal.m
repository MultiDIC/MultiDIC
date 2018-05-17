function [varargout]=patchNormal(F,V)

% [N,Vn,Nv]=patchNormal(F,V)
% ------------------------------------------------------------------------
% Normals are derived based on cross product of triangle edge vectors. Each
% triangle is constructed using the two points of its face and the mean of
% the face. 
%
%
% To do: Check for co-linear edges ?
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2014/06/02 %Updated general lay-out and commenting
% 2015/09/22 %Fixed to allow for 2D patch data
% 2016/11/15 %Added handling of non-triangular faces
%------------------------------------------------------------------------

%%

%Deal with 2D patch data
if size(V,2)==2
    V(:,3)=0; 
end

%Get vertices
Vp=V; %The input vertices (the original V is used later on)

%Get mean face coordinates for normal vectors
X=Vp(:,1); Y=Vp(:,2); Z=Vp(:,3);
if size(F,1)==1
    Vn=[mean(X(F)',2) mean(Y(F)',2) mean(Z(F)',2)];
else
    Vn=[mean(X(F),2) mean(Y(F),2) mean(Z(F),2)];
end
sizV=size(Vp);
Vp=[Vp; Vn]; %Add mean patch points to vertex list

%Vertex indices for triangle edge sets
indVertex=[1:size(F,2); 2:size(F,2)+1]';
indVertex(end,end)=1;

%Constructing triangle faces matrix
Ft=nan(size(F,1),3); %initialse

%Derive face normals for all triangles composing the patch face
N=nan(size(F,1),3,size(F,2));% 1 face normal per subtriangle
for q=1:1:size(F,2)
    
    Ft(:,1:2)=F(:,[indVertex(q,1) indVertex(q,2)]);
    Ft(:,3)=(sizV(1)+1):size(Vp,1); %Set last points as face mean
    
    %Getting triangle surface normal (cross product of two edge vectors)
    vec1=[Vp(Ft(:,2),1)-Vp(Ft(:,1),1)  Vp(Ft(:,2),2)-Vp(Ft(:,1),2)  Vp(Ft(:,2),3)-Vp(Ft(:,1),3)]; %First edge vector
    vec2=[Vp(Ft(:,3),1)-Vp(Ft(:,1),1)  Vp(Ft(:,3),2)-Vp(Ft(:,1),2)  Vp(Ft(:,3),3)-Vp(Ft(:,1),3)]; %Second edge vector
    
    %Derive triangle face normal using cross product
    Nq=cross(vec1,vec2,2); 
    Nq=Nq./(sqrt(sum(Nq.^2,2))*ones(1,size(Nq,2))); %Normalizing vector length
    
    N(:,:,q)=Nq; %Add layer
end

N=mean(N,3);
N=N./(sqrt(sum(N.^2,2))*ones(1,size(N,2))); %Normalizing vector length

%% Collect output

switch nargout
    case 1 %Only face normals
        varargout{1}=N;
    case 2 %Face normals and face centres
        varargout{1}=N;
        varargout{2}=Vn;
    case 3 %Face normals, face centres and vertex normals                
        varargout{1}=N;
        varargout{2}=Vn;
        [Nv]=faceToVertexMeasure(F,V,N);
        varargout{3}=Nv;
end

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

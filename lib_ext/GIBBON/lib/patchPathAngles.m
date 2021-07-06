function [indAngles]=patchPathAngles(F,V,ind,isClosedLoop)

ind=ind(:);

%Get non-unique edge set
E=patchEdges(F,0);
Es=sort(E,2); 

[NE,~]=edgeNormal(F,V);
sizVirt=size(V,1)*ones(1,2); % ind_Es=sub2indn(sizVirt,Es);

%Check if path order is consistent with edges
if ~any(all(E(:,1)==ind(1) & E(:,2)==ind(2),2))
    orderFlipped=1;
    ind=flip(ind); %Flip if not consistent 
else
    orderFlipped=0;
end

if isClosedLoop==1
    es=[ind(1:end-1) ind(2:end); ind(end) ind(1); ]; %Closed loop
else
    es=[ind(1:end-1) ind(2:end);]; %Open segment
end
ind_es=sub2indn(sizVirt,sort(es,2));

ne=zeros(size(es,1),size(NE,2));
for q=1:1:size(NE,2)    
    S=sparse(Es(:,1),Es(:,2),NE(:,q),sizVirt(1),sizVirt(2),size(NE,1));% sparse(i,j,v,m,n,nz)
    ne(:,q)=S(ind_es);
end

[~,~,NF]=patchNormal(F,V); 

if isClosedLoop==1
    ne_1=ne;
    ne_2=[ne(end,:); ne(1:end-1,:)];
    nf=NF(ind,:);    
else
    ne_1=ne(2:end,:);
    ne_2=ne(1:end-1,:);
    nf=NF(ind(2:end-1),:);
end

%Compute angles
A=pi-real(acos(dot(ne_1,ne_2,2)));

%Make sure angle conforms to mesh side
X=vecnormalize(cross(ne_1,ne_2,2));
logic1=dot(nf,X,2)<0;
logicAngle= (A<(pi) & logic1);
A(logicAngle)=(2*pi)-A(logicAngle);

if isClosedLoop==1
    indAngles=A;
else
    indAngles=nan(size(ind));
    indAngles(2:end-1)=A;
end

if orderFlipped==1
    indAngles=flip(indAngles); %Flip data back
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

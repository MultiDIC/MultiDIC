function [C]=patchConnectivity(varargin)

% function [C]=patchConnectivity(F,V)
% -----------------------------------------------------------------------
% This functions creates connectivity matrices for the input patch data
% defined by the faces F and the vertices V. The output is a structure
% containing the connectivity matrices:
%
% C.vertex.vertex
% C.vertex.face
% C.vertex.edge
% 
% C.edge.face
% C.edge.vertex
% C.edge.edge
% 
% C.face.vertex
% C.face.face
% C.face.edge
% 
% Change log: 
% 2018/08/22 
% -----------------------------------------------------------------------

%% Parse input
switch nargin
    case 1
        F=varargin{1};
        V=[]; 
    case 2
        F=varargin{1};
        V=varargin{2};
end

if isempty(V)
    numVertices=max(F(:)); %Assume all points are used in F
else
    numVertices=size(V,1);
end
numFaces=size(F,1);
numFaceVertices=size(F,2);

%% Face-edge connectivity
E=patchEdges(F,0); %The non-unique edge set
E_sort=sort(E,2); %Sorted in column dir so 1 2 looks the same as 2 1
indEdges=sub2indn(numVertices*ones(1,2),E_sort); %Create "virtual" indices
[~,ind1,indEdges_2]=unique(indEdges); %Get indices for unique edges
E_uni=E(ind1,:); %Get unique edges
ind_F_E=reshape(indEdges_2,numFaceVertices,numFaces)'; %Reshape to get results
numEdges=size(E_uni,1);
numEdgeVertices=size(E_uni,2);

%% Edge-face connectivity
ind=(1:1:numFaces)'; %Indices for all faces
ind=ind(:,ones(1,numFaceVertices)); %Indices copied over so it is the size of F
ind=ind(:); %Force as column
ind_E_F=sparse(ind_F_E(:),ind,ind,numEdges,numFaces,numel(ind)); %Create sparse form of connectivity matrix
ind_E_F=sort(ind_E_F,2,'descend'); %Sort the sparse array
ind_E_F=full(ind_E_F(:,[1 2])); %Keep relevant columns, convert to full array

%% Vertex-face connectivity
ind=(1:1:numFaces)'; 
ind=ind(:,ones(1,numFaceVertices));
ind=ind(:); 
ind_V_F=sparse(F(:),ind,ind,numVertices,numFaces);
ind_V_F=sort(ind_V_F,2,'descend');
[~,J,~] = find(ind_V_F);
ind_V_F=full(ind_V_F(:,1:max(J)));

%% Vertex-edge connectivity
ind=(1:1:numEdges)';
ind=ind(:,ones(1,numEdgeVertices)); 
ind=ind(:);
ind_V_E=sparse(E_uni(:),ind,ind,numVertices,numEdges);
ind_V_E=sort(ind_V_E,2,'descend');
[~,J,~] = find(ind_V_E);
ind_V_E=full(ind_V_E(:,1:max(J)));

%% Vertex-vertex connectivity
EV=[E_uni;fliplr(E_uni)];
ind_V_V=sparse(EV(:,1),EV(:,2),EV(:,2),numVertices,numVertices);
ind_V_V=sort(ind_V_V,2,'descend');
[~,J,~] = find(ind_V_V);
ind_V_V=full(ind_V_V(:,1:max(J)));

%% Face-face connectivity 
A=ind_E_F(ind_F_E(:),:);
ind_F_F=reshape(A,numFaces,numel(A)/numFaces);
ind=(1:1:numFaces)';
ind=ind(:,ones(1,size(ind_F_F,2)));
ind=ind(:);
logicValid=ind_F_F(:)>0;
ind_F_F=ind_F_F(logicValid);
ind=ind(logicValid);
ind_F_F=sparse(ind,ind_F_F(:),ind_F_F(:),numFaces,numFaces,numel(ind));
ind_F_F(inddiag(ind_F_F))=0;
ind_F_F=sort(ind_F_F,2,'descend');
ind_F_F=full(ind_F_F(:,1:numFaceVertices));

%% Edge-edge connectivity
A=ind_V_E(E_uni(:),:);
ind_E_E=reshape(A,numEdges,numel(A)./numEdges);
ind=(1:1:numEdges)';
ind=ind(:,ones(1,size(ind_E_E,2)));
ind=ind(:);
logicValid=ind_E_E(:)>0;
ind=ind(logicValid);
ind_E_E=ind_E_E(logicValid);
ind_E_E=sparse(ind(:),ind_E_E(:),ind_E_E(:),numEdges,numEdges,numel(ind));
ind_E_E(inddiag(ind_E_E))=0;
ind_E_E=sort(ind_E_E,2,'descend');
[~,J,~] = find(ind_E_E);
ind_E_E=full(ind_E_E(:,1:max(J)));

%% Collect output in structure

C.vertex.vertex=ind_V_V;
C.vertex.face=ind_V_F;
C.vertex.edge=ind_V_E;

C.edge.face=ind_E_F;
C.edge.vertex=E_uni;
C.edge.edge=ind_E_E;

C.face.vertex=F;
C.face.face=ind_F_F;
C.face.edge=ind_F_E;

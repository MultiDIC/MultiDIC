function indVertexBowtied=findIndVertexBowtied(F,V)
Eb=patchBoundary(F,V);

% Get patch connectivity 
C=patchConnectivity(F,V);
E=C.edge.vertex;
vertexEdgeConnectivity=C.vertex.edge;

% Work out membership to boundary set
sizVirt=size(V,1)*ones(1,2); 
ind_uniqueEdges=sub2indn(sizVirt,sort(E,2)); 
ind_boundaryEdges=sub2indn(sizVirt,sort(Eb,2)); 
logicIsBoundaryVertex=ismember(ind_uniqueEdges,ind_boundaryEdges);

vertexBoundaryEdgeConnectivity=vertexEdgeConnectivity; 
vertexBoundaryEdgeConnectivity(vertexBoundaryEdgeConnectivity>0)=logicIsBoundaryVertex(vertexBoundaryEdgeConnectivity(vertexBoundaryEdgeConnectivity>0));

indBoundaryVertices=unique(Eb(:)); 

logicVertexBowtied=sum(vertexBoundaryEdgeConnectivity,2)>2;

indVertexBowtied=find(logicVertexBowtied);



end
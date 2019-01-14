function [F1,V1,ind1]=triSurfSelfTriangulateBoundary(F1,V1,ind1,angleThreshold,isClosedLoop)

if ind1(1)==ind1(end)
    ind1=ind1(1:end-1);
    cropNeeded=1;
else
    cropNeeded=0;
end

%Get non-unique edge set
E=patchEdges(F1,0);

%Check if path order is consistent with edges
if ~any(all(E(:,1)==ind1(1) & E(:,2)==ind1(2),2))
    orderFlipped=1;
    ind1=flip(ind1); %Flip if not consistent    
else
    orderFlipped=0;
end

while 1
    
    [A]=patchPathAngles(F1,V1,ind1,isClosedLoop);
    
    indSharp=find(A<=angleThreshold);
    if isempty(indSharp)
        break
    end
    A_sharp=A(indSharp);
    
    if isClosedLoop==1
        es=[ind1(1:end-1) ind1(2:end); ind1(end) ind1(1); ]; %Closed loop
    else
        es=[ind1(1:end-1) ind1(2:end);]; %Open segment
    end
    
    [~,indNow]=min(A_sharp); %Sharpest first
    indVertexNow=ind1(indSharp(indNow));
    indEdges=(es(any(es==indVertexNow,2),:));
    f=fliplr(ind1(ismember(ind1,indEdges(:)))');
    
    F1=[F1;f];    
    
    ind1=ind1(ind1~=indVertexNow);
    ind1=ind1(:);    
    
end

if orderFlipped
    ind1=flip(ind1);
end

if cropNeeded
    ind1(end+1)=ind1(1);
end
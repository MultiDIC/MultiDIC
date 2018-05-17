function [V1_trace]=triSurfRaySetIntersect(V1,N1,F2,V2,optStruct)

numSteps=size(V1,1);
V1_trace=nan(numSteps,3);
c=1;
% hw=waitbar(c/numSteps,['Ray tracing...',num2str(round(100.*c/numSteps)),'%']);
L1=false(numSteps,1);
logicShot=false(size(F2,1),1);
for q=1:1:numSteps
    v1=V1(q,:);
    n1=N1(q,:);
    [V_intersect,logicIntersect,~] = triangleRayIntersection(v1(ones(size(F2,1),1),:),n1(ones(size(F2,1),1),:),V2,F2,optStruct);
    logicShot=logicShot | logicIntersect;
    V_intersect=V_intersect(logicIntersect,:);
        
    if nnz(logicIntersect)>0
%         [~,indMin]=minDist(v1,V_intersect);                
        d=sqrt(sum(v1(ones(size(V_intersect,1),1),:)-V_intersect,2));
        [~,indMin]=min(d);
        V1_trace(q,:)=V_intersect(indMin,:);        
        L1(q)=1;
    end
%     waitbar(c/numSteps,hw,['Ray tracing...',num2str(round(100.*c/numSteps)),'%']);
    c=c+1;
end
% close(hw);


end
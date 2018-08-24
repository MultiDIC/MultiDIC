function [CT1,CT2] = removeOverlapSurface_temp_2018_07_26(F1,F2,V1,V2,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trying to speed up this function by removing groups of faces at a time
% and not one by one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% [CT1,CT2] = removeOverlapSurface(F1,F2,V1,V2)
% [CT1,CT2] = removeOverlapSurface(F1,F2,V1,V2,Q1,Q2)
% [CT1,CT2] = removeOverlapSurface(F1,F2,V1,V2,Q1,Q2,minGap)
% [CT1,CT2] = removeOverlapSurface(F1,F2,V1,V2,[],[],minGap)

%%
% Inputs:
% F1,F2: triangular faces (nF-by-3)
% V1,V2: vertices 3D corrdinates (nV-by-3)
% Q1,Q2: quality measures(larger is worse) (nV-by-3)
% minGap: minimum distance between surfaces (scalar)

% output:
% CT1,CT2: logical, which faces to keep (nF-by-3)
%%
nargin=numel(varargin);
switch nargin
    case 0
        % if no quality measure is given, use distances between the surfaces as quality measure
        Q1=triSurfSetDist(F1,V1,F2,V2,'ray');
        Q2=triSurfSetDist(F2,V2,F1,V1,'ray');
        nOverlap=max([sum(~isnan(Q1)) sum(~isnan(Q2))]);
        minGap=.6*min([nanmean(patchEdgeLengths(F1,V1)) nanmean(patchEdgeLengths(F2,V2))]);
    case 2
        Q1=varargin{1};
        Q2=varargin{2};
        nOverlap=sum(~isnan(triSurfSetDist(F1,V1,F2,V2,'ray')));
        minGap=.6*min([nanmean(patchEdgeLengths(F1,V1)) nanmean(patchEdgeLengths(F2,V2))]);
    case 3
        if ~isempty(varargin{1}) && ~isempty(varargin{2})
            Q1=varargin{1};
            Q2=varargin{2};
            nOverlap=sum(~isnan(triSurfSetDist(F1,V1,F2,V2,'ray')));
        else
            Q1=triSurfSetDist(F1,V1,F2,V2,'ray');
            Q2=triSurfSetDist(F2,V2,F1,V1,'ray');
            nOverlap=max([sum(~isnan(Q1)) sum(~isnan(Q2))]);
        end
        if ~isempty(varargin{3})
            minGap=varargin{3};
        else
            minGap=.6*min([nanmean(patchEdgeLengths(F1,V1)) nanmean(patchEdgeLengths(F2,V2))]);
        end
    otherwise
        error('Wrong number of input arguments');
end

%% Parse input
%Ray trace 1 onto 2
optStruct.eps      = 1e-6;
optStruct.triangle = 'two sided';
optStruct.ray      = 'ray';
optStruct.border   = 'normal';

%% Get vertex surface normal vectors
[~,~,N1]=patchNormal(F1,V1);
[~,~,N2]=patchNormal(F2,V2);

%% remove overlapping surfaces

CT1=true(size(F1,1),1); % faces to keep from surface 1
CT2=true(size(F2,1),1); % faces to keep from surface 2

logicNoHit1=false(size(V1,1),1); % vertices from surface 1 which do not overlap with surface 2
logicNoHit2=false(size(V2,1),1); % vertices from surface 2 which do not overlap with surface 1

%remove overlapping surfaces
hw=waitbar(0,'Resolving Overlap...');
meanEdgeLength=nanmean([patchEdgeLengths(F1,V1) ; patchEdgeLengths(F2,V2)]);
icount=0;
while 1
    icount=icount+1;
    waitbar(icount/nOverlap);
    % Get boundary edges indices
    [Eb1]=patchBoundary(F1(CT1,:),V1);
    [Eb2]=patchBoundary(F2(CT2,:),V2);
    % unique vertices of boundary edges
    indBoundary1=unique(Eb1(:));
    indBoundary2=unique(Eb2(:));
    % remove boundary vertices which were not hit by ray tracing previously
    indBoundary1=indBoundary1(~ismember(indBoundary1,find(logicNoHit1)));
    indBoundary2=indBoundary2(~ismember(indBoundary2,find(logicNoHit2)));
    
    % if both boundaries are not empty, find which one has the worst point and remove its face
    if ~isempty(indBoundary1) && ~isempty(indBoundary2)
        % ray tracing from vertices of boundary 1 to faces of 2 (V1_trace is the size of indBoundary1
        [V1_trace]=triSurfRaySetIntersect(V1(indBoundary1,:),N1(indBoundary1,:),F2(CT2,:),V2,optStruct);
        V1_trace_diff=V1(indBoundary1,:)-V1_trace;
        V1_trace_dist=sqrt(sum(V1_trace_diff.^2,2));
        logicOverlapClose1=V1_trace_dist<2*meanEdgeLength;
%         logicOverlap1=~any(isnan(V1_trace),2); % true if the boundary vertices of 1 hit the other surface but is not too far
        c1=logicOverlapClose1; % true if the boundary vertices of 1 hit the other surface
        logicV1in=false(size(V1,1),1); % color all vertices false
        logicV1in(indBoundary1)=c1; % true if vertices are on the boundary and hit the other surface
        indV1in=find(logicV1in);
        logicNoHit1(indBoundary1(~c1))=1; % add vertices to boundary no hit list
        [worstQ1,indBoundaryworstQ1]=max(Q1(logicV1in)); % vertex with worst score
        
        [V2_trace]=triSurfRaySetIntersect(V2(indBoundary2,:),N2(indBoundary2,:),F1(CT1,:),V1,optStruct);
        V2_trace_diff=V2(indBoundary2,:)-V2_trace;
        V2_trace_dist=sqrt(sum(V2_trace_diff.^2,2));
        logicOverlapClose2=V2_trace_dist<2*meanEdgeLength;
%         logicOverlap2=~any(isnan(V2_trace),2); % true if the boundary vertices of 1 hit the other surface
        c2=logicOverlapClose2; % true if the boundary vertices of 1 hit the other surface
        logicV2in=false(size(V2,1),1); % color all vertices false
        logicV2in(indBoundary2)=c2; % true if vertices are on the boundary and hit the other surface
        indV2in=find(logicV2in);
        logicNoHit2(indBoundary2(~c2))=1; % add vertices to no hit list
        [worstQ2,indBoundaryworstQ2]=max(Q2(logicV2in)); % vertex with worst score
        
        if worstQ1>worstQ2
            logicV1Remove=false(size(V1,1),1); % color all vertices false
            logicV1Remove(indV1in(indBoundaryworstQ1))=true(1); % true if vertices are on the boundary and hit the other surface
            logicF1Remove=~any(logicV1Remove(F1),2); % false for faces for which none of the vertices is in logicV1Remove
            CT1 = CT1 & logicF1Remove; % in CT1 remove (turn to false) boundary face
        else
            logicV2Remove=false(size(V2,1),1); % color all vertices false
            logicV2Remove(indV2in(indBoundaryworstQ2))=true(1); % true if vertices are on the boundary and hit the other surface
            logicF2Remove=~any(logicV2Remove(F2),2); % false for faces for which none of the vertices is in logicV1Remove
            CT2 = CT2 & logicF2Remove; % in CT2 remove (turn to false) boundary face
        end
        
    elseif ~isempty(indBoundary1) && isempty(indBoundary2)
        % ray tracing from vertices of boundary 1 to faces of 2 (V1_trace is the size of indBoundary1
        [V1_trace]=triSurfRaySetIntersect(V1(indBoundary1,:),N1(indBoundary1,:),F2(CT2,:),V2,optStruct);
        V1_trace_diff=V1(indBoundary1,:)-V1_trace;
        V1_trace_dist=sqrt(sum(V1_trace_diff.^2,2));
        logicOverlapClose1=V1_trace_dist<2*meanEdgeLength;
%         logicOverlap1=~any(isnan(V1_trace),2); % true if the boundary vertices of 1 hit the other surface
        c1=logicOverlapClose1; % true if the boundary vertices of 1 hit the other surface
        logicV1in=false(size(V1,1),1); % color all vertices false
        logicV1in(indBoundary1)=c1; % true if vertices are on the boundary and hit the other surface
        indV1in=find(logicV1in);
        logicNoHit1(indBoundary1(~c1))=1; % add vertices to no hit list
        [~,indBoundaryworstQ1]=max(Q1(logicV1in)); % vertex with worst score
        
        logicV1Remove=false(size(V1,1),1); % color all vertices false
        logicV1Remove(indV1in(indBoundaryworstQ1))=true(1); % true if vertices are on the boundary and hit the other surface
        logicF1Remove=~any(logicV1Remove(F1),2); % false for faces for which none of the vertices is in logicV1Remove
        CT1 = CT1 & logicF1Remove; % in CT1 remove (turn to false) boundary face
        
    elseif isempty(indBoundary1) && ~isempty(indBoundary2)
        [V2_trace]=triSurfRaySetIntersect(V2(indBoundary2,:),N2(indBoundary2,:),F1(CT1,:),V1,optStruct);
                V2_trace_diff=V2(indBoundary2,:)-V2_trace;
        V2_trace_dist=sqrt(sum(V2_trace_diff.^2,2));
        logicOverlapClose2=V2_trace_dist<2*meanEdgeLength;
%         logicOverlap2=~any(isnan(V2_trace),2); % true if the boundary vertices of 1 hit the other surface
        c2=logicOverlapClose2; % true if the boundary vertices of 1 hit the other surface
        logicV2in=false(size(V2,1),1); % color all vertices false
        logicV2in(indBoundary2)=c2; % true if vertices are on the boundary and hit the other surface
        indV2in=find(logicV2in);
        logicNoHit2(indBoundary2(~c2))=1; % add vertices to no hit list
        [~,indBoundaryworstQ2]=max(Q2(logicV2in)); % vertex with worst score
        
        logicV2Remove=false(size(V2,1),1); % color all vertices false
        logicV2Remove(indV2in(indBoundaryworstQ2))=true(1); % true if vertices are on the boundary and hit the other surface
        logicF2Remove=~any(logicV2Remove(F2),2); % false for faces for which none of the vertices is in logicV1Remove
        CT2 = CT2 & logicF2Remove; % in CT2 remove (turn to false) boundary face
        
    end
    
    if isempty(indBoundary1) && isempty(indBoundary2)
        break
    end
    
end

%% then remove surfaces that are too close (not overlapping but too close)
logicNoHit1=false(size(V1,1),1); % vertices from surface 1 which do not overlap with surface 2
logicNoHit2=false(size(V2,1),1); % vertices from surface 2 which do not overlap with surface 1

X1=V1(:,1);
Y1=V1(:,2);
Z1=V1(:,3);
X2=V2(:,1);
Y2=V2(:,2);
Z2=V2(:,3);

while 1
    icount=icount+1;
    waitbar(icount/nOverlap);
    
    % Get boundary edges indices
    [Eb1]=patchBoundary(F1(CT1,:),V1);
    [Eb2]=patchBoundary(F2(CT2,:),V2);
    % unique vertices of boundary edges
    indBoundary1=unique(Eb1(:));
    indBoundary2=unique(Eb2(:));
    % remove boundary vertices which were not hit by ray tracing previously
    indBoundary1=indBoundary1(~ismember(indBoundary1,find(logicNoHit1)));
    indBoundary2=indBoundary2(~ismember(indBoundary2,find(logicNoHit2)));
        
    % if both boundaries are not empty, find which one has the worst point and remove its face
    if ~isempty(indBoundary1) && ~isempty(indBoundary2)
        Eb1_crop=Eb1(all(ismember(Eb1,indBoundary1),2),:);
        V1m=[mean(X1(Eb1_crop),2) mean(Y1(Eb1_crop),2) mean(Z1(Eb1_crop),2)];
        Eb2_crop=Eb2(all(ismember(Eb2,indBoundary2),2),:);
        V2m=[mean(X2(Eb2_crop),2) mean(Y2(Eb2_crop),2) mean(Z2(Eb2_crop),2)];
        
        % distance vertices of boundary 1 to vertices of 2
        D1=minDist(V1(indBoundary1,:),[V2(indBoundary2,:); V2m]);
        [logicDist1]=D1<minGap; % true if small distance
        c1=logicDist1; % true if the boundary vertices of 1 very close
        logicV1in=false(size(V1,1),1); % color all vertices false
        logicV1in(indBoundary1)=c1; % true if vertices are on the boundary and hit the other surface
        indV1in=find(logicV1in);
        [worstQ1,indBoundaryworstQ1]=max(Q1(logicV1in)); % vertex with worst score
        
        % distance vertices of boundary 1 to vertices of 2
        D2=minDist(V2(indBoundary2,:),[V1(indBoundary1,:); V1m]);
        [logicDist2]=D2<minGap; % true if small distance
        c2=logicDist2; % true if the boundary vertices of 1 hit the other surface or very close
        logicV2in=false(size(V2,1),1); % color all vertices false
        logicV2in(indBoundary2)=c2; % true if vertices are on the boundary and hit the other surface
        indV2in=find(logicV2in);
        [worstQ2,indBoundaryworstQ2]=max(Q2(logicV2in)); % vertex with worst score
        
        if worstQ1>worstQ2
            logicV1Remove=false(size(V1,1),1); % color all vertices false
            logicV1Remove(indV1in(indBoundaryworstQ1))=true(1); % true if vertices are on the boundary and hit the other surface
            logicF1Remove=~any(logicV1Remove(F1),2); % false for faces for which none of the vertices is in logicV1Remove
            CT1 = CT1 & logicF1Remove; % in CT1 remove (turn to false) boundary face
            logicNoHit1(indBoundary1(~c1))=1; % add vertices to boundary no hit list
        else
            logicV2Remove=false(size(V2,1),1); % color all vertices false
            logicV2Remove(indV2in(indBoundaryworstQ2))=true(1); % true if vertices are on the boundary and hit the other surface
            logicF2Remove=~any(logicV2Remove(F2),2); % false for faces for which none of the vertices is in logicV1Remove
            CT2 = CT2 & logicF2Remove; % in CT2 remove (turn to false) boundary face
            logicNoHit2(indBoundary2(~c2))=1; % add vertices to no hit list
        end
    
    else
        break
    end
    
end
close(hw)
% 

end
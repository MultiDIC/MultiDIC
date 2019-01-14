function [CT1,CT2] = removeOverlapSurface_temp_2018_10_30(F1,F2,V1,V2,varargin)

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
        minGap=.5*nanmean([patchEdgeLengths(F1,V1) ; patchEdgeLengths(F2,V2)]);
    case 2
        Q1=varargin{1};
        Q2=varargin{2};
        nOverlap=sum(~isnan(triSurfSetDist(F1,V1,F2,V2,'ray')));
        minGap=.5*nanmean([patchEdgeLengths(F1,V1) ; patchEdgeLengths(F2,V2)]);
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
            minGap=.5*nanmean([patchEdgeLengths(F1,V1) ; patchEdgeLengths(F2,V2)]);
        end
    otherwise
        error('Wrong number of input arguments');
end

%% plot on?
plotOn=0;

%% Parse input
%Ray trace 1 onto 2
optStruct.eps      = 1e-6;
optStruct.triangle = 'two sided';
optStruct.ray      = 'ray';
optStruct.border   = 'normal';

%% Get vertex surface normal vectors
[~,~,N1]=patchNormal(F1,V1);
[~,~,N2]=patchNormal(F2,V2);

%%

CT1=true(size(F1,1),1); % faces to keep from surface 1
CT2=true(size(F2,1),1); % faces to keep from surface 2

% Get boundary edges indices
[Eb1]=patchBoundary(F1(CT1,:),V1);
[Eb2]=patchBoundary(F2(CT2,:),V2);
% unique vertices of boundary edges
indBoundary1=unique(Eb1(:));
indBoundary2=unique(Eb2(:));

logicNoHit1=false(size(V1,1),1); % vertices from surface 1 which do not overlap with surface 2
logicNoHit2=false(size(V2,1),1); % vertices from surface 2 which do not overlap with surface 1

% first remove overlapping surfaces
hw=waitbar(0,'Resolving Overlap...');
meanEdgeLength=nanmean([patchEdgeLengths(F1,V1) ; patchEdgeLengths(F2,V2)]);
icount=0;

if plotOn
    hf=cFigure; hold on
    % subplot(1,2,1); hold on
    hp1=gpatch(F1,V1,'none','k',.5);
    hs1=scatterV(V1(indBoundary1,:),25,Q1(indBoundary1),'Filled');
    axisGeom;
    hc=colorbar;
    title(hc,'distance');
    % subplot(1,2,2); hold on
    hp2=gpatch(F2,V2,'none','r',.5);
    hs2=scatterV(V2(indBoundary2,:),25,Q2(indBoundary2),'Filled');
    axisGeom;
    caxis([0 max([Q1(indBoundary1); Q2(indBoundary2)])]);
    drawnow
end

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
            
            if plotOn
                hp1.Faces=F1(CT1,:);
                hs1.XData=V1(indBoundary1,1);
                hs1.YData=V1(indBoundary1,2);
                hs1.ZData=V1(indBoundary1,3);
                hs1.CData=Q1(indBoundary1);
            end
        else
            logicV2Remove=false(size(V2,1),1); % color all vertices false
            logicV2Remove(indV2in(indBoundaryworstQ2))=true(1); % true if vertices are on the boundary and hit the other surface
            logicF2Remove=~any(logicV2Remove(F2),2); % false for faces for which none of the vertices is in logicV1Remove
            CT2 = CT2 & logicF2Remove; % in CT2 remove (turn to false) boundary face
            
            if plotOn
                hp2.Faces=F2(CT2,:);
                hs2.XData=V2(indBoundary2,1);
                hs2.YData=V2(indBoundary2,2);
                hs2.ZData=V2(indBoundary2,3);
                hs2.CData=Q2(indBoundary2);
            end
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
        
        if plotOn
            hp1.Faces=F1(CT1,:);
            hs1.XData=V1(indBoundary1,1);
            hs1.YData=V1(indBoundary1,2);
            hs1.ZData=V1(indBoundary1,3);
            hs1.CData=Q1(indBoundary1);
        end
        
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
        
        if plotOn
            hp2.Faces=F2(CT2,:);
            hs2.XData=V2(indBoundary2,1);
            hs2.YData=V2(indBoundary2,2);
            hs2.ZData=V2(indBoundary2,3);
            hs2.CData=Q2(indBoundary2);
        end
    end
    
    if isempty(indBoundary1) && isempty(indBoundary2)
        break
    end
    drawnow
end

hs1.XData=[]; hs1.YData=[]; hs1.ZData=[]; hs1.CData=[];
hs2.XData=[]; hs2.YData=[]; hs2.ZData=[]; hs2.CData=[];

%% then remove surfaces that are too close
% logicNoHit1=false(size(V1,1),1); % vertices from surface 1 which do not overlap with surface 2
% logicNoHit2=false(size(V2,1),1); % vertices from surface 2 which do not overlap with surface 1
 if ~isempty(F1(CT1,:)) && ~isempty(F2(CT2,:))
    
X1=V1(:,1);
Y1=V1(:,2);
Z1=V1(:,3);
X2=V2(:,1);
Y2=V2(:,2);
Z2=V2(:,3);

% Get boundary edges indices
[Eb1]=patchBoundary(F1(CT1,:),V1);
[Eb2]=patchBoundary(F2(CT2,:),V2);
% unique vertices of boundary edges
indBoundary1=unique(Eb1(:));
indBoundary2=unique(Eb2(:));

 Eb1_crop=Eb1(all(ismember(Eb1,indBoundary1),2),:);
    V1m=[mean(X1(Eb1_crop),2) mean(Y1(Eb1_crop),2) mean(Z1(Eb1_crop),2)];
    V1m1=[X1(Eb1_crop(:,1))/3+2*X1(Eb1_crop(:,2))/3 Y1(Eb1_crop(:,1))/3+2*Y1(Eb1_crop(:,2))/3 Z1(Eb1_crop(:,1))/3+2*Z1(Eb1_crop(:,2))/3];
    V1m2=[X1(Eb1_crop(:,2))/3+2*X1(Eb1_crop(:,1))/3 Y1(Eb1_crop(:,2))/3+2*Y1(Eb1_crop(:,1))/3 Z1(Eb1_crop(:,2))/3+2*Z1(Eb1_crop(:,1))/3];
    
    Eb2_crop=Eb2(all(ismember(Eb2,indBoundary2),2),:);
    V2m=[mean(X2(Eb2_crop),2) mean(Y2(Eb2_crop),2) mean(Z2(Eb2_crop),2)];
    V2m1=[X2(Eb2_crop(:,1))/3+2*X2(Eb2_crop(:,2))/3 Y2(Eb2_crop(:,1))/3+2*Y2(Eb2_crop(:,2))/3 Z2(Eb2_crop(:,1))/3+2*Z2(Eb2_crop(:,2))/3];
    V2m2=[X2(Eb2_crop(:,2))/3+2*X2(Eb2_crop(:,1))/3 Y2(Eb2_crop(:,2))/3+2*Y2(Eb2_crop(:,1))/3 Z2(Eb2_crop(:,2))/3+2*Z2(Eb2_crop(:,1))/3];
    
    % distance vertices of boundary 1 to vertices of 2
    D1=minDist(V1(indBoundary1,:),[V2(indBoundary2,:); V2m; V2m1; V2m2]);
    D2=minDist(V2(indBoundary2,:),[V1(indBoundary1,:); V1m; V1m1; V1m2]);
    

                hp1.Faces=F1(CT1,:);
        hs1.XData=V1(indBoundary1,1);
        hs1.YData=V1(indBoundary1,2);
        hs1.ZData=V1(indBoundary1,3);
        hs1.CData=D1;
        hp2.Faces=F2(CT2,:);
        hs2.XData=V2(indBoundary2,1);
        hs2.YData=V2(indBoundary2,2);
        hs2.ZData=V2(indBoundary2,3);
        hs2.CData=D2;
        drawnow
        caxis([0 minGap]);
    
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
    %     indBoundary1=indBoundary1(~ismember(indBoundary1,find(logicNoHit1)));
    %     indBoundary2=indBoundary2(~ismember(indBoundary2,find(logicNoHit2)));
    
    % if both boundaries are not empty, find which one has the worst point and remove its face
    %     if ~isempty(indBoundary1) && ~isempty(indBoundary2)
    % find boundary edges and points in the middle of the edges
    Eb1_crop=Eb1(all(ismember(Eb1,indBoundary1),2),:);
    V1m=[mean(X1(Eb1_crop),2) mean(Y1(Eb1_crop),2) mean(Z1(Eb1_crop),2)];
    V1m1=[X1(Eb1_crop(:,1))/4+3*X1(Eb1_crop(:,2))/4 Y1(Eb1_crop(:,1))/4+3*Y1(Eb1_crop(:,2))/4 Z1(Eb1_crop(:,1))/4+3*Z1(Eb1_crop(:,2))/4];
    V1m2=[X1(Eb1_crop(:,2))/4+3*X1(Eb1_crop(:,1))/4 Y1(Eb1_crop(:,2))/4+3*Y1(Eb1_crop(:,1))/4 Z1(Eb1_crop(:,2))/4+3*Z1(Eb1_crop(:,1))/4];
    
    Eb2_crop=Eb2(all(ismember(Eb2,indBoundary2),2),:);
    V2m=[mean(X2(Eb2_crop),2) mean(Y2(Eb2_crop),2) mean(Z2(Eb2_crop),2)];
    V2m1=[X2(Eb2_crop(:,1))/4+3*X2(Eb2_crop(:,2))/4 Y2(Eb2_crop(:,1))/4+3*Y2(Eb2_crop(:,2))/4 Z2(Eb2_crop(:,1))/4+3*Z2(Eb2_crop(:,2))/4];
    V2m2=[X2(Eb2_crop(:,2))/4+3*X2(Eb2_crop(:,1))/4 Y2(Eb2_crop(:,2))/4+3*Y2(Eb2_crop(:,1))/4 Z2(Eb2_crop(:,2))/4+3*Z2(Eb2_crop(:,1))/4];
    
    % distance vertices of boundary 1 to vertices of 2
    D1=minDist(V1(indBoundary1,:),[V2(indBoundary2,:); V2m; V2m1; V2m2]);
    [logicDist1]=D1<minGap; % true if small distance
    c1=logicDist1; % true if the boundary vertices of 1 very close
    logicV1in=false(size(V1,1),1); % color all vertices false
    logicV1in(indBoundary1)=c1; % true if vertices are on the boundary and are close the other surface
    indV1in=find(logicV1in);
    [worstQ1,indBoundaryworstQ1]=max(Q1(logicV1in)); % vertex with worst score
    
    % distance vertices of boundary 1 to vertices of 2
    D2=minDist(V2(indBoundary2,:),[V1(indBoundary1,:); V1m; V1m1; V1m2]);
    [logicDist2]=D2<minGap; % true if small distance
    c2=logicDist2; % true if the boundary vertices of 1 hit the other surface or very close
    logicV2in=false(size(V2,1),1); % color all vertices false
    logicV2in(indBoundary2)=c2; % true if vertices are on the boundary and hit the other surface
    indV2in=find(logicV2in);
    [worstQ2,indBoundaryworstQ2]=max(Q2(logicV2in)); % vertex with worst score
    
    if sum(logicDist1)>0 && sum(logicDist2)>0
        if worstQ1>worstQ2
            logicV1Remove=false(size(V1,1),1); % color all vertices false
            logicV1Remove(indV1in(indBoundaryworstQ1))=true(1); % true if vertices are on the boundary and hit the other surface
            logicF1Remove=~any(logicV1Remove(F1),2); % false for faces for which none of the vertices is in logicV1Remove
            CT1 = CT1 & logicF1Remove; % in CT1 remove (turn to false) boundary face
            %             logicNoHit1(indBoundary1(~c1))=1; % add vertices to boundary no hit list
        else
            logicV2Remove=false(size(V2,1),1); % color all vertices false
            logicV2Remove(indV2in(indBoundaryworstQ2))=true(1); % true if vertices are on the boundary and hit the other surface
            logicF2Remove=~any(logicV2Remove(F2),2); % false for faces for which none of the vertices is in logicV1Remove
            CT2 = CT2 & logicF2Remove; % in CT2 remove (turn to false) boundary face
            %             logicNoHit2(indBoundary2(~c2))=1; % add vertices to no hit list
        end
    elseif  sum(logicDist1)>0
         logicV1Remove=false(size(V1,1),1); % color all vertices false
            logicV1Remove(indV1in(indBoundaryworstQ1))=true(1); % true if vertices are on the boundary and hit the other surface
            logicF1Remove=~any(logicV1Remove(F1),2); % false for faces for which none of the vertices is in logicV1Remove
            CT1 = CT1 & logicF1Remove; % in CT1 remove (turn to false) boundary face
    elseif sum(logicDist2)>0
         logicV2Remove=false(size(V2,1),1); % color all vertices false
            logicV2Remove(indV2in(indBoundaryworstQ2))=true(1); % true if vertices are on the boundary and hit the other surface
            logicF2Remove=~any(logicV2Remove(F2),2); % false for faces for which none of the vertices is in logicV1Remove
            CT2 = CT2 & logicF2Remove; % in CT2 remove (turn to false) boundary face
    else
        break
    end
    
    if plotOn
        hp1.Faces=F1(CT1,:);
        hs1.XData=V1(indBoundary1,1);
        hs1.YData=V1(indBoundary1,2);
        hs1.ZData=V1(indBoundary1,3);
        hs1.CData=D1;
        hp2.Faces=F2(CT2,:);
        hs2.XData=V2(indBoundary2,1);
        hs2.YData=V2(indBoundary2,2);
        hs2.ZData=V2(indBoundary2,3);
        hs2.CData=D2;
        drawnow
    end
end

 end
 
if plotOn

%     close(hf)
end

close(hw)

end
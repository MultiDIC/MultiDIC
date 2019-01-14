function [F,V,C]=delaunayZip(F1,V1,F2,V2,inputStruct)

% function [Fn]=delaunayZip(F1,V1,F2,V2,inputStruct)
%-------------------------------------------------------------------------
%
%
% Change log:
% 2018/05/09: Added growing of region size if error occurs
% 2018/05/09: Added stepwise removal of last triangle since it does not
% include next points.
% 2018/11/02: Added self-triangulation as an option (default on) 
% 2018/11/02: Added walk-back for paths so that one will not be too long
% 2018/11/02: Fixed bug in relation to one of the paths "running out of
% points". Remainder is triangulated properly. 
%-------------------------------------------------------------------------

%%

maxD=max([patchEdgeLengths(F1,V1);patchEdgeLengths(F2,V2)]);

defaultInputStruct.ind1=[];
defaultInputStruct.ind2=[];
defaultInputStruct.distLocal=2*maxD;
defaultInputStruct.startInd=[];
defaultInputStruct.plotOn=0;
defaultInputStruct.selfTriangulate=1; %Option to self-triangulate surfaces first
defaultInputStruct.angleThreshold=(60/180)*pi; %Angular threshold for self-triangulation
[inputStruct]=structComplete(inputStruct,defaultInputStruct,0); %Complement provided with default if missing or empty

ind1=inputStruct.ind1;
ind2=inputStruct.ind2;
distLocal=inputStruct.distLocal;
startInd=inputStruct.startInd;
selfTriangulate=inputStruct.selfTriangulate;
angleThreshold=inputStruct.angleThreshold;
plotOn=inputStruct.plotOn;

%% Self-triangulate if needed

if selfTriangulate==1
    if ind1(1)==ind1(end)
        isClosedLoop=1;
    else
        isClosedLoop=0;
    end
    numFacesInitial_1=size(F1,1);
    [F1,V1,ind1]=triSurfSelfTriangulateBoundary(F1,V1,ind1,angleThreshold,isClosedLoop);
    numFaces_1=size(F1,1);
    numFacesInitial_2=size(F2,1);
    [F2,V2,ind2]=triSurfSelfTriangulateBoundary(F2,V2,ind2,angleThreshold,isClosedLoop);
    numFaces_2=size(F2,1);  
    
    L1=false(size(F1,1),1);
    if numFaces_1>numFacesInitial_1
        L1(end-(numFaces_1-numFacesInitial_1-1):end)=1;
    end
    L2=false(size(F2,1),1);
    if numFaces_2>numFacesInitial_2
        L2(end-(numFaces_2-numFacesInitial_2-1):end)=1;
    end
else
    L1=zeros(size(F1,1),1);
    L2=zeros(size(F2,1),1);
end

%%
[F,V,C]=joinElementSets({F1,F2},{V1,V2});
[~,~,Nv]=patchNormal(F,V);

%%

ind2=ind2+size(V1,1);
startInd(2)=startInd(2)+size(V1,1);

%%
%Create edges list
E=[ind1(1:end-1) ind1(2:end); ind2(1:end-1) ind2(2:end)];

%%

logicNotUsed=ismember((1:1:size(V,1))',E); %Logic to keep track of points that are used

%Initialize other parameters
Fn=[]; %Faces
Cn=[]; %"Colors"=step count for face group
c=1; %While loop counter variable

indGroup=startInd(:); %Initialize current group
numGroup=numel(indGroup); %Initialize current number of members of the current group
numGroupPrevious=0;

if plotOn==1
    h1=[]; %Initiate empty plot handle
    markerSize=25;
    
    hf=cFigure; hold on;
    gpatch(F1(~L1,:),V1,'r','none',0.2);
    gpatch(F1(L1,:),V1,'rw','k',1);
    gpatch(F2(~L2,:),V2,'b','none',0.2);
    gpatch(F2(L2,:),V2,'bw','k',1);
    plotV(V(ind1,:),'r.-','LineWidth',1,'MarkerSize',10);
    plotV(V(ind2,:),'b.-','LineWidth',1,'MarkerSize',10);
    
    axisGeom;
    camlight headlight;
    colormap(gjet(250));
    colorbar;
    drawnow;
end

%Turn off Delaunay warning as this case is handled properly
warning('off','MATLAB:delaunayTriangulation:ConsConsSplitWarnId');

%%
while 1
    
    numGroupStep=1;
    lastTry=0;
    while 1
        try
            if plotOn==1 %%Plot if plotting is on
                delete(h1); h1=[];
            end
            
            if plotOn==1 %%Plot if plotting is on
                figure(hf);
                h1(end+1)=plotV(V(startInd,:),'y.','MarkerSize',markerSize);
                drawnow
            end
            
            %Grow the current region
            while 1 %Loop to form local group (stuff attached to current group within a given distance)
                logicMember=any(ismember(E,indGroup),2); %Logic for all edges touching the current groupt
                E_sub=E(logicMember,:); %The subset of touching edges
                
                indGroup=unique([indGroup; E_sub(:)]); %Grow group with point indices in the edges that are touching
                
                d1=sqrt(sum((V(indGroup,:)-V(startInd(1)*ones(numel(indGroup),1),:)).^2,2));
                d2=sqrt(sum((V(indGroup,:)-V(startInd(2)*ones(numel(indGroup),1),:)).^2,2));
                D=min(d1,d2);
                
                logicKeep= (D<distLocal) & (logicNotUsed(indGroup));
                indGroup=indGroup(logicKeep); %Remove points that are too far
                if numGroup==numel(indGroup) %Compare current group size to previous step
                    break %break while loop if the group is no longer growing
                end
                numGroup=numel(indGroup); %Get new current group size
            end
            
            %Get curve start and end points
            logicMember=all(ismember(E,indGroup),2); %Logic for all edges touching the current group
            E_sub=E(logicMember,:); %The subset of touching edges
            
            [~,~,~,vCount]=cunique(E_sub); %Get vertex occurance counts
            indEndPoints=unique(E_sub(vCount==1));
            indEndPoints1=indEndPoints(ismember(indEndPoints,ind1));
            indEndPoints2=indEndPoints(ismember(indEndPoints,ind2));
            
            %Compose sub-curves
            logicSub1=all(ismember(E_sub,ind1),2);
            %             if ~any(logicSub1)
            %                 warning('None of the current edges nodes are a member of ind1');
            %             end
            E_sub1=E_sub(logicSub1,:);
            [indListSub1]=edgeListToCurve(E_sub1);
            indListSub1=indListSub1(:);
            
            logicSub2=all(ismember(E_sub,ind2),2);
            %             if ~any(logicSub2)
            %                 warning('None of the current edges nodes are a member of ind2');
            %             end
            E_sub2=E_sub(logicSub2,:);
            [indListSub2]=edgeListToCurve(E_sub2);
            indListSub2=indListSub2(:);
            
            %Alter segments so they can be merged into closed loop curve
            if isempty(Fn) %First iteration
                d=sqrt(sum((V([indListSub2(1) indListSub2(end)],:)-V(indListSub1(1)*ones(1,2),:)).^2,2));
                [~,indClosest]=min(d);
                if indClosest~=2 %If the start is closest to the other start
                    indListSub2=flipud(indListSub2);
                end
            else
                indStart1=find(ismember(indListSub1,startInd));
                if indStart1~=1
                    indListSub1=flipud(indListSub1);
                end
                indStart2=find(ismember(indListSub2,startInd));
                if indStart2==1
                    indListSub2=flipud(indListSub2);
                end
            end
            
            if any(ismember(indGroup,ind1))==0
                startInd
                ind1
                indGroup
            end
            
            %             if any(ismember(indGroup,ind1)) && any(ismember(indGroup,ind2))
            %             else
            %                 %One of the two sets contains only 1 member, triangulate
            %                 %remaining points manually
            %
            %                 %-----------------------------------
            %                 if nnz(ismember(indGroup,ind1))==0
            %                     indLeft=flipud(indListSub2);
            %                     sKeep=startInd(ismember(startInd,ind1));
            %                 else%if nnz(ismember(indGroup,ind2))==0
            %                     indLeft=indListSub1;
            %                     sKeep=startInd(ismember(startInd,ind2));
            %                 end
            %
            %                 f=[];
            %                 for qStep=1:1:numel(indLeft)-1
            %                     sVary=indLeft(1);
            %                     f=[f; [sKeep sVary indLeft(2)]];
            %                     indLeft=indLeft(2:end); %First is part of startInd
            %                 end
            %             end
            
            
            %Shorten curves if they stick out to much wrt end of other curve
            d=sqrt(sum((V(indListSub2,:)-V(indListSub1(end)*ones(1,numel(indListSub2)),:)).^2,2));
            [~,indClosest]=min(d);
            if indClosest~=1 && indClosest~=numel(indListSub2)
                indListSub2=indListSub2(indClosest:end);
            end
            
            d=sqrt(sum((V(indListSub1,:)-V(indListSub2(1)*ones(1,numel(indListSub1)),:)).^2,2));
            [~,indClosest]=min(d);
            if indClosest~=numel(indListSub1) && indClosest~=1
                indListSub1=indListSub1(1:indClosest);
            end
            
            %Create closed curve
            if nnz(ismember(indGroup,ind1))==0
                indListSub=[startInd(ismember(startInd,ind1)); indListSub2; startInd(ismember(startInd,ind1))]; %Closed loop
            elseif nnz(ismember(indGroup,ind2))==0
                indListSub=[startInd(ismember(startInd,ind2)); indListSub1; startInd(ismember(startInd,ind2))]; %Closed loop
            else
                indListSub=[indListSub1;indListSub2;indListSub1(1)]; %Closed loop
            end
            
            indGroup=indListSub(1:end-1);
            
            if plotOn==1 %%Plot if plotting is on
                figure(hf);
                if ~isempty(indEndPoints1)
                    h1(end+1)=plotV(V(indEndPoints1,:),'r.','MarkerSize',markerSize);
                end
                if ~isempty(indEndPoints2)
                    h1(end+1)=plotV(V(indEndPoints2,:),'b.','MarkerSize',markerSize);
                end
                if ~isempty(indListSub1)
                    h1(end+1)=plotV(V(indListSub1,:),'r-','LineWidth',3);
                end
                if ~isempty(indListSub2)
                    h1(end+1)=plotV(V(indListSub2,:),'b-','LineWidth',3);
                end
                if ~isempty(indListSub)
                    h1(end+1)=plotV(V(indListSub,:),'k-','LineWidth',2);
                end
                drawnow
            end
            
            %Rotate current point set
            V_now=V(indListSub(1:end-1),:); %Current closed curve coordinate set
            [R]=pointSetPrincipalDir(V_now); %Fit local coordinate system with 3rd direction pointing outward of local planar-ish region
            V_now_R=V_now*R; %Rotate coordinate set to prepare for 2D Delaunay based triangulation
            
            %Do 2D Delaunay triangulation
            DT_contraints=[(1:numel(indGroup))' ([2:numel(indGroup) 1])']; %Constraints are edges forming boundary
            DT = delaunayTriangulation(V_now_R(:,[1 2]),DT_contraints); % Initial triangulation
            f=DT.ConnectivityList; %Get faces set
            L = isInterior(DT); %Remove faces not inside region
            f=f(L,:); %Faces excluding external faces (only keep those inside constraint edges)
            f=indGroup(f); %Change indices to overall system
            
            %Remove first and last triangles if there are more than 3 triangles, this
            %will improve triangulation quality as these triangles were not fully
            %embedded in the point set.
            %                 if size(f,1)>2
            %                     logicKeep=~(any(ismember(f,indEndPoints1(~ismember(indEndPoints1,Fn))),2)...
            %                         & any(ismember(f,indEndPoints2(~ismember(indEndPoints2,Fn))),2))...
            %                         | any(ismember(f,startInd),2);
            %                     f=f(logicKeep,:);
            %
            %                     [~,~,IND_FF]=tesIND(f,V);
            %                     logicNeighbours=sum(IND_FF>1,2)>1;
            %                     logicNeighbours(any(ismember(f,startInd),2))=1;
            %                     f=f(logicNeighbours,:);
            %                 end
            
            break
            
        catch ME
            if lastTry==1
                rethrow(ME)
            end
            
            if numGroupStep==numGroupPrevious
                lastTry=1;
            end
            distLocal=distLocal+maxD;
            warning(['Failed using current step size, increasing to: ',num2str(distLocal)]);
            numGroupPrevious=numGroupStep;
        end
    end
    distLocal=inputStruct.distLocal; %reset
    
    if dot(mean(patchNormal(f,V),1),mean(Nv(indListSub(1:end-1),:),1))<0
        f=fliplr(f);
    end
    
    if plotOn==1 %%Plot if plotting is on
        figure(hf);
        V_now_mean=mean(V_now); %Mean of coordinate set
        h1(end+1)=plotV(V_now_mean,'kx','MarkerSize',markerSize);
        gpatch(f,V,c*ones(size(f,1),1),'k',1);
        drawnow
    end
    
    %Collect faces and color data
    Fn=[Fn;f];
    Cn=[Cn;c*ones(size(f,1),1)];
    
    E_Fn=sort(patchEdges(Fn,0),2);
    ind_E_Fn= reshape(sub2indn(size(V,1)*ones(1,2),E_Fn),size(Fn));
    
    Es=sort(E,2);
    ind_E= sub2indn(size(V,1)*ones(1,2),Es);
    
    logicEdgesUsed=ismember(ind_E,ind_E_Fn);
    
    E_sub=E(~logicEdgesUsed,:);
    indUnusedPoints=unique(E_sub);
    logicNotUsed=false(size(logicNotUsed));
    logicNotUsed(indUnusedPoints)=1;
    
    if any(logicNotUsed)==0
        break
    end
    
    %Get curve start and end points
    [~,~,~,vCount]=cunique(E_sub); %Get vertex occurance counts
    indEndPoints=unique(E_sub(vCount==1));
    indEndPoints1=indEndPoints(ismember(indEndPoints,ind1));
    indEndPoints2=indEndPoints(ismember(indEndPoints,ind2));
    
    e=patchEdges(Fn,1);
    edgeEndPoints=e(any(ismember(e,indEndPoints1),2)&any(ismember(e,indEndPoints2),2),:);
    logicMemberOfLast=any(ismember(edgeEndPoints,f),2);
    edgeEndPoints=edgeEndPoints(logicMemberOfLast,:);
    if ~isempty(edgeEndPoints)
        startInd=edgeEndPoints(1,:);
    else
        startInd=[indListSub1(end) indListSub2(1)];
    end
    
    c=c+1;
end

F=[F;Fn];
C=[C; max(C(:))+Cn];

%Turn Delaunay warning back on
warning('on','MATLAB:delaunayTriangulation:ConsConsSplitWarnId');



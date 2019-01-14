function [DIC3DStitched]= DIC3DsurfaceStitch(DIC3DAllPairsResults,pairIndList)

% pairIndList = lists of surfaces indices to stitch

if ~isempty(pairIndList)
    %% assign first surfaces for stitching
    
    nPairs=numel(DIC3DAllPairsResults);
    nPairsForStitching=size(pairIndList,2);
    
    pairInd1=pairIndList(1);
    colors=gjet(nPairs+nPairsForStitching-1);
    
    % surfaces
    s1=DIC3DAllPairsResults{pairInd1};
    DIC3DStitched=struct;
    DIC3DStitched.Faces=s1.Faces;
    DIC3DStitched.Points3D=s1.Points3D;
    DIC3DStitched.corrComb=s1.corrComb;
    DIC3DStitched.FaceColors=s1.FaceColors;
    DIC3DStitched.FacePairInds=pairInd1*ones(size(s1.Faces,1),1);
    DIC3DStitched.pairIndices=zeros(nPairs,2);
    DIC3DStitched.pairIndices(pairInd1,:)=s1.cameraPairInd;
    DIC3DStitched.PointPairInds=pairInd1*ones(size(s1.Points3D{1},1),1);
    
    numFreeBytes=freeMemory;
    for imesh=1:nPairsForStitching-1
        
        pairInd2=pairIndList(imesh+1);
        
        % surfaces
        s1=DIC3DStitched;
        s2=DIC3DAllPairsResults{pairInd2};
        
        %faces
        F1=s1.Faces;
        F2=s2.Faces;
        
        %vertices
        V1=s1.Points3D{1};
        V2=s2.Points3D{1};
        
        %face Colors as pair index
        FC1=s1.FacePairInds;
        FC2=(pairInd2)*ones(size(F2,1),1);
        
        %face Colors as texture
        FCT1=s1.FaceColors;
        FCT2=s2.FaceColors;
        
        % remove faces with NaN
        F1NaNLogic=any(isnan(V1(F1)),2);
        F2NaNLogic=any(isnan(V2(F2)),2);
        F1(F1NaNLogic,:)=[];
        F2(F2NaNLogic,:)=[];
        FC1(F1NaNLogic,:)=[];
        FC2(F2NaNLogic,:)=[];
        FCT1(F1NaNLogic,:)=[];
        FCT2(F2NaNLogic,:)=[];
        
        %% plot original surfaces
        hf(imesh)=cFigure; hold all
        gtitle(['Stitching [' num2str(pairIndList(1):pairIndList(imesh)) '] to [' num2str(pairInd2) ']'],20);
        
        subplot(1,3,1); hold all
        gpatch(F1,V1,colors(FC1,:),'k',.5);
        gpatch(F2,V2,colors(FC2,:),'k',.5);
        axisGeom
        %         legend({num2str(pairIndList(1):pairIndList(imesh)), num2str(pairInd2)});
        
        %%
        % reduced surfaces
        minGap=.4*nanmean([patchEdgeLengths(F1,V1) ; patchEdgeLengths(F2,V2)]);
        minDistValue=5*minGap;
        
        % Q1=
        % Q2=
        % [CT1,CT2] = removeOverlapSurface(F1,F2,V1,V2,Q1,Q2,minGap);
        %         [CT1,CT2] = removeOverlapSurface(F1,F2,V1,V2,[],[],minGap);
        [CT1,CT2] = removeOverlapSurface_temp_2018_10_30(F1,F2,V1,V2,[],[],minGap);
        
        %reduced faces
        F1r=F1(CT1,:);
        F2r=F2(CT2,:);
        FC1r=FC1(CT1,:);
        FC2r=FC2(CT2,:);
        FCT1r=FCT1(CT1,:); % colors that represent the pair index
        FCT2r=FCT2(CT2,:);
        % remove NaN faces
        F1rNaNLogic=any(isnan(V1(F1r)),2);
        F2rNaNLogic=any(isnan(V2(F2r)),2);
        F1r(F1rNaNLogic,:)=[];
        F2r(F2rNaNLogic,:)=[];
        FC1r(F1rNaNLogic,:)=[];
        FC2r(F2rNaNLogic,:)=[];
        FCT1r(F1rNaNLogic,:)=[];
        FCT2r(F2rNaNLogic,:)=[];
        
        % assign NaN Vertices
        V1r=V1;
        V1rInds=unique(F1r(:));
        V1rlogic=false(size(V1,1),1);
        V1rlogic(V1rInds)=true;
        V1r(~V1rlogic,:)=NaN;
        V2r=V2;
        V2rInds=unique(F2r(:));
        V2rlogic=false(size(V2,1),1);
        V2rlogic(V2rInds)=true;
        V2r(~V2rlogic,:)=NaN;
        
        % boundaries of reduced surfaces
        E1=patchBoundary(F1r,V1r);
        E2=patchBoundary(F2r,V2r);
        
        %% remove faces for which all 3 edges are on the border
        indAllEdges1=[];
        for ii=1:size(F1r,1)
            
            if (ismember([F1r(ii,1) F1r(ii,2)],E1,'rows') || ismember([F1r(ii,2) F1r(ii,1)],E1,'rows')) && (ismember([F1r(ii,2) F1r(ii,3)],E1,'rows') || ismember([F1r(ii,3) F1r(ii,2)],E1,'rows')) && (ismember([F1r(ii,1) F1r(ii,3)],E1,'rows') || ismember([F1r(ii,3) F1r(ii,1)],E1,'rows'))
                indAllEdges1=[indAllEdges1 ii];
            end
            
        end
        F1r(indAllEdges1,:)=[];
        FC1r(indAllEdges1,:)=[];
        FCT1r(indAllEdges1,:)=[];
        V1r=V1;
        V1rInds=unique(F1r(:));
        V1rlogic=false(size(V1,1),1);
        V1rlogic(V1rInds)=true;
        V1r(~V1rlogic,:)=NaN;
        
        indAllEdges2=[];
        for ii=1:size(F2r,1)
            
            if (ismember([F2r(ii,1) F2r(ii,2)],E2,'rows') || ismember([F2r(ii,2) F2r(ii,1)],E2,'rows')) && (ismember([F2r(ii,2) F2r(ii,3)],E2,'rows') || ismember([F2r(ii,3) F2r(ii,2)],E2,'rows')) && (ismember([F2r(ii,1) F2r(ii,3)],E2,'rows') || ismember([F2r(ii,3) F2r(ii,1)],E2,'rows'))
                indAllEdges2=[indAllEdges2 ii];
            end
            
        end
        F2r(indAllEdges2,:)=[];
        FC2r(indAllEdges2,:)=[];
        FCT2r(indAllEdges2,:)=[];
        V2r=V2;
        V2rInds=unique(F2r(:));
        V2rlogic=false(size(V2,1),1);
        V2rlogic(V2rInds)=true;
        V2r(~V2rlogic,:)=NaN;
        
        % new boundaries
        E1=patchBoundary(F1r,V1r);
        E2=patchBoundary(F2r,V2r);
        
        %% plot reduced surfaces
        set(0, 'currentfigure', hf(imesh));
        subplot(1,3,2); hold all
        gpatch(F1r,V1r,colors(FC1r,:),'k',.5);
        gpatch(F2r,V2r,colors(FC2r,:),'k',.5);
        axisGeom
        
        % % plot boundary edges
        % he1=gpatch(E1,V1r,'none','c',.5); he1.LineWidth=2;
        % he2=gpatch(E2,V2r,'none','m',.5); he2.LineWidth=2;
        
        % plot each group seperately
        Fcombined=[F1r;F2r+size(V1r,1)];
        Ccombined=[FC1r; FC2r];
        FCTcombined=[FCT1r; FCT2r];
        
        if ~isempty(F1r) && ~isempty(F2r)
        %% check if the boundaries form more than 1 closed group
        indVertexBowtied1=findIndVertexBowtied(F1r,V1r);
        indVertexBowtied2=findIndVertexBowtied(F2r,V2r);
        
        optionStruct=struct;
        optionStruct.indExclude=indVertexBowtied1;
        G1=tesgroup(E1,optionStruct);
        optionStruct=struct;
        optionStruct.indExclude=indVertexBowtied2;
        G2=tesgroup(E2,optionStruct);

%         G1=tesgroup(E1);
%         G2=tesgroup(E2);
        nBorders1=size(G1,2);
        nBorders2=size(G2,2);
        
        for ib1=1:nBorders1
            for ib2=1:nBorders2

                E1now=E1(G1(:,ib1),:);
                E2now=E2(G2(:,ib2),:);
                
                hp1=gpatch(E1now,V1r,'c','c',.5); hp1.LineWidth=2;
                hp2=gpatch(E2now,V2r,'g','g',.5); hp2.LineWidth=2;
                
                %% find closest points to start marching from (from 1 to 2)
                c1=edgeListToCurve(E1now)';
                c1(1)=[];
                
                c2=edgeListToCurve(E2now)';
                c2(1)=[];
                
                % for each point in the boundary of 1, find the closest point in the boundary of 2
                [D12,~]=minDist(V1r(c1,:),V2r(c2,:),[],0,numFreeBytes);
                
                %% find the edges which are close and group them
                
                D12V=zeros(size(V1r,1),1);
                D12V(c1)=D12;
                LogicD12=all(D12V(E1now)<minDistValue,2);
                E1close=E1now(LogicD12,:);
                if ~isempty(E1close)
                    Groups1=tesgroup(E1close);
                    % delete groups that have only one memeber
                    Groups1(:,sum(Groups1)==1)=[];
                    numGroups1=size(Groups1,2);
                    
                    Fz=cell(numGroups1,1);
                    for ii=1:numGroups1
                        
                        curve1=edgeListToCurve(E1close(Groups1(:,ii),:))';
                        curve1=unique(curve1,'stable');
                        plotV(V1r(curve1,:),'sk','MarkerFaceColor','r');
                        
                        % find the vertices on the other mesh that are close to this group
                        [~,indD12]=minDist(V1r(curve1,:),V2r(c2,:),[],0,numFreeBytes);
                        
                        indD12Diff=diff(indD12);
                        if sum(sign(indD12Diff))<0 % if series is going down, flip it
                            indD12=flip(indD12);
                        end
                        
                        % take all numbers between extremes
                        if indD12(1)<=indD12(end)
                            %                                 indD12=indD12(1):indD12(end);
                            indD12=min(indD12):max(indD12);
                        else
                            %                                 indD12=[indD12(1):length(c2) 1:indD12(end)];
                            [~,indMin]=min(indD12);
                            min1=min(indD12(1:indMin-1));
                            max1=max(indD12(indMin:end));
                            indD12=[min1:length(c2) 1:max1];
                        end
                        
                        curve2=c2(indD12);
                        plotV(V2r(curve2,:),'ok','MarkerFaceColor','y');
                        
                        % flip curve 2 if V2(curve2(1),:) is closer to V1(curve1(end),:) than V1(curve1(1),:)
                        distV21V11=pdist([V1(curve1(1),:);V2(curve2(1),:)],'euclidean');
                        distV21V1end=pdist([V1(curve1(end),:);V2(curve2(1),:)],'euclidean');
                        if distV21V1end<distV21V11
                            curve2=flip(curve2);
                        end
                        
                        % retract vertices from the beginning and end of the lists if they are further
                        while length(curve1)>=2 && length(curve2)>=2
                            if pdist2(V1r(curve1(1),:),V2r(curve2(1),:))>pdist2(V1r(curve1(2),:),V2r(curve2(1),:))
                                curve1(1)=[];
                                %             plotV(V1r(curve1(1),:),'sk','MarkerSize',10,'MarkerFaceColor','b'); % plot closest point on mesh 2
                            elseif pdist2(V1r(curve1(1),:),V2r(curve2(1),:))>pdist2(V1r(curve1(1),:),V2r(curve2(2),:))
                                curve2(1)=[];
                                %             plotV(V2r(curve2(1),:),'sk','MarkerSize',10,'MarkerFaceColor','r'); % plot closest point on mesh 2
                            else
                                break
                            end
                        end
                        while length(curve1)>=2 && length(curve2)>=2
                            if pdist2(V1r(curve1(end),:),V2r(curve2(end),:))>pdist2(V1r(curve1(end-1),:),V2r(curve2(end),:))
                                curve1(end)=[];
                                %             plotV(V1r(curve1(end),:),'sk','MarkerSize',10,'MarkerFaceColor','b'); % plot closest point on mesh 2
                            elseif pdist2(V1r(curve1(end),:),V2r(curve2(end),:))>pdist2(V1r(curve1(end),:),V2r(curve2(end-1),:))
                                curve2(end)=[];
                                %             plotV(V2r(curve2(end),:),'sk','MarkerSize',10,'MarkerFaceColor','r'); % plot closest point on mesh 2
                            else
                                break
                            end
                        end
                        
                        if length(curve1)>=2 && length(curve2)>=2
                            plotV(V1r(curve1(1),:),'sk','MarkerSize',10,'MarkerFaceColor','k');
                            plotV(V2r(curve2(1),:),'sk','MarkerSize',10,'MarkerFaceColor','w');
                            plotV(V1r(curve1(end),:),'^k','MarkerSize',10,'MarkerFaceColor','k');
                            plotV(V2r(curve2(end),:),'^k','MarkerSize',10,'MarkerFaceColor','w');
                            
                            %Create edges list
                            ind1=curve1;
                            ind2=(size(V1r,1)+curve2);
                            ind=[ind1;flipud(ind2);ind1(1)];
                            V=[V1r;V2r];
                            E=[ind(1:end-1) ind(2:end)];
                            
                            set(0, 'currentfigure', hf(imesh));
                            hp=gpatch(E,V,'none','m',.5); hp.LineWidth=2;
                            
                            startInd=[curve1(1) curve2(1)]';
                            
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % stitching using delaunayZip
                            %Create input structure
                            inputStruct=struct;
                            inputStruct.plotOn=1;
                            inputStruct.distLocal=4*max([max(patchEdgeLengths(F1r,V1r)) max(patchEdgeLengths(F2r,V2r))]);
                            inputStruct.startInd=startInd;
                            inputStruct.ind1=curve1;
                            inputStruct.ind2=curve2;
                            
                            try
                                % zip
                                [FzTemp,~,CzTemp]=delaunayZip(F1r,V1r,F2r,V2r,inputStruct);
                                Fz{ii}=FzTemp(CzTemp>=3,:);
                            catch
                                warning(['Problem zipping [' num2str(pairIndList(1):pairIndList(imesh)) '] to [' num2str(pairInd2) '] . Continuing without zipping']);
                                % take the zipped part (without the original surfaces)
                                Fz{ii}=[]; 
%                                 [num2str(imesh) '_' num2str(ib1) '_' num2str(ib2)  '_' num2str(ii)]
%                                 save(['C:\Users\Dana\Dropbox (Personal)\Research\DANA_KEVIN\MATLAB\Projects\surfaceStitching\zipInput_' num2str(imesh) '_' num2str(ib1) '_' num2str(ib2)  '_' num2str(ii)],'F1r','V1r','F2r','V2r','inputStruct'); 
                            end
                            
                            set(0, 'currentfigure', hf(imesh));
                            %                     gpatch(Fz{ii},[V1r;V2r],CzTemp(CzTemp>=3,:),'k',.5);
                            
                            Fcombined=[Fcombined; Fz{ii}];
                            Ccombined=[Ccombined; (nPairs+imesh)*ones(size(Fz{ii},1),1)];
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            % sample face colors from neighboring faces into the stitched faces
                            FCTz=zeros(size(Fz{ii},1),1);
                            for iif=1:size(FCTz,1)
                                % find touching faces
                                Ftemp=Fz{ii};
                                
                                Etemp=[Ftemp(:,1) Ftemp(:,2); Ftemp(:,2) Ftemp(:,3); Ftemp(:,3) Ftemp(:,1)];
                                logicTouch=any(reshape(sum(ismember(Etemp,Ftemp(iif,:)),2)==2,size(Ftemp,1),3),2);
                                logicTouch(iif)=0;
                                FCTz(iif)=mean(FCTcombined(logicTouch));
                                
                            end
                            FCTcombined=[FCTcombined; FCTz];
                        end
                    end
                    
                end
            end
            
        end
        
        Vcombined=[V1r;V2r];
        DIC3DStitched.PointPairInds=[DIC3DStitched.PointPairInds; pairInd2*ones(size(V2r,1),1)];
        
        % plot stitched surface
        set(0, 'currentfigure', hf(imesh));
        hs3=subplot(1,3,3); hold all
        gpatch(Fcombined,Vcombined,Ccombined,'k',1);
        axisGeom;
        colormap gjet
        icolorbar([1 nPairs+nPairsForStitching-1])
        
        nT=numel(DIC3DStitched.Points3D);
        
        for it=1:nT
            Points3D{it}=[s1.Points3D{it}; s2.Points3D{it}];
            Points3D{it}(any(isnan(Vcombined),2),:)=NaN;
            corrComb{it}=[s1.corrComb{it}; s2.corrComb{it}];
        end
        end 
        %% update DIC3DStitched
        
        DIC3DStitched.Faces=Fcombined;
        DIC3DStitched.Points3D=Points3D;
        DIC3DStitched.corrComb=corrComb;
        DIC3DStitched.FaceColors=FCTcombined;
        DIC3DStitched.FacePairInds=Ccombined;
        DIC3DStitched.pairIndices(pairInd2,:)=s2.cameraPairInd;
        
    end
    
    %% Fill holes - by groups of boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F=DIC3DStitched.Faces;
    V=DIC3DStitched.Points3D{1};
    [Eb,E,indBoundary]=patchBoundary(F,V); % find boundaries
    
    %     gpatch(Eb,V,'none','m',1,5);
    
    [G,G_iter]=tesgroup(Eb); % group boundaries
    for ii=1:size(G,2) % loop over all groups of closed boundaries
        EbNow=Eb(G(:,ii),:);
        if size(EbNow,1)==3 % if the hole is the size of one triangle
            F(end+1,:)=unique(EbNow)' ;% add this face to face list
            [IND_F,IND_V,IND_FF]=tesIND(F,V); % connectivity
            % find the faces that are connected to both vertex 1 and 2
            indF12=intersect(nonzeros(IND_F(F(end,1),:)),nonzeros(IND_F(F(end,2),:)));
            indF12(indF12==size(F,1),:)=[];
            % find the faces that are connected to both vertex 1 and 3
            indF13=intersect(nonzeros(IND_F(F(end,1),:)),nonzeros(IND_F(F(end,3),:)));
            indF13(indF13==size(F,1),:)=[];
            % find the faces that are connected to both vertex 2 and 3
            indF23=intersect(nonzeros(IND_F(F(end,3),:)),nonzeros(IND_F(F(end,2),:)));
            indF23(indF23==size(F,1),:)=[];
            
            % assign the median pairInd
            DIC3DStitched.FacePairInds(end+1)=median(DIC3DStitched.FacePairInds([indF12 indF13 indF23]));
            % assign the mean color
            DIC3DStitched.FaceColors(end+1)=mean(DIC3DStitched.FaceColors([indF12 indF13 indF23]));
            
            DIC3DStitched.Faces=F;
            
        end
    end
    
    
    %% append surfaces that are not part of the stitching
    % surfaces not stitched
    pairInds=1:nPairs;
    pairsIndsNotStitched=find(~ismember(pairInds,pairIndList));
    for ipair=pairsIndsNotStitched
        DIC3DStitched.Faces=[DIC3DStitched.Faces; DIC3DAllPairsResults{ipair}.Faces+size(DIC3DStitched.Points3D{1},1)];
        DIC3DStitched.FaceColors=[DIC3DStitched.FaceColors; DIC3DAllPairsResults{ipair}.FaceColors];
        
        for it=1:nT
            DIC3DStitched.Points3D{it}=[DIC3DStitched.Points3D{it}; DIC3DAllPairsResults{ipair}.Points3D{it}];
            DIC3DStitched.corrComb{it}=[DIC3DStitched.corrComb{it}; DIC3DAllPairsResults{ipair}.corrComb{it}];
        end
        DIC3DStitched.FacePairInds=[DIC3DStitched.FacePairInds; (ipair)*ones(size(DIC3DAllPairsResults{ipair}.Faces,1),1)];
        DIC3DStitched.PointPairInds=[DIC3DStitched.PointPairInds; (ipair)*ones(size(DIC3DAllPairsResults{ipair}.Points3D{it},1),1)];
        DIC3DStitched.pairIndices(ipair,:)=DIC3DAllPairsResults{ipair}.cameraPairInd;
    end
    
    %% close stitching figures?
    closeButton = questdlg('Close stitching figures?', 'Close stitching figures?', 'Yes', 'No', 'Yes');
    switch closeButton
        case 'Yes'
            for imesh=1:nPairsForStitching-1
                close(hf(imesh));
            end
        case 'No'
    end
    
else % if no stitching required
    
    nPairs=numel(DIC3DAllPairsResults);
    pairInds=1:nPairs;
    nT=numel(DIC3DAllPairsResults{1}.Points3D);
    s1=DIC3DAllPairsResults{1};
    DIC3DStitched=struct;
    DIC3DStitched.Faces=s1.Faces;
    DIC3DStitched.Points3D=s1.Points3D;
    DIC3DStitched.corrComb=s1.corrComb;
    DIC3DStitched.FaceColors=s1.FaceColors;
    DIC3DStitched.FacePairInds=ones(size(s1.Faces,1),1);
    DIC3DStitched.pairIndices=zeros(nPairs,2);
    DIC3DStitched.pairIndices(1,:)=s1.cameraPairInd;
    DIC3DStitched.PointPairInds=ones(size(s1.Points3D{1},1),1);
    
    for ipair=pairInds(2:end)
        DIC3DStitched.Faces=[DIC3DStitched.Faces; DIC3DAllPairsResults{ipair}.Faces+size(DIC3DStitched.Points3D{1},1)];
        DIC3DStitched.FaceColors=[DIC3DStitched.FaceColors; DIC3DAllPairsResults{ipair}.FaceColors];
        
        for it=1:nT
            DIC3DStitched.Points3D{it}=[DIC3DStitched.Points3D{it}; DIC3DAllPairsResults{ipair}.Points3D{it}];
            DIC3DStitched.corrComb{it}=[DIC3DStitched.corrComb{it}; DIC3DAllPairsResults{ipair}.corrComb{it}];
        end
        DIC3DStitched.FacePairInds=[DIC3DStitched.FacePairInds; (ipair)*ones(size(DIC3DAllPairsResults{ipair}.Faces,1),1)];
        DIC3DStitched.PointPairInds=[DIC3DStitched.PointPairInds; (ipair)*ones(size(DIC3DAllPairsResults{ipair}.Points3D{it},1),1)];
        DIC3DStitched.pairIndices(ipair,:)=DIC3DAllPairsResults{ipair}.cameraPairInd;
    end
    
end

%% plot stitched surface
optStruct.lineColor='k';
anim8_DIC3DPP_faceMeasure(DIC3DStitched,'FacePairInds',0,optStruct);
anim8_DIC3DPP_pointMeasure(DIC3DStitched,'pairInd',0,optStruct);

end

function [deformationStruct] = triSurfaceDeformation(F,Vref,Vdef)
%% function for computing deformation on triangular elements
%
% INPUTS:
% * F: nFaces-by-3 array representing the list of vertices for triangular faces
% * Vref: nVertices-by-3 array representing the 3d positions of the vertices of F in the reference positions
% * Vdef: a cell array where each cell represents a deformed time frame 
% and contains a nVertices-by-3 array representing the 3d positions of the
% vertices of F in the reference positions
%
% OUTPUTS: 
% * deformationStruct: a structure with the following fields:
%   - Fmat % deformation gradient tensor
%   - Cmat % deformation gradient tensor
%   - Lamda1 % deformation gradient tensor
%   - Lamda2 % deformation gradient tensor
%   - J % dilatation (volume change. in this case area change)
%   - Emat % Lagrangian strain tensor
%   - emat % Euler-Almansi strain tensor
%   - Emgn % Lagrangian strain magnitude
%   - emgn % Almansi strain magnitude
%   - d3 % Normal to the faces.
%   - Epc1 % smallest planar Lagrangian principal strain
%   - Epc2 % largest planar Lagrangian principal strain
%   - Epc1vec % 1st planar Lagrangian principal strain direction (corresponds to Epc1) in the reference conf.
%   - Epc1vecCur % 2nd planar Lagrangian principal strain direction (corresponds to Epc1) in the reference conf.
%   - Epc2vec % 1st planar Lagrangian principal strain direction (corresponds to Epc1) transformed into the deformed (current) conf. so it is planar in the current state
%   - Epc2vecCur % 2nd planar Lagrangian principal strain direction (corresponds to Epc1) transformed into the deformed (current) conf. so it is planar in the current state
%   - epc1 % smallest planar Almansi principal strain
%   - epc2 % largest planar Almansi principal strain
%   - epc1vec % 1st planar Almansi principal strain direction (corresponds to Epc1). it is planar in the current conf.
%   - epc2vec % 2nd planar Almansi principal strain direction (corresponds to epc2). it is planar in the current conf.
%
% The deformation calculation is based on the Triangular Cosserat Point Theory (TCPE). 
% Solav, Dana, et al. "Bone pose estimation in the presence of soft tissue artifact using triangular cosserat point elements." Annals of biomedical engineering 44.4 (2016): 1181-1190.
% Solav, Dana, M. B. Rubin, and Alon Wolf. "Soft tissue artifact compensation using triangular cosserat point elements (TCPEs)." International Journal of Engineering Science 85 (2014): 1-9.

%%
logicVcell = iscell(Vdef);
if ~logicVcell
    error('the third input must be a cell array with one cell for each frame (time) containing a vertices 3D matrix).');
end

deformationStruct=struct;

% number of frames
nFrames=numel(Vdef);
% number of triangular faces
nFaces=size(F,1);

% preallocate cells
D1=cell(nFrames,1); 
D2=cell(nFrames,1); 
D3=cell(nFrames,1); 
d1=cell(nFrames,1); 
d2=cell(nFrames,1); 
d3=cell(nFrames,1);
Drec1=cell(nFrames,1); 
Drec2=cell(nFrames,1); 
Dnorm=cell(nFrames,1);
Fmat=cell(nFrames,1);
Cmat=cell(nFrames,1);
Lamda1=cell(nFrames,1);
Lamda2=cell(nFrames,1);
E=cell(nFrames,1); 
e=cell(nFrames,1); 
J=cell(nFrames,1); 
Emgn=cell(nFrames,1); 
emgn=cell(nFrames,1);
Epc1=cell(nFrames,1); 
Epc1vec=cell(nFrames,1); 
Epc1vecCur=cell(nFrames,1);
Epc2=cell(nFrames,1); 
Epc2vec=cell(nFrames,1); 
Epc2vecCur=cell(nFrames,1);
epc1=cell(nFrames,1); 
epc1vec=cell(nFrames,1); 
epc2=cell(nFrames,1); 
epc2vec=cell(nFrames,1);
EShearMax=cell(nFrames,1);
eShearMax=cell(nFrames,1);
EShearMaxVec1=cell(nFrames,1);
EShearMaxVec2=cell(nFrames,1);
EShearMaxVecCur1=cell(nFrames,1);
EShearMaxVecCur2=cell(nFrames,1);
eShearMaxVec1=cell(nFrames,1);
eShearMaxVec2=cell(nFrames,1);
Eeq=cell(nFrames,1);
eeq=cell(nFrames,1);
Area=cell(nFrames,1);

hw = waitbar(0,'Calculating deformations and strains');

for itime=1:nFrames
    waitbar(itime/(nFrames));
    
    % preallocation for speed
    D1{itime}=zeros(size(F,1),3); 
    D2{itime}=zeros(size(F,1),3); 
    D3{itime}=zeros(size(F,1),3);
    d1{itime}=zeros(size(F,1),3); 
    d2{itime}=zeros(size(F,1),3); 
    d3{itime}=zeros(size(F,1),3);
    Drec1{itime}=zeros(size(F,1),3); 
    Drec2{itime}=zeros(size(F,1),3); 
    Dnorm{itime}=zeros(size(F,1),1);
    Fmat{itime}=zeros(3,3,size(F,1));
    Cmat{itime}=zeros(3,3,size(F,1));
    Lamda1{itime}=zeros(size(F,1),1);
    Lamda2{itime}=zeros(size(F,1),1);
    E{itime}=zeros(3,3,size(F,1)); 
    e{itime}=zeros(3,3,size(F,1));
    J{itime}=zeros(size(F,1),1); 
    Emgn{itime}=zeros(size(F,1),1); 
    emgn{itime}=zeros(size(F,1),1);
    Epc1{itime}=zeros(size(F,1),1); 
    Epc1vec{itime}=zeros(size(F,1),3); 
    Epc1vecCur{itime}=zeros(size(F,1),3);
    Epc2{itime}=zeros(size(F,1),1); 
    Epc2vec{itime}=zeros(size(F,1),3); 
    Epc2vecCur{itime}=zeros(size(F,1),3);
    epc1{itime}=zeros(size(F,1),1); 
    epc1vec{itime}=zeros(size(F,1),3);
    epc2{itime}=zeros(size(F,1),1); 
    epc2vec{itime}=zeros(size(F,1),3);
    EShearMax{itime}=zeros(size(F,1),1);
    eShearMax{itime}=zeros(size(F,1),1);
    EShearMaxVec1{itime}=zeros(size(F,1),3);
    EShearMaxVec2{itime}=zeros(size(F,1),3);
    EShearMaxVecCur1{itime}=zeros(size(F,1),3);
    EShearMaxVecCur2{itime}=zeros(size(F,1),3);
    eShearMaxVec1{itime}=zeros(size(F,1),3);
    eShearMaxVec2{itime}=zeros(size(F,1),3);
    Eeq{itime}=zeros(size(F,1),1);
    eeq{itime}=zeros(size(F,1),1);
    Area{itime}=zeros(size(F,1),1);
    
    for itri=1:nFaces
        
        % reference director vectors
        D1{itime}(itri,:)=Vref(F(itri,2),:)-Vref(F(itri,1),:);
        D2{itime}(itri,:)=Vref(F(itri,3),:)-Vref(F(itri,1),:);
        D3{itime}(itri,:)=cross(D1{itime}(itri,:),D2{itime}(itri,:))/norm(cross(D1{itime}(itri,:),D2{itime}(itri,:)));
        
        % current director vectors
        d1{itime}(itri,:)=Vdef{itime}(F(itri,2),:)-Vdef{itime}(F(itri,1),:);
        d2{itime}(itri,:)=Vdef{itime}(F(itri,3),:)-Vdef{itime}(F(itri,1),:);
        d3{itime}(itri,:)=cross(d1{itime}(itri,:),d2{itime}(itri,:))/norm(cross(d1{itime}(itri,:),d2{itime}(itri,:)));
        
        % reciprocal vectors
        Dnorm{itime}(itri)=cross(D1{itime}(itri,:),D2{itime}(itri,:))*D3{itime}(itri,:)';
        Drec1{itime}(itri,:) = cross(D2{itime}(itri,:),D3{itime}(itri,:))/Dnorm{itime}(itri);
        Drec2{itime}(itri,:) = cross(D3{itime}(itri,:),D1{itime}(itri,:))/Dnorm{itime}(itri);
        
        % area
        Area{itime}(itri)=0.5*Dnorm{itime}(itri);
        
        % deformation gradient tensor
        for ii=1:3
            for jj=1:3
                Fmat{itime}(ii,jj,itri)=d1{itime}(itri,ii)*Drec1{itime}(itri,jj)+d2{itime}(itri,ii)*Drec2{itime}(itri,jj)+d3{itime}(itri,ii)*D3{itime}(itri,jj);
            end
        end
        
        if sum(sum(isnan(Fmat{itime}(:,:,itri))))==0 % if F doesn't have NaNs in it
            
            % Cauchy-Green deformation tensor
            Cmat{itime}(:,:,itri)=Fmat{itime}(:,:,itri)'*Fmat{itime}(:,:,itri);
            [~,eigValC]=eig(Cmat{itime}(:,:,itri));
            % principal stretches
            Lamdas=sqrt([eigValC(1,1) eigValC(2,2) eigValC(3,3)]);
            [~,ind1]=min(abs(Lamdas-1)); % find the index of the stretch that equals 1
            Lamdas(ind1)=[]; % delete the stretch that equals 1
            Lamda1{itime}(itri)=Lamdas(1);
            Lamda2{itime}(itri)=Lamdas(2);
        
            % dilitation
            J{itime}(itri)=det(Fmat{itime}(:,:,itri));
           
            % Lagrangian finite strain tensor
            E{itime}(:,:,itri)=0.5*(Fmat{itime}(:,:,itri)'*Fmat{itime}(:,:,itri)-eye(3));
            
            % Eulerian-Almansi finite strain tensor
            e{itime}(:,:,itri)=0.5*(eye(3)-inv(Fmat{itime}(:,:,itri)*Fmat{itime}(:,:,itri)'));
            
            % Lagrangian Strain magnitude
            Emgn{itime}(itri) = norm(E{itime}(:,:,itri),'fro');
            
            % Eulerian Strain magnitude
            emgn{itime}(itri) = norm(e{itime}(:,:,itri),'fro');
            
            % Lagrangian Strain eigenvalues and eigenvactors (principal strains and their directions)
            [eigVecE,eigValE]=eig(E{itime}(:,:,itri));
            
            % find the 2 eigenvectors of E that are on the plane of the triangle (remove D3) because the strain in the D3 direction is zero anyway
            [~,D3Ind]=max(abs(eigVecE'*D3{itime}(itri,:)')); % by taking the dot of every eigen vector with D3 we find which eigenvector is D3 (the only one that's not zero)
            eigVecPlanInd=[1 2 3];
            eigVecPlanInd(eigVecPlanInd==D3Ind)=[]; % the indeces of the eigenvectors on the plane of the triangle
            [Epc1{itime}(itri),Epc1Ind]=min([eigValE(eigVecPlanInd(1),eigVecPlanInd(1)) eigValE(eigVecPlanInd(2),eigVecPlanInd(2))]); % smallest principal strain that is no zero
            [Epc2{itime}(itri),Epc2Ind]=max([eigValE(eigVecPlanInd(1),eigVecPlanInd(1)) eigValE(eigVecPlanInd(2),eigVecPlanInd(2))]); % largest principal strain that is no zero
            Epc1vec{itime}(itri,:)=eigVecE(:,eigVecPlanInd(Epc1Ind)); % direction of EPrinc1 in reference configuration
            Epc1vecCur{itime}(itri,:)=Fmat{itime}(:,:,itri)*Epc1vec{itime}(itri,:)'; % direction of EPrinc1 transformed into the deformed configuration
            Epc1vecCur{itime}(itri,:)=Epc1vecCur{itime}(itri,:)/norm(Epc1vecCur{itime}(itri,:)); % direction of EPrinc1 in deformed configuration- normalized
            Epc2vec{itime}(itri,:)=eigVecE(:,eigVecPlanInd(Epc2Ind)); % same for EPrinc2
            Epc2vecCur{itime}(itri,:)=Fmat{itime}(:,:,itri)*Epc2vec{itime}(itri,:)';
            Epc2vecCur{itime}(itri,:)=Epc2vecCur{itime}(itri,:)/norm(Epc2vecCur{itime}(itri,:));
            
            % Eulerian Strain eigenvalues and eigenvactors
            [eigVece,eigVale]=eig(e{itime}(:,:,itri));
            % find the 2 eigenvectors of E that are on the plane of the triangle (remove D3) because the strain in the D3 direction is zero anyway
            [~,d3Ind]=max(abs(eigVece'*d3{itime}(itri,:)')); % by taking the dot of every eigen vector with D3 we find which eigenvector is D3 (the only one that's not zero)
            eigVecPlanInd=[1 2 3];
            eigVecPlanInd(eigVecPlanInd==d3Ind)=[]; % the indeces of the eigenvectors on the plane of the triangle
            [epc1{itime}(itri),epc1Ind]=min([eigVale(eigVecPlanInd(1),eigVecPlanInd(1)) eigVale(eigVecPlanInd(2),eigVecPlanInd(2))]); % smallest principal strain that is no zero
            [epc2{itime}(itri),epc2Ind]=max([eigVale(eigVecPlanInd(1),eigVecPlanInd(1)) eigVale(eigVecPlanInd(2),eigVecPlanInd(2))]); % largest principal strain that is no zero
            epc1vec{itime}(itri,:)=eigVece(:,eigVecPlanInd(epc1Ind)); % direction of ePrinc1 in deformed configuration
            epc2vec{itime}(itri,:)=eigVece(:,eigVecPlanInd(epc2Ind)); % same for EPrinc2
                 
            % Max shear strain and its direction
            % Lagrangian
            EShearMax{itime}(itri)=.5*(Epc2{itime}(itri)-Epc1{itime}(itri));
            EShearMaxVec1{itime}(itri,:)=(1/sqrt(2))*(Epc1vec{itime}(itri,:)+Epc2vec{itime}(itri,:));
            EShearMaxVec2{itime}(itri,:)=(1/sqrt(2))*(Epc2vec{itime}(itri,:)-Epc1vec{itime}(itri,:));
            EShearMaxVecCur1{itime}(itri,:)=(1/sqrt(2))*(Epc1vecCur{itime}(itri,:)+Epc2vecCur{itime}(itri,:));
            EShearMaxVecCur2{itime}(itri,:)=(1/sqrt(2))*(Epc2vecCur{itime}(itri,:)-Epc1vecCur{itime}(itri,:));
            % Eulerian
            eShearMax{itime}(itri)=.5*(epc2{itime}(itri)-epc1{itime}(itri));
            eShearMaxVec1{itime}(itri,:)=(1/sqrt(2))*(epc1vec{itime}(itri,:)+epc2vec{itime}(itri,:));
            eShearMaxVec2{itime}(itri,:)=(1/sqrt(2))*(epc2vec{itime}(itri,:)-epc1vec{itime}(itri,:));
            
            % Equivalent strain (von-mises)
            Edev=E{itime}(:,:,itri)-(1/3)*trace(E{itime}(:,:,itri))*eye(3);
            Eeq{itime}(itri)=sqrt((2/3)*sum(sum(Edev.*Edev)));
            edev=e{itime}(:,:,itri)-(1/3)*trace(e{itime}(:,:,itri))*eye(3);
            eeq{itime}(itri)=sqrt((2/3)*sum(sum(edev.*edev)));
            
        else % if F has NaNs, all the resulting measures are also NaN
            Cmat{itime}(:,:,itri)=[NaN NaN NaN; NaN NaN NaN; NaN NaN NaN];
            Lamda1{itime}(itri)=NaN;
            Lamda2{itime}(itri)=NaN;
            E{itime}(:,:,itri)=[NaN NaN NaN; NaN NaN NaN; NaN NaN NaN];
            e{itime}(:,:,itri)=[NaN NaN NaN; NaN NaN NaN; NaN NaN NaN];
            Emgn{itime}(itri) = NaN; 
            emgn{itime}(itri) = NaN; 
            J{itime}(itri)=NaN;
            Epc1{itime}(itri)=NaN; 
            Epc1vec{itime}(itri,:)=[NaN NaN NaN]; 
            Epc1vecCur{itime}(itri,:)=[NaN NaN NaN];
            Epc2{itime}(itri)=NaN; 
            Epc2vec{itime}(itri,:)=[NaN NaN NaN]; 
            Epc2vecCur{itime}(itri,:)=[NaN NaN NaN];
            epc1{itime}(itri)=NaN; 
            epc1vec{itime}(itri,:)=[NaN NaN NaN];
            epc2{itime}(itri)=NaN; 
            epc2vec{itime}(itri,:)=[NaN NaN NaN];
            EShearMax{itime}(itri)=NaN;
            eShearMax{itime}(itri)=NaN;
            Eeq{itime}(itri)=NaN;
            eeq{itime}(itri)=NaN;
            Area{itime}(itri)=NaN;
        end
        
    end
    
end

delete(hw);

% save all results in the structure

deformationStruct.Fmat=Fmat; % deformation gradient tensor
deformationStruct.Cmat=Cmat; % deformation gradient tensor
deformationStruct.Lamda1=Lamda1; % deformation gradient tensor
deformationStruct.Lamda2=Lamda2; % deformation gradient tensor
deformationStruct.J=J; % dilatation (volume change. in this case area change)
deformationStruct.Emat=E; % Lagrangian strain tensor
deformationStruct.emat=e; % Euler-Almansi strain tensor
deformationStruct.Emgn=Emgn; % Lagrangian strain magnitude
deformationStruct.emgn=emgn; % Almansi strain magnitude
deformationStruct.d3=d3; % Normal to the faces.
deformationStruct.Epc1=Epc1; % smallest planar Lagrangian principal strain
deformationStruct.Epc2=Epc2; % largest planar Lagrangian principal strain
deformationStruct.Epc1vec=Epc1vec; % 1st planar Lagrangian principal strain direction (corresponds to Epc1) in the reference conf.
deformationStruct.Epc1vecCur=Epc1vecCur; % 2nd planar Lagrangian principal strain direction (corresponds to Epc1) in the reference conf.
deformationStruct.Epc2vec=Epc2vec; % 1st planar Lagrangian principal strain direction (corresponds to Epc1) transformed into the deformed (current) conf. so it is planar in the current state
deformationStruct.Epc2vecCur=Epc2vecCur; % 2nd planar Lagrangian principal strain direction (corresponds to Epc1) transformed into the deformed (current) conf. so it is planar in the current state
deformationStruct.epc1=epc1; % smallest planar Almansi principal strain
deformationStruct.epc2=epc2; % largest planar Almansi principal strain
deformationStruct.epc1vec=epc1vec; % 1st planar Almansi principal strain direction (corresponds to Epc1). it is planar in the current conf.
deformationStruct.epc2vec=epc2vec; % 2nd planar Almansi principal strain direction (corresponds to epc2). it is planar in the current conf.
deformationStruct.EShearMax=EShearMax; % Max shear strain (Lagrangian)
deformationStruct.eShearMax=eShearMax; % Max shear strain (Eulerian)
deformationStruct.EShearMaxVec1=EShearMaxVec1;
deformationStruct.EShearMaxVec2=EShearMaxVec2;
deformationStruct.EShearMaxVecCur1=EShearMaxVecCur1;
deformationStruct.EShearMaxVecCur2=EShearMaxVecCur2;
deformationStruct.eShearMaxVec2=eShearMaxVec2;
deformationStruct.Eeq=Eeq; % Equivalent strain (Lagrangian)
deformationStruct.eeq=eeq; % Equivalent strain (Eulerian)
deformationStruct.Area=Area;

end



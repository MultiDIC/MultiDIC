%% STEP 2 Run Ncorr analysis on sets of images from a pair of 2 cameras
% 3D reconstruction of points correlated with Ncorr using DLT 
% for running this step, you need to have DIC2DpairResults structures of
% the pairs you want to reconstruct from STEP2, as well as the DLT
% parameters of the cameras from STEP1.


%%
clearvars; close all

%% CHOOSE PATHS OPTIONS
% select DIC2DpairResults structures
PathInitial='C:\Users\Dana\Dropbox (Personal)\MIT\research\DIC\360_DIC\LimbTest\data';
structPaths = uipickfiles('FilterSpec',PathInitial,'Prompt','Select one or multiple 2D-DIC results structures');
nPairs=numel(structPaths);
indPairs=zeros(nPairs,2);
for ip=1:nPairs
    DIC2D{ip}=load(structPaths{ip});
    DIC2D{ip}=DIC2D{ip}.DIC2DpairResults;
    indPairs(ip,:)=[DIC2D{ip}.nCamRef DIC2D{ip}.nCamDef];
end
indCams=unique(indPairs);
nCams=numel(indCams);
% select the folder where the DLT parameter files are stored
folderPathInitial='C:\Users\Dana\Dropbox (Personal)\MIT\research\DIC\360_DIC\LimbTest\data\';
DLTparameterFolder=uigetdir(folderPathInitial,['Select the folder containing DLT parameters for cameras ' num2str(indCams')]);
for ic=1:nCams
    DLTparStruct=load([DLTparameterFolder '\DLTpar_cam' num2str(indCams(ic),'%02u')]);
    DLTpar{ic}=DLTparStruct.L;
end

% save 3D-DIC results? choose save path and overwrite options
[save3DDIClogic,savePath]=Qsave3DDICresults(structPaths);


%% 3D reconstruction using Direct Linear Transformation

DIC3DAllPairsResults=cell(nPairs,1);
for ip=1:nPairs
    
    DIC3DpairResults=struct;
    
    nCamRef=indPairs(ip,1);
    nCamDef=indPairs(ip,2);
    
    DIC3DpairResults.cameraPairInd=[nCamRef nCamDef];
    
    L1=DLTpar{find(indCams==nCamRef)}; % DLT parameters for 1st camera
    L2=DLTpar{find(indCams==nCamDef)}; % DLT parameters for 2nd camera
        
    nImages=DIC2D{ip}.nImages;
    Points=DIC2D{ip}.Points;
    CorCoeff=DIC2D{ip}.CorCoeffVec;
    F=DIC2D{ip}.Faces;
    FC=DIC2D{ip}.FaceColors;

    DIC3DpairResults.Faces=F;
    DIC3DpairResults.FaceColors=FC;
    
    % allocate result variables
    DIC3DpairResults.Points3D=cell(nImages,1);
    DIC3DpairResults.Disp.DispVec=cell(nImages,1);
    DIC3DpairResults.Disp.DispMgn=cell(nImages,1);
    DIC3DpairResults.FaceCentroids=cell(nImages,1);
    DIC3DpairResults.corrComb=cell(nImages,1);
        
    for ii=1:nImages
        % correlated points
        P1=Points{ii};
        P2=Points{ii+nImages};
        
%         % load the cropped rectangle from the original image  
%         rect1=load([imagePath '\rect_' num2str(nCamRef) '_' num2str(ii,'%03u') '.mat']); rect1=round(rect1.rect);
%         rect2=load([imagePath '\rect_' num2str(nCamDef) '_' num2str(ii,'%03u') '.mat']); rect2=round(rect2.rect);
%         
        % modify the points to the full image
        P1m=P1; 
%         P1m=P1m+[rect1(1)-1 rect1(2)-1];
        P2m=P2; 
%         P2m=P2m+[rect2(1)-1 rect2(2)-1];

        % Solve the DLT system
        P3D=DLT11Reconstruction(P1m,P2m,L1,L2);
        DIC3DpairResults.Points3D{ii}=P3D;
        
        % Combined (worst) correlation coefficients
        DIC3DpairResults.corrComb{ii}=max([CorCoeff{ii} CorCoeff{ii+nImages}],[],2);
        
        % compute face centroids
        for iface=1:size(F,1)
            DIC3DpairResults.FaceCentroids{ii}(iface,:)=mean(P3D(F(iface,:),:));
        end

        % Compute displacements between frames (per point)
        DispVec=DIC3DpairResults.Points3D{ii}-DIC3DpairResults.Points3D{1};
        DIC3DpairResults.Disp.DispVec{ii}=DispVec;
        DIC3DpairResults.Disp.DispMgn{ii}=sqrt(DispVec(:,1).^2+DispVec(:,2).^2+DispVec(:,3).^2);

    end
    
    % compute deformation and strains (per triangular face)
    deformationStruct=triSurfaceDeformation(F,DIC3DpairResults.Points3D{1},DIC3DpairResults.Points3D);
    DIC3DpairResults.Deform=deformationStruct;
    
    % compute triangle regularity (isotropy index) based on reference image
    [FisoInd]=faceIsotropyIndex(F,DIC3DpairResults.Points3D{1});
    DIC3DpairResults.FaceIsoInd=FisoInd;
    
    DIC3DAllPairsResults{ip}=DIC3DpairResults;
            
     if save3DDIClogic
        save([savePath '\DIC3DpairResults_C' num2str(nCamRef,'%02u') '_C' num2str(nCamDef,'%02u')] ,'DIC3DpairResults');
     end
    
end

 if save3DDIClogic
        save([savePath '\DIC3DAllPairsResults'] ,'DIC3DAllPairsResults');
 end

 %% Plot pair results?
 
plotButton = questdlg('Plot 3D-DIC results?', 'Plot?', 'Yes', 'No', 'Yes');
switch plotButton
    case 'Yes'
        plotMultiDICPairResults(DIC3DAllPairsResults);
    case 'No'
end
  
 %%
return

%% ADD LATER (AFTER MERGING)
%% Compute rigid body transformation between point clouds
REC_pair_transformed=cell(nPairs,nSets);
RotMat=cell(nPairs,nSets);
TransVec=cell(nPairs,nSets);
for iset=indSets
    for ipair=1:nPairs
        TruP=REC_pair{ipair,1};
        RecP=REC_pair{ipair,iset};
        [RotMat{ipair,iset},TransVec{ipair,iset},REC_pair_transformed{ipair,iset}]=rigidTransformation(RecP,TruP);        
    end
end

%% Compute displacements between sets - after transformation
REC_TransDispVec=cell(nPairs,nSets); REC_TransDispTot=cell(nPairs,nSets);
for iset=indSets
    for ipair=1:nPairs
        REC_TransDispVec{ipair,iset}=REC_pair_transformed{ipair,iset}-REC_pair_transformed{ipair,1};
        REC_TransDispTot{ipair,iset}=sqrt(REC_TransDispVec{ipair,iset}(:,1).^2+REC_TransDispVec{ipair,iset}(:,2).^2+REC_TransDispVec{ipair,iset}(:,3).^2);
    end
end

if saveOn==1
    save([savePath '\3Dreconstruction_allVar.mat']);
end




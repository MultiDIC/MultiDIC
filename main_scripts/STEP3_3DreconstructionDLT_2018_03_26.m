%% STEP 3: 3D reconstruction
% 3D reconstruction of points correlated with Ncorr using DLT 
% for running this step, you need to have DIC2DpairResults structures of
% the pairs you want to reconstruct from STEP2, as well as the DLT
% parameters of the cameras from STEP1.

%%
clearvars; close all

fs=get(0, 'DefaultUIControlFontSize');
set(0, 'DefaultUIControlFontSize', 10);

%% CHOOSE PATHS OPTIONS

% select DIC2DpairResults structures
PathInitial=pwd;
structPaths = uipickfiles('FilterSpec',PathInitial,'Prompt','Select one or multiple 2D-DIC results structures');
nPairs=numel(structPaths);
indPairs=zeros(nPairs,2);
DIC2D=cell(nPairs,1);
for ip=1:nPairs
    DIC2D{ip}=load(structPaths{ip});
    DIC2D{ip}=DIC2D{ip}.DIC2DpairResults;
    indPairs(ip,:)=[DIC2D{ip}.nCamRef DIC2D{ip}.nCamDef];
end
indCams=unique(indPairs(:));
nCams=numel(indCams);

% select the folder where the DLT parameter files are stored
folderPathInitial=pwd;
DLTparameterFolder = uipickfiles('FilterSpec',folderPathInitial,'Prompt',['Select the folder containing DLT parameters for cameras ' num2str(indCams')]);
DLTstructAllCams=cell(nCams,1);
for ic=1:nCams
    DLTstructTemp=load([DLTparameterFolder{1} '\DLTstruct_cam_' num2str(indCams(ic))]);
    DLTstructAllCams{ic}=DLTstructTemp.DLTstructCam;
end

% remove distortion? If yes, choose the folder where the parameters exist
distortionRemovalButton = questdlg('Remove distortion from calibration points?', 'Remove distortion?', 'Yes', 'No', 'Yes');
switch distortionRemovalButton
    case 'Yes'
        distortionRemovalLogic=true(1);
        distortionParametersPath=uigetdir([],'Select a folder where the cameraCBparameters are located');
    case 'No'
        distortionRemovalLogic=false(1);
end

% save 3D-DIC results? choose save path and overwrite options
[save3DDIClogic,savePath]=Qsave3DDICresults(structPaths);

%% if distortion removal is required, calculate the new DLT parameters
switch distortionRemovalLogic
    case 0
        
    case 1
        distortionPar=cell(nCams,1);
        for ic=1:nCams
            icam=indCams(ic);
            P2Dtemp=DLTstructAllCams{ic}.imageCentroids;
            columns=DLTstructAllCams{ic}.columns;
            C3Dtrue=DLTstructAllCams{ic}.C3Dtrue;
            % undistort image points
            distortionParTemp=load([distortionParametersPath '\cameraCBparameters_cam_' num2str(icam)]);
            distortionPar{ic}=distortionParTemp.cameraCBparameters.cameraParameters;
            [P2Dtemp] = undistortPoints(P2Dtemp,distortionPar{ic});
            % calculate new DLT parameters
            C3Dtemp=reshape(C3Dtrue(:,columns,:),size(C3Dtrue,1)*numel(columns),3);
            L=DLT11Calibration(P2Dtemp,C3Dtemp);
            % assign to DLTstructPairs
            DLTstructAllCams{ic}.DLTparams=L;
        end
    
end

%% 3D reconstruction using Direct Linear Transformation
DIC3DAllPairsResults=cell(nPairs,1);

for ip=1:nPairs % loop over stereo pairs
    hw = waitbar(0,['Please wait while reconstructing points and calculating deformations for pair ' num2str(ip) '...']);
    
    % create 3D-DIc results struct
    DIC3DpairResults=struct;
    
    % camera indices of current pair
    nCamRef=indPairs(ip,1);
    nCamDef=indPairs(ip,2);
    
    iCamRef=find(indCams==nCamRef);
    iCamDef=find(indCams==nCamDef);
    
    DIC3DpairResults.cameraPairInd=[nCamRef nCamDef];
    
    % DLT parameters
    L1=DLTstructAllCams{iCamRef}.DLTparams; % DLT parameters for 1st camera
    L2=DLTstructAllCams{iCamDef}.DLTparams; % DLT parameters for 2nd camera
    
    % extract information from 2d-dic results
    nImages=DIC2D{ip}.nImages;
    CorCoeff=DIC2D{ip}.CorCoeffVec;
    F=DIC2D{ip}.Faces;
    FC=DIC2D{ip}.FaceColors;

    DIC3DpairResults.Faces=F;
    DIC3DpairResults.FaceColors=FC;
    
    switch distortionRemovalLogic
        case 0 % if no distortion correction, use original points and set distortion model to 'none'
            Points=DIC2D{ip}.Points;
            DIC3DpairResults.distortionModel{1}='none';
            DIC3DpairResults.distortionModel{2}='none';
        case 1 % if yes distortion correction, correct points and set distortion model
            Points=DIC2D{ip}.Points;
            DIC3DpairResults.distortionModel{1}=distortionPar{iCamRef};
            DIC3DpairResults.distortionModel{2}=distortionPar{iCamDef};
            step=min([40,size(Points{1},1)]);
            % remove distortion from points
            % first camera
            for ii=1:nImages 
                waitbar(ii/(6*nImages));
                
                % remove nan points
                PointsNaNlogic=any(isnan(Points{ii}),2);
                PointNoNan=Points{ii}(~PointsNaNlogic,:);
                distortionParNow=distortionPar{iCamRef};
                PointsTemp=zeros(size(PointNoNan));
                times=1:step:size(PointNoNan,1);
                % correct points in small sets because large sets are slow
                for jj=times(1:end-1)
                    PointsTemp(jj:jj+step,:) = undistortPoints(PointNoNan(jj:jj+step,:),distortionParNow);
                end
                PointsTemp(times(end):end,:) = undistortPoints(PointNoNan(times(end):end,:),distortionParNow);
                Points{ii}(~PointsNaNlogic,:) = PointsTemp;                
            end
            % repeat for second camera
            for ii=nImages+1:2*nImages 
                waitbar(ii/(6*nImages));
                
                PointsNaNlogic=any(isnan(Points{ii}),2);
                PointNoNan=Points{ii}(~PointsNaNlogic,:);
                distortionParNow=distortionPar{iCamDef};
                PointsTemp=zeros(size(PointNoNan));
                times=1:step:size(PointNoNan,1);
                for jj=times(1:end-1)
                    PointsTemp(jj:jj+step,:) = undistortPoints(PointNoNan(jj:jj+step,:),distortionParNow);
                end
                PointsTemp(times(end):end,:) = undistortPoints(PointNoNan(times(end):end,:),distortionParNow);
                Points{ii}(~PointsNaNlogic,:) = PointsTemp;
            end 
            
    end

    % pre-allocate 3D-DIC result variables
    DIC3DpairResults.Points3D=cell(nImages,1);
    DIC3DpairResults.Points3Dtransformed=cell(nImages,1);
    DIC3DpairResults.Disp.DispVec=cell(nImages,1);
    DIC3DpairResults.Disp.DispMgn=cell(nImages,1);
    DIC3DpairResults.Disp.DispVecTransformed=cell(nImages,1);
    DIC3DpairResults.Disp.DispMgnTransformed=cell(nImages,1);
    DIC3DpairResults.FaceCentroids=cell(nImages,1);
    DIC3DpairResults.corrComb=cell(nImages,1);
    DIC3DpairResults.FaceCorrComb=cell(nImages,1);
    DIC3DpairResults.FaceIsoInd=cell(nImages,1);
    DIC3DpairResults.RBM.RotMat=cell(nImages,1);
    DIC3DpairResults.RBM.TransVec=cell(nImages,1);
    
    for ii=1:nImages % loop over images (time frames)
        waitbar(1/3+ii/(3*nImages));
        
        % correlated points from 2 cameras
        P1=Points{ii};
        P2=Points{ii+nImages};
        
        % Solve the DLT system
        P3D=DLT11Reconstruction(P1,P2,L1,L2);
        DIC3DpairResults.Points3D{ii}=P3D;
        
        % Combined (worst) correlation coefficients
        DIC3DpairResults.corrComb{ii}=max([CorCoeff{ii} CorCoeff{ii+nImages}],[],2);
        % Face correlation coefficient (worst)
        DIC3DpairResults.FaceCorrComb{ii}=max(DIC3DpairResults.corrComb{ii}(F),[],2);
        
        % compute face centroids
        for iface=1:size(F,1)
            DIC3DpairResults.FaceCentroids{ii}(iface,:)=mean(P3D(F(iface,:),:));
        end

        % Compute displacements between frames (per point)
        DispVec=DIC3DpairResults.Points3D{ii}-DIC3DpairResults.Points3D{1};
        DIC3DpairResults.Disp.DispVec{ii}=DispVec;
        DIC3DpairResults.Disp.DispMgn{ii}=sqrt(DispVec(:,1).^2+DispVec(:,2).^2+DispVec(:,3).^2);
        
        % Compute rigid body transformation between point clouds
        [RotMat,TransVec,Points3Dtransformed]=rigidTransformation(DIC3DpairResults.Points3D{ii},DIC3DpairResults.Points3D{1});
        DIC3DpairResults.RBM.RotMat{ii}=RotMat;
        DIC3DpairResults.RBM.TransVec{ii}=TransVec;
        DIC3DpairResults.Points3Dtransformed{ii}=Points3Dtransformed;
        
        % Compute displacements between sets - after transformation
        DispVec=DIC3DpairResults.Points3Dtransformed{ii}-DIC3DpairResults.Points3Dtransformed{1};
        DIC3DpairResults.Disp.DispVecTransformed{ii}=DispVec;
        DIC3DpairResults.Disp.DispMgnTransformed{ii}=sqrt(DispVec(:,1).^2+DispVec(:,2).^2+DispVec(:,3).^2);
        
        % compute face centroids - after transformation
        for iface=1:size(F,1)
            DIC3DpairResults.FaceCentroidsTransformed{ii}(iface,:)=mean(Points3Dtransformed(F(iface,:),:));
        end
    end
    
    % compute deformation and strains (per triangular face)
    deformationStruct=triSurfaceDeformation(F,DIC3DpairResults.Points3D{1},DIC3DpairResults.Points3D);
    DIC3DpairResults.Deform=deformationStruct;
    
    % compute triangle regularity (isotropy index) based on reference image
    for ii=1:nImages
        waitbar(2/3+ii/(3*nImages));
        [FisoInd]=faceIsotropyIndex(F,DIC3DpairResults.Points3D{ii});
        DIC3DpairResults.FaceIsoInd{ii}=FisoInd;
    end
    
    DIC3DAllPairsResults{ip}=DIC3DpairResults;
            
     if save3DDIClogic
        save([savePath '\DIC3DpairResults_C_' num2str(nCamRef) '_C_' num2str(nCamDef)] ,'DIC3DpairResults');
     end
    
     delete(hw) 
end

%% save all pairs

if save3DDIClogic
    saveName=[savePath '\DIC3DAll' num2str(nPairs) 'PairsResults.mat'];
    icount=1;
    while exist(saveName,'file')
        saveName=[savePath '\DIC3DAll' num2str(nPairs) 'PairsResults(' num2str(icount) ').mat'];
        icount=icount+1;
    end
    save(saveName,'DIC3DAllPairsResults');
end


 %% Plot pair results?
 
plotButton = questdlg('Plot 3D-DIC results from camera pairs?', 'Plot?', 'Yes', 'No', 'Yes');
switch plotButton
    case 'Yes'
        optStruct=struct;
        optStruct.zDirection=-1;
        optStruct.FaceAlpha=1;
        plotMulti3DPairResults(DIC3DAllPairsResults,optStruct); 
%         plotMulti3DPairResultsRBM(DIC3DAllPairsResults,optStruct);

    case 'No'
end

%% finish
h=msgbox('STEP3 is completed');
h.CurrentAxes.Children.FontSize=11;

set(0, 'DefaultUIControlFontSize', fs);

return

%% ADD LATER (AFTER MERGING)
% %% Compute rigid body transformation between point clouds
% REC_pair_transformed=cell(nPairs,nSets);
% RotMat=cell(nPairs,nSets);
% TransVec=cell(nPairs,nSets);
% for iset=indSets
%     for ipair=1:nPairs
%         TruP=REC_pair{ipair,1};
%         RecP=REC_pair{ipair,iset};
%         [RotMat{ipair,iset},TransVec{ipair,iset},REC_pair_transformed{ipair,iset}]=rigidTransformation(RecP,TruP);        
%     end
% end
% 
% %% Compute displacements between sets - after transformation
% REC_TransDispVec=cell(nPairs,nSets); REC_TransDispTot=cell(nPairs,nSets);
% for iset=indSets
%     for ipair=1:nPairs
%         REC_TransDispVec{ipair,iset}=REC_pair_transformed{ipair,iset}-REC_pair_transformed{ipair,1};
%         REC_TransDispTot{ipair,iset}=sqrt(REC_TransDispVec{ipair,iset}(:,1).^2+REC_TransDispVec{ipair,iset}(:,2).^2+REC_TransDispVec{ipair,iset}(:,3).^2);
%     end
% end
% 
% if saveOn==1
%     save([savePath '\3Dreconstruction_allVar.mat']);
% end


%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>

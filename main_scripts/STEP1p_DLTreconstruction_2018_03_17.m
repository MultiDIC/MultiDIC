%% STEP 1p : calculate and plot reconstruction errors for camera stereo pairs
% This main scripts reconstruct 3D points based on pairs of cameras and
% their DLT structures, and calculates the errors betweent the
% reconstructed points and the true points.
% reconstruction errors are plotted and saved

%%
clearvars; close all

fs=get(0, 'DefaultUIControlFontSize');
set(0, 'DefaultUIControlFontSize', 10);

%% CHOOSE PATHS and OPTIONS

% initial path
PathInitial=pwd; % change as needed
% select at least one pair (at least two cameras) DLT structures forreconstruction
filePaths = uipickfiles('FilterSpec',PathInitial,'Prompt','Select DLT structures from at least one stereo camera pair');
if ~iscell(filePaths)
    return
end

% save camera parameters? choose save path and overwrite options
[saveDLTpairsLogic,savePath]=QsaveDLTpairStruct(filePaths);
Ncam=numel(filePaths);
indCams=zeros(1,Ncam);
for ic=1:Ncam % camera indeces from image paths
    DLTstructTemp=load(filePaths{ic});
    DLTcamStructs{ic}=DLTstructTemp.DLTstructCam;
    indCams(ic)=DLTcamStructs{ic}.indCam;
    C3Dtrue{ic}=DLTcamStructs{ic}.C3Dtrue;
end

% check that all the structures are based on the same calibration target
for ic=2:Ncam
    if ~isequal(C3Dtrue{1},C3Dtrue{ic})
        error('All DLT structures must use the same calibration target 3D true points');
    end
end
NcTotal=size(C3Dtrue{1},2);
Nr=size(C3Dtrue{1},1);

% figures path
figuresPath=[savePath '\calibrationFigures'];
warning('off','MATLAB:MKDIR:DirectoryExists');
mkdir(figuresPath);

% select camera pairs indices
for ic=1:Ncam-1 % initial guess for camera pairs are consecutive indices
    indPairsInitial(ic,:)=[indCams(ic) indCams(ic+1)];
end
answer = inputdlg({['Enter the indices of pairs of cameras [' num2str(indCams) '] as semicolon seperated pairs:']},'Input',[1,50],{mat2str(indPairsInitial)});
indPairs=str2num(answer{1});

% remove distortion? If yes, choose the folder where the parameters are stored
distortionRemovalButton = questdlg('Remove distortion from calibration points?', 'Remove distortion?', 'Yes', 'No', 'Yes');
switch distortionRemovalButton
    case 'Yes'
        distortionRemovalLogic=true(1);
        distortionParametersPath=uigetdir([],'Select a folder where the cameraCBparameters are located');
    case 'No'
        distortionRemovalLogic=false(1);
end

%% 3D recronstruction of overlapping points for image pairs 

% create DLTstructPairs structures
DLTstructPairs=struct;
DLTstructPairs.indCams=indCams;
DLTstructPairs.indPairs=indPairs;
DLTstructPairs.truePoints=C3Dtrue{1};

% if distortion correction is required, correct the points and calculate new DLT parameters
switch distortionRemovalLogic
    case 0 % use original points and DLT parameters
        for ic=1:Ncam
            DLTstructPairs.imageCentroids{ic}=DLTcamStructs{ic}.imageCentroids;
            DLTstructPairs.columns{ic}=DLTcamStructs{ic}.columns;
            DLTstructPairs.DLTparams{ic}=DLTcamStructs{ic}.DLTparams;
            DLTstructPairs.distortionModel{ic}='none';
        end
    case 1 % correct points and DLT parameters
        for ic=1:Ncam
            icam=indCams(ic);
            P2Dtemp=DLTcamStructs{ic}.imageCentroids;
            columns=DLTcamStructs{ic}.columns;
            % undistort image points
            distortionPar=load([distortionParametersPath '\cameraCBparameters_cam_' num2str(icam)]);
            distortionPar=distortionPar.cameraCBparameters.cameraParameters;
            [P2Dtemp] = undistortPoints(P2Dtemp,distortionPar);
            % calculate new DLT parameters
            C3Dtemp=reshape(C3Dtrue{1}(:,columns,:),Nr*numel(columns),3);
            L=DLT11Calibration(P2Dtemp,C3Dtemp);
            % assign to DLTstructPairs
            DLTstructPairs.imageCentroids{ic}=P2Dtemp;
            DLTstructPairs.columns{ic}=columns;
            DLTstructPairs.DLTparams{ic}=L;
            DLTstructPairs.distortionModel{ic}=distortionPar;
        end
end

% calculate recronstructed points and recronstruction errors
DLTreconstructPoints=cell(1,size(DLTstructPairs.indPairs,1));
DLTreconstructErr=cell(1,size(DLTstructPairs.indPairs,1));
for ipair=1:size(DLTstructPairs.indPairs,1) % loop over pairs
    
    % camera indeces of current pair
    icam1=DLTstructPairs.indPairs(ipair,1);
    icount1=find(DLTstructPairs.indCams==icam1);
    icam2=DLTstructPairs.indPairs(ipair,2);
    icount2=find(DLTstructPairs.indCams==icam2);
    
    % columns for camera 1 and 2
    columns1=DLTstructPairs.columns{icount1};
    columns2=DLTstructPairs.columns{icount2};
    
    % overlapping columns for current pair
    [columnsMutual,ind1,ind2]=intersect(columns1,columns2,'stable');
    NpPair=numel(columnsMutual)*Nr;
    
    % overlapping image points from camera 1 and 2
    P2D1=DLTstructPairs.imageCentroids{icount1}((ind1(1)-1)*Nr+1:ind1(end)*Nr,:);
    P2D2=DLTstructPairs.imageCentroids{icount2}((ind2(1)-1)*Nr+1:ind2(end)*Nr,:);
    
    % DLT parameters for camera 1 and 2
    L1=DLTstructPairs.DLTparams{icount1}; % DLT parameters for 1st camera
    L2=DLTstructPairs.DLTparams{icount2}; % DLT parameters for 2nd camera
     
    % recronstruction errors between true points and reconstructed points
    % Solve the linear system with the LS method for the unknown X Y and Z
    P3D=DLT11Reconstruction(P2D1,P2D2,L1,L2);
    % save into cell array
    DLTreconstructPoints{ipair}=P3D;
    
    % find the true points for this pair
    C3DtruePair=C3Dtrue{1}(:,columnsMutual,:);
    P3Dtrue=reshape(C3DtruePair,NpPair,3);

    %recronstruction error vectors
    DLTreconstructErr{ipair}=P3D-P3Dtrue;
    
end
DLTstructPairs.reconstructPoints=DLTreconstructPoints;
DLTstructPairs.reconstructErrors=DLTreconstructErr;

%% save DLTstruct

if saveDLTpairsLogic  
    saveName=[savePath '\DLTstructPairs.mat'];
    icount=1;
    while exist(saveName,'file')
        saveName=[savePath '\DLTstructPairs(' num2str(icount) ').mat'];
        icount=icount+1;
    end
    save(saveName,'DLTstructPairs');
end

%% Plot recronstructed points and recronstruction errors
% this section can be executed also by loading an existing DLTstruct

plotCalibrationResults(DLTstructPairs,saveDLTpairsLogic,figuresPath);

%% finish

h = msgbox('STEP1p is completed');
h.CurrentAxes.Children.FontSize=11;

set(0, 'DefaultUIControlFontSize', fs);

%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>

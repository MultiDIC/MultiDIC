%% STEP 0: Calculate Distortion Parameters
% This is a main script to perform the following steps:
% 1) Estimate camera parameters from checkerboard images
% 2) Use these parameters to correct the distortion
% 3) Plot the camera parameters before and after the correction
% 4) Save the parameters to correct the distortion from image points in Step 1p and Step 3.

%%
clear all; close all

fs=get(0, 'DefaultUIControlFontSize');
set(0, 'DefaultUIControlFontSize', 10);

%% NEW or REPEAT?
% select repeat if you already ran analysis with these same images, but want to repeat with a different distortion model
% repeat is faster because the point detection is skipped (which is the slowest part in this code)

repeatButton = questdlg('New analysis of repeat calibration with a different model?', 'New analysis of repeat calibration with a different model?', 'New', 'Repeat', 'New');
switch repeatButton
    case 'Repeat'
        repeatLogic=true(1);
    case 'New'
        repeatLogic=false(1);
end

%% CHOOSE PATHS OPTIONS

if repeatLogic % if a repeated analysis, request to point to the existing parameters
    initialPath=pwd;
    folder_paths = uipickfiles('FilterSpec',initialPath,'Prompt','Select one or multiple cameraCBparameters files');
else
    % initial image path. Change this path if you want the UI to start with a specific path. Otherwise leave [] or pwd.
    imagePath=pwd;
    % select the folder containing the checkerboard images (if imagePath=[] then the initial path is the current path)
    folder_paths = uipickfiles('FilterSpec',imagePath,'Prompt','Select one or multiple folders containing checkerboard images for analysis');
end

% save camera parameters? y/n. if y, choose save path and overwrite options
[saveCameraParametersLogic,savePath]=QsaveCameraParameters(folder_paths);

% save undistorted images? y/n. if yes, warn for overwriting. undistorted images will be saved in savePath
[saveUndistortedImagesLogic,overWriteUDimagesLogic]=QsaveUndistortedImages;

% figures path to save the plotted results
figuresPath=fullfile(savePath, 'figures');
warning('off','MATLAB:MKDIR:DirectoryExists');
mkdir(figuresPath);

%% SELECT CHECKERBOARD PARAMETERS

if repeatLogic % if repeated analysis, extract the checkerboard parameters from the files
    cameraCBparameters1=load(folder_paths{1});
    Nrows=cameraCBparameters1.cameraCBparameters.boardSize(1); 
    Ncols=cameraCBparameters1.cameraCBparameters.boardSize(2); 
    squareSize=cameraCBparameters1.cameraCBparameters.squareSize;
else % new analysis. ask for checkerboard parameters
    %initial parameters: these are the default parameters which appear in the dialog box. Change them if you want other numbers to appear instead
    Nrows=15; % Number of black rows (should be uneven)
    Ncols=20; % Number of black columns (should be even)
    squareSize=10; %[mm]
    % dialog box
    answer = inputdlg({'Enter number of rows (odd):','Enter number of columns (even):','Enter square size [mm]:'},'Input',[1,50],{num2str(Nrows),num2str(Ncols),num2str(squareSize)});
    % extract answers
    Nrows=str2double(answer{1}); 
    Ncols=str2double(answer{2}); 
    squareSize=str2double(answer{3});
end

%% SELECT DISTORTION MODEL

% dialog box for selecting the distortion model. 
% The default is the full model (3 radial parameters, 2 tangential (1 stands for true), and a skew parameter (1 = true).
answer = inputdlg({'Enter number of radial distortion coefficients (2 or 3):',...
    'Estimate tangential distortion? (1 or 0 for yes/no):',...
    'Estimate skew? (1 or 0 for yes/no)'},...
    'Input',[1,70],{'3','1','1'});
optStruct=struct;
optStruct.NumRadialDistortionCoefficients=str2num(answer{1});
optStruct.EstimateTangentialDistortion=logical(str2num(answer{2}));
optStruct.EstimateSkew=logical(str2num(answer{3}));

%% compute and plot camera parameters for each camera

Ncam=numel(folder_paths); % number of cameras in this analysis
cameraCBparametersAllCams=cell(Ncam,1); % assign cell array for all camera parmaters

for ic=1:Ncam % loop over all cameras
    
    % Extract images and info
    if repeatLogic % if repeated analysis, extract cameraCBparameters and image info from files
        cameraCBparameters=load(folder_paths{ic});
        cameraCBparameters=cameraCBparameters.cameraCBparameters;
        CBimagesInfo=cameraCBparameters.imagesInfo;
        I1=imread(CBimagesInfo.imageFileNames{1}); % read first image from path
        I=zeros([size(I1) CBimagesInfo.Nimages],'uint8'); % reallocate image var
        for ip=1:CBimagesInfo.Nimages
            I(:,:,:,ip)=imread(CBimagesInfo.imageFileNames{ip}); % load all images into I
        end
        CBimagesInfo.I=I;
        % plot all images in one figure
        plotAllCameraImages(CBimagesInfo);
        % re-calculate the distortion parameters with the selected model 
        hm=msgbox(['Please wait while computing distortion parameters for camera ' num2str(CBimagesInfo.icam)]);
        cameraCBparameters=RecalculateCBcalibrationParameters(cameraCBparameters,optStruct);
        delete(hm);
    else
        % if New, extract only image info
        CBimagesInfo=extractImagesInfo(folder_paths{ic});
        % plot all images in one figure
        plotAllCameraImages(CBimagesInfo);
        % Extract images, Detect the checkerboard points, calculate camera paramters, and save a structure containing all necessary parameters
        set(0, 'DefaultUIControlFontSize', 11);
        hm=msgbox(['Please wait while computing distortion parameters for camera ' num2str(CBimagesInfo.icam)]);
        cameraCBparameters=calculateCBcalibrationParameters(CBimagesInfo,squareSize,optStruct);
        delete(hm);
            % check if detected boardsize matches entered values
            if (cameraCBparameters.boardSize(1)~=Nrows) || (cameraCBparameters.boardSize(2)~=Ncols)
                error('Detected number of columns or rows does not match entered values');
            end
    end
    
    % plot camera parameters and reorojection errors before and after distortion correction
    plot_camera_parameters_2tabs(cameraCBparameters);
        if saveCameraParametersLogic
            savefig(fullfile(figuresPath,[ 'params_cam' num2str(CBimagesInfo.icam)]));
        end
    
    % plot reprojected points vs. true points and straight lines on each image
    plot_reprojectVSreal_points(CBimagesInfo,cameraCBparameters);
    
    % undistort images and save if required
    hm=msgbox(['Please wait while correcting distortion from images of ' num2str(CBimagesInfo.icam)]);
    undistortImagesSavePlot(CBimagesInfo,cameraCBparameters,saveUndistortedImagesLogic,overWriteUDimagesLogic,savePath)
    delete(hm);
    
    % remove images with high errors
    answer=questdlg('Do you want to remove some of the images with higher errors?','Remove images?','Yes','No','No');
    switch answer
        case 'Yes'
            % reprojection errors per image
            Error=cameraCBparameters.cameraParameters.ReprojectionErrors;
            ErrorMgn=squeeze(sqrt(Error(:,1,:).^2+Error(:,2,:).^2));
            meanError=squeeze(mean(ErrorMgn,1));
            prctile80Error=prctile(meanError,80);
            % select threshold
            answer=inputdlg('Select threshold (mean error in pixels) for removing images','select threshold',1,{num2str(prctile80Error)});
            threshold=str2num(answer{1});
            indImgHighError=cameraCBparameters.imagesUsed(find(meanError>threshold));
            if length(indImgHighError)>0
                msgbox([ num2str(length(indImgHighError)) ' images will be removed from the calibration parameters calculation']);
                fileNamesImgHighError=cameraCBparameters.imagesInfo.imageFileNames(indImgHighError)';
                % delete high error images
                CBimagesInfo.I(:,:,:,indImgHighError)=[];
                CBimagesInfo.imageFileNames(indImgHighError)=[];
                CBimagesInfo.Nimages=numel(CBimagesInfo.imageFileNames);
                % plot all images in one figure
                plotAllCameraImages(CBimagesInfo);
                % Extract images, Detect the checkerboard points, calculate camera paramters, and save a structure containing all necessary parameters
                hm=msgbox(['Please wait while re-computing distortion parameters for camera ' num2str(CBimagesInfo.icam)]);
                cameraCBparameters=calculateCBcalibrationParameters(CBimagesInfo,squareSize,optStruct);
                delete(hm);
                % plot camera parameters and reorojection errors before and after distortion correction
                plot_camera_parameters_2tabs(cameraCBparameters);
                if saveCameraParametersLogic
                    savefig(fullfile(figuresPath,[ 'params_cam', num2str(CBimagesInfo.icam)]));
                end
            end
    end
    
    % save camera parameters into the cell array of all cameras
    cameraCBparametersAllCams{ic}=cameraCBparameters;
    % save parameters into savePath
    if saveCameraParametersLogic
        save(fullfile(savePath, ['cameraCBparameters_cam_', num2str(cameraCBparameters.icam)]),'cameraCBparameters');
    end
    
end

% save cell array containing the camera parameters for all cameras in this analysis
if saveCameraParametersLogic
    save(fullfile(savePath, 'cameraCBparametersAllCams'),'cameraCBparametersAllCams');
end

%% plot camera instrinsic statistics for all cameras (if more than 1 camera)

Ncam=numel(cameraCBparametersAllCams);
if Ncam>1
    plotButton = questdlg('Plot intrinsic parameters statistics for all cameras?', 'Plot?', 'Yes', 'No', 'Yes');    
    if strcmp(plotButton,'Yes')         
        plotIntrinsicStatsAll(cameraCBparametersAllCams);
                if saveCameraParametersLogic   
                    savefig(fullfile(figuresPath, 'IntrinsicStats')); 
                end        
    end
end

%% finish

h=msgbox('STEP0 is completed');
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
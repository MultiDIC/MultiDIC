%% STEP0 Calculate Distortion Parameters
%Estimate camera parameters from checkerboard images and undistort them

%%
clear all; close all

%% CHOOSE PATHS OPTIONS
% initial image path
imagePath='C:\Users\Dana\Dropbox (Personal)\MIT\research\DIC\360_DIC\LimbTest\RasPi_Intrinsic_Calibration\CB_images';
% select the folder containing the checkerboard images (if imagePath=[] then the initial path is the current path)
image_folder_paths = uipickfiles('FilterSpec',imagePath,'Prompt','Select one or multiple folders containing checkerboard images for analysis');
% save camera parameters? choose save path and overwrite options
[saveCameraParameters,savePath]=QsaveCameraParameters(image_folder_paths);
% save undistorted images? warn for overwriting
[saveUndistortedImages,overWriteUDimages]=QsaveUndistortedImages;
% figures path
figuresPath=[savePath '\figures'];
mkdir(figuresPath);

%% CHECKERBOARD PARAMETERS
%initial parameters
Nrows=15; % Number of black rows (should be uneven)
Ncols=20; % Number of black columns (should be even)
squareSize=0.01; %[m]
% dialog box
answer = inputdlg({'Enter number of rows (uneven):','Enter number of columns (even):','Enter square size in meters:'},...
    'Input',[1,50],{num2str(Nrows),num2str(Ncols),num2str(squareSize)});
Nrows=str2double(answer{1}); Ncols=str2double(answer{2}); squareSize=str2double(answer{3});

%% CALIBRATION MODEL PARAMETERS
% dialog box
answer = inputdlg({'Enter number of radial distortion coefficients (2 or 3):',...
    'Estimate tangential distortion? (1 or 0 for yes/no):',...
    'Estimate skew? (1 or 0 for yes/no)'},...
    'Input',[1,70],{num2str(3),'0','0'});
optStruct=struct;
optStruct.NumRadialDistortionCoefficients=str2num(answer{1});
optStruct.EstimateTangentialDistortion=logical(str2num(answer{2}));
optStruct.EstimateSkew=logical(str2num(answer{3}));

%% compute camera parameters and display original images with straight lines
Ncam=numel(image_folder_paths);
if exist([savePath '\cameraCBparametersAllCams.mat'],'file')
    load([savePath '\cameraCBparametersAllCams.mat']);
else
    cameraCBparametersAllCams=cell(Ncam,1);
    cameraCBparametersJAllCams=cell(Ncam,1);
end

for ic=1:Ncam % loop over all cameras
    
    % Extract images and info;
    [CBimagesInfo]=extractImagesInfo(image_folder_paths{ic});
    
    % plot all images in one figure
    plotAllCameraImages(CBimagesInfo);
    
    % Extract images, Detect the checkerboard points, calculate camera paramters, and save a structure containing all necessary parameters
    
    cameraCBparameters=calculateCBcalibrationParameters(CBimagesInfo,squareSize,optStruct);
    
    % check if detected boardsize matches entered values
    if (cameraCBparameters.boardSize(1)~=Nrows) || (cameraCBparameters.boardSize(2)~=Ncols)
        error('Detected number of columns or rows does not match entered values');
    end
    
    % plot camera reprojection errors and Extrinsics
    plot_camera_parameters(cameraCBparameters);
    if saveCameraParameters
        savefig([figuresPath '\params_cam' num2str(CBimagesInfo.icam)]);
    end
    
    % plot camera reprojection errors and Extrinsics after "undistortion"
    plot_camera_parameters_AUD(cameraCBparameters);
    if saveCameraParameters
        savefig([figuresPath '\paramsAUD_cam' num2str(CBimagesInfo.icam)]);
    end
    
    % plot reprojected points vs. true points and straight lines on each image
    plot_reprojectVSreal_points(CBimagesInfo,cameraCBparameters);
%     if saveCameraParameters
%         savefig([figuresPath '\reprojections_cam' num2str(CBimagesInfo.icam)]);
%     end
    
    % undistort images and save if required
    undistortImagesSavePlot(CBimagesInfo,cameraCBparameters,overWriteUDimages,overWriteUDimages,image_folder_paths{ic})
%     if saveCameraParameters
%         savefig([figuresPath '\reprojectionsAUD_cam' num2str(CBimagesInfo.icam)]);
%     end
    
    % save camera parameters into the cell array of all cameras
    cameraCBparametersAllCams{ic}=cameraCBparameters;
    % save parameters into savePath
    if saveCameraParameters
        save([savePath '\cameraCBparameters_cam' num2str(cameraCBparameters.icam,'%02i')],'cameraCBparameters');
    end
    
end

if saveCameraParameters
    save([savePath '\cameraCBparametersAllCams'],'cameraCBparametersAllCams');
end

%% plot camera instrinsic statistics
plotIntrinsicStats(cameraCBparametersAllCams,0);
if saveCameraParameters,     savefig([figuresPath '\IntrinsicStats']); end

plotIntrinsicStats(cameraCBparametersAllCams,1);
if saveCameraParameters,     savefig([figuresPath '\IntrinsicStatsAUD']); end

plotIntrinsicStatsBA(cameraCBparametersAllCams);
if saveCameraParameters,     savefig([figuresPath '\IntrinsicStatsBAUD']); end

msgbox('STEP0 is finished');
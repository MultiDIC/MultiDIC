function []=plot_camera_parameters_2tabs(cameraCBparameters)
%% function for plotting the camera parameters in STEP0
%
% INPUTS:
% * cameraCBparameters:  a structure containing all the calibration parameters created in STEP0
%
% OUTPUTS:
% * a figure plotting the intrinsic, extrinsic, and reprojection errors for
% one camera, before and after distortion correction.

%%

f = figure('name',['Camera ' num2str(cameraCBparameters.icam) ' Parameters'],'units','normalized','outerposition',[.05 .05 .9 .9]);

tabgp = uitabgroup(f);

%% before correction
tab(1) = uitab(tabgp,'Title','Before Distortion Correction');
axes('Parent',tab(1));

%% top row
% show Intrinsic parameters and their errors
subplot(6,8,[1 2 3 9 10 11 17 18 19])
showIntrinsics(cameraCBparameters);

% visualize the distortion model
%     subplot(6,7,[3 4 10 11 17 18])
%     visulaizeDistortionModel(cameraCBparameters);

% extrinsic visualization (using MATLAB function)
subplot(6,8,[4:8 12:16 20:24])
showExtrinsics(cameraCBparameters.cameraParameters);

%% bottom row
% show reprojection errors statistics in boxplot
subplot(6,8,[33 34 41 42])
showReprojectionErrorsBox(cameraCBparameters.cameraParameters);

% All reprojection error per image (point distribution)
subplot(6,8,[35 36 37 43 44 45])
showReprojectionPointsErrors(cameraCBparameters.cameraParameters);

% mean reprojection error per image (using MATLAB function)
subplot(6,8,[38 39 40 46 47 48])
showReprojectionErrors(cameraCBparameters.cameraParameters);

drawnow

%% after correction

tab(2) = uitab(tabgp,'Title','After Distortion Correction');
axes('Parent',tab(2));

cameraCBparametersMod=cameraCBparameters;
cameraCBparametersMod.cameraParameters=cameraCBparameters.cameraParametersAUD;
cameraCBparametersMod.estimationErrors=cameraCBparameters.estimationErrorsAUD;


%% top row
% show Intrinsic parameters and their errors
subplot(6,8,[1 2 3 9 10 11 17 18 19])
showIntrinsics(cameraCBparametersMod);

% visualize the distortion model
%     subplot(6,7,[3 4 10 11 17 18])
%     visulaizeDistortionModel(cameraCBparametersMod);

% extrinsic visualization
subplot(6,8,[4:8 12:16 20:24])
showExtrinsics(cameraCBparametersMod.cameraParametersAUD);

%% bottom row
% show reprojection errors statistics in boxplot
subplot(6,8,[33 34 41 42])
showReprojectionErrorsBox(cameraCBparametersMod.cameraParametersAUD);

% All reprojection error per image (point distribution)
subplot(6,8,[35 36 37 43 44 45])
showReprojectionPointsErrors(cameraCBparametersMod.cameraParametersAUD);

% mean reprojection error per image
subplot(6,8,[38 39 40 46 47 48])
showReprojectionErrors(cameraCBparametersMod.cameraParametersAUD);


end

 
%% 
% MultiDIC: a MATLAB Toolbox for Multi-View 3D Digital Image Correlation
% 
% License: <https://github.com/MultiDIC/MultiDIC/blob/master/LICENSE.txt>
% 
% Copyright (C) 2018  Dana Solav
% 
% If you use the toolbox/function for your research, please cite our paper:
% <https://engrxiv.org/fv47e>
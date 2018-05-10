function [cameraCBparameters] = RecalculateCBcalibrationParameters(cameraCBparameters,optStruct)
%% function for re-calculating the distortion parameters in STEP0.
% This function is called only in case the user selected a repeated
% analysis.
%
% INPUTS:
% * cameraCBparameters: structure previously created in STEP0
% * optStruct: a structure containing the new parameters of the distortion
%   model selected in the repeated analysis.
%
% OUTPUTS:
% * a figure of all the checkerboard images

%%

%  extract the points already detected in previous analysis
imagePoints=cameraCBparameters.imagePoints;
worldPoints=cameraCBparameters.cameraParameters.WorldPoints;

% Re-calculate camera parameters with the new distortion model
[params,~,estimationErrors] = estimateCameraParameters(imagePoints,worldPoints,...
    'NumRadialDistortionCoefficients',optStruct.NumRadialDistortionCoefficients,'EstimateTangentialDistortion',optStruct.EstimateTangentialDistortion,'EstimateSkew',optStruct.EstimateSkew);

% feed new results into cameraCBparameters
cameraCBparameters.cameraParameters=params;
cameraCBparameters.estimationErrors=estimationErrors;

%parameters after undistortion
imagePointsUndistorted=zeros(size(cameraCBparameters.imagePoints));
for ii=1:size(cameraCBparameters.imagePoints,3)
    imagePointNow=cameraCBparameters.imagePoints(:,:,ii);
    [imagePointsUndistorted(:,:,ii)] = undistortPoints(imagePointNow,cameraCBparameters.cameraParameters);
end
[paramsJ,~,estimationErrorsJ] = estimateCameraParameters(imagePointsUndistorted,worldPoints,...
    'NumRadialDistortionCoefficients',optStruct.NumRadialDistortionCoefficients,'EstimateTangentialDistortion',optStruct.EstimateTangentialDistortion,'EstimateSkew',optStruct.EstimateSkew);

% feed results after undistortion into cameraCBparameters and return
cameraCBparameters.cameraParametersAUD=paramsJ;
cameraCBparameters.estimationErrorsAUD=estimationErrorsJ;
cameraCBparameters.imagePointsAUD=imagePointsUndistorted;

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
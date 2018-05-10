function [cameraCBparameters] = calculateCBcalibrationParameters(CBimagesInfo,squareSize,varargin)
%% function for calculating the distortion parameters in STEP0.
% This function is called only in case the user selected a repeated
% analysis.
%
% INPUTS:
% * CBimagesInfo: structure previously created in STEP0
% * squareSize: scalar, in meters.
% * optStruct (optional input): a structure containing the distortion model. If not given, the default (full) model is used
%
% OUTPUTS:
% * cameraCBparameters: a structure containing all the calibration parameters, with the following fields:
% - icam
% - cameraParameters
% - imagesUsed
% - estimationErrors
% - boardSize
% - squareSize
% - imagePoints
% - imagesInfo
% - cameraParametersAUD
% - estimationErrorsAUD
% - imagePointsAUD

%%

Narg=numel(varargin);
switch Narg
    case 0
        optStruct=struct;
    case 1
        optStruct=varargin{1};
    otherwise
        error('wrong number of input arguments');
end

% fill in missing fields
if ~isfield(optStruct, 'NumRadialDistortionCoefficients')
    optStruct.NumRadialDistortionCoefficients=3;
end
if ~isfield(optStruct, 'EstimateTangentialDistortion')
    optStruct.EstimateTangentialDistortion=true;
end
if ~isfield(optStruct, 'EstimateSkew')
    optStruct.EstimateSkew=true;
end

% Detect the checkerboard points
[imagePoints, boardSize, imagesUsedCB] = detectCheckerboardPoints(CBimagesInfo.I);


% generate the real checkerboard points from the known square size and board size (the 1000 factor is to convert from meters to mm)
worldPoints = generateCheckerboardPoints(boardSize,squareSize*1000);

% calculate camera parameters
[params,imagesUsedECP,estimationErrors] = estimateCameraParameters(imagePoints,worldPoints,...
    'NumRadialDistortionCoefficients',optStruct.NumRadialDistortionCoefficients,'EstimateTangentialDistortion',optStruct.EstimateTangentialDistortion,'EstimateSkew',optStruct.EstimateSkew);

imagesUsedFinal=find(imagesUsedCB)';
imagesUsedFinal(find(~imagesUsedECP))=[];

% create cameraCBparameters structure
cameraCBparameters=struct;
% feed results
cameraCBparameters.icam=CBimagesInfo.icam;
cameraCBparameters.cameraParameters=params;
cameraCBparameters.imagesUsed=imagesUsedFinal;
cameraCBparameters.estimationErrors=estimationErrors;
cameraCBparameters.boardSize=boardSize;
cameraCBparameters.squareSize=squareSize;
cameraCBparameters.imagePoints=imagePoints;

imagesInfo=rmfield(CBimagesInfo,'I');
cameraCBparameters.imagesInfo=imagesInfo;

%parameters after undistortion
imagePointsUndistorted=zeros(size(cameraCBparameters.imagePoints));
for ii=1:size(cameraCBparameters.imagePoints,3)
    imagePointNow=cameraCBparameters.imagePoints(:,:,ii);
    [imagePointsUndistorted(:,:,ii)] = undistortPoints(imagePointNow,cameraCBparameters.cameraParameters);
end
[paramsJ,~,estimationErrorsJ] = estimateCameraParameters(imagePointsUndistorted,worldPoints,...
    'NumRadialDistortionCoefficients',optStruct.NumRadialDistortionCoefficients,'EstimateTangentialDistortion',optStruct.EstimateTangentialDistortion,'EstimateSkew',optStruct.EstimateSkew);

% feed results after undistortion
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
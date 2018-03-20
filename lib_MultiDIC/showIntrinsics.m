function []=showIntrinsics(cameraCBparameters)
%% function for plotting the radial distortion, tangential distortion, and skew and their uncenterties, as text, in STEP0
% The function is used inside the function plot_camera_parameters_2tabs 
%
% INPUTS:
% * cameraCBparameters:  a structure containing all the calibration parameters created in STEP0
%
%%
% plot radial distortion, tangential distortion, and skew and their uncenterties

% focal length
fx=cameraCBparameters.cameraParameters.FocalLength(1);
fy=cameraCBparameters.cameraParameters.FocalLength(2);
fxError=cameraCBparameters.estimationErrors.IntrinsicsErrors.FocalLengthError(1);
fyError=cameraCBparameters.estimationErrors.IntrinsicsErrors.FocalLengthError(2);
FocalLengthChar=[ '[' num2str(fx,5) ' \pm ' num2str(fxError,2), ' ,  '   num2str(fy,5) ' \pm ' num2str(fyError,2) ']' ];

% principal point
px=cameraCBparameters.cameraParameters.PrincipalPoint(1);
py=cameraCBparameters.cameraParameters.PrincipalPoint(2);
pxError=cameraCBparameters.estimationErrors.IntrinsicsErrors.PrincipalPointError(1);
pyError=cameraCBparameters.estimationErrors.IntrinsicsErrors.PrincipalPointError(2);;
PrincipalPointChar=[ '[' num2str(px,5) ' \pm ' num2str(pxError,2), ' ,  '   num2str(py,5) ' \pm ' num2str(pyError,2) ']' ];

% skew
if ~cameraCBparameters.cameraParameters.EstimateSkew
    SkewChar=['0 (not estimated)'];
else
    skew=cameraCBparameters.cameraParameters.Skew;
    skewError=cameraCBparameters.estimationErrors.IntrinsicsErrors.SkewError;
    skewDeg=atan(skew/fy);
    skewDegError=atan(skewError/fy);
    SkewChar=[num2str(skew,4) ' \pm ' num2str(skewError,4) '  (' num2str(skewDeg,4)  '\circ \pm ' num2str(skewDegError,4) '\circ )'];
end

% tangential distortion
if ~cameraCBparameters.cameraParameters.EstimateTangentialDistortion
    TangentialDistortionChar=['[0 0] (not estimated)'];
else
    p1=cameraCBparameters.cameraParameters.TangentialDistortion(1);
    p2=cameraCBparameters.cameraParameters.TangentialDistortion(2);
    p1Error=cameraCBparameters.estimationErrors.IntrinsicsErrors.TangentialDistortionError(1);
    p2Error=cameraCBparameters.estimationErrors.IntrinsicsErrors.TangentialDistortionError(2);
    TangentialDistortionChar=[ '[' num2str(p1,4) ' \pm ' num2str(p1Error,4),' ,  '   num2str(p2,4) ' \pm ' num2str(p2Error,4) ']' ];
end

% radial distortion
k=cameraCBparameters.cameraParameters.RadialDistortion;
kError=cameraCBparameters.estimationErrors.IntrinsicsErrors.RadialDistortionError;
if cameraCBparameters.cameraParameters.NumRadialDistortionCoefficients==3
    RadialDistortionChar=[ '[' num2str(k(1),4) ' \pm ' num2str(kError(1),4),' ,  '   num2str(k(2),4) ' \pm ' num2str(kError(2),4),' ,  '   num2str(k(3),4) ' \pm ' num2str(kError(3),4) ']' ];
else
   RadialDistortionChar=[ '[' num2str(k(1),4) ' \pm ' num2str(kError(1),4), ' ,  '   num2str(k(2),4) ' \pm ' num2str(kError(2),4) '] (k3 not estimated)'];
end


text(0.01,1,{'Instrinsic Parameters'},'fontsize',14,'fontweight','bold');
text(0.01,0.81,{'\bf Focal length [fx,fy]'; ['\rm' FocalLengthChar]},'fontsize',12);
text(0.01,0.61,{'\bf Principal Point'; ['\rm' PrincipalPointChar]},'fontsize',12);
text(0.01,0.41,{'\bf RadialDistortion [k1,k2,k3]'; ['\rm' RadialDistortionChar]},'fontsize',12);
text(0.01,0.21,{'\bf Tangential Distortion [p1,p2]'; ['\rm' TangentialDistortionChar]},'fontsize',12);
text(0.01,0.01,{'\bf Skew [s]'; ['\rm' SkewChar]},'fontsize',12);
axis off

% title('Intrinsic Parameters');

end
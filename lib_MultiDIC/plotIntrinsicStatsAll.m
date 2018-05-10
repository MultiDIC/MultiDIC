function [] = plotIntrinsicStatsAll(cameraCBparametersAllCams)
%% function for plotting intrinsics statistics for all cameras in STEP0
%
% INPUTS:
% * cameraCBparametersAllCams: a cell array containing cameraCBparameters structures for all cameras
%
% OUTPUTS:
% * a figure plotting intrinsics statistics for all cameras
% tab 1-  plot original parameters
% tab 2- plot parameters after distortion
% tab 3- original and corrected together

%%
hf=figure('units','normalized','outerposition',[.1 .1 .8 .8]);

Ncam=numel(cameraCBparametersAllCams);

tabgp = uitabgroup(hf);

%% BEFORE CORRECTION

tab(1) = uitab(tabgp,'Title','Before distortion correction');
axes('Parent',tab(1));


for ic=1:Ncam
    RadialDistortion(ic,:)=cameraCBparametersAllCams{ic}.cameraParameters.RadialDistortion;
    RadialDistortionError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrors.IntrinsicsErrors.RadialDistortionError;
    TangentialDistortion(ic,:)=cameraCBparametersAllCams{ic}.cameraParameters.TangentialDistortion;
    TangentialDistortionError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrors.IntrinsicsErrors.TangentialDistortionError;
    Skew(ic,:)=cameraCBparametersAllCams{ic}.cameraParameters.Skew;
    SkewError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrors.IntrinsicsErrors.SkewError;
    FocalLength(ic,:)=cameraCBparametersAllCams{ic}.cameraParameters.FocalLength;
    FocalLengthError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrors.IntrinsicsErrors.FocalLengthError;
    PrincipalPoint(ic,:)=cameraCBparametersAllCams{ic}.cameraParameters.PrincipalPoint;
    PrincipalPointError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrors.IntrinsicsErrors.PrincipalPointError;
    imageSizeAll(ic,:)=cameraCBparametersAllCams{ic}.imagesInfo.imageSize;
    camInds(ic)=cameraCBparametersAllCams{ic}.icam;
end

hf.Name=['camera intrinsic parameters statistics for cameras [' num2str(camInds) ']'];

subplot(2,5,1)
if size(RadialDistortion,2)==3
    boxplot(RadialDistortion,'Labels',{'k1','k2','k3'});
else
    boxplot(RadialDistortion,'Labels',{'k1','k2'});
end
ylim([-max(max(abs(RadialDistortion))) max(max(abs(RadialDistortion)))]);
hline = refline(0,0); hline.Color = 'g'; hline.LineStyle='--';
title('Radial Distorsion');
subplot(2,5,2,'align')
boxplot(TangentialDistortion,'Labels',{'p1','p2'});
%     ylim([-max(max(abs(TangentialDistortion))) max(max(abs(TangentialDistortion)))]);
hline =refline(0,0); hline.Color = 'g'; hline.LineStyle='--';
title('Tangential Distorsion');
subplot(2,5,3,'align')
boxplot(Skew,'Labels',{'s'});
hline = refline(0,0); hline.Color = 'g'; hline.LineStyle='--';
title('Skew');
subplot(2,5,4,'align')
boxplot(FocalLength,'Labels',{'fx','fy'});
title('Focal Length [pixel]');
subplot(2,5,5,'align')
boxplot(PrincipalPoint,'Labels',{'px','py'});
% plot principal points refline only if all cameras has same size
if sum(sum(imageSizeAll-mean(imageSizeAll)))==0
    hline =refline(0,imageSizeAll(1,1)/2); hline.Color = 'g'; hline.LineStyle='--';
    hline =refline(0,imageSizeAll(1,2)/2); hline.Color = 'g'; hline.LineStyle='--';
end
title('Principal Point [pixel]');


subplot(2,5,6,'align')
if size(RadialDistortionError,2)==3
    boxplot(RadialDistortionError,'Labels',{'k1','k2','k3'});
else
    boxplot(RadialDistortionError,'Labels',{'k1','k2'});
end
title('Radial Distorsion Error');
subplot(2,5,7,'align')
boxplot(TangentialDistortionError,'Labels',{'p1','p2'});
title('Tangential Distorsion Error');
subplot(2,5,8,'align')
boxplot(SkewError,'Labels',{'s'});
title('Skew Error');
subplot(2,5,9,'align')
boxplot(FocalLengthError,'Labels',{'fx','fy'});
title('Focal Length Error [pixel]');
subplot(2,5,10,'align')
boxplot(PrincipalPointError,'Labels',{'px','py'});
title('Principal Point Error [pixel]');


%% AFTER CORRECTION

tab(2) = uitab(tabgp,'Title','After distortion correction');
axes('Parent',tab(2));

for ic=1:Ncam
    RadialDistortion(ic,:)=cameraCBparametersAllCams{ic}.cameraParametersAUD.RadialDistortion;
    RadialDistortionError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrorsAUD.IntrinsicsErrors.RadialDistortionError;
    TangentialDistortion(ic,:)=cameraCBparametersAllCams{ic}.cameraParametersAUD.TangentialDistortion;
    TangentialDistortionError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrorsAUD.IntrinsicsErrors.TangentialDistortionError;
    Skew(ic,:)=cameraCBparametersAllCams{ic}.cameraParametersAUD.Skew;
    SkewError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrorsAUD.IntrinsicsErrors.SkewError;
    FocalLength(ic,:)=cameraCBparametersAllCams{ic}.cameraParametersAUD.FocalLength;
    FocalLengthError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrorsAUD.IntrinsicsErrors.FocalLengthError;
    PrincipalPoint(ic,:)=cameraCBparametersAllCams{ic}.cameraParametersAUD.PrincipalPoint;
    PrincipalPointError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrorsAUD.IntrinsicsErrors.PrincipalPointError;
    imageSizeAll(ic,:)=cameraCBparametersAllCams{ic}.imagesInfo.imageSize;
    camInds(ic)=cameraCBparametersAllCams{ic}.icam;
end

subplot(2,5,1)
if size(RadialDistortion,2)==3
    boxplot(RadialDistortion,'Labels',{'k1','k2','k3'});
else
    boxplot(RadialDistortion,'Labels',{'k1','k2'});
end
ylim([-max(max(abs(RadialDistortion))) max(max(abs(RadialDistortion)))]);
hline = refline(0,0); hline.Color = 'g'; hline.LineStyle='--';
title('Radial Distorsion');
subplot(2,5,2,'align')
boxplot(TangentialDistortion,'Labels',{'p1','p2'});
%     ylim([-max(max(abs(TangentialDistortion))) max(max(abs(TangentialDistortion)))]);
hline =refline(0,0); hline.Color = 'g'; hline.LineStyle='--';
title('Tangential Distorsion');
subplot(2,5,3,'align')
boxplot(Skew,'Labels',{'s'});
hline = refline(0,0); hline.Color = 'g'; hline.LineStyle='--';
title('Skew [deg]');
subplot(2,5,4,'align')
boxplot(FocalLength,'Labels',{'fx','fy'});
title('Focal Length [pixel]');
subplot(2,5,5,'align')
boxplot(PrincipalPoint,'Labels',{'px','py'});
% plot principal points refline only if all cameras has same size
if sum(sum(imageSizeAll-mean(imageSizeAll)))==0
    hline =refline(0,imageSizeAll(1,1)/2); hline.Color = 'g'; hline.LineStyle='--';
    hline =refline(0,imageSizeAll(1,2)/2); hline.Color = 'g'; hline.LineStyle='--';
end
title('Principal Point [pixel]');


subplot(2,5,6,'align')
if size(RadialDistortionError,2)==3
    boxplot(RadialDistortionError,'Labels',{'k1','k2','k3'});
else
    boxplot(RadialDistortionError,'Labels',{'k1','k2'});
end
title('Radial Distorsion Error');
subplot(2,5,7,'align')
boxplot(TangentialDistortionError,'Labels',{'p1','p2'});
title('Tangential Distorsion Error');
subplot(2,5,8,'align')
boxplot(SkewError,'Labels',{'s'});
title('Skew Error [deg]');
subplot(2,5,9,'align')
boxplot(FocalLengthError,'Labels',{'fx','fy'});
title('Focal Length Error [pixel]');
subplot(2,5,10,'align')
boxplot(PrincipalPointError,'Labels',{'px','py'});
title('Principal Point Error [pixel]');

%% BEFORE AND AFTYER

tab(3) = uitab(tabgp,'Title','Before (b) and after (a) distortion correction');
axes('Parent',tab(3));

for ic=1:Ncam
    RadialDistortion(ic,:)=cameraCBparametersAllCams{ic}.cameraParameters.RadialDistortion;
    RadialDistortionError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrors.IntrinsicsErrors.RadialDistortionError;
    TangentialDistortion(ic,:)=cameraCBparametersAllCams{ic}.cameraParameters.TangentialDistortion;
    TangentialDistortionError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrors.IntrinsicsErrors.TangentialDistortionError;
    Skew(ic,:)=cameraCBparametersAllCams{ic}.cameraParameters.Skew;
    SkewError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrors.IntrinsicsErrors.SkewError;
    FocalLength(ic,:)=cameraCBparametersAllCams{ic}.cameraParameters.FocalLength;
    FocalLengthError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrors.IntrinsicsErrors.FocalLengthError;
    PrincipalPoint(ic,:)=cameraCBparametersAllCams{ic}.cameraParameters.PrincipalPoint;
    PrincipalPointError(ic,:)=cameraCBparametersAllCams{ic}.estimationErrors.IntrinsicsErrors.PrincipalPointError;
    imageSizeAll(ic,:)=cameraCBparametersAllCams{ic}.imagesInfo.imageSize;
    
    RadialDistortionAUD(ic,:)=cameraCBparametersAllCams{ic}.cameraParametersAUD.RadialDistortion;
    RadialDistortionErrorAUD(ic,:)=cameraCBparametersAllCams{ic}.estimationErrorsAUD.IntrinsicsErrors.RadialDistortionError;
    TangentialDistortionAUD(ic,:)=cameraCBparametersAllCams{ic}.cameraParametersAUD.TangentialDistortion;
    TangentialDistortionErrorAUD(ic,:)=cameraCBparametersAllCams{ic}.estimationErrorsAUD.IntrinsicsErrors.TangentialDistortionError;
    SkewAUD(ic,:)=cameraCBparametersAllCams{ic}.cameraParametersAUD.Skew;
    SkewErrorAUD(ic,:)=cameraCBparametersAllCams{ic}.estimationErrorsAUD.IntrinsicsErrors.SkewError;
    FocalLengthAUD(ic,:)=cameraCBparametersAllCams{ic}.cameraParametersAUD.FocalLength;
    FocalLengthErrorAUD(ic,:)=cameraCBparametersAllCams{ic}.estimationErrorsAUD.IntrinsicsErrors.FocalLengthError;
    PrincipalPointAUD(ic,:)=cameraCBparametersAllCams{ic}.cameraParametersAUD.PrincipalPoint;
    PrincipalPointErrorAUD(ic,:)=cameraCBparametersAllCams{ic}.estimationErrorsAUD.IntrinsicsErrors.PrincipalPointError;
    
    camInds(ic)=cameraCBparametersAllCams{ic}.icam;
end

subplot(2,14,[1 2 3])
if size(RadialDistortion,2)==3
    boxplot([RadialDistortion(:,1) RadialDistortionAUD(:,1) RadialDistortion(:,2) RadialDistortionAUD(:,2) RadialDistortion(:,3) RadialDistortionAUD(:,3)],'Labels',{'k1(b)','k1(a)','k2(b)','k2(a)','k3(b)','k3(a)'});
    ylim([-max(max(abs([RadialDistortion RadialDistortionAUD]))) max(max(abs([RadialDistortion RadialDistortionAUD])))]);
    ylimits=ylim;
    line([2.5 2.5],[ylimits(1) ylimits(2)],'color','k');
    line([4.5 4.5],[ylimits(1) ylimits(2)],'color','k');
else
    boxplot([RadialDistortion(:,1) RadialDistortionAUD(:,1) RadialDistortion(:,2) RadialDistortionAUD(:,2)],'Labels',{'k1(b)','k1(a)','k2(b)','k2(a)'});
    ylim([-max(max(abs([RadialDistortion RadialDistortionAUD]))) max(max(abs([RadialDistortion RadialDistortionAUD])))]);
    ylimits=ylim;
    line([2.5 2.5],[ylimits(1) ylimits(2)],'color','k');
end
hline = refline(0,0); hline.Color = 'g'; hline.LineStyle='--';
title('Radial Distorsion Parameters');

subplot(2,14,[5 6],'align')
boxplot([TangentialDistortion(:,1) TangentialDistortionAUD(:,1) TangentialDistortion(:,2) TangentialDistortionAUD(:,2)],'Labels',{'p1(b)','p1(a)','p2(b)','p2(a)'});
ylimits=ylim;
line([2.5 2.5],[ylimits(1) ylimits(2)],'color','k');
%     ylim([-max(max(abs(TangentialDistortion))) max(max(abs(TangentialDistortion)))]);
hline =refline(0,0); hline.Color = 'g'; hline.LineStyle='--';
title('Tangential Distorsion Parameters');

subplot(2,14,8,'align')
boxplot([Skew SkewAUD],'Labels',{'s(b)','s(a)'});
hline = refline(0,0); hline.Color = 'g'; hline.LineStyle='--';
title('Skew Parameter');

subplot(2,14,[10 11],'align')
boxplot([FocalLength(:,1) FocalLengthAUD(:,1) FocalLength(:,2) FocalLengthAUD(:,2)],'Labels',{'fx(b)','fx(a)','fy(b)','fy(a)'});
ylimits=ylim;
line([2.5 2.5],[ylimits(1) ylimits(2)],'color','k');
title('Focal Length [pixel]');

subplot(2,14,[13 14],'align')
boxplot([PrincipalPoint(:,1) PrincipalPointAUD(:,1) PrincipalPoint(:,2) PrincipalPointAUD(:,2)],'Labels',{'px(b)','px(a)','py(b)','py(a)'});
ylimits=ylim;
line([2.5 2.5],[ylimits(1) ylimits(2)],'color','k');

% plot principal points refline only if all cameras has same size
if sum(sum(imageSizeAll-mean(imageSizeAll)))==0
    hline =refline(0,imageSizeAll(1,1)/2); hline.Color = 'g'; hline.LineStyle='--';
    hline =refline(0,imageSizeAll(1,2)/2); hline.Color = 'g'; hline.LineStyle='--';
end
title('Principal Point [pixel]');

subplot(2,14,[15 16 17],'align')
if size(RadialDistortionError,2)==3
    boxplot([RadialDistortionError(:,1) RadialDistortionErrorAUD(:,1) RadialDistortionError(:,2) RadialDistortionErrorAUD(:,2) RadialDistortionError(:,3) RadialDistortionErrorAUD(:,3)],'Labels',{'k1(b)','k1(a)','k2(b)','k2(a)','k3(b)','k3(a)'})
    ylimits=ylim;
    line([2.5 2.5],[ylimits(1) ylimits(2)],'color','k');
    line([4.5 4.5],[ylimits(1) ylimits(2)],'color','k');
else
    boxplot([RadialDistortionError(:,1) RadialDistortionErrorAUD(:,1) RadialDistortionError(:,2) RadialDistortionErrorAUD(:,2)],'Labels',{'k1(b)','k1(a)','k2(b)','k2(a)'})
    ylimits=ylim;
    line([2.5 2.5],[ylimits(1) ylimits(2)],'color','k');
end
title('Radial Distorsion Errors');

subplot(2,14,[19 20],'align')
boxplot([TangentialDistortionError(:,1) TangentialDistortionErrorAUD(:,1) TangentialDistortionError(:,2) TangentialDistortionErrorAUD(:,2)],'Labels',{'p1(b)','p1(a)','p2(b)','p2(a)'});
ylimits=ylim;
line([2.5 2.5],[ylimits(1) ylimits(2)],'color','k');
title('Tangential Distorsion Errors');

subplot(2,14,22,'align')
boxplot([SkewError SkewErrorAUD],'Labels',{'s(b)','s(a)'});
title('Skew Error');

subplot(2,14,[24 25],'align')
boxplot([FocalLengthError(:,1) FocalLengthErrorAUD(:,1) FocalLengthError(:,2) FocalLengthErrorAUD(:,2)],'Labels',{'fx(b)','fx(a)','fy(b)','fy(a)'});
ylimits=ylim;
line([2.5 2.5],[ylimits(1) ylimits(2)],'color','k');
title('Focal Length Errors [pixel]');

subplot(2,14,[27 28],'align')
boxplot([PrincipalPointError(:,1) PrincipalPointErrorAUD(:,1) PrincipalPointError(:,2) PrincipalPointErrorAUD(:,2)],'Labels',{'px(b)','px(a)','py(b)','py(a)'});
ylimits=ylim;
line([2.5 2.5],[ylimits(1) ylimits(2)],'color','k');
title('Principal Point Errors [pixel]');



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
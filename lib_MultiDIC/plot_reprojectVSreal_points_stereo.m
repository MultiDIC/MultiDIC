function []=plot_reprojectVSreal_points_stereo(I1,I2,imagePoints,stereoParams,IndpairsUsed)
%% function for plot each image with reprojected points vs. real points and straight lines, and reprojection error statistics in STEP0
%
% INPUTS:
% * CBimagesInfo
% * cameraCBparameters:  a structure containing all the calibration parameters created in STEP0
%
%%

% plot each image with reprojected points vs. real points and straight lines
Nimages=length(IndpairsUsed);
reprojectedPoints1=stereoParams.CameraParameters1.ReprojectedPoints;
reprojectedPoints2=stereoParams.CameraParameters2.ReprojectedPoints;
reprojectionErrors1=stereoParams.CameraParameters1.ReprojectionErrors;
reprojectionErrors2=stereoParams.CameraParameters2.ReprojectionErrors;

if (stereoParams.CameraParameters1.ImageSize(1)>=stereoParams.CameraParameters1.ImageSize(2))
    n=1; m=8; p1=1:3; p2=4; p3=5:7; p4=8;
else
   n=2; m=6; p1=1:5; p2=6; p3=7:11; p4=12;
end
    

% PLOT
% define tabs
f = figure('name','Reprojected points on all images. Scroll between tabs to view the different images','units','normalized','outerposition',[.1 .1 .8 .8]);
tabgp = uitabgroup(f);

icount=0;
for iplot=IndpairsUsed'
    icount=icount+1;
    %plot
    tab(icount) = uitab(tabgp,'Title',['IM' num2str(iplot)]);
    
    axes('Parent',tab(icount)); % somewhere to plot
    
    subplot(n,m,p1)
    imshow(I1(:,:,:,iplot)); hold all;
    plot(imagePoints(:,1,icount,1), imagePoints(:,2,icount,1),'go','linewidth',1.5);
    plot(reprojectedPoints1(:,1,icount),reprojectedPoints1(:,2,icount),'r+','linewidth',1.5);
    title(['Camera left, Image ' num2str(iplot)]);
    drawnow
    legend('Detected Points','Reprojected Points');
    hold off;
    
    reprojectionErrorsNow1=reprojectionErrors1(:,:,icount);
    reprojectionErrorsMgnNow1=sqrt(sum(reprojectionErrorsNow1.^2,2));
    reprojectionErrorsNow1=[reprojectionErrorsNow1 reprojectionErrorsMgnNow1];
    
    subplot(n,m,p2)
    boxplot(reprojectionErrorsNow1,'Labels',{'X','Y','Mgn'});
    ylim([-max(abs(reprojectionErrorsNow1(:))) max(abs(reprojectionErrorsNow1(:)))]);
    title({'Reprojection error'; 'statistics [pix]'});
    
    subplot(n,m,p3)
    imshow(I2(:,:,:,iplot)); hold all;
    plot(imagePoints(:,1,icount,2), imagePoints(:,2,icount,2),'go','linewidth',1.5);
    plot(reprojectedPoints2(:,1,icount),reprojectedPoints2(:,2,icount),'r+','linewidth',1.5);
    title(['Camera right, Image ' num2str(iplot)]);
    drawnow
    legend('Detected Points','Reprojected Points');
    hold off;
    
    reprojectionErrorsNow2=reprojectionErrors2(:,:,icount);
    reprojectionErrorsMgnNow2=sqrt(sum(reprojectionErrorsNow2.^2,2));
    reprojectionErrorsNow2=[reprojectionErrorsNow2 reprojectionErrorsMgnNow2];
    
    subplot(n,m,p4)
    boxplot(reprojectionErrorsNow2,'Labels',{'X','Y','Mgn'});
    ylim([-max(abs(reprojectionErrorsNow2(:))) max(abs(reprojectionErrorsNow2(:)))]);
    title({'Reprojection error'; 'statistics [pix]'});
    
end


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
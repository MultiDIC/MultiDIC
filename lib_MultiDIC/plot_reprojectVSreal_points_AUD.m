function []=plot_reprojectVSreal_points_AUD(CBimagesInfo,cameraCBparameters)
%% function for plotting each image with reprojected points vs. real points, straight lines, and reprojection error statistics, in STEP0
%
% INPUTS:
% * I: matrix containing all CB images (Height x Width x nColors x Nimages)
% * cameraCBparameters: structure containing the following fields:
%       - icam, cameraParameters, imagesUsed, estimationErrors, boardSize, imagePoints

%%
imageIndecesToDisplay=1:numel(cameraCBparameters.imagesUsed);

% plot each image with reprojected points vs. real points and straight lines
I=CBimagesInfo.I;
Nimages=numel(imageIndecesToDisplay);
imagePoints=cameraCBparameters.imagePointsAUD(:,:,imageIndecesToDisplay);
reprojectedPoints=cameraCBparameters.cameraParametersAUD.ReprojectedPoints(:,:,imageIndecesToDisplay);
reprojectionErrors=cameraCBparameters.cameraParametersAUD.ReprojectionErrors(:,:,imageIndecesToDisplay);
Nrows=cameraCBparameters.boardSize(1);
Ncols=cameraCBparameters.boardSize(2);
icam=cameraCBparameters.icam;

reprojectedPointsMat=reshape(reprojectedPoints,Nrows-1,Ncols-1,2,Nimages);
reprojectedPointsHorizonal=reprojectedPointsMat([1 size(reprojectedPointsMat,1)],:,:,:);
reprojectedPointsVertical=reprojectedPointsMat(:,[1 size(reprojectedPointsMat,2)],:,:);

% PLOT
% define tabs
f = figure('name',['Reprojected points on all images taken by camera ' num2str(icam) ' After distortion correction. Scroll between tabs to view the different images'],'units','normalized','outerposition',[.1 .1 .8 .8]);;
tabgp = uitabgroup(f);

icount=0;
for iplot=cameraCBparameters.imagesUsed(imageIndecesToDisplay)
    icount=icount+1;
    %plot
    tab(icount) = uitab(tabgp,'Title',['IM' num2str(iplot)]);
    
    axes('Parent',tab(icount)); % somewhere to plot
    
    subplot(1,6,1:5)
    
    imshow(I(:,:,:,iplot)); hold all;
    plot(imagePoints(:,1,icount), imagePoints(:,2,icount),'go','linewidth',1.5);
    plot(reprojectedPoints(:,1,icount),reprojectedPoints(:,2,icount),'r+','linewidth',1.5);
    title(['Camera ' num2str(icam) ' Image ' num2str(iplot)]);
    drawnow
    % plot straight lines
    plot(squeeze(reprojectedPointsHorizonal(:,1,1,icount)),squeeze(reprojectedPointsHorizonal(:,1,2,icount)),'-c','linewidth',1.5);
    plot(squeeze(reprojectedPointsVertical(1,:,1,icount)),squeeze(reprojectedPointsVertical(1,:,2,icount)),'-m','linewidth',1.5);
    for icol=2:Ncols-1
        plot(squeeze(reprojectedPointsHorizonal(:,icol,1,icount)),squeeze(reprojectedPointsHorizonal(:,icol,2,icount)),'-c','linewidth',1.5);
    end
    for irow=1:Nrows-1
        plot(squeeze(reprojectedPointsVertical(irow,:,1,icount)),squeeze(reprojectedPointsVertical(irow,:,2,icount)),'-m','linewidth',1.5);
    end
    legend('Detected Points','Reprojected Points','straight horizontal lines','straight vertical lines');
    hold off;
    
    reprojectionErrorsNow=reprojectionErrors(:,:,icount);
    reprojectionErrorsMgnNow=sqrt(sum(reprojectionErrorsNow.^2,2));
    reprojectionErrorsNow=[reprojectionErrorsNow reprojectionErrorsMgnNow];
    
    subplot(1,6,6)
    boxplot(reprojectionErrorsNow,'Labels',{'X','Y','Mgn'});
    ylim([-max(abs(reprojectionErrorsNow(:))) max(abs(reprojectionErrorsNow(:)))]);
    title({'Reprojection error'; 'statistics [pix]'});
    
    drawnow
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
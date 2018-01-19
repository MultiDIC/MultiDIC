% for each camera, undistort the the calibration image based on the
% parameters found in STEP0
% then mask it and turn it into BW
% Then, find the 2D centroids on each image and save to Centroids_unsorted.mat
% Then, Sort the 2D positions of the centroids according to number of row and columns and save to 'Centroids_sorted.mat'
% Then, calculate the DLT calibration parameters based on the known 3D positions of the points

clearvars; close all

%% CHOOSE PATHS and OPTIONS

% initial image path
imagePathInitial='C:\Users\Dana\Dropbox (Personal)\MIT\research\DIC\360_DIC\LimbTest\data';
% select the folder containing the calibration images (if imagePathInitial=[] then the initial path is the current path)
imagePaths = uipickfiles('FilterSpec',imagePathInitial,'Prompt','Select one or multiple images of 3D calibration target for analysis');
% save camera parameters? choose save path and overwrite options
[saveDLTParametersLogic,savePath]=QsaveDLTParameters(imagePaths);
Ncam=numel(imagePaths);
indCams=zeros(1,Ncam);
for ic=1:Ncam % camera indeces from image paths
    [~,name,fileType] = fileparts(imagePaths{ic});
    indCams(ic)=str2num(name(end-1:end));
end
% save processed images? warn for overwriting
[saveProcessedImagesLogic,overWriteImages]=QsaveDLTImages;
% select the file containing the 3D world points of the calibration target
DefaultName='C:\Users\Dana\Dropbox (Personal)\MIT\research\DIC\360_DIC\LimbTest\data\test\CylinderCoor.mat';
[FileName,PathName,~] = uigetfile('','select the file containing the 3D world points of the calibration target',DefaultName);
C3DtrueStruct=load([PathName FileName]);
fieldNames=fieldnames(C3DtrueStruct);
C3Dtrue=getfield(C3DtrueStruct,fieldNames{1});
NcTotal=size(C3Dtrue,2);
Nr=size(C3Dtrue,1);
% figures path
figuresPath=[savePath '\figures'];
mkdir(figuresPath);
% remove distortion? If yes, choose the folder where the parameters exist
distortionRemovalButton = questdlg('Remove distortion from images?', 'Remove distortion?', 'Yes', 'No', 'Yes');
switch distortionRemovalButton
    case 'Yes'
        distortionRemovalLogic=true(1);
        distortionParametersPath=uigetdir([],'Select a folder where the cameraCBparameters are located');
    case 'No'
        distortionRemovalLogic=false(1);
end
% select camera pairs
for ic=1:Ncam-1
    indPairsInitial(ic,:)=[indCams(ic) indCams(ic+1)];
end
answer = inputdlg({['Enter the pairs of cameras [' num2str(indCams) '] as semicolon seperated pairs:']},'Input',[1,50],{mat2str(indPairsInitial)});
indPairs=str2num(answer{1});
Npairs=size(indPairs,1);

%% load images, undistort and mask them, and find centroids

colFirst=zeros(Ncam,1); colLast=zeros(Ncam,1); NcCam=zeros(Ncam,1);
Centroids_sorted=cell(1,Ncam);

icount=0;
for icam=indCams
    icount=icount+1;
    
    IM = imread(imagePaths{icount});
    
    % rotate image
%     IM = imrotate(IM,imageRotationAngle);
    
    % turn to gray if RGB
    if size(IM,3)==3
        IMgr = rgb2gray(IM);
    else
        IMgr = IM;
    end
    
    % undistort image based on camera parameters obtained by checkerboard calibration
    if distortionRemovalLogic
        CBparams=load([distortionParametersPath '\cameraCBparameters_cam' num2str(icam,'%02i')]);
        CBparams=CBparams.cameraCBparameters;
        [IMgrUD,~] = undistortImage(IMgr,CBparams.cameraParameters);
    else
        IMgrUD=IMgr;
    end
    
    if saveProcessedImagesLogic
        imName=[savePath '\IMgrUD_' num2str(icam,'%02i') fileType];
        if exist(imName,'file') && ~overWriteImages
                waitfor(warndlg({'Processed image '; imName ;' will be overwritten'}));
        end
        imwrite(IMgrUD,imName);
    end
    
    %plot original and gray undistorted
    fh=figure;
    fh.Units='normalized'; fh.Position=[.02 .1 .9 .8];
    subplot(1,2,1)
    imshow(IM);
    title('Original Image');
    subplot(1,2,2)
    imshow(IMgrUD);
    if distortionRemovalLogic
        title('Grayscale undistorted Image');
    else
        title('Grayscale Image');
    end
    suptitle({['Camera ' num2str(icam)]; 'Click on the figure to continue'});
    waitforbuttonpress
    
    % increase image resolution for better centroid finding
    resFactor=3;
    IMgrUDhr = imresize(IMgrUD,resFactor);
    
    centroidsOK='No'; % check if cerntroids are correct
    while ~strcmp(centroidsOK,'Yes')
        % crop polygon and turn to white everything out of polygon
        polyOK='No';
        while ~strcmp(polyOK,'Yes')
            figure;
            title({'Draw polygon around ROI'; 'Polygon vertices can be dragged after polygon is closed, Double click on the polygon to finish'}); hold on
            [PG,xi,yi] = roipoly(IMgrUDhr);
            % turn to white everything outside the polygon
            IMgrcr=IMgrUDhr;
            IMgrcr(~PG)=255;
            % plot masked image
            figure;
            warning('off','images:initSize:adjustingMag');
            imshow(zeros(size(IMgrUDhr))); hold on;
            h=imshow(IMgrUDhr); hold on; drawnow
            maskAlpha=double(PG);
            maskAlpha(PG==0)=.4;
            h.AlphaData=maskAlpha;
            title('Masked image. Click on the figure to continue');
            axis off
            drawnow
            waitforbuttonpress
            % query if mask is correct
            polyOK = questdlg('polygon mask OK? (if you select "No", this GUI will start over)', 'polygon mask OK?', 'Yes', 'No', 'Yes');
        end
        
        answer = inputdlg({'First column number:';'Last column number:'},'Input',[1,30]);
        colFirst(icount)=str2num(answer{1});
        colLast(icount)=str2num(answer{2});

        % calculate number of columns for this camera
        if colLast(icount)>colFirst(icount)
            NcCam(icount)=colLast(icount)-colFirst(icount)+1;
        else
            NcCam(icount)=colLast(icount)-colFirst(icount)+NcTotal+1;
        end
        
        % number of points on this image
        Np=NcCam(icount)*Nr;
                  
        [C] = plot_centroids_threshold(IMgrcr,Np);
        C=C/resFactor;
        
        centroidsOK=questdlg('centroids OK? (if you select "No", it will start over)', 'centroids OK?', 'Yes', 'No', 'Yes');
    end
    
    % Sort the  centroids     
    Ptemp=C;
    NpNew=Np; % New number of points (1 column is deleted every loop)
    PS=zeros(size(C)); % matrix for sorted points 
    for irow=1:Nr
        % select the Nc points with largest y (lowest in the image, corresponds to smallest z)
        [PtempYsort,PtempYsortInd]=sort(Ptemp(:,2),1,'descend');
        PtempNow=Ptemp(PtempYsortInd(1:NcCam(icount)),:);
        
        % sort them by x (left to right)
        [PtempNowXSsort,PtempNowXSsortInd]=sort(PtempNow(:,1));
        PStemp=PtempNow(PtempNowXSsortInd,:);     
        PS(irow+((1:NcCam(icount))-1)*Nr,:)=PStemp; % save this column into PS
        Ptemp(PtempYsortInd(1:NcCam(icount)),:)=[]; % delete these points from p3
        NpNew=NpNew-NcCam(icount); % update number of remaining points
                
    end
    
    close all
    
    % plot sorted centroids
    fh=figure; hold all
    imshow(IMgrUD); hold on
    fh.Units='normalized'; fh.Position=[.02 .1 .9 .8];
    title(['Sorted points, camera ' num2str(icam) '. Click on the figure to continue']);
    plot(PS(:,1),PS(:,2),'+y','markersize',3);
    pause(.5)
    plot(PS(:,1),PS(:,2),'b','linewidth',1);
    pause(.5)
    numCams=arrayfun(@num2str,1:Np,'unif',0);
    text(PS(:,1),PS(:,2),numCams,'color','r','fontsize',12);
    pause(.5)
    gap=.2;
    xRange=max(PS(:,1))-min(PS(:,1));
    yRange=max(PS(:,2))-min(PS(:,2));
    xlim([min(PS(:,1))-gap*xRange max(PS(:,1))+gap*xRange]);
    ylim([min(PS(:,2))-gap*yRange max(PS(:,2))+gap*yRange]);
    waitforbuttonpress
    
     % save centroids in cell
    Centroids_sorted{icount}=PS;
    
end

%% Direct linear transformation from the 2D image to 3D target points (save DLT transformation matrix L for each camera)

C3D=cell(Ncam,1); % 3D world true points visible by each camera
DLTparAllCams=cell(Ncam,2); % cell array to store the DLT parameters first column is camera number and second column is the DLT parameters

icount=0;
for icam=indCams
    icount=icount+1;
    
    Np=NcCam(icount)*Nr;
    
    % 2D coordinates from image
    C2D=Centroids_sorted{icount};
    
    % 3D true points coordinates 
    % column indeces
    if colLast(icount)>colFirst(icount)
        cInds=colFirst(icount):colLast(icount);
    else
        cInds=[colFirst(icount):NcTotal 1:colLast(icount)];
    end    
    C3Dtemp=C3Dtrue(:,cInds,:);
    C3D{icount}=reshape(C3Dtemp,Np,3);
    
    % Solve the linear system with the LS method for the unknown 11 DLT parameters that reflect the relationships between the world 3D and the image 2D (vector L)
    L=DLT11Calibration(C2D,C3D{icount});
    
    % save L in cell array
    DLTparAllCams{icount}=L;
    
    if saveDLTParametersLogic
        save([savePath '\DLTpar_cam' num2str(icam,'%02i')],'L');
    end    

end

DLTstruct=struct;
DLTstruct.indCams=indCams;
DLTstruct.indPairs=indPairs;
DLTstruct.DLTparams=DLTparAllCams;
DLTstruct.columns=[colFirst colLast];
DLTstruct.imageCentroids=Centroids_sorted;

%% 3D recronstruction of overlapping points for image pairs 
% calculating reprojected points and reprojection errors

DLTreprojectPoints=cell(size(DLTstruct.indPairs,1),1);
DLTreprojectErr=cell(size(DLTstruct.indPairs,1),1);
for ipair=1:size(DLTstruct.indPairs,1)
    
    icamL=DLTstruct.indPairs(ipair,1);
    icountL=find(DLTstruct.indCams==icamL);
    icamR=DLTstruct.indPairs(ipair,2);
    icountR=find(DLTstruct.indCams==icamR);
    
    % overlapping columns for current pair
    cFirst=DLTstruct.columns(icountL,1);
    cLast=DLTstruct.columns(icountR,2);
    if cLast>cFirst
        NpPair=(cLast-cFirst+1)*Nr;
    else
        NpPair=(cLast-cFirst+NcTotal+1)*Nr;
    end
    P2DL=DLTstruct.imageCentroids{icountL}(1:NpPair,:);
    P2DR=DLTstruct.imageCentroids{icountR}(size(DLTstruct.imageCentroids{icountR},1)-NpPair+1:end,:);
 
    LR=DLTstruct.DLTparams{icountR}; % DLT parameters for 1st camera
    LL=DLTstruct.DLTparams{icountL}; % DLT parameters for 2nd camera

    % reprojection errors between true points and reconstructed points
    % Solve the linear system with the LS method for the unknown X Y and Z
    P3D=DLT11Reconstruction(P2DR,P2DL,LR,LL);
    % save into cell array
    DLTreprojectPoints{ipair}=P3D;
    
    % find the true points for this pair
    if cLast>cFirst
        cInds=cFirst:cLast;
    else
        cInds=[cFirst:NcTotal 1:cLast];
    end 
    P3Dtrue=reshape(C3Dtrue(:,cInds,:),NpPair,3);

    %reprojection error vectors
    DLTreprojectErr{ipair}=P3D-P3Dtrue;
    
end
DLTstruct.reprojectPoints=DLTreprojectPoints;
DLTstruct.reprojectErr=DLTreprojectErr;

%% save DLTstruct
if saveDLTParametersLogic
    save([savePath '\DLTstruct'],'DLTstruct');
end

%% Plot reprojected points and reprojection errors
% this section can be executed also by loading an existing DLTstruct

Colors=gjet(Ncam+1);

% plot reprojected points  
hf1=cFigure;
hf1.Name='Reprojected points vs. true points';
legendstrings=cell(Npairs,1);
for ipair=1:Npairs    
    plotV(DLTstruct.reprojectPoints{ipair},'x','Color',Colors(ipair,:)); hold on
    legendstrings{ipair}=['Pair ' num2str(ipair)];
end
P3DtrueArray=reshape(C3Dtrue,NcTotal*Nr,3);
plotV(P3DtrueArray,'s','Color','k'); hold on
legendstrings{ipair+1}='True points';
legend(legendstrings);
axisGeom
title('reprojected points [mm]');

if saveDLTParametersLogic
    savefig([figuresPath '\reprojectPoints']);
end

% plot reprojected errors and their statisticss
hf2=cFigure;
reprojectErrArray=cell2mat(DLTstruct.reprojectErr(:));
hf2.Name='Reprojection errors';
errMaxTotal=0;
ax1=subplot(1,2,1);
for ipair=1:Npairs
    NpPair=size(DLTstruct.reprojectErr{ipair},1);
    quiver3(zeros(NpPair,1),zeros(NpPair,1),zeros(NpPair,1),DLTstruct.reprojectErr{ipair}(:,1),DLTstruct.reprojectErr{ipair}(:,2),DLTstruct.reprojectErr{ipair}(:,3),'Color',Colors(ipair,:),'AutoScale','off'); hold on
    errMax=max(abs(DLTstruct.reprojectErr{ipair}(:)));
    errMaxTotal=max([errMaxTotal errMax]); 
    legendstrings{ipair}=['Pair ' num2str(ipair)];
end
legendstrings{ipair+1}='mean';
quiver3(0,0,0,mean(reprojectErrArray(:,1)),mean(reprojectErrArray(:,2)),mean(reprojectErrArray(:,3)),'Color','k','AutoScale','off','linewidth',4); hold on
legend(legendstrings);
axisGeom
xlim([-errMaxTotal errMaxTotal]); ylim([-errMaxTotal errMaxTotal]); zlim([-errMaxTotal errMaxTotal]);
title('Reprojection errors per pair [mm]');

ax2=subplot(1,2,2);
hb=boxplot(reprojectErrArray,{'x','y','z'});
subplot(ax2)
ylim([-errMaxTotal errMaxTotal]);
hrf=refline(0,0); hrf.Color='g'; hrf.LineStyle='--';
title('Reprojection error statistics for all pairs [mm]');

if saveDLTParametersLogic
    savefig([figuresPath '\reprojectErrors']);
end

msgbox('STEP1 is finished');

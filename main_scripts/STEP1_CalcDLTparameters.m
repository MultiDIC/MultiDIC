%% STEP 1: Calculate DLT parameters (Stereo Calibration)
% This main scripts calculates the DLT parameters for each camera based on images of the calibration object and the known positins of the
% control points on the calibration object. This script is designed to work with a cylindrical object.
% For each camera, the image is maskef and turned it into BW to find the 2D centroids
% Then, the centroids are sorted according to number of row and columns
% Then, the DLT calibration parameters are calculated and saved

%%
clearvars; close all

fs=get(0, 'DefaultUIControlFontSize');
set(0, 'DefaultUIControlFontSize', 10);

%% CHOOSE PATHS and OPTIONS

% initial image path. Edit this path if you want the code to open a different folder than the current path.
imagePathInitial=pwd;

% select the calibration images (if imagePathInitial=[] then the initial path is the current path)
% the calibration images must end with the camera number after '_', or be just the camera number
imagePaths = uipickfiles('FilterSpec',imagePathInitial,'Prompt','Select stereo calibration images from one camera or multiple cameras');
if ~iscell(imagePaths)
    return
end

% save camera parameters? choose save path and overwrite options
[saveDLTParametersLogic,savePath]=QsaveDLTParameters(imagePaths);
Ncam=numel(imagePaths);
indCams=zeros(1,Ncam);
for ic=1:Ncam % camera indeces from image paths
    [~,name,fileType] = fileparts(imagePaths{ic});
    indCamCell=strsplit(name,'_');
    indCams(ic)=str2num(indCamCell{end});
end

% select the MAT file containing the 3D world points of the calibration target
calibObjectPathInitial=[];
[FileName,PathName,~] = uigetfile('','Select the file containing the 3D world points of the calibration target',calibObjectPathInitial);
if ~ischar(FileName)
    return
end
C3DtrueStruct=load([PathName FileName]);
fieldNames=fieldnames(C3DtrueStruct);
% true 3d coordinates of the calibration target (nRows x nColumns x xyz)
C3Dtrue=getfield(C3DtrueStruct,fieldNames{1});
NcTotal=size(C3Dtrue,2); % total number of columns
Nr=size(C3Dtrue,1); % number of rows

%% load images, mask them, and find centroids

colFirst=zeros(Ncam,1); colLast=zeros(Ncam,1); NcCam=zeros(Ncam,1);
Centroids_sorted=cell(1,Ncam);

icount=0;
for icam=indCams % loop over all cameras
    icount=icount+1;
    
    % load image
    IM = imread(imagePaths{icount});
        
    % turn to gray if RGB
    if size(IM,3)==3
        IMgr = rgb2gray(IM);
    else
        IMgr = IM;
    end
    
%     %plot original and gray 
%     fh=figure;
%     fh.Units='normalized'; fh.Position=[.05 .05 .9 .8];
%     subplot(1,2,1)
%     imshow(IM);
%     title('Original Image');
%     subplot(1,2,2)
%     imshow(IMgr);
%     title('Grayscale Image'); 
%     suptitle({['Camera ' num2str(icam)]; 'Click on the figure to continue'});
%     waitforbuttonpress
    
    % temporarily increase image resolution for better centroid finding
    resFactor=3;
    IMgrHR = imresize(IMgr,resFactor); % the default is cubic interpolation
    
    centroidsOK='No'; % check if cerntroids are correct
    while ~strcmp(centroidsOK,'Yes')
        % crop polygon and turn to white everything out of polygon
        polyOK='No';
        
        while ~strcmp(polyOK,'Yes') % loop as long as polygon is not OK
            % plot the calibration object image and let user draw polygon mask
            fh=figure;
            fh.Name=['Camera ' num2str(icam)];
            fh.Units='normalized'; fh.Position=[.05 .05 .9 .8];
            title({'Draw polygon around ROI'; 'Polygon vertices can be dragged after polygon is closed, Double click on the polygon to finish'}); hold on
            [PG,xi,yi] = roipolyFit(IMgrHR);
 
            % turn to white everything outside the polygon
            IMgrcr=IMgrHR;
            IMgrcr(~PG)=255;
            
            % plot masked image
            warning('off','images:initSize:adjustingMag');
            imshow(zeros(size(IMgrHR)), 'InitialMagnification', 'fit'); hold on;
            h=imshow(IMgrHR, 'InitialMagnification', 'fit'); hold on; drawnow
            maskAlpha=double(PG);
            maskAlpha(PG==0)=.4;
            h.AlphaData=maskAlpha;
            title('Masked image. Click on the figure to continue');
            axis off
            drawnow
            waitforbuttonpress
            
            % query if mask is correct
            polyOK = questdlg('Is the polygon mask correct? (if you select "No", this GUI will start over)', 'polygon mask OK?', 'Yes', 'No', 'Yes');
        end
        
        % ask for column numbers
        options.Resize = 'on';
        answer = inputdlg({'First column number:';'Last column number:'},'Input',[1 40],{'',''},options);
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
        
        % plot the centroids with the calculated threshold and let the user change it if necessary
        [C] = plot_centroids_threshold(IMgrcr,Np);
        % scale back to real size
        C=C/resFactor;
        
        %chexk if centroids are ok
        centroidsOK=questdlg('centroids OK? (if you select "No", it will start over)', 'centroids OK?', 'Yes', 'No', 'Yes');
    end
    
    % Sort the  centroids     
    Ptemp=C;
    NpNew=Np; % New number of points (1 column is deleted every loop)
    PS=zeros(size(C)); % matrix for sorted points 
    for irow=1:Nr          % for each row
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
    
    % plot the sorted centroids
    fh=figure; hold all
    imshow(IMgr); hold on
    fh.Units='normalized'; fh.Position=[.05 .05 .9 .8];
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
close(fh);

%% calculate Direct linear transformation from the 2D image and 3D target points (save DLT parameter matrix L for each camera)

C3D=cell(Ncam,1); % 3D world true points visible by each camera

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
    
    % Solve the linear system for the unknown 11 DLT parameters that reflect the relationships between the world 3D and the image 2D 
    L=DLT11Calibration(C2D,C3D{icount});
       
    % create DLTstructCam structure withh all the data
    DLTstructCam=struct;
    DLTstructCam.indCam=icam;
    DLTstructCam.columns=cInds;
    DLTstructCam.imageCentroids=C2D;
    DLTstructCam.DLTparams=L;
    DLTstructCam.C3Dtrue=C3Dtrue;

    % save the DLT parameters structure for each camera
    if saveDLTParametersLogic
        save([savePath '\DLTstruct_cam_' num2str(icam)],'DLTstructCam');
    end    

end

%% finish

h = msgbox('STEP1 is completed');
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
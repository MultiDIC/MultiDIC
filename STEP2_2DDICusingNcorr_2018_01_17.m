%% STEP 2 Run Ncorr analysis on sets of images from a pair of 2 cameras
% The complete set of 2n images includes 2 sets of images taken simultaneously
% from 2 cameras (2 views). The first n images are from the "reference"
% camera and the last n are the "deformed" camera. It's not really
% important which one is defined as Ref and Def, as long as it's
% consistant. The 1st image from the 1st camera is defined as the reference
% image.

clearvars; close all

%% CHOOSE PATHS OPTIONS
% initial image path
folderPathInitial='C:\Users\Dana\Dropbox (Personal)\MIT\research\DIC\360_DIC\LimbTest\data\HH_2017_10_02\HH_2017_10_02_t01\dynamic\testHH1';
% select the folder containing the analysis images (if imagePathInitial=[] then the initial path is the current path)
folderPathRef=uigetdir(folderPathInitial,'Select the folders containing images from the reference camera');
folderPathDef=uigetdir(folderPathInitial,'Select the folders containing images from the second camera');
folderPaths=cell(1,2);
folderPaths{1}=folderPathRef;
folderPaths{2}=folderPathDef;
% camera indeces for current analysis    
nCamRef=str2double(folderPaths{1}(end-1:end));
nCamDef=str2double(folderPaths{2}(end-1:end));
% save 2D-DIC results? choose save path and overwrite options
[save2DDIClogic,savePath]=Qsave2DDICresults(folderPaths);
% remove distortion? If yes, choose the folder where the parameters exist
distortionRemovalButton = questdlg('Remove distortion from images?', 'Remove distortion?', 'Yes', 'No', 'Yes');
switch distortionRemovalButton
    case 'Yes'
        distortionRemovalLogic=true(1);
        distortionParametersPath=uigetdir([],'Select a folder where the cameraCBparameters are located');
    case 'No'
        distortionRemovalLogic=false(1);
        distortionParametersPath=[];
end
%% create structure for saving the 2DDIC results
DIC2DpairResults = struct;

DIC2DpairResults.nCamRef=nCamRef;
DIC2DpairResults.nCamDef=nCamDef;

%%  load images from the paths, convert to gray and undistort, and create IMset cell for Ncorr
[ImPaths,ImSet]=createImageSet(folderPaths,distortionRemovalLogic,distortionParametersPath);
DIC2DpairResults.nImages=numel(ImPaths)/2;
DIC2DpairResults.ImPaths=ImPaths;

%% animate the 2 sets of images to be correlated with Ncorr
anim8_DIC_images(ImPaths,ImSet);
pause
%% choose ROI
% This is a GUI for choosing the ROI instead of choosing the ROI in the
% NCorr softwhere (too small). It also allows the assistance of SIFT
% matches (it helps locating the overlapping region, but is time costly)

chooseMaskButton = questdlg('Create new mask for correlation, use saved mask, or use Ncorr to draw mask?', 'mask options?', 'New', 'Saved','Ncorr', 'New'); % existing mask should be in savePath
switch chooseMaskButton
    case 'New'
        nROI=1;
        logSIFT=0;
        SIFTthresh=1.5;
        ROImask = selectROI(ImPaths,nROI,logSIFT,SIFTthresh);
        % save image mask
        % The format is ROIMask_C01_C02, where 01 is the reference camera of the pair, and 02 is the "deformed" camera of the pair.
        save([savePath '\ROIMask' '_C' num2str(nCamRef,'%02u') '_C' num2str(nCamDef,'%02u')],'ROImask');
        DIC2DpairResults.ROImask=ROImask;
    case 'Saved'
        load([savePath '\ROIMask' '_C' num2str(nCamRef,'%02u') '_C' num2str(nCamDef,'%02u')]);
        DIC2DpairResults.ROImask=ROImask;
    case 'Ncorr'
end

msgbox('Wait while initializing Ncorr');

%% Start Ncorr 2D analysis
% open Ncorr
handles_ncorr = ncorr;
% set reference image
handles_ncorr.set_ref(ImSet{1});
% set current image
handles_ncorr.set_cur(ImSet);
% set ROI (skip this step if you want to select the ROI in Ncorr)
if ~strcmp(chooseMaskButton,'Ncorr')
    handles_ncorr.set_roi_ref(ROImask);
end

% Set analysis in Ncorr and wait
disp('press enter in the command window when Ncorr analysis is finished');
pause

%% Extract results from Ncorr and calculate correlated image points, correlation coefficients, faces and face colors
[Points,CorCoeffVec,F,CF] = extractNcorrResults(handles_ncorr,ImSet{1});
DIC2DpairResults.Points=Points;
DIC2DpairResults.CorCoeffVec=CorCoeffVec;
DIC2DpairResults.Faces=F;
DIC2DpairResults.FaceColors=CF;
if ~strcmp(chooseMaskButton,'Ncorr')
    DIC2DpairResults.ROImask=handles_ncorr.reference.roi.mask;
end

%% save important variables for further analysis (write text files of correlated 2D points, their cirrelation coefficients, triangular faces, and face colors
if save2DDIClogic
    save([savePath '\DIC2DpairResults_C' num2str(nCamRef,'%02u') '_C' num2str(nCamDef,'%02u')] ,'DIC2DpairResults');
end

%% plot?
plotButton = questdlg('Plot correlated points on images?', 'Plot?', 'Yes', 'No', 'Yes');
switch plotButton
    case 'Yes'
        plotNcorrPairResults(DIC2DpairResults);
    case 'No'
end

%% close Ncorr figure
close(handles_ncorr.handles_gui.figure);

msgbox('STEP2 is finished');

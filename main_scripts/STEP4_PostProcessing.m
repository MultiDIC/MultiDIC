%% STEP 4: Post Processing
% Upload 3D reconstructed points (either as seperate pairs or a stitched surface)
% and calculate diplacements, deformations, and strains.

%%
 clearvars; close all

fs=get(0, 'DefaultUIControlFontSize');
set(0, 'DefaultUIControlFontSize', 10);

%% CHOOSE PATHS OPTIONS

% select DIC3DpairResults structures
PathInitial=pwd;
[file,path] = uigetfile(PathInitial,'Select a 3D-DIC results structure (stitched or unstitched)');

DIC3D=load([path file]);
DIC3Dname=fieldnames(DIC3D);
DIC3D=DIC3D.(DIC3Dname{1});

% save 3D-DIC post processing results? choose save path and overwrite options
[save3DDIClogic,savePath]=Qsave3DDICPPresults(path);


%% 3D reconstruction using Direct Linear Transformation

nImages= numel(DIC3D.Points3D);

% pre-allocate 3D-DIC result variables
DIC3D.Points3D_ARBM=cell(1,nImages);
DIC3D.Disp.DispVec=cell(1,nImages);
DIC3D.Disp.DispMgn=cell(1,nImages);
DIC3D.Disp.DispVec_ARBM=cell(1,nImages);
DIC3D.Disp.DispMgn_ARBM=cell(1,nImages);
DIC3D.FaceCentroids=cell(1,nImages);
DIC3D.FaceCentroids_ARBM=cell(1,nImages);
DIC3D.FaceCorrComb=cell(1,nImages);
DIC3D.FaceIsoInd=cell(1,nImages);
DIC3D.RBM.RotMat=cell(1,nImages);
DIC3D.RBM.TransVec=cell(1,nImages);

F=DIC3D.Faces;

hw = waitbar(0,'Calculating displacements and rigid body motion');
for ii=1:nImages % loop over images (time frames)
    waitbar(ii/(nImages));
    
    % Face correlation coefficient (worst)
    DIC3D.FaceCorrComb{ii}=max(DIC3D.corrComb{ii}(F),[],2);
    
    % compute face centroids
    for iface=1:size(F,1)
        DIC3D.FaceCentroids{ii}(iface,:)=mean(DIC3D.Points3D{ii}(F(iface,:),:));
    end
    
    % Compute displacements between frames (per point)
    DispVec=DIC3D.Points3D{ii}-DIC3D.Points3D{1};
    DIC3D.Disp.DispVec{ii}=DispVec;
    DIC3D.Disp.DispMgn{ii}=sqrt(DispVec(:,1).^2+DispVec(:,2).^2+DispVec(:,3).^2);
    
    % Compute rigid body transformation between point clouds
    [RotMat,TransVec,Points3D_ARBM]=rigidTransformation(DIC3D.Points3D{ii},DIC3D.Points3D{1});
    DIC3D.RBM.RotMat{ii}=RotMat;
    DIC3D.RBM.TransVec{ii}=TransVec;
    DIC3D.Points3D_ARBM{ii}=Points3D_ARBM;
    
    % Compute displacements between sets - after RBM
    DispVec=DIC3D.Points3D_ARBM{ii}-DIC3D.Points3D_ARBM{1};
    DIC3D.Disp.DispVec_ARBM{ii}=DispVec;
    DIC3D.Disp.DispMgn_ARBM{ii}=sqrt(DispVec(:,1).^2+DispVec(:,2).^2+DispVec(:,3).^2);
    
    % compute face centroids - after transformation
    for iface=1:size(F,1)
        DIC3D.FaceCentroids_ARBM{ii}(iface,:)=mean(Points3D_ARBM(F(iface,:),:));
    end
    
end
delete(hw);

% compute deformation and strains (per triangular face)
deformationStruct=triSurfaceDeformation(F,DIC3D.Points3D{1},DIC3D.Points3D);
DIC3D.Deform=deformationStruct;

deformationStruct_ARBM=triSurfaceDeformation(F,DIC3D.Points3D_ARBM{1},DIC3D.Points3D_ARBM);
DIC3D.Deform_ARBM=deformationStruct_ARBM;

% compute triangle regularity (isotropy index)
for ii=1:nImages
    [FisoInd]=faceIsotropyIndex(F,DIC3D.Points3D{ii});
    DIC3D.FaceIsoInd{ii}=FisoInd;
end

DIC3DPPresults=DIC3D;


%% save results
nPairs=size(DIC3DPPresults.pairIndices,1);
if save3DDIClogic
    saveName=fullfile(savePath, ['DIC3DPPresults_' num2str(nPairs) 'Pairs.mat']);
    icount=1;
    while exist(saveName,'file')
        saveName=fullfile(savePath, ['DIC3DPPresults_' num2str(nPairs) 'Pairs(' num2str(icount) ').mat']);
        icount=icount+1;
    end
    save(saveName,'DIC3DPPresults','-v7.3');
end

%% Plot results?

plotButton = questdlg('Plot 3D-DIC post-processing results?', 'Plot?', 'Yes', 'No', 'Yes');
switch plotButton
    case 'Yes'
        plotMoreLogic=true;
        while plotMoreLogic
            optStruct=struct;
            optStruct.zDirection=1;
            optStruct.FaceAlpha=1;
            optStruct.smoothLogic=1;
            % PLOT
            plot3DDICPPresults(DIC3DPPresults,optStruct);
            
            plotMoreButton = questdlg('Plot more results?', 'Plot?', 'Yes', 'No', 'Yes');
            switch plotMoreButton
                case 'Yes'
                    plotMoreLogic=true;
                case 'No'
                    plotMoreLogic=false;
            end
        end
    case 'No'
end

%% finish
h=msgbox('STEP4 is completed');
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

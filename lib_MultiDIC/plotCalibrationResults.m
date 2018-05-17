function []=plotCalibrationResults(varargin)
%% function for plot calibration (reprojection) results in step 1p
% INPUT options:
% plotCalibrationResults
% plotCalibrationResults(DLTstruct)
% plotCalibrationResults(DLTstruct,saveFiguresLogic,figuresPath)
% 

%%
nargin=numel(varargin);
switch nargin
    case 0
        [FileName,PathName,~]=uigetfile('','Select DLTstruct file');
        DLTstruct=load([PathName FileName]);
        DLTstruct=DLTstruct.DLTstruct;
        saveButton = questdlg('Save calibration restults figures?', 'Save?', 'Yes', 'No', 'Yes');
        switch saveButton
            case 'Yes'
                saveFiguresLogic=true(1);
                figuresPath = uigetdir('','Select a folder for saving the figures');
            case 'No'
                saveFiguresLogic=false(1);
                figuresPath=[];
        end
    case 1
        DLTstruct=varargin{1};
        saveButton = questdlg('Save calibration restults figures?', 'Save?', 'Yes', 'No', 'Yes');
        switch saveButton
            case 'Yes'
                saveFiguresLogic=true(1);
                figuresPath = uigetdir('','Select a folder for saving the figures');
            case 'No'
                saveFiguresLogic=false(1);
                figuresPath=[];
        end
    case 3
        DLTstruct=varargin{1};
        saveFiguresLogic=varargin{2};
        figuresPath=varargin{3};
    otherwise
        error('wrong number of input arguments');
end

%%
Npairs=size(DLTstruct.indPairs,1);
Colors=gjet(Npairs+1);

%% plot reconstructed points
hf1=cFigure;
hf1.Units='normalized'; hf1.Position=[.05 .05 .9 .8];

hf1.Name='Reconstructed points vs. true points';
legendstrings=cell(Npairs,1);
for ipair=1:Npairs
    plotV(DLTstruct.reconstructPoints{ipair},'x','Color',Colors(ipair,:)); hold on
    legendstrings{ipair}=['Pair ' num2str(ipair) ' (cameras ' num2str(DLTstruct.indPairs(ipair,1)) '+' num2str(DLTstruct.indPairs(ipair,2)) ')'];
end
NcTotal=size(DLTstruct.truePoints,2);
Nr=size(DLTstruct.truePoints,1);
P3DtrueArray=reshape(DLTstruct.truePoints,NcTotal*Nr,3);
plotV(P3DtrueArray,'s','Color','k'); hold on
legendstrings{ipair+1}='True points';
legend(legendstrings);
axisGeom
title('reconstructed points [mm]');

%% save figure 1
if saveFiguresLogic
        saveName=[figuresPath '\reconstructedPoints.fig'];
    icount=1;
    while exist(saveName,'file')
        saveName=[figuresPath '\reconstructedPoints(' num2str(icount) ').fig'];
        icount=icount+1;
    end
    savefig(saveName);
end

%% plot reconstruction errors and their statisticss
hf2=cFigure;
hf2.Units='normalized'; hf2.Position=[.05 .05 .9 .8];

reconstructErrArray=cell2mat(DLTstruct.reconstructErrors(:));
hf2.Name='Reconstruction errors';
errMaxTotal=0;
ax1=subplot(1,3,[1 2]);
for ipair=1:Npairs
    NpPair=size(DLTstruct.reconstructErrors{ipair},1);
    quiver3(zeros(NpPair,1),zeros(NpPair,1),zeros(NpPair,1),DLTstruct.reconstructErrors{ipair}(:,1),DLTstruct.reconstructErrors{ipair}(:,2),DLTstruct.reconstructErrors{ipair}(:,3),'Color',Colors(ipair,:),'AutoScale','off'); hold on
    errMax=max(abs(DLTstruct.reconstructErrors{ipair}(:)));
    errMaxTotal=max([errMaxTotal errMax]);
    RMSE(ipair)=rms(sqrt(sum(DLTstruct.reconstructErrors{ipair}.^2,2)));
    legendstrings{ipair}=['Pair ' num2str(ipair) ', RMSE = ' num2str(RMSE(ipair),3)];
end
quiver3(0,0,0,mean(reconstructErrArray(:,1)),mean(reconstructErrArray(:,2)),mean(reconstructErrArray(:,3)),'Color','k','AutoScale','off','linewidth',4); hold on
legendstrings{ipair+1}=['Mean, |Mean| = ' num2str(norm(mean(reconstructErrArray(:,:))),3)];
legend(legendstrings);
axisGeom
xlim([-errMaxTotal errMaxTotal]); ylim([-errMaxTotal errMaxTotal]); zlim([-errMaxTotal errMaxTotal]);
RMSEtot=rms(sqrt(sum(reconstructErrArray.^2,2)));
title({'Reconstruction errors per camera pair [mm]';['RMSE=' num2str(RMSEtot,3)]},'fontsize',16);

ax2=subplot(1,3,3);
reprojectErrArrayMgn=[reconstructErrArray sqrt(sum(reconstructErrArray.^2,2))];
hb=boxplot(reprojectErrArrayMgn,{'x','y','z','Mgn'});
set(gca,'FontSize',16);
subplot(ax2);
ylim([-max(max(abs(reprojectErrArrayMgn))) max(max(abs(reprojectErrArrayMgn)))]);
hrf1=refline(0,0); hrf1.Color='g'; hrf1.LineStyle='--';
hrf2=line([3.5 4.5],[RMSEtot RMSEtot]); hrf2.Color='m'; hrf2.LineStyle='--';
hl=line([3.5 3.5],ylim); hl.Color='k';
title('Reconstruction error statistics for all pairs [mm]','fontsize',16);
legend([hrf1 hrf2],{'Zero','RMSE(Mgn)'});

%% save figure 2
if saveFiguresLogic
        saveName=[figuresPath '\reconstructionErrors.fig'];
    icount=1;
    while exist(saveName,'file')
        saveName=[figuresPath '\reconstructionErrors(' num2str(icount) ').fig'];
        icount=icount+1;
    end
    savefig(saveName);
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
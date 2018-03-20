function []=plotNcorrPairResults(varargin)
%% function for plotting 2D-DIC results imported from Ncorr in step 2
% plotNcorrPairResults;
% plotNcorrPairResults(DIC2DpairResults);

%% 
switch nargin
    case 0
        % ask user to load ncorr results for current camera pair
        [FileName,PathName,~] = uigetfile('','Select a DIC2DpairResults file from a pair of cameras to visualize results');
        load([PathName FileName]);
    case 1
        % use given struct
        DIC2DpairResults=varargin{1};
end

%% select plot options
answer = inputdlg({'Enter maximum correlation coefficient in the colorbar (leave blank for max)',...
                   'Enter maximum correlation coefficient to keep point (leave blank for keeping all points)'},...
                   'Input',[1,50]);
               
CorCoeffDispMax=str2double(answer{1}); % maximal correlation coefficient for display in colorbar (use [] for default which is max)
CorCoeffCutOff=str2double(answer{2}); % maximal correlation coefficient to display point (use [] for default which is keep all points)

%% load images and animate
ImPaths=DIC2DpairResults.ImPaths;
ImSet=cell(length(ImPaths),1);
for ii=1:length(ImPaths)
    ImSet{ii}=imread(ImPaths{ii});
end

% anim8_DIC_images(ImPaths,ImSet);

%% plot correlated points - ANIMATION
% plot Ref image on left and all current on right
anim8_DIC_images_corr_points_1_2n(ImSet,DIC2DpairResults,CorCoeffCutOff,CorCoeffDispMax);
% plot ref camera on left and Def camera on right
anim8_DIC_images_corr_points_n_n(ImSet,DIC2DpairResults,CorCoeffCutOff,CorCoeffDispMax);
% plot correlated points with colors on the faces
% plot ref camera on left and Def camera on right
anim8_DIC_images_corr_faces_n_n(ImSet,DIC2DpairResults,CorCoeffCutOff,CorCoeffDispMax);

end
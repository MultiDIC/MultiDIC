function []=undistortImagesSavePlot(CBimagesInfo,cameraCBparameters,saveUndistortedImagesLogic,overWriteUDimagesLogic,savePath)
%% function for undistorting images and save them if required, in STEP0
%
% INPUTS:
% * CBimagesInfo
% * cameraCBparameters:  a structure containing all the calibration parameters created in STEP0
% * saveUndistortedImagesLogic: save (true) or not (false)
% * overWriteUDimagesLogic: overwrite existing undistorted images (true) or not (false)
% * savePath: path for saving the undistorted images

%%
% undistort images and save if required
CBimagesInfoJ=CBimagesInfo;
J=zeros(size(CBimagesInfoJ.I),'uint8');
for ip=1:CBimagesInfoJ.Nimages
    [J(:,:,:,ip),~] = undistortImage(CBimagesInfo.I(:,:,:,ip),cameraCBparameters.cameraParameters,'OutputView','same');
end
if saveUndistortedImagesLogic
    warning('off','MATLAB:MKDIR:DirectoryExists');
    mkdir(fullfile(savePath, 'undistortedImages', num2str(CBimagesInfo.icam)));
    for ip=1:CBimagesInfo.Nimages
        [~,name,~]=fileparts(CBimagesInfo.imageFileNames{ip});
        % save undistorted image as J_001_originalNAme.ext
        imName=fullfile(savePath, 'undistortedImages', num2str(CBimagesInfo.icam), ['J_' num2str(ip,'%03i') '_' name CBimagesInfo.imageType]);;
        % warn about overwriting if already exists
        if exist(imName,'file')
            if overWriteUDimagesLogic
                if strcmp(CBimagesInfo.imageType,'.png')
                    imwrite(J(:,:,:,ip),imName,'Compression','none');
                else
                    imwrite(J(:,:,:,ip),imName,'Quality',100);
                end
            else
                waitfor(warndlg({'Undistorted image'; name ;' already exist so it will not be overwritten'}));
            end
        else
                if strcmp(CBimagesInfo.imageType,'.png')
                    imwrite(J(:,:,:,ip),imName,'Compression','none');
                else
                    imwrite(J(:,:,:,ip),imName,'Quality',100);
                end
        end
    end
end
CBimagesInfoJ.I=J;

% plot all undistorted images with reprojected points after correction
plot_reprojectVSreal_points_AUD(CBimagesInfoJ,cameraCBparameters);

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
function [saveLogic,savePath]=QsaveCameraParameters(image_folder_paths)
%% function for asking the user for saving options for the camera parmaters in STEP0.
% Including overwriting options
%
% INPUTS:
% * image_folder_paths: 1-by-NumberOfCameras cell array, containing in each
%   cell the path to the folder contatining the checkerboard images. The
%   camera index is retrieved from the folder name (_Ncam)
%
% OUTPUTS:
% * saveLogic: logic true=save, false=not save
% * savePath: a path were to save the results. [] if saveLogic=false.

%%

UIControlFontSize = get(0, 'DefaultUIControlFontSize');
set(0, 'DefaultUIControlFontSize', 10);
saveButton = questdlg('Save camera parameters?', 'Save?', 'Yes', 'No', 'Yes');

switch saveButton
    case 'Yes'
        saveLogic=true(1);
        continueLog=false(1);
        
        while ~continueLog
            % Path where to save the camera parameters
            initalPath=fileparts(image_folder_paths{1});
            savePath = uigetdir(initalPath,'Select a folder for saving the results');
            
            % check if camera parameters for these cameras already exist 
            filenames=cell(1,numel(image_folder_paths));
            fileExistAll=false(1,numel(image_folder_paths));
            for ic=1:numel(image_folder_paths)
                icam=str2num(image_folder_paths{ic}(end-1:end));
                filenames{ic}=[savePath '\cameraCBparameters_cam' num2str(icam,'%02i') '.mat'];
                if exist(filenames{ic},'file')
                    fileExistAll(1,ic)=true(1);
                end
            end
            
            % if camera parameters already exist, prompt and ask if to
            % overwrite or not. If not, ask for an alternative folder.
            if any(fileExistAll)
                fileExistTrue=find(fileExistAll);
                
                promptMessage=cell(numel(fileExistTrue)+2,1);
                promptMessage{1}='The following output files already exists: ';
                promptMessage{numel(fileExistTrue)+4}='Do you want to overwrite?';
                for ii=1:numel(fileExistTrue)
                    [~,fileName,ext]=fileparts(filenames{fileExistTrue(ii)});
                    promptMessage{ii+2}=[fileName ext];
                end
                
                titleBarCaption = 'Overwrite?';
                buttonText = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
                switch buttonText
                    case 'Yes'
                        % User wants to overwrite. % Set flag to continue
                        continueLog=true(1);
                    case 'No'
                        % User does not want to overwrite. % Set flag to not do the write.
                        waitfor(warndlg('Choose another folder to save the results','OK'));
                        % let user choose a different folder to save the file
                end
            else
                % if file doesn't esixt, cotinue.
                continueLog=true(1);
            end
            
        end
        
        
    case 'No'
        saveLogic=false(1);
        savePath=[];
end

set(0, 'DefaultUIControlFontSize', UIControlFontSize);


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
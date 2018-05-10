function [saveLogic,savePath]=QsaveDLTpairStruct(filePaths)
%% function for asking the user for saving options for the reprojection results in STEP1p.
% Including overwriting options
%
% INPUTS:
% * filePaths: 1-by-NumberOfPairs cell array, containing in each
%   cell the path to the DLTparameters_cam# file
%
% OUTPUTS:
% * saveLogic: logic true=save, false=not save
% * savePath: a path were to save the results. [] if saveLogic=false.

%%
saveButton = questdlg('Save DLT pairs reprojection structure?', 'Save?', 'Yes', 'No', 'Yes');

switch saveButton
    case 'Yes'
        saveLogic=true(1);
        continueLog=false(1);
        
        while ~continueLog
            % Path where to save the camera parameters
            savePath = uigetdir(fileparts(filePaths{1}),'Select a folder for saving the results');
            
            % check if DLT parameters in this folder already exist
            fileExistAll=false(1,numel(filePaths));
            for ic=1:numel(filePaths)
                icam=str2num(filePaths{ic}(end-1:end));
                filenames{ic}=[savePath '\DLTparameters_cam' num2str(icam,'%02i') '.mat'];
                if exist(filenames{ic},'file')
                    fileExistAll(1,ic)=true(1);
                end
            end
            if any(fileExistAll)
                fileExistTrue=find(fileExistAll);
                
                promptMessage=cell(numel(fileExistTrue)+2,1);
                promptMessage{1}='These output files already exists: ';
                promptMessage{numel(fileExistTrue)+2}='Do you want to overwrite?';
                for ii=1:numel(fileExistTrue)
                    promptMessage{ii+1}=filenames{fileExistTrue(ii)};
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
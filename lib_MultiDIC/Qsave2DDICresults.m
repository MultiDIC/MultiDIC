function [saveLogic,savePath]=Qsave2DDICresults(imagePaths)
%% function for asking the user for saving options for the 2D-DIC results in STEP2.
% Including overwriting options
%
% INPUTS:
% * imagePaths: 1-by-2 cell array, containing in each
%   cell the path to the speckle images. The camera indices are retrieved
%   from the folder names
%
% OUTPUTS:
% * saveLogic: logic true=save, false=not save
% * savePath: a path were to save the results. [] if saveLogic=false.

%%
saveButton = questdlg('Save 2D-DIC results?', 'Save?', 'Yes', 'No', 'Yes');

switch saveButton
    case 'Yes'
        saveLogic=true(1);
        continueLog=false(1);
        
        while ~continueLog
            % Path where to save the camera parameters
            savePath = uigetdir(fileparts(imagePaths{1}),'Select a folder for saving the results');
            
            % check if DIC2DpairResults in this folder already exist
            icam1=imagePaths{1}(end-1:end);
            icam2=imagePaths{2}(end-1:end);
            filename=fullfile(savePath, ['DIC2DpairResults_C' icam1 'C' icam2]);
            if exist(filename,'file')
                
                promptMessage={'This results file already exists: '; filename; 'Do you want to overwrite?'};
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
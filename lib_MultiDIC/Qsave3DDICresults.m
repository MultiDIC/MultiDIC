function [saveLogic,savePath]=Qsave3DDICresults(structPaths)
%% function for asking the user for saving options for the 3D-DIC results in STEP3.
% Including overwriting options
%
% INPUTS:
% * structPaths: 1-by-Npairs cell array, containing in each
%   cell the path to the DIC2DpairResults_C_#1_C_#2 structure file
%
% OUTPUTS:
% * saveLogic: logic true=save, false=not save
% * savePath: a path were to save the results. [] if saveLogic=false.

%%
saveButton = questdlg('Save 3D-DIC results?', 'Save?', 'Yes', 'No', 'Yes');

switch saveButton
    case 'Yes'
        saveLogic=true(1);
        continueLog=false(1);
        
        while ~continueLog
            % Path where to save the camera parameters
            savePath = uigetdir(fileparts(structPaths{1}),'Select a folder for saving the results');
            
            % check if DIC3DpairResults in this folder already exist
            filenames=cell(1,numel(structPaths));
            fileExistAll=false(1,numel(structPaths));
            for ic=1:numel(structPaths)
                [~,name,~]=fileparts(structPaths{ic});
                filenames{ic}=[savePath '\DIC3DpairResults_' name(end-6:end) '.mat'];
                if exist(filenames{ic},'file')
                    fileExistAll(1,ic)=true(1);
                end
            end
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
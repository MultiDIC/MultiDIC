function [saveLogic,savePath]=QsaveDLTParameters(imagePaths)

saveButton = questdlg('Save DLT parameters?', 'Save?', 'Yes', 'No', 'Yes');

switch saveButton
    case 'Yes'
        saveLogic=true(1);
        continueLog=false(1);
        
        while ~continueLog
            % Path where to save the camera parameters
            savePath = uigetdir(fileparts(imagePaths{1}),'Select a folder for saving the results');
            
            % check if DLT parameters in this folder already exist
            fileExistAll=false(1,numel(imagePaths));
            for ic=1:numel(imagePaths)
                icam=str2num(imagePaths{ic}(end-1:end));
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
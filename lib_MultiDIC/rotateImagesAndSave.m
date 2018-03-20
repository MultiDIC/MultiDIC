function []=rotateImagesAndSave(varargin)
%% rotateImagesAndSave
%  rotateImagesAndSave()
%  rotateImagesAndSave(startPath)
%  rotateImagesAndSave(startPath,OverWriteLogic)
%  rotateImagesAndSave(startPath,OverWriteLogic,rotateValue)

%%
nArg=numel(varargin);
switch nArg
    
    case 0
        startPath=pwd;
        imagePaths = uipickfiles('FilterSpec',startPath,'Prompt','Select images to rotate');
        OverWriteButton = questdlg('Overwrite or save in a new folder?', 'Remove distortion?', 'Overwrite', 'Select a new folder for saving', 'Overwrite');
        switch OverWriteButton
            case 'Select a new folder for saving'
                folderName = uigetdir(startPath,'Select a folder for saving the images');
                OverWriteLogic=false(1);
            case 'Overwrite'
                OverWriteLogic=true(1);
        end
        rotateButton = inputdlg(sprintf('Rotate Images: \n insert a value [deg] in the counterclockwise direction'), 'rotate Images?', 1,{'0'});
        rotateValue = str2double(rotateButton{1});

    case 1
        startPath=varargin{1};
        if isempty(startPath)
            startPath=pwd;
        end
        imagePaths = uipickfiles('FilterSpec',startPath,'Prompt','Select images to rotate');
        OverWriteButton = questdlg('Overwrite or save in a new folder?', 'Remove distortion?', 'Overwrite', 'Select a new folder for saving', 'Overwrite');
        switch OverWriteButton
            case 'Select a new folder for saving'
                folderName = uigetdir(startPath,'Select a folder for saving the images');
                OverWriteLogic=false(1);
            case 'Overwrite'
                OverWriteLogic=true(1);
        end
        rotateButton = inputdlg(sprintf('Rotate Images: \n insert a value [deg] in the counterclockwise direction'), 'rotate Images?', 1,{'0'});
        rotateValue = str2double(rotateButton{1});
    
    case 2
        startPath=varargin{1};
        if isempty(startPath)
            startPath=pwd;
        end
        imagePaths = uipickfiles('FilterSpec',startPath,'Prompt','Select images to rotate');
        OverWriteLogic=varargin{2};
        rotateButton = inputdlg(sprintf('Rotate Images: \n insert a value [deg] in the counterclockwise direction'), 'rotate Images?', 1,{'0'});
        rotateValue = str2double(rotateButton{1});
    
    case 3
        startPath=varargin{1};
        if isempty(startPath)
            startPath=pwd;
        end
        imagePaths = uipickfiles('FilterSpec',startPath,'Prompt','Select images to rotate');
        OverWriteLogic=varargin{2};
        if isempty(OverWriteLogic)
            OverWriteButton = questdlg('Overwrite or save in a new folder?', 'Remove distortion?', 'Overwrite', 'Select a new folder for saving', 'Overwrite');
            switch OverWriteButton
                case 'Select a new folder for saving'
                    folderName = uigetdir(startPath,'Select a folder for saving the images');
                    OverWriteLogic=false(1);
                case 'Overwrite'
                    OverWriteLogic=true(1);
            end
        end
        rotateValue = varargin{3};

    otherwise
        error('wrong number of input arguments');
end


for ii=1:numel(imagePaths)
    IM=imread(imagePaths{ii});
    IMr = imrotate(IM,rotateValue);
    switch OverWriteLogic
        case 0
            [~,imName,ext]=fileparts(imagePaths{ii});
            imwrite(IMr,[folderName '\' imName ext]);
        case 1
            imwrite(IMr,imagePaths{ii});
    end
end


end
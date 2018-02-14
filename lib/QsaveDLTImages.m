
function [saveLogic,overWriteLogic]=QsaveDLTImages()

saveButton = questdlg('Save DLT processed Images?', 'Save?', 'Yes', 'No', 'Yes');
switch saveButton
    case 'Yes'
        saveLogic=true(1);
        overWriteButton = questdlg('Overwrite existing images if they exist?', 'overwrite?', 'Yes', 'No', 'Yes');
        switch overWriteButton
            case 'Yes'
                overWriteLogic=true(1);
            case 'No'
                overWriteLogic=false(1);
        end
    case 'No'
        saveLogic=false(1);
        overWriteLogic=false(1);
end


end
function [saveLogic,overWriteLogic]=QsaveUndistortedImages()
%% function for asking the user for saving options for the undistorted images in STEP0.
% Including overwriting options
%
% OUTPUTS:
% * saveLogic: logic true=save, false=not save
% * overWriteLogic

%%
UIControlFontSize = get(0, 'DefaultUIControlFontSize');
set(0, 'DefaultUIControlFontSize', 11);

saveButton = questdlg('Save Undistorted Images?', 'Save?', 'Yes', 'No', 'Yes');
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
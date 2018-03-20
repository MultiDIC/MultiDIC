function outstate = ncorr_util_properimgfmt(img,name)
% This function determines if the input images (either a numeric array or a 
% string input of the filename with path) have the correct format. The
% string allows us to test for saved images without explicitly loading them, 
% while the array allows us to test images which were loaded through the 
% matlab workspace and may not be explicitly saved.
%
% Inputs -----------------------------------------------------------------%
%   img - numerical array or string; 
%   name - string; used to display errordlg if array is incorrect format
%
% Outputs ----------------------------------------------------------------%
%   outstate - integer; returns either out.cancelled, out.failed, or out.success.
%
% Limiting size of 5x5 is due to the fact that the b-spline coefficients
% need a size of at least 5 in each direction to be convolved with the 
% kernel in ncorr_class_img.form_bcoef(). This is subject to change.
% Also note that existance of the image must be checked (in the case a 
% filename is passed) before calling this function.
    
    % Initialize outputs
    outstate = out.failed;
    
    if (isnumeric(img))
        % Test array directly
        if ((isa(img,'uint8') || isa(img,'uint16') || isa(img,'double')) && ...
            (length(size(img)) == 2 || length(size(img)) == 3) && (size(img,3) == 1 || size(img,3) == 3) && ...
            size(img,1) >= 5 && size(img,2) >= 5)
            % Image format is correct
            outstate = out.success;
        end
    elseif (ischar(img))
        % Test using imfinfo
        info_img = imfinfo(img);        
        if ((strcmp(info_img.ColorType,'truecolor') || strcmp(info_img.ColorType,'grayscale') || strcmp(info_img.ColorType,'indexed')) && ...
            info_img.Height >= 5 && info_img.Width >= 5)
            % Image format is correct
            outstate = out.success;
        end
    end
    
    % Display error dialogue if image format isnt correct
    if (outstate == out.failed)        
        h_error = errordlg([name ' must be of class uint8, uint16, or double, have size of mxn or mxnx3, and must have a height and width greater than or equal to 5.'],'Error','modal');
        uiwait(h_error);
    end
end

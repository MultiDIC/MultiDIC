function [img,outstate] = ncorr_util_loadsavedimg(img_save)
% This function loads a saved image using info in img_save and returns the 
% corresponding ncorr_class_img.
%
% Inputs -----------------------------------------------------------------%
%   img_save - struct; struct('type',{},'gs',{},'name',{},'path',{},'roi',{})
%   note that 'type' is either 'load', 'lazy', or 'file'. For 'lazy', gs is
%   empty;
%
% Outputs ----------------------------------------------------------------%
%   outstate - integer; returns either out.cancelled, out.failed, or out.success.
% 
% Note that if image is 'lazy' or 'file', it will be loaded based on its 
% filename and path. When the image is not found on the path specified, 
% this function will also check the current directory. This allows saved 
% data to be used on other computers. 
% There is currently no checking done for the possibility of the user to 
% change the size of saved image.
    
    % Initialize outputs
    outstate = out.failed;
    img = ncorr_class_img.empty;

    img_prelim = ncorr_class_img;
    if (strcmp(img_save.type,'load'))
        % 'load' image; gs is specified
        img_prelim.set_img(img_save.type,struct('img',img_save.gs,'name',img_save.name,'path',img_save.path)); 
        
        outstate = out.success;
    else
        % Must find image since it is either 'lazy' or 'file'
        if (exist(fullfile(img_save.path,img_save.name),'file'))
            % File exists in path specified
            if (strcmp(img_save.type,'lazy'))
                % This is a 'lazy' image
                img_prelim.set_img(img_save.type,struct('name',img_save.name,'path',img_save.path));
        
                outstate = out.success;
            else
                % This is a 'file' image
                img_prelim.set_img(img_save.type,struct('img',imread(fullfile(img_save.path,img_save.name)),'name',img_save.name,'path',img_save.path));
        
                outstate = out.success;
            end
        else
            % Try finding image in current directory
            if (exist(fullfile(pwd,img_save.name),'file'))
                % Filename exists in current directory, use 'pwd' as the path
                if (strcmp(img_save.type,'lazy'))
                    % This is a 'lazy' image
                    img_prelim.set_img(img_save.type,struct('name',img_save.name,'path',pwd));
        
                    outstate = out.success;
                else
                    % This is a 'file' image
                    img_prelim.set_img(img_save.type,struct('img',imread(fullfile(pwd,img_save.name)),'name',img_save.name,'path',pwd));
        
                    outstate = out.success;
                end                       
            end
        end
    end
    
    % Assign output
    if (outstate == out.success)
        img = img_prelim;
    end
end


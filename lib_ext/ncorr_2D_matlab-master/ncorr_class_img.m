classdef ncorr_class_img < handle
% This is the class definition for an image.
    
    % Properties ---------------------------------------------------------%
    properties(Access = private)
        % These properties require a get* function call
        img;            % uint8
        gs;             % double
        bcoef;          % double
    end
    
    properties(SetAccess = private)
        type;           % string
        height;         % double
        width;          % double
        max_gs;         % double
        min_gs;         % double
        border_bcoef;   % int
        name;           % string
        path;           % string
    end
   
    % Methods ------------------------------------------------------------%
    methods(Access = public)
        % Constructor
        function obj = ncorr_class_img
            obj.img = uint8([]);          
            obj.gs = [];       
            obj.bcoef = [];        
            obj.type = '';           
            obj.height = [];        
            obj.width = [];
            obj.max_gs = [];         
            obj.min_gs = [];    
            obj.border_bcoef = [];     
            obj.name = '';      
            obj.path = '';           
        end
        
        function set_img(obj,type_i,data_i)
        % This function will set the image and other useful properties.
        %
        % Inputs ---------------------------------------------------------%
        %   type_i - string; describes how the image was loaded. Supported
        %   types shown below:
        %       'load' - full image loaded manually through workspace.
        %       data_i is struct('img',{},'name',{},'path,{}); path
        %       will be empty.
        %       'file' - full image loaded through image file. data_i is 
        %       struct('img',{},'name',{},'path,{}).
        %       'lazy' - full image that doesn't save image data and loads
        %       it on demand. data_i is struct('name',{},'path,{}).
        %       'reduced' - reduced image data is saved directly; this is 
        %       usually only used for display purposes as a buffer. 
        %       data_i is struct('img',{},'name',{},'path,{}).
        %   data_i - struct;
        %
        % Outputs ------------------------------------------------------%
        %   None
        % 
        % Note, returns error if type is incorrect. Proper image data format
        % needs to be determined with ncorr_util_properimgfmt() first before
        % setting the image. The name and path also need to be checked
        % first to make sure they are valid before calling this method.
               
            if (~strcmp(type_i,'file') && ~strcmp(type_i,'load') && ~strcmp(type_i,'lazy') && ~strcmp(type_i,'reduced'))                
                error('Incorrect type provided.');
            end
            
            % Input image array should be a monochrome 16-bit image 
            % (uint16). 8-bit (8uint) RGB-color and 8-bit monochrome 
            % images are also supported. Image size must either be mxn or
            % mxnx3. RGB-color images are converted to monochrome for 
            % computational purposes. 

            % Both grayscale values and mxnx3 RGB are stored. This is 
            % because when overlaying images, the colormap will not 
            % apply to it if it is an mxnx3 sized array. Thus, the "img" 
            % property is only used for display purposes when being 
            % overlaid with another plot and a colormap is being used.

            % .jpg will be either mxnx1 or mxnx3
            % .png will be either mxnx1 or mxnx3
            % .bmp will be either mxnx1 or mxnx3
            % .tiff will be either mxnx1, mxnx3, or  mxnx4. NOTE:
            % mxnx4 uses a CMYK colorspace which is not
            % supported

            % Only set gs and image for 'file', 'load', and 'reduced'
            if (strcmp(type_i,'file') || strcmp(type_i,'load') || strcmp(type_i,'reduced'))
                if (size(data_i.img,3) == 3)
                    % Image is RGB-color
                    % Get image:
                    img_prelim = im2uint8(data_i.img);
                    % Get gs:
                    gs_prelim = im2double(rgb2gray(data_i.img));
                else
                    % Image is monochrome
                    % Get image:            
                    img_prelim = uint8([]); % Must set type and initialize          
                    img_prelim(:,:,1) = im2uint8(data_i.img);     
                    img_prelim(:,:,2) = im2uint8(data_i.img);     
                    img_prelim(:,:,3) = im2uint8(data_i.img);
                    % Get gs:
                    gs_prelim = im2double(data_i.img);
                end

                % Set image:
                obj.img = img_prelim;
                % Set gs:
                obj.gs = gs_prelim;
            end                       

            % Find b-coefficients - only do this for 'file' and 'load'. 
            % 'reduced' images are only used for display purposes and
            % not computation. 'lazy' bcoefficients are set as they 
            % are used and then discarded afterward
            border_bcoef_prelim = 20; % MUST BE GREATER THAN OR EQUAL TO 2!!!!
            if (strcmp(type_i,'file') || strcmp(type_i,'load'))                    
                % Pad data first
                obj.bcoef = zeros(size(obj.gs,1)+2*border_bcoef_prelim,size(obj.gs,2)+2*border_bcoef_prelim);
                obj.bcoef(border_bcoef_prelim+1:end-border_bcoef_prelim,border_bcoef_prelim+1:end-border_bcoef_prelim) = obj.gs;
                % Fill corners
                obj.bcoef(1:border_bcoef_prelim,1:border_bcoef_prelim) = obj.bcoef(border_bcoef_prelim+1,border_bcoef_prelim+1); % topleft
                obj.bcoef(1:border_bcoef_prelim,end-(border_bcoef_prelim-1):end) = obj.bcoef(border_bcoef_prelim+1,end-border_bcoef_prelim); % topright
                obj.bcoef(end-(border_bcoef_prelim-1):end,end-(border_bcoef_prelim-1):end) = obj.bcoef(end-border_bcoef_prelim,end-border_bcoef_prelim); % bottomright
                obj.bcoef(end-(border_bcoef_prelim-1):end,1:border_bcoef_prelim) = obj.bcoef(end-border_bcoef_prelim,border_bcoef_prelim+1); % bottomleft
                % Fill sides
                obj.bcoef(1:border_bcoef_prelim,border_bcoef_prelim+1:end-border_bcoef_prelim) = repmat(obj.bcoef(border_bcoef_prelim+1,border_bcoef_prelim+1:end-border_bcoef_prelim),border_bcoef_prelim,1); % top
                obj.bcoef(border_bcoef_prelim+1:end-border_bcoef_prelim,end-(border_bcoef_prelim-1):end) = repmat(obj.bcoef(border_bcoef_prelim+1:end-border_bcoef_prelim,end-border_bcoef_prelim),1,border_bcoef_prelim); % right
                obj.bcoef(end-(border_bcoef_prelim-1):end,border_bcoef_prelim+1:end-border_bcoef_prelim) = repmat(obj.bcoef(end-border_bcoef_prelim,border_bcoef_prelim+1:end-border_bcoef_prelim),border_bcoef_prelim,1); % bottom
                obj.bcoef(border_bcoef_prelim+1:end-border_bcoef_prelim,1:border_bcoef_prelim) = repmat(obj.bcoef(border_bcoef_prelim+1:end-border_bcoef_prelim,border_bcoef_prelim+1),1,border_bcoef_prelim); % left

                % Form/set b-coefficients
                obj.bcoef = ncorr_class_img.form_bcoef(obj.bcoef);
            end       

            % Set type:
            obj.type = type_i;

            % For height and width for lazy, use imfinfo
            if (strcmp(type_i,'lazy'))               
                info_img = imfinfo(fullfile(data_i.path,data_i.name));

                % Set height:
                obj.height = info_img.Height;
                % Set width:
                obj.width = info_img.Width;
            else
                % Set height:
                obj.height = size(gs_prelim,1);
                % Set width:
                obj.width = size(gs_prelim,2);
            end

            % Set max_gs: set to 1 for now
            obj.max_gs = 1;
            % Set min_gs: set to 0 for now
            obj.min_gs = 0;        
            % Set b-ceofficient border:                
            % This parameter determines the border padding around
            % the image. MUST BE GREATER THAN OR EQUAL TO 2!!!
            obj.border_bcoef = border_bcoef_prelim;                 
            % Set name:
            obj.name = data_i.name;    
            % Set path:
            obj.path = data_i.path;
        end
        
        function img_reduced = reduce(obj,spacing)
        % This function will return a reduced image for display purposes.
        % The image is filtered before downsizing, although the filtering
        % isnt optimized - it only uses a basic heuristic.
        %
        % Inputs ---------------------------------------------------------%
        %   spacing - integer; spacing parameter used to reduce the image
        %
        % Outputs --------------------------------------------------------%
        %   img_reduced - ncorr_class_img; reduced image
        % 
        % Returns error if image has not been set yet.
        
            if (isempty(obj.type))     
                error('Image has not been set yet');
            end
            
            % Initialize output
            img_reduced = ncorr_class_img;

            % Reduce the image if spacing is greater than zero
            img_reduced_buffer = obj.get_img();
            if (spacing > 0)
                % Low pass filter image first:
                for i = 0:size(img_reduced_buffer,3)-1
                    % NOTE: imfilter works on uint8 and uint16; output
                    % will be consistent with input class
                    img_reduced_buffer(:,:,i+1) = imfilter(img_reduced_buffer(:,:,i+1),fspecial('gaussian',[(spacing+1) (spacing+1)],(spacing+1)/2),'same');
                end

                % Resample
                img_reduced_buffer = img_reduced_buffer(1:spacing+1:end,1:spacing+1:end,:);
            end

            % Set image - set the same name and path
            img_reduced.set_img('reduced',struct('img',img_reduced_buffer,'name',obj.name,'path',obj.path));
        end
        
        function img_f = formatted(obj)
        % This function returns the formatted image, which can be inputted to a
        % mex function. Returns as a struct since private properties can't
        % be passed. Mex functions can receive ncorr_class_img as either
        % a class or structure.
        %
        % Inputs ---------------------------------------------------------%
        %   none;
        %
        % Outputs --------------------------------------------------------%
        %   img_f - ncorr_class_img; formatted image
        %
        % Returns error if image has not been set yet.
        
            if (isempty(obj.type))
                error('Image has not been set yet');
            end
            
            % Must send img and gs arrays directly since MEX functions
            % currently do not have the get* functions implemented.
            img_f.img = obj.get_img();
            img_f.gs = obj.get_gs();

            img_f.bcoef = [];
            if (~strcmp(obj.type,'reduced'))
                % Get b-spline coefficients if this isnt a reduced
                % image.
                img_f.bcoef = obj.get_bcoef();
            end

            img_f.type = obj.type;
            img_f.height = obj.height;
            img_f.width = obj.width;
            img_f.max_gs = obj.max_gs;
            img_f.min_gs = obj.min_gs;
            img_f.border_bcoef = int32(obj.border_bcoef());  
            img_f.name = obj.name;
            img_f.path = obj.path;
        end
        
        % ----------------------------------------------------------------%
        % Get functions --------------------------------------------------%
        % ----------------------------------------------------------------%
        
        function img_o = get_img(obj)   
        % This function returns the img array (mxnx3). For 'lazy' loading,
        % it will read the image and then return it on demand.
        %
        % Inputs ---------------------------------------------------------%
        %   none;
        %
        % Outputs --------------------------------------------------------%
        %   img_o - double arrray; mxnx3   
        %
        % Returns error if image has not been set yet. If image is moved
        % during analysis, it will return an mxnx3 array of zeros. Images
        % need to be prechecked to make sure they have a valid format.
        % There is currently no handling for the possibility of the user to
        % change the image size during processing.
        
            if (isempty(obj.type))
                error('Image has not been set yet');
            end
            
            % Initialize
            img_o = uint8([]); % Must set type to uint8

            % For 'file' and 'load' and 'reduced' just return the image. For 'lazy'
            % you must read the image file
            if (strcmp(obj.type,'file') || strcmp(obj.type,'load') || strcmp(obj.type,'reduced'))
                img_o = obj.img;
            else
                % Check to make sure image exists first
                if (exist(fullfile(obj.path,obj.name),'file'))
                    % Read file
                    img_buf = imread(fullfile(obj.path,obj.name));

                    % Images must be mxnx3
                    if (size(img_buf,3) == 3)
                        % Image is RGB-color
                        img_o = im2uint8(img_buf);
                    else
                        % Image is monochrome              
                        img_o(:,:,1) = im2uint8(img_buf);     
                        img_o(:,:,2) = im2uint8(img_buf);     
                        img_o(:,:,3) = im2uint8(img_buf);
                    end
                else
                    % File doesn't exist. 
                    msg_error{1} = 'Image could not be located. Make sure it was not moved from the locations listed here:';
                    msg_error{2} = fullfile(obj.path,obj.name);
                    h_error = errordlg(msg_error,'Error','modal');   
                    uiwait(h_error);  

                    % Return a zero array of the correct size
                    img_o(1:obj.height,1:obj.width,1:3) = 0;
                end                
            end 
        end
        
        function gs_o = get_gs(obj)    
        % This function returns the gs array (mxn). For 'lazy' loading,
        % it will read the image and then return it on demand.
        %
        % Inputs ---------------------------------------------------------%
        %   none;
        %
        % Outputs --------------------------------------------------------%
        %   gs_o - double arrray; mxn   
        %
        % Returns error if image has not been set yet. If image is moved
        % during analysis, it will return an mxn array of zeros. Images
        % need to be prechecked to make sure they have a valid format.
        % There is currently no handling for the possibility of the user to
        % change the image size during processing.
        
            if (isempty(obj.type))
                error('Image has not been set yet');
            end
            
            % Initialize
            gs_o = []; 

            % For 'file' and 'load' and 'reduced' just return the image. For 'lazy'
            % you must read the image file
            if (strcmp(obj.type,'file') || strcmp(obj.type,'load') || strcmp(obj.type,'reduced'))
                gs_o = obj.gs;
            elseif (strcmp(obj.type,'lazy'))
                % Check to make sure image exists first
                if (exist(fullfile(obj.path,obj.name),'file'))
                    % Read file and return it
                    img_buf = imread(fullfile(obj.path,obj.name));  

                    % Convert to gs
                    % For a mxnx3, use rgb2gray then im2double. If it is mxn, 
                    % use im2double directly
                    if (size(img_buf,3) == 3)
                        % Image is RGB-color
                        gs_o = im2double(rgb2gray(img_buf));
                    else
                        % Image is monochrome
                        gs_o = im2double(img_buf);
                    end
                else
                    % File doesn't exist.
                    msg_error{1} = 'Image could not be located. Make it was not moved from the locations listed here:';
                    msg_error{2} = fullfile(obj.path,obj.name);
                    h_error = errordlg(msg_error,'Error','modal');   
                    uiwait(h_error);

                    % Return a zero array of the correct size
                    gs_o(1:obj.height,1:obj.width) = 0;
                end                
            end     
        end
                
        function bcoef_o = get_bcoef(obj) 
        % This function returns the bcoef array (mxn). For 'lazy' loading,
        % it will read the image and then return it on demand.
        %
        % Inputs ---------------------------------------------------------%
        %   none;
        %
        % Outputs --------------------------------------------------------%
        %   bcoef_o - double arrray; mxn           
        %
        % Returns error if image has not been set yet. Also returns error
        % if trying to get bcoefficients for a reduced image. If image is 
        % moved during analysis, it will return an mxn array of zeros. 
        % Images need to be prechecked to make sure they have a valid format.
        % There is currently no handling for the possibility of the user to
        % change the image size during processing.
        
            if (isempty(obj.type))                
                error('Image has not been set yet');
            elseif (strcmp(obj.type,'reduced'))
                error('B-coefficients should never be obtained for reduced images.');
            end
            
            % Initialize
            bcoef_o = [];   

            % For 'file' and 'load' just return them. For 'lazy', load the
            % image and then get the bcoefficients. For 'reduced', return
            % an error since they should never be used with reduced images.            
            if (strcmp(obj.type,'file') || strcmp(obj.type,'load'))
                bcoef_o = obj.bcoef;
            elseif (strcmp(obj.type,'lazy'))                
                % Check to make sure image exists first
                if (exist(fullfile(obj.path,obj.name),'file'))
                    % Read file and return it
                    img_buf = imread(fullfile(obj.path,obj.name));  

                    % Convert to gs
                    % For a mxnx3, use rgb2gray then im2double. If it is mxn, 
                    % use im2double directly
                    if (size(img_buf,3) == 3)
                        % Image is RGB-color
                        gs_buf = im2double(rgb2gray(img_buf));
                    else
                        % Image is monochrome
                        gs_buf = im2double(img_buf);
                    end

                    % Form bcoefficent values for gs
                    % Pad data first
                    bcoef_o = zeros(obj.height+2*obj.border_bcoef,obj.width+2*obj.border_bcoef);
                    bcoef_o(obj.border_bcoef+1:end-obj.border_bcoef,obj.border_bcoef+1:end-obj.border_bcoef) = gs_buf;
                    % Fill corners
                    bcoef_o(1:obj.border_bcoef,1:obj.border_bcoef) = bcoef_o(obj.border_bcoef+1,obj.border_bcoef+1); % topleft
                    bcoef_o(1:obj.border_bcoef,end-(obj.border_bcoef-1):end) = bcoef_o(obj.border_bcoef+1,end-obj.border_bcoef); % topright
                    bcoef_o(end-(obj.border_bcoef-1):end,end-(obj.border_bcoef-1):end) = bcoef_o(end-obj.border_bcoef,end-obj.border_bcoef); % bottomright
                    bcoef_o(end-(obj.border_bcoef-1):end,1:obj.border_bcoef) = bcoef_o(end-obj.border_bcoef,obj.border_bcoef+1); % bottomleft
                    % Fill sides
                    bcoef_o(1:obj.border_bcoef,obj.border_bcoef+1:end-obj.border_bcoef) = repmat(bcoef_o(obj.border_bcoef+1,obj.border_bcoef+1:end-obj.border_bcoef),obj.border_bcoef,1); % top
                    bcoef_o(obj.border_bcoef+1:end-obj.border_bcoef,end-(obj.border_bcoef-1):end) = repmat(bcoef_o(obj.border_bcoef+1:end-obj.border_bcoef,end-obj.border_bcoef),1,obj.border_bcoef); % right
                    bcoef_o(end-(obj.border_bcoef-1):end,obj.border_bcoef+1:end-obj.border_bcoef) = repmat(bcoef_o(end-obj.border_bcoef,obj.border_bcoef+1:end-obj.border_bcoef),obj.border_bcoef,1); % bottom
                    bcoef_o(obj.border_bcoef+1:end-obj.border_bcoef,1:obj.border_bcoef) = repmat(bcoef_o(obj.border_bcoef+1:end-obj.border_bcoef,obj.border_bcoef+1),1,obj.border_bcoef); % left

                    % Form/set b-coefficients
                    bcoef_o = ncorr_class_img.form_bcoef(bcoef_o);
                else
                    % File doesn't exist.
                    msg_error{1} = 'Image could not be located. Make sure it was not moved from the locations listed here:';
                    msg_error{2} = fullfile(obj.path,obj.name);
                    h_error = errordlg(msg_error,'Error','modal');   
                    uiwait(h_error);                        

                    % Return a zero array of the correct size
                    bcoef_o = zeros(obj.height+2*obj.border_bcoef,obj.width+2*obj.border_bcoef);
                end      
            end
        end
    end
   
    methods(Static)
        function plot_bcoef = form_bcoef(plot_gs)    
        % This function returns the b-spline coefficients of the input
        % array data. Note that the size returned will be the same as the
        % input array.
        %
        % Inputs -----------------------------------------------------%
        %   plot_gs - double array; data whose b-spline coefficients need
        %   to be determined
        %
        % Outputs ----------------------------------------------------%
        %   plot_bcoef - double array; biquintic b-spline coefficients
        %   of data
        %
        % Note, input array needs to be at least 5 (the length of biquintic 
        % b-spline kernel) in size for each dimension **OR** be empty, or else an 
        % error will be thrown. The empty array is supported because for 
        % regions which are skipped or become empty, they are left an empty 
        % place holder which can be sent to this function. For conciseness,
        % they are not checked for their size. Typically, even if an array 
        % is smaller than 5x5, it will either be padded first by border_interp 
        % (which must be greater than or equal to 2) or left empty in the 
        % ncorr_alg_extrapdata function, so all arrays will meet the size 
        % requirements if used properly.
        
            if (~isempty(plot_gs) && (size(plot_gs,1) < 5 || size(plot_gs,2) < 5))
                error('Array for obtaining b-spline coefficients must be greater then or equal to 5x5 OR be empty.');
            end
        
            % Initialize
            plot_bcoef = zeros(size(plot_gs));

            % Initialize b coefficient kernel
            kernel_b = [1/120 13/60 11/20 13/60 1/120];

            % FFT across rows first --------------------------------------%
            % Form kernel vector
            kernel_b_x = zeros(1,size(plot_gs,2));
            kernel_b_x(1:3) = kernel_b(3:end);
            kernel_b_x(end-1:end) = kernel_b(1:2);
            kernel_b_x = fft(kernel_b_x);

            % FFT across rows
            for i = 0:size(plot_gs,1)-1
                plot_bcoef(i+1,:) = ifft(fft(plot_gs(i+1,:))./kernel_b_x);
            end

            % Now FFT across columns -------------------------------------%
            % Form kernel vector
            kernel_b_y = zeros(size(plot_gs,1),1);
            kernel_b_y(1:3) = kernel_b(3:end);
            kernel_b_y(end-1:end) = kernel_b(1:2);
            kernel_b_y = fft(kernel_b_y);

            % FFT across columns
            for i = 0:size(plot_gs,2)-1
                plot_bcoef(:,i+1) = ifft(fft(plot_bcoef(:,i+1))./kernel_b_y);
            end
        end
    end
end

classdef ncorr < handle
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%
% Programmed by: Justin Blaber
% Principle Investigator: Antonia Antoniou
%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%
% Start Ncorr with:
%   
%   handles_ncorr = ncorr;
%
%   NOTE: Images/masks can be set using the GUI menu or by using
%   "handles_ncorr.set_ref(data)", "handles_ncorr.set_cur(data)", 
%   "handles_ncorr.set_roi_ref(data)", or "handles_ncorr.set_roi_cur(data)". 
%   These will manually set the reference image, current image(s), or 
%   region of interest if these have been calculated or
%   modified using MATLAB and exist in the main workspace. Sometimes the
%   reference and current images can be obtained using an average of a
%   series of images. The region of interest can also be calculated using
%   various thresholding and edge detecting algorithms. Displacement and 
%   strain data can be obtained through handles_ncorr.data_dic in case the 
%   user is interested in doing further calculations on the deformation 
%   data. In addition, if the menu freezes for some reason due to an error, 
%   it can be unfrozen with "handles_ncorr.refresh()". Lastly, calling 
%   "delete(handles_ncorr)" will force close Ncorr. 
%
%
% MAIN CONSIDERATIONS:
%
%   1) All program data and appdata in the ncorr object are SET by callbacks directly. 
%   2) "Downstream" program data (i.e. dependent data) and appdata are 
%   CLEARED with the "clear_downstream()" function before setting the 
%   new data.
%   3) UI objects are modified by the update* functions through the wrapper 
%   function "util_wrapcallbackall" and not within callback directly.
%
%   The point of doing this was to make each callback as clear and
%   uncluttered as possible. The clearing of "downstream" data and updating
%   of UI controls are mainly for upkeep.
%
%
% MAIN FLOW OF CALLBACK:
%
%   1) Update UI controls. Generally, if a menu option updates information, 
%      it will freeze the top menu to prevent users from calling other 
%      options while the first one is running. 
%   2) Check if data will be overwritten with: "overwrite = clear_downstream('callback')"
%      This returns a logical value which is true if the user chooses
%      to proceed and also clears the downstream data. 
%   3) Perform calculations
%   4) If calculations succeed, store the data directly in the callback
%   5) Update UI controls. This unfreezes the top menu and updates it, 
%      as well as sets any images (such as the reference/current image) 
%      if they were loaded.
%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%
% NCORR_GUI* FUNCTIONS:
%
%   These functions contain GUIs which are called by the callbacks within
%   this "ncorr" object. These are, unlike "ncorr," functions that
%   utilize "uiwait" and nested functions instead of an object. Thus, when 
%   the function returns, all workspace data within those functions get 
%   cleared. Some of these ncorr_gui* functions have an input called 
%   "params_init" which contains initialization data if the function has 
%   been called before. Most of these functions use a modal window because 
%   they must be set before proceeding. These figure are generally not 
%   resizeable, and their size is checked with "ncorr_util_formatwindow" to 
%   ensure the top right edges of the window are within the screen.
%
%   The initialization of the ncorr_gui_functions generally proceeds as:
%
%       1) Initialize outputs
%       2) Create GUI handles
%       3) Run the constructor function
%  
%   The callbacks have the general flow:
%
%       1) Possibly freeze menu if it exists
%       2) Get data with "getappdata()"
%       3) Perform calculations
%       4) Set data with "setappdata()"
%       5) Possibly unfreeze menu if it exists and update UI controls
%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%
% NCORR_ALG* FUNCTIONS:
%
%   These functions contain source-code for algorithms used with the ncorr
%   object callbacks and ncorr_gui* functions. Most of these algorithms are
%   mex files written in C++ and must be compiled with "compile_func_cpp_mex()" 
%   or mex before using ncorr.
%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%
% NCORR_UTIL* FUNCTIONS:
%
%   These functions contain certain utilities used by ncorr, such as
%   checking whether the window is within the screen, formatting a colobar,
%   etc.
%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

    properties (SetAccess = private)   
        % GUI handles:
        handles_gui;
        
        % DIC data:
        reference; 
        current;    
        data_dic;
        
        % Installation data:
        support_openmp;
        total_cores;
    end
    
    % Constructor    
    methods(Access = public)
        function obj = ncorr
            % Set default UI control font size here before UI opens. Larger
            % font sizes can make the GUI too large. This seems to be an
            % error in earlier versions of Matlab
            set(0,'defaultuicontrolfontsize',8);
            
            % Initialize GUI and get GUI handles
            obj.handles_gui = init_gui(obj);     
            
            % Run c-tor 
            feval(ncorr_util_wrapcallbacktrycatch(@obj.constructor,obj.handles_gui.figure));
        end
        
        % Constructor ----------------------------------------------------%
        function constructor(obj)      
            % Check MATLAB version and toolbox compatibility -------------%
            if (obj.check_matlabcompat() ~= out.success)
                % Force close
                obj.callback_close_function([],[],true); 
                return;
            end
            
            % Check MEX installation and load/set ncorr_installinfo.txt --%
            if (obj.check_mexinstall() ~= out.success)
                % Force close
                obj.callback_close_function([],[],true); 
                return;
            end
            
            % Set path ---------------------------------------------------%
            % Path to ncorr install directory must be set or else program 
            % can freeze when loading images from a different directory - 
            % not sure why this happens.
            % Note that if ncorr.m is not in the current directory, then
            % path has already been set.
            listing = dir;            
            if (isempty(strfind(lower(path),lower(pwd))) && any(strcmp('ncorr.m',{listing.name})))  %#ok<STREMP>
                % Ask user if its okay to add path
                contbutton = questdlg('Path has not been set. Path must be set before running program. Press yes to add path.','Continue Operation','Yes','No','Yes');     
                if (strcmp(contbutton,'Yes'))
                    % Add Path
                    path(path,pwd);     
                else           
                    % Force close, ncorr may not to work if path is not set
                    obj.callback_close_function([],[],true);
                    return;     
                end
            end
                        
            % Initialize opengl ------------------------------------------%
            % In earlier versions of Matlab this will fix inverted plots
            % Plotting tools are also run based on opengl
            if (ispc)
                data_opengl = opengl('data');
                if (~data_opengl.Software)
                    % Set opengl to software
                    opengl('software'); 
                end
            end                
            
            % Start timer to fetch name of the handle that points to ncorr.
            % Use a timer because, to my knowledge, there isn't a callback
            % option if the handle pointing to ncorr gets cleared.
            timer_ncorr  = timer('executionMode', 'fixedRate', ...
                                 'Period', 1/5, ...
                                 'TimerFcn', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_update_title(obj,hObject,eventdata),obj.handles_gui.figure));
            start(timer_ncorr);
          
            % Program Data -----------------------------------------------%            
            % Reference image info
            obj.reference = struct('imginfo',{},'roi',{});  
            % Current image(s) info
            obj.current = struct('imginfo',{},'roi',{});    
            % DIC data            
            obj.data_dic = struct('displacements',struct('plot_u_dic',{}, ...                                                % U plot after DIC (can be regular or eulerian)
                                                         'plot_v_dic',{}, ...                                                % V plot after DIC (can be regular or eulerian)
                                                         'plot_corrcoef_dic',{}, ...                                         % Correlaton coef plot after DIC (can be regular or eulerian)
                                                         'roi_dic',{}, ...                                                   % ROI after DIC
                                                         'plot_u_ref_formatted',{}, ...                                      % formatted U plot WRT reference image used for displaying displacements
                                                         'plot_v_ref_formatted',{}, ...                                      % formatted V plot WRT reference image used for displaying displacements
                                                         'roi_ref_formatted',{}, ...                                         % ROI after formatting used for plotting displacements
                                                         'plot_u_cur_formatted',{}, ...                                      % formatted U plot WRT current image used for displaying displacements
                                                         'plot_v_cur_formatted',{}, ...                                      % formatted V plot WRT current image used for displaying displacements
                                                         'roi_cur_formatted',{}), ...                                        % ROI after formatting used for plotting displacements
                                  'dispinfo',struct('type',{}, ...                                                           % Type of DIC: either regular or backward
                                                    'radius',{}, ...                                                         % Radius for DIC
                                                    'spacing',{}, ...                                                        % Spacing for DIC                                                    
                                                    'cutoff_diffnorm',{}, ...                                                % Cutoff for norm of the difference vector
                                                    'cutoff_iteration',{}, ...                                               % Cutoff for the number of iterations
                                                    'total_threads',{}, ...                                                  % Number of threads for computation
                                                    'stepanalysis',struct('enabled',{},'type',{},'auto',{},'step',{}), ...   % Indicates whether or not to implement step analysis for high strain
                                                    'subsettrunc',{}, ...                                                    % Indicates whether or not to implement subset truncation for DIC analysis                
                                                    'imgcorr',struct('idx_ref',{},'idx_cur',{}), ...                         % Image correspondences
                                                    'pixtounits',{}, ...                                                     % Ratio of "units" to pixels. Assumes pixels are square
                                                    'units',{}, ...                                                          % String to display units
                                                    'cutoff_corrcoef',{}, ...                                                % Correlation coefficient cutoff for each formatted displacement plot
                                                    'lenscoef',{}), ...                                                      % Radial lens distortion coefficient
                                  'strains',struct('plot_exx_ref_formatted',{}, ...                                          % Exx Green-Lagragian strain plot
                                                   'plot_exy_ref_formatted',{}, ...                                          % Exy Green-Lagragian strain plot
                                                   'plot_eyy_ref_formatted',{}, ...                                          % Exy Green-Lagragian strain plot
                                                   'roi_ref_formatted',{}, ...                                               % ROI used for plotting strains 
                                                   'plot_exx_cur_formatted',{}, ...                                          % Exx Eulerian-Almansi strain plot 
                                                   'plot_exy_cur_formatted',{}, ...                                          % Exy Eulerian-Almansi strain plot 
                                                   'plot_eyy_cur_formatted',{}, ...                                          % Exy Eulerian-Almansi strain plot 
                                                   'roi_cur_formatted',{}), ...                                              % ROI used for plotting strains                          
                                  'straininfo',struct('radius',{}, ...                                                       % Strain radius used for calculating strains
                                                      'subsettrunc',{}));                                                    % Indicates whether or not to implement subset truncation for strain analysis                
                                       
            % Appdata - setting these to '[]' technically does nothing, but
            % it is included here to show what appdata are/will be present in the
            % figure
            % handles_plot are the handles for the data plots:
            setappdata(obj.handles_gui.figure,'handles_plot',[]);   
            % num_cur is the current image which is currently being displayed:
            setappdata(obj.handles_gui.figure,'num_cur',[]);        
            % timer_ncorr is the timer for this instantiation of ncorr;
            % it is used for fetching the name of the handle points to
            % ncorr in the current workspace
            setappdata(obj.handles_gui.figure,'timer',timer_ncorr);       
            
            % Set visible
            set(obj.handles_gui.figure,'Visible','on');
        end
        
        % Destructor -----------------------------------------------------%
        % This ensures data is deleted if delete is called directly on the
        % ncorr handle.
        function delete(obj) 
            if (ishandle(obj.handles_gui.figure))
                % Setting 3rd argument to true results in a force close.
                obj.callback_close_function([],[],true);  
            end
        end  
        
        % Public methods -------------------------------------------------%
        % These are used for setting the reference image, current image, or
        % ROI manually. 
        function set_ref(obj,ref_prelim)
            % Form wrapped function
            handle_set_ref = obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(ref_prelim)set_ref_source(obj,ref_prelim),obj.handles_gui.figure));
            
            % Call wrapped function
            handle_set_ref(ref_prelim);
        end
                
        function set_cur(obj,cur_prelim)
            % Form wrapped function
            handle_set_cur = obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(cur_prelim)set_cur_source(obj,cur_prelim),obj.handles_gui.figure));
            
            % Call wrapped function
            handle_set_cur(cur_prelim);
        end
                    
        function set_roi_ref(obj,mask_prelim)
            % Form wrapped function
            handle_set_roi_ref = obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(mask_prelim)set_roi_ref_source(obj,mask_prelim),obj.handles_gui.figure));
            
            % Call wrapped function
            handle_set_roi_ref(mask_prelim);
        end    
        
        function set_roi_cur(obj,mask_prelim)
            % Form wrapped function
            handle_set_roi_cur = obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(mask_prelim)set_roi_cur_source(obj,mask_prelim),obj.handles_gui.figure));
            
            % Call wrapped function
            handle_set_roi_cur(mask_prelim);
        end    
        
        function refresh(obj)
            % Form wrapped function
            handle_refresh = obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@()refresh_source(obj),obj.handles_gui.figure));
            
            % Call wrapped function
            handle_refresh();
        end
    end   
    
    methods(Access = private)
        % ----------------------------------------------------------------%
        % Local Utility Function(s) --------------------------------------%
        % ----------------------------------------------------------------%
        
        function handle_wrapcallback = util_wrapcallbackall(obj,handle_callback)
        % This function wraps the callback with the update* UI functions. It
        % also checks to see if an error was generated, and if it was, it
        % will rethrow it and recomplete the finishing update* UI functions
        % to ensure the menu unfreezes. Note that the error will not be
        % rethrown if the figure is closed.
        %
        % Inputs ---------------------------------------------------------%
        %   handle_callback - function handle;
        %
        % Outputs --------------------------------------------------------%
        %   handle_wrapcallback - function handle;

            handle_wrapcallback = @updatestatecallback;

            function updatestatecallback(varargin)
                try
                    % Beginning of callback - freeze menu;
                    obj.freeze_menu();  
                    
                    % Run callback
                    handle_callback(varargin{:});
                    
                    % End of callback; 
                    % Update axes
                    obj.update_axes('set');                        
                    % Update UI menu
                    obj.update_topmenu();
                    % Update state text
                    obj.update_dicstatetext();
                    % Unfreeze top layer of menu
                    obj.unfreeze_menu(); 
                catch err
                    % If figure still exists, update UI and then rethrow
                    % error
                    if (isvalid(obj) && ishandle(obj.handles_gui.figure))
                        % Error is generated in callback
                        obj.update_axes('set');                        
                        % Update UI menu
                        obj.update_topmenu();
                        % Update state text
                        obj.update_dicstatetext();
                        % Unfreeze top layer of menu
                        obj.unfreeze_menu(); 
                        
                        % Rethrow error
                        rethrow(err);
                    end
                end
            end
        end
        
        function name = util_get_name(obj)
        % This function will obtain the name of the variable in the main 
        % workspace that points to this ncorr object
            % Initialize output
            name = '';
            
            % Get all variables in base workspace
            basevars = evalin('base','whos');
            
            % Grab variables of class ncorr
            ncorrvars = basevars(strcmp({basevars.class},class(obj)));  
            
            % Cycle through variables and find which one is equivalent to
            % this object
            for i = 1:length(ncorrvars)
                if (eq(evalin('base',ncorrvars(i).name),obj))
                    % Its possible there are more than one handle that
                    % point to this object, if that's the case just return
                    % the first one
                    name = ncorrvars(i).name;
                    break;
                end
            end
        end
                
        %-----------------------------------------------------------------%
        % Handle Functions Source ----------------------------------------%
        %-----------------------------------------------------------------%      
        
        function set_ref_source(obj,ref_prelim) 
        % This function allows the manual uploading of a reference image;
        % only a single image is permitted. 
            % Check to make sure input is of the correct form
            if (ncorr_util_properimgfmt(ref_prelim,'Reference image') == out.success)
                % See if data will be overwritten 
                overwrite = obj.clear_downstream('set_ref');
                if (overwrite)                      
                    % Set data
                    obj.reference(1).imginfo = ncorr_class_img;
                    obj.reference(1).imginfo.set_img('load',struct('img',ref_prelim,'name','reference','path',''));
                    obj.reference(1).roi = ncorr_class_roi.empty;
                end
            end
        end
        
        function set_cur_source(obj,cur_prelim)
        % This function allows the manual uploading of a current image -
        % if a single image is uploaded it must be a double array. If multiple
        % images are uploaded, they must be in a cell array. 
            % Check to make sure input is a cell           
            if (~iscell(cur_prelim))
                % Wrap cur_prelim in a cell
                cur_prelim = {cur_prelim};
            end
            
            % Make sure each image has the correct format
            cursameformat = true;
            for i = 0:length(cur_prelim)-1
                if (ncorr_util_properimgfmt(cur_prelim{i+1},'Current image') ~= out.success)
                    cursameformat = false;
                    break;
                end
            end
            
            if (cursameformat)
                % See if data will be overwritten
                overwrite = obj.clear_downstream('set_cur');                
                if (overwrite)
                    % It's possible to run out of memory here, so put this
                    % in a try-catch block
                    try
                        % Set data
                        for i = 0:length(cur_prelim)-1
                            obj.current(i+1).imginfo = ncorr_class_img;
                            obj.current(i+1).imginfo.set_img('load',struct('img',cur_prelim{i+1},'name',['current_' num2str(i+1)],'path',''));
                            obj.current(i+1).roi = ncorr_class_roi.empty;
                        end
                        setappdata(obj.handles_gui.figure,'num_cur',length(cur_prelim)-1);
                    catch %#ok<CTCH>
                        h_error = errordlg('Loading current images failed, most likely because Ncorr ran out of memory.','Error','modal');  
                        uiwait(h_error);
                        
                        % Clear all current images because exception could have
                        % been thrown midway through setting current images
                        obj.current(:) = [];
                        setappdata(obj.handles_gui.figure,'num_cur',[]);
                    end
                end
            else
                h_error = errordlg('Loading current images failed because the image data format was incorrect.','Error','modal');  
                uiwait(h_error);
            end
        end
        
        function set_roi_ref_source(obj,mask_prelim) 
        % This function allows the manual uploading of a region of
        % interest.
            % Check to make sure reference image has been loaded first -
            % for menu callbacks this isn't necessary since these options
            % will only become callable if "upstream" data has been loaded
            if (~isempty(obj.reference))
                % Make sure ROI is the same size as the reference image and
                % is of class logical. 
                if (isequal(size(mask_prelim),[obj.reference.imginfo.height obj.reference.imginfo.width]) && islogical(mask_prelim))
                    % Form roi_prelim
                    roi_prelim = ncorr_class_roi;
                    roi_prelim.set_roi('load',struct('mask',mask_prelim,'cutoff',2000));
                    % Make sure ROI is not empty
                    if (roi_prelim.get_fullregions > 0)
                        % At this point, roi_prelim fits the critera for a ROI. Show ROI
                        % overlayed and ask user if ROI is appropriate:        
                        outstate = ncorr_gui_loadroi(obj.reference.imginfo, ...
                                                     roi_prelim, ...
                                                     get(obj.handles_gui.figure,'OuterPosition'));
                                                    
                        if (outstate == out.success)
                            % See if data will be overwritten
                            overwrite = obj.clear_downstream('set_roi');
                            if (overwrite)
                                % At this point, ROI is acceptable, has been
                                % approved by the user, and will now be set.                             
                                % Package data:
                                dispinfo_template.type = 'regular';
                                dispinfo_template.radius = []; 
                                dispinfo_template.spacing = [];
                                dispinfo_template.cutoff_diffnorm = [];
                                dispinfo_template.cutoff_iteration = [];
                                dispinfo_template.total_threads = [];
                                dispinfo_template.stepanalysis = struct('enabled',{},'type',{},'auto',{},'step',{});
                                dispinfo_template.subsettrunc = [];
                                dispinfo_template.imgcorr = struct('idx_ref',{},'idx_cur',{});
                                dispinfo_template.pixtounits = [];
                                dispinfo_template.units = '';
                                dispinfo_template.cutoff_corrcoef = [];
                                dispinfo_template.lenscoef = [];
                                
                                % Store the data:
                                obj.data_dic.dispinfo(1) = dispinfo_template;                    
                                obj.reference(1).roi = roi_prelim;
                            end
                        end
                    else
                        h_error = errordlg('ROI must contain a large contiguous region.','Error','modal');
                        uiwait(h_error);
                    end 
                else
                    h_error = errordlg('Input must be of class logical and the same size as the reference image','Error','modal');
                    uiwait(h_error);
                end
            else
                h_error = errordlg('Reference image has not been loaded yet','Error','modal');
                uiwait(h_error);
            end
        end     
        
        function set_roi_cur_source(obj,mask_prelim) 
        % This function allows the manual uploading of a region of
        % interest. This ROI corresponds to the LAST current image.
            % Check to make sure current image(s) have been loaded first -
            % for menu callbacks this isn't necessary since these options
            % will only become callable if "upstream" data has been loaded
            if (~isempty(obj.current))
                % Make sure ROI is the same size as the last current image
                if (isequal(size(mask_prelim),[obj.current(end).imginfo.height obj.current(end).imginfo.width]) && islogical(mask_prelim))
                    % Form roi_prelim
                    roi_prelim = ncorr_class_roi;
                    roi_prelim.set_roi('load',struct('mask',mask_prelim,'cutoff',2000)); 
                    % Make sure ROI is not empty
                    if (roi_prelim.get_fullregions > 0)
                        % At this point, roi_prelim fits the critera for a ROI. Show ROI
                        % overlayed and ask user if ROI is appropriate:        
                        outstate = ncorr_gui_loadroi(obj.current(end).imginfo, ...
                                                     roi_prelim, ...
                                                     get(obj.handles_gui.figure,'OuterPosition'));
                                                    
                        if (outstate == out.success)
                            % See if data will be overwritten
                            overwrite = obj.clear_downstream('set_roi');
                            if (overwrite)
                                % At this point, ROI is acceptable, has been
                                % approved by the user, and will now be set.
                                % Package data:
                                dispinfo_template.type = 'backward';
                                dispinfo_template.radius = []; 
                                dispinfo_template.spacing = [];
                                dispinfo_template.cutoff_diffnorm = [];
                                dispinfo_template.cutoff_iteration = [];
                                dispinfo_template.total_threads = [];
                                dispinfo_template.stepanalysis = struct('enabled',{},'type',{},'auto',{},'step',{});
                                dispinfo_template.subsettrunc = [];
                                dispinfo_template.imgcorr = struct('idx_ref',{},'idx_cur',{});
                                dispinfo_template.pixtounits = [];
                                dispinfo_template.units = '';
                                dispinfo_template.cutoff_corrcoef = [];
                                dispinfo_template.lenscoef = [];
                                
                                % Store the data:
                                obj.data_dic.dispinfo(1) = dispinfo_template;                    
                                obj.current(end).roi = roi_prelim;
                            end
                        end
                    else
                        h_error = errordlg('ROI must contain a large contiguous region.','Error','modal');
                        uiwait(h_error);
                    end 
                else
                    h_error = errordlg('Input must be of class logical and the same size as the last current image','Error','modal');
                    uiwait(h_error);
                end
            else
                h_error = errordlg('Current image(s) have not been loaded yet.','Error','modal');
                uiwait(h_error);
            end
        end       
        
        function refresh_source(obj) %#ok<MANU>
        % Call this function if there's an error and the menu needs to be
        % unfrozen.
        end
        
        %-----------------------------------------------------------------%
        % Menu Callbacks -------------------------------------------------%
        %-----------------------------------------------------------------%
        function callback_topmenu_loadref(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for loading a reference image from the
        % GUI.
            % false means lazy loading isnt used
            [ref_prelim,outstate] = ncorr_util_loadimgs(false);
            if (outstate == out.success && length(ref_prelim) == 1)     
            % See if data will be overwritten
                overwrite = obj.clear_downstream('set_ref');
                if (overwrite)
                    % Set Image
                    obj.reference(1).imginfo = ref_prelim;
                    obj.reference(1).roi = ncorr_class_roi.empty;
                end
            elseif (outstate == out.success && length(ref_prelim) > 1)
                h_error = errordlg('Please select only one reference image.','Error','modal');  
                uiwait(h_error);
            end
        end

        function callback_topmenu_loadcur(obj,hObject,eventdata,lazyparam) %#ok<INUSL>
        % This is the callback for loading current image(s) from the GUI.
            [cur_prelim,outstate] = ncorr_util_loadimgs(lazyparam);
            if (outstate == out.success)     
                overwrite = obj.clear_downstream('set_cur');
                if (overwrite)
                    % Set Image
                    for i = 0:length(cur_prelim)-1
                        obj.current(i+1).imginfo = cur_prelim(i+1);
                        obj.current(i+1).roi = ncorr_class_roi.empty;
                    end
                    % Display last current image
                    setappdata(obj.handles_gui.figure,'num_cur',length(cur_prelim)-1);
                end
            end
        end
        
        function callback_topmenu_loaddata(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for loading data from a previous analysis. 
            % Get data filename
            [filename,pathname] = uigetfile({'*.mat'},'Select the previous DIC data (must be .dat)');
            if (~isequal(filename,0) && ~isequal(pathname,0))
                try
                    % Load data
                    struct_load = load(fullfile(pathname,filename));
                    loadsuccess = true;
                catch %#ok<CTCH>
                    % Probably out of memory
                    loadsuccess = false;
                end
                
                if (loadsuccess)
                    if (isfield(struct_load,'reference_save') && isfield(struct_load,'current_save') && isfield(struct_load,'data_dic_save'))
                        % The data has the correct fields
                        if (isa(struct_load.reference_save,'struct') && ...
                            isa(struct_load.current_save,'struct') && ... 
                            isa(struct_load.data_dic_save,'struct'))
                            % Make sure data_dic has the correct
                            % fields. Easiest way to do this is to assign
                            % the saved data to a structure with the
                            % correct fields, an error will result if the
                            % fields are wrong
                            data_dic_prelim = obj.data_dic; % This will copy fields - only used temporarily
                            try
                                data_dic_prelim(1) = struct_load.data_dic_save; %#ok<NASGU>
                                loaddatasuccess = true;
                            catch %#ok<CTCH>
                                % Some fields are not correct                                                     
                                loaddatasuccess = false;  
                            end
                            
                            if (loaddatasuccess)                                
                                % See if data will be overwritten
                                overwrite = obj.clear_downstream('all');      
                                if (overwrite)     
                                    % Put in try block - its possible to run
                                    % out of memory here.
                                    try               
                                        % Do reference image first
                                        [ref_prelim,outstate_ref] = ncorr_util_loadsavedimg(struct_load.reference_save);
                                        
                                        % Do current images next
                                        cur_prelim = ncorr_class_img.empty;
                                        for i = 0:length(struct_load.current_save)-1
                                            [cur_buffer,outstate_cur] = ncorr_util_loadsavedimg(struct_load.current_save(i+1));
                                            if (outstate_cur == out.success)
                                                cur_prelim(i+1) = cur_buffer;
                                            else
                                                break;
                                            end
                                        end       

                                        if (outstate_ref == out.success && outstate_cur == out.success)        
                                            % Store data:
                                            % Set reference:
                                            obj.reference(1).imginfo = ref_prelim;
                                            obj.reference(1).roi = struct_load.reference_save.roi;

                                            % Set current:
                                            for i = 0:length(struct_load.current_save)-1
                                                obj.current(i+1).imginfo = cur_prelim(i+1);
                                                obj.current(i+1).roi = struct_load.current_save(i+1).roi;
                                            end                
                                            setappdata(obj.handles_gui.figure,'num_cur',length(struct_load.current_save)-1);

                                            % Set dic data:
                                            obj.data_dic(1) = struct_load.data_dic_save;
                                        else
                                            % Images could not be found
                                            % Form error message
                                            msg_error{1} = 'Images could not be located. Make sure they are in the current directory or were not moved from the locations listed here:';
                                            msg_error{2} = fullfile(struct_load.reference_save.path,struct_load.reference_save.name);
                                            max_display = 10;
                                            for i = 0:min(length(struct_load.current_save),max_display)-1
                                                msg_error{end+1} = fullfile(struct_load.current_save(i+1).path,struct_load.current_save(i+1).name); %#ok<AGROW>
                                            end                         
                                            if (length(struct_load.current_save) > max_display)
                                                msg_error{end+1} = '...';
                                            end
                                            h_error = errordlg(msg_error,'Error','modal');   
                                            uiwait(h_error);
                                        end
                                    catch %#ok<CTCH>
                                        h_error = errordlg('Loading failed, most likely because Ncorr ran out of memory.','Error','modal');  
                                        uiwait(h_error);

                                        % Clear everything - its possible
                                        % exception was caught while
                                        % loading data midway. 
                                        obj.reference(:) = [];
                                        obj.current(:) = [];
                                        setappdata(obj.handles_gui.figure,'num_cur',[]);
                                        fields_data_dic = fieldnames(obj.data_dic); 
                                        for i = 0:length(fields_data_dic)-1
                                            obj.data_dic.(fields_data_dic{i+1})(:) = [];
                                        end  
                                    end
                                end
                            else                                   
                                h_error = errordlg('The data is not valid, please load data saved specifically from ncorr.','Error','modal');  
                                uiwait(h_error);
                            end
                        else         
                            h_error = errordlg('The data is not valid, please load data saved specifically from ncorr.','Error','modal');   
                            uiwait(h_error);
                        end
                    else
                        h_error = errordlg('The data is not valid, please load data saved specifically from ncorr.','Error','modal');
                        uiwait(h_error);
                    end
                else
                    h_error = errordlg('Loading failed, most likely because Ncorr ran out of memory.','Error','modal');  
                    uiwait(h_error);
                end
            end            
        end

        function callback_topmenu_savedata(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for saving data. 
            % Get save data filename
            [filename,pathname] = uiputfile({'*.mat'},'Save DIC data (must be .dat)');
            if (~isequal(filename,0) && ~isequal(pathname,0))
                % Check if data will be overwritten in directory
                overwrite = true;
                if (exist(fullfile(pathname,filename),'file'))
                    contbutton = questdlg('File already exists. Do you want to overwrite?','Continue Operation','Yes','No','Yes');
                    if (strcmp(contbutton,'No'))
                        overwrite = false;
                    end
                end
                
                if (overwrite)  
                    % Must deal with image data, ROI data, and data_dic.                
                    %   1) Images: Save the file names if they are 'file' or 'lazy'
                    %   If they were set through the workspace (i.e. 'load') then save the
                    %   image as double precision.       
                    %   2) ROIs: For the ROIs, save each directly.
                    %   3) Data: Save data_dic directly.
                    reference_save = struct('type',{},'gs',{},'name',{},'path',{},'roi',{});
                    current_save = struct('type',{},'gs',{},'name',{},'path',{},'roi',{});
                
                    % Images ---------------------------------------------%
                    % Reference image:
                    reference_save(1).type = obj.reference.imginfo.type;
                    reference_save.gs = [];
                    if (strcmp(obj.reference.imginfo.type,'load'))
                        % Must save gs data directly
                        reference_save.gs = obj.reference.imginfo.get_gs();
                    end
                    reference_save.name = obj.reference.imginfo.name;
                    reference_save.path = obj.reference.imginfo.path;
                    reference_save.roi = obj.reference.roi; %#ok<STRNU>
                    
                    % Current image(s):
                    for i = 0:length(obj.current)-1
                        current_save(i+1).type = obj.current(i+1).imginfo.type;
                        current_save(i+1).gs = [];
                        if (strcmp(obj.current(i+1).imginfo.type,'load'))
                            % Must save gs data directly        
                            current_save(i+1).gs = obj.current(i+1).imginfo.get_gs();
                        end       
                        current_save(i+1).name = obj.current(i+1).imginfo.name;
                        current_save(i+1).path = obj.current(i+1).imginfo.path;    
                        current_save(i+1).roi = obj.current(i+1).roi;
                    end
                    
                    % Data -----------------------------------------------%
                    data_dic_save = obj.data_dic; %#ok<NASGU>
                    
                    % Save -----------------------------------------------%
                    try
                        % The v7.3 flag helps save larger files
                        save(fullfile(pathname,filename),'reference_save','current_save','data_dic_save','-v7.3');
                    catch %#ok<CTCH>
                        % If saving fails, its generally because the files
                        % are too large.
                        h_error = errordlg('Saving failed, probably because the amount of data being saved is too large.','Error','modal');
                        uiwait(h_error);
                    end
                end
            end
        end

        function callback_topmenu_clear(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for the clearing all the data.
            obj.clear_downstream('all');    
        end     
        
        function callback_topmenu_sethandle(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for setting the handle to ncorr in the case
        % that it gets cleared from the workspace
            % Get Name
            name = util_get_name(obj);
            if (isempty(name))
                % Prompt user for what handle he/she wants to use
                [handle_name,outstate] = gui_sethandle(get(obj.handles_gui.figure,'OuterPosition'));
                if (outstate == out.success)
                    % Assign handles to base workspace
                    assignin('base',handle_name,obj);
                end
            else
                % Handle already exists
                e = msgbox(['Handle already exists and is called: ' name],'WindowStyle','modal'); 
                uiwait(e);
            end    
        end
        
        function callback_topmenu_reinstall(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for reinstalling the mex files.
            % Check if ncorr.m is in current directory            
            listing = dir;            
            if (any(strcmp('ncorr.m',{listing.name})))        
                % Ask user if its okay to reinstall
                contbutton = questdlg('Are you sure you want to reinstall? All data will be lost if proceeding.','Continue Operation','Yes','No','Yes');  
                if (strcmp(contbutton,'Yes'))
                    % Clear ncorr_installinfo.txt file - this is the standard
                    % way to reinstall
                    filename = 'ncorr_installinfo.txt';
                    fid = fopen(filename,'w'); % This will clear file  
                    if (fid ~= -1)    
                        fclose(fid); % Close file

                        % Force close 
                        obj.callback_close_function([],[],true);  

                        % Reopen ncorr - this will cause reinstallation since
                        % ncorr_installinfo.txt has been cleared
                        handles_ncorr = ncorr;

                        % Get name
                        name = util_get_name(obj);
                        if (~isempty(name))
                            % If there is a name then assign the name to the
                            % ncorr object.
                            assignin('base',name,handles_ncorr);
                        end
                    else
                        % Error
                        h_error = errordlg('For some reason ncorr was not able to clear "ncorr_installinfo.txt" file, reinstall cannot take place.','Error','modal');                    
                        uiwait(h_error);
                    end   
                end
            else
                % Error
                h_error = errordlg('Please navigate to folder containing ncorr.m first before reinstalling.','Error','modal');                    
                uiwait(h_error);
            end
        end

        function callback_exit_callback(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for exiting the program.
        % 3rd argument is set to false, which queries user about
        % closing if DIC data has been computed.
            obj.callback_close_function([],[],false);  
        end

        function callback_topmenu_set_roi_ref(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for loading the region of interest from the 
        % GUI
            % Get ROIs ---------------------------------------------------%
            % Check if ROI has been set before - must be from the regular
            % analysis, as backward DIC analysis will also set a reference ROI
            params_init = ncorr_class_roi.empty;
            if (~isempty(obj.reference) && ~isempty(obj.reference.roi) && ...
                ~isempty(obj.data_dic.dispinfo) && strcmp(obj.data_dic.dispinfo.type,'regular'))
                params_init(1) = obj.reference.roi; 
            end
            
            [roi_prelim,outstate] = ncorr_gui_setrois(obj.reference.imginfo, ...
                                                      get(obj.handles_gui.figure,'OuterPosition'), ...
                                                      params_init);                                                                                         
            
            if (outstate == out.success)
                % See if data will be overwritten
                overwrite = obj.clear_downstream('set_roi');                
                if (overwrite) 
                    % Package data:
                    dispinfo_template.type = 'regular';
                    dispinfo_template.radius = []; 
                    dispinfo_template.spacing = [];
                    dispinfo_template.cutoff_diffnorm = [];
                    dispinfo_template.cutoff_iteration = [];
                    dispinfo_template.total_threads = [];
                    dispinfo_template.stepanalysis = struct('enabled',{},'type',{},'auto',{},'step',{});
                    dispinfo_template.subsettrunc = [];
                    dispinfo_template.imgcorr = struct('idx_ref',{},'idx_cur',{});
                    dispinfo_template.pixtounits = [];
                    dispinfo_template.units = '';
                    dispinfo_template.cutoff_corrcoef = [];
                    dispinfo_template.lenscoef = [];
                    
                    % Store the data:
                    obj.data_dic.dispinfo(1) = dispinfo_template;
                    obj.reference(1).roi = roi_prelim;
                end  
            end
        end
        
        function callback_topmenu_set_roi_cur(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for loading the region of interest from
        % the GUI
            % Get ROIs ---------------------------------------------------%
            % Check if ROI has been set before - must be from the backward
            % analysis, as backward dic analysis will also set a current ROI
            params_init = ncorr_class_roi.empty;
            if (~isempty(obj.current) && ~isempty(obj.current(end).roi) && ...
                ~isempty(obj.data_dic.dispinfo) && strcmp(obj.data_dic.dispinfo.type,'backward'))
                params_init(1) = obj.current(end).roi; 
            end
            
            [roi_prelim,outstate] = ncorr_gui_setrois(obj.current(end).imginfo, ...
                                                      get(obj.handles_gui.figure,'OuterPosition'), ...
                                                      params_init);                                                                                         
            
            if (outstate == out.success)
                % See if data will be overwritten
                overwrite = obj.clear_downstream('set_roi');                
                if (overwrite)     
                    % Package data:
                    dispinfo_template.type = 'backward';
                    dispinfo_template.radius = []; 
                    dispinfo_template.spacing = [];
                    dispinfo_template.cutoff_diffnorm = [];
                    dispinfo_template.cutoff_iteration = [];
                    dispinfo_template.total_threads = [];
                    dispinfo_template.stepanalysis = struct('enabled',{},'type',{},'auto',{},'step',{});
                    dispinfo_template.subsettrunc = [];
                    dispinfo_template.imgcorr = struct('idx_ref',{},'idx_cur',{});
                    dispinfo_template.pixtounits = [];
                    dispinfo_template.units = '';
                    dispinfo_template.cutoff_corrcoef = [];
                    dispinfo_template.lenscoef = [];
                    
                    % Store the data:
                    obj.data_dic.dispinfo(1) = dispinfo_template;
                    obj.current(end).roi = roi_prelim;
                end  
            end
        end
        
        function callback_topmenu_setdicparameters(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for setting DIC parameters
            % Check if DIC parameters were set before - if dic parameters
            % are being set, the data_dic.dispinfo field is gauranteed
            % nonempty, so you don't need to check this first.
            params_init = {};
            if (~isempty(obj.data_dic.dispinfo.radius) &&  ...
                ~isempty(obj.data_dic.dispinfo.spacing) && ...
                ~isempty(obj.data_dic.dispinfo.cutoff_diffnorm) &&  ...
                ~isempty(obj.data_dic.dispinfo.cutoff_iteration) && ...
                ~isempty(obj.data_dic.dispinfo.total_threads) && ...
                ~isempty(obj.data_dic.dispinfo.stepanalysis) && ...
                ~isempty(obj.data_dic.dispinfo.subsettrunc))
                params_init = {obj.data_dic.dispinfo.radius,  ...
                               obj.data_dic.dispinfo.spacing, ...
                               obj.data_dic.dispinfo.cutoff_diffnorm, ...
                               obj.data_dic.dispinfo.cutoff_iteration, ...
                               obj.data_dic.dispinfo.total_threads, ...
                               obj.data_dic.dispinfo.stepanalysis, ...
                               obj.data_dic.dispinfo.subsettrunc}; 
            end
            
            % Get Params -------------------------------------------------%
            if (strcmp(obj.data_dic.dispinfo.type,'regular'))
                % Show reference image for regular analysis
                [radius_prelim,spacing_prelim,cutoff_diffnorm_prelim,cutoff_iteration_prelim,total_threads_prelim,stepanalysis_prelim,subsettrunc_prelim,outstate] = ncorr_gui_setdicparams(obj.reference.imginfo, ...
                                                                                                                                                                                            obj.reference.roi, ...
                                                                                                                                                                                            obj.support_openmp, ...
                                                                                                                                                                                            obj.total_cores, ...
                                                                                                                                                                                            get(obj.handles_gui.figure,'OuterPosition'), ...
                                                                                                                                                                                            params_init);   
            else
                % Show last current image for backward analysis
                [radius_prelim,spacing_prelim,cutoff_diffnorm_prelim,cutoff_iteration_prelim,total_threads_prelim,stepanalysis_prelim,subsettrunc_prelim,outstate] = ncorr_gui_setdicparams(obj.current(end).imginfo, ...
                                                                                                                                                                                            obj.current(end).roi, ...
                                                                                                                                                                                            obj.support_openmp, ...
                                                                                                                                                                                            obj.total_cores, ...
                                                                                                                                                                                            get(obj.handles_gui.figure,'OuterPosition'), ...
                                                                                                                                                                                            params_init);   
            end
            
            if (outstate == out.success)
                % See if data will be overwritten
                overwrite = obj.clear_downstream('set_dicparameters');                
                if (overwrite)      
                    % Package data:
                    dispinfo_template.type = obj.data_dic.dispinfo.type;
                    dispinfo_template.radius = radius_prelim; 
                    dispinfo_template.spacing = spacing_prelim;
                    dispinfo_template.cutoff_diffnorm = cutoff_diffnorm_prelim;
                    dispinfo_template.cutoff_iteration = cutoff_iteration_prelim;
                    dispinfo_template.total_threads = total_threads_prelim;
                    dispinfo_template.stepanalysis = stepanalysis_prelim;
                    dispinfo_template.subsettrunc = subsettrunc_prelim;                    
                    dispinfo_template.imgcorr = struct('idx_ref',{},'idx_cur',{});
                    dispinfo_template.pixtounits = [];
                    dispinfo_template.units = '';
                    dispinfo_template.cutoff_corrcoef = [];
                    dispinfo_template.lenscoef = [];
                    
                    % Store the data:
                    obj.data_dic.dispinfo(1) = dispinfo_template;
                end  
            end
        end
        
        function callback_topmenu_dic(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for performing DIC.
            % Merge images. This makes it convenient for the idx
            % tracker if step analysis is used
            imgs = [obj.reference obj.current];  
            
            % Initialize 
            rois_dic_prelim = ncorr_class_roi.empty;
            displacements_prelim = struct('plot_u',{},'plot_v',{},'plot_corrcoef',{},'plot_validpoints',{});    
            imgcorr_prelim = struct('idx_ref',{},'idx_cur',{});
            seedinfo_buffer = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{});
            if (strcmp(obj.data_dic.dispinfo.type,'regular'))  
                % Imgcorr starts out empty so this guarantees at least one
                % iteration. Make sure to exit when all images have been
                % analyzed.
                while (isempty(imgcorr_prelim) || imgcorr_prelim(end).idx_cur < length(imgs)-1)                    
                    % Update idx_ref
                    if (isempty(imgcorr_prelim))
                        imgcorr_prelim(1).idx_ref = 0;
                    else
                        % Use previous last current image as the new
                        % reference image
                        imgcorr_prelim(end+1).idx_ref = imgcorr_prelim(end).idx_cur; %#ok<AGROW>
                    end
                    
                    % Update idx_cur
                    if (obj.data_dic.dispinfo.stepanalysis.enabled && strcmp(obj.data_dic.dispinfo.stepanalysis.type,'leapfrog'))
                        % If leapfrog is used, then set current image to a
                        % preset increment
                        imgcorr_prelim(end).idx_cur = min(imgcorr_prelim(end).idx_ref+obj.data_dic.dispinfo.stepanalysis.step,length(imgs)-1);
                    else
                        % If this is not leapfrog, set the current image to
                        % the last image. This number will get updated
                        % depending on how many images are successfully
                        % analyzed and seeded
                        imgcorr_prelim(end).idx_cur = length(imgs)-1;
                    end
                    
                    % Do DIC analysis - only send initial parameters if
                    % stepanalysis is enabled with automatic seed
                    % propagation and an interation has already been done.
                    params_init = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{});
                    if (obj.data_dic.dispinfo.stepanalysis.enabled && ...
                        obj.data_dic.dispinfo.stepanalysis.auto && ...
                        ~isempty(seedinfo_buffer))
                        % Get previous seeds as initial parameters if step 
                        % analysis and automatic seed propagation are enabled. 
                        % If these are provided then the seeds are automatically 
                        % updated
                        params_init = seedinfo_buffer;
                    end
                    
                    % Returns failed if exception is thrown; for both 
                    % failed and cancel just return.
                    [displacements_buffer,rois_dic_buffer,seedinfo_buffer,outstate_dic] = ncorr_alg_dicanalysis(imgs(imgcorr_prelim(end).idx_ref+1:imgcorr_prelim(end).idx_cur+1), ...
                                                                                                                obj.data_dic.dispinfo.radius, ...
                                                                                                                obj.data_dic.dispinfo.spacing, ...
                                                                                                                obj.data_dic.dispinfo.cutoff_diffnorm, ...
                                                                                                                obj.data_dic.dispinfo.cutoff_iteration, ...
                                                                                                                obj.data_dic.dispinfo.total_threads, ...
                                                                                                                obj.data_dic.dispinfo.stepanalysis.enabled, ...
                                                                                                                obj.data_dic.dispinfo.subsettrunc, ...
                                                                                                                imgcorr_prelim(end).idx_ref, ...
                                                                                                                length(obj.current), ...
                                                                                                                get(obj.handles_gui.figure,'OuterPosition'), ...
                                                                                                                params_init);
                    if (outstate_dic ~= out.success)  
                        % DIC analysis was cancelled or failed by throwing
                        % an exception, for both cases return
                        return;
                    end
                    
                    % Update ROIs
                    for i = 0:length(displacements_buffer)-1   
                        imgs(i+imgcorr_prelim(end).idx_ref+2).roi = imgs(imgcorr_prelim(end).idx_ref+1).roi.update_roi(displacements_buffer(i+1).plot_u, ...
                                                                                                                       displacements_buffer(i+1).plot_v, ...
                                                                                                                       rois_dic_buffer(i+1), ...
                                                                                                                       [imgs(i+imgcorr_prelim(end).idx_ref+2).imginfo.height imgs(i+imgcorr_prelim(end).idx_ref+2).imginfo.width], ...
                                                                                                                       obj.data_dic.dispinfo.spacing, ...
                                                                                                                       obj.data_dic.dispinfo.radius);
                    end                     

                    % Store outputs
                    for i = 0:length(displacements_buffer)-1
                        displacements_prelim(i+imgcorr_prelim(end).idx_ref+1) = displacements_buffer(i+1);
                        rois_dic_prelim(i+imgcorr_prelim(end).idx_ref+1) = rois_dic_buffer(i+1);
                    end

                    % Update imgcorr_prelim if stepanalysis is enabled 
                    % - this is necessary since not all images are 
                    % necessarily analyzed.
                    if (obj.data_dic.dispinfo.stepanalysis.enabled)
                        imgcorr_prelim(end).idx_cur = imgcorr_prelim(end).idx_ref + length(displacements_buffer);
                    end                        
                end
            else  
                % Imgcorr starts out empty so this guarantees at least one
                % iteration. Make sure to exit when all images have been
                % analyzed.
                while (isempty(imgcorr_prelim) || imgcorr_prelim(end).idx_cur > 0)                          
                    % Update idx_ref
                    if (isempty(imgcorr_prelim))
                        imgcorr_prelim(1).idx_ref = length(imgs)-1;
                    else
                        % Use previous last current image as the new
                        % reference image
                        imgcorr_prelim(end+1).idx_ref = imgcorr_prelim(end).idx_cur; %#ok<AGROW>
                    end
                    
                    % Update idx_cur
                    if (obj.data_dic.dispinfo.stepanalysis.enabled && strcmp(obj.data_dic.dispinfo.stepanalysis.type,'leapfrog'))
                        % If leapfrog is used, then set current image to a
                        % preset increment
                        imgcorr_prelim(end).idx_cur = max(imgcorr_prelim(end).idx_ref-obj.data_dic.dispinfo.stepanalysis.step,0);
                    else
                        % If this is not leapfrog, set the current image to
                        % the last image. This number will get updated
                        % depending on how many images are successfully
                        % analyzed and seeded
                        imgcorr_prelim(end).idx_cur = 0;
                    end   
                    
                    % Do DIC analysis - only send initial parameters if
                    % stepanalysis is enabled with automatic seed
                    % propagation and an interation has already been done.
                    params_init = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{});
                    if (obj.data_dic.dispinfo.stepanalysis.enabled && ...
                        obj.data_dic.dispinfo.stepanalysis.auto && ...
                        ~isempty(seedinfo_buffer))
                        % Get previous seeds as initial parameters if step 
                        % analysis and automatic seed propagation are enabled. 
                        % If these are provided then the seeds are automatically 
                        % updated
                        params_init = seedinfo_buffer;
                    end
                                        
                    % Returns failed if exception is thrown; for both 
                    % failed and cancel just return.
                    [displacements_buffer,rois_dic_buffer,seedinfo_buffer,outstate_dic] = ncorr_alg_dicanalysis(imgs(imgcorr_prelim(end).idx_ref+1:-1:imgcorr_prelim(end).idx_cur+1), ...
                                                                                                                obj.data_dic.dispinfo.radius, ...
                                                                                                                obj.data_dic.dispinfo.spacing, ...
                                                                                                                obj.data_dic.dispinfo.cutoff_diffnorm, ...
                                                                                                                obj.data_dic.dispinfo.cutoff_iteration, ...
                                                                                                                obj.data_dic.dispinfo.total_threads, ...
                                                                                                                obj.data_dic.dispinfo.stepanalysis.enabled, ...
                                                                                                                obj.data_dic.dispinfo.subsettrunc, ...
                                                                                                                length(obj.current)-imgcorr_prelim(end).idx_ref, ...
                                                                                                                length(obj.current), ...
                                                                                                                get(obj.handles_gui.figure,'OuterPosition'), ...
                                                                                                                params_init);
                    if (outstate_dic ~= out.success)  
                        % DIC analysis was cancelled or failed by throwing
                        % an exception, for both cases return
                        return;
                    end
                    
                    % Update ROIs
                    for i = 0:length(displacements_buffer)-1                  
                        imgs(imgcorr_prelim(end).idx_ref-i).roi = imgs(imgcorr_prelim(end).idx_ref+1).roi.update_roi(displacements_buffer(i+1).plot_u, ...
                                                                                                                     displacements_buffer(i+1).plot_v, ...
                                                                                                                     rois_dic_buffer(i+1), ...
                                                                                                                     [imgs(imgcorr_prelim(end).idx_ref-i).imginfo.height imgs(imgcorr_prelim(end).idx_ref-i).imginfo.width], ...
                                                                                                                     obj.data_dic.dispinfo.spacing, ...
                                                                                                                     obj.data_dic.dispinfo.radius);
                    end

                    % Update displacements                              
                    for i = 0:length(displacements_buffer)-1  
                        % Displacements must be updated if more than two 
                        % images were analyzed. Also, leave the last 
                        % displacement field in the buffer alone since 
                        % it's already WRT the correct configuration
                        if (i > 0)                                        
                            % Update displacement fields    
                            params_init = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{});
                            if (obj.data_dic.dispinfo.stepanalysis.auto)
                                params_init = seedinfo_buffer(:,:,end-i);
                            end
                            
                            % Returns failed if exception is thrown; for both 
                            % failed and cancel just return.
                            [displacement_buffer_update,roi_dic_buffer_update,seedinfo_buffer_update,outstate_dic] = ncorr_alg_dicanalysis([imgs(imgcorr_prelim(end).idx_ref-length(displacements_buffer)+i+1) imgs(imgcorr_prelim(end).idx_ref-length(displacements_buffer)+1)], ...
                                                                                                                                           obj.data_dic.dispinfo.radius, ...
                                                                                                                                           obj.data_dic.dispinfo.spacing, ...
                                                                                                                                           obj.data_dic.dispinfo.cutoff_diffnorm, ...
                                                                                                                                           obj.data_dic.dispinfo.cutoff_iteration, ...
                                                                                                                                           obj.data_dic.dispinfo.total_threads, ...
                                                                                                                                           obj.data_dic.dispinfo.stepanalysis.enabled, ...
                                                                                                                                           obj.data_dic.dispinfo.subsettrunc, ...
                                                                                                                                           length(obj.current)-imgcorr_prelim(end).idx_ref+i, ...
                                                                                                                                           length(obj.current), ...
                                                                                                                                           get(obj.handles_gui.figure,'OuterPosition'), ...
                                                                                                                                           params_init); %#ok<ASGLU>
                            if (outstate_dic ~= out.success)  
                                % DIC analysis was cancelled or failed by throwing
                                % an exception, for both cases return
                                return;
                            end

                            % Store Displacements
                            displacements_prelim(imgcorr_prelim(end).idx_ref-length(displacements_buffer)+i) = displacement_buffer_update;
                            rois_dic_prelim(imgcorr_prelim(end).idx_ref-length(displacements_buffer)+i) = roi_dic_buffer_update;                                                                                                                    
                        else
                            % Directly store displacements
                            displacements_prelim(imgcorr_prelim(end).idx_ref) = displacements_buffer(end);
                            rois_dic_prelim(imgcorr_prelim(end).idx_ref) = rois_dic_buffer(end);
                        end
                    end

                    % Update imgcorr_prelim - this is necessary since
                    % not all images are necessarily analyzed.
                    if (obj.data_dic.dispinfo.stepanalysis.enabled)
                        imgcorr_prelim(end).idx_cur = imgcorr_prelim(end).idx_ref - length(displacements_buffer);
                    end                        
                end
            end    
            
            % Save the data            
            % See if data will be overwritten
            overwrite = obj.clear_downstream('dicanalysis');
            if (overwrite)
                % Package the data:
                for i = 0:length(obj.current)-1
                    displacements_template(i+1).plot_u_dic = displacements_prelim(i+1).plot_u; %#ok<AGROW>
                    displacements_template(i+1).plot_v_dic = displacements_prelim(i+1).plot_v; %#ok<AGROW>
                    displacements_template(i+1).plot_corrcoef_dic = displacements_prelim(i+1).plot_corrcoef; %#ok<AGROW>
                    displacements_template(i+1).roi_dic = rois_dic_prelim(i+1); %#ok<AGROW>
                    displacements_template(i+1).plot_u_ref_formatted = []; %#ok<AGROW>
                    displacements_template(i+1).plot_v_ref_formatted = []; %#ok<AGROW>
                    displacements_template(i+1).roi_ref_formatted = []; %#ok<AGROW>
                    displacements_template(i+1).plot_u_cur_formatted = []; %#ok<AGROW>
                    displacements_template(i+1).plot_v_cur_formatted = []; %#ok<AGROW>
                    displacements_template(i+1).roi_cur_formatted = []; %#ok<AGROW>
                end
                dispinfo_template.type = obj.data_dic.dispinfo.type;
                dispinfo_template.radius = obj.data_dic.dispinfo.radius; 
                dispinfo_template.spacing = obj.data_dic.dispinfo.spacing;
                dispinfo_template.cutoff_diffnorm = obj.data_dic.dispinfo.cutoff_diffnorm;
                dispinfo_template.cutoff_iteration = obj.data_dic.dispinfo.cutoff_iteration;
                dispinfo_template.total_threads = obj.data_dic.dispinfo.total_threads;
                dispinfo_template.stepanalysis = obj.data_dic.dispinfo.stepanalysis;
                dispinfo_template.subsettrunc = obj.data_dic.dispinfo.subsettrunc;
                dispinfo_template.imgcorr = imgcorr_prelim;
                dispinfo_template.pixtounits = [];
                dispinfo_template.units = '';
                dispinfo_template.cutoff_corrcoef = [];
                dispinfo_template.lenscoef = [];

                % Store the data:           
                % Updated ROIs
                obj.reference.roi = imgs(1).roi;
                for i = 0:length(obj.current)-1
                    obj.current(i+1).roi = imgs(i+2).roi;
                end
                % Displacement info                    
                obj.data_dic.dispinfo(1) = dispinfo_template;
                for i = 0:length(obj.current)-1
                    % Store displacements
                    obj.data_dic.displacements(i+1) = displacements_template(i+1);   
                end

                % Tell user analysis is done
                msgbox('DIC Analysis completed successfully. Press ok to finish.','WindowStyle','modal');
            end    
        end
        
        function callback_topmenu_formatdisp(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for formatting the displacements.
            % Merge images. This makes it convenient for the idx
            % tracker
            imgs = [obj.reference.imginfo obj.current.imginfo];  
            
            % Check if format has been done before
            params_init = cell(0);
            if (~isempty(obj.data_dic.dispinfo.pixtounits) && ...
                ~isempty(obj.data_dic.dispinfo.units) && ...
                ~isempty(obj.data_dic.dispinfo.cutoff_corrcoef) && ...
                ~isempty(obj.data_dic.dispinfo.lenscoef))
                params_init = {obj.data_dic.dispinfo.pixtounits, ...
                               obj.data_dic.dispinfo.units, ...
                               obj.data_dic.dispinfo.cutoff_corrcoef, ... 
                               obj.data_dic.dispinfo.lenscoef};
            end
            
            % Initialize
            plots_disp_f_prelim = struct('plot_u_ref_formatted',{},'plot_v_ref_formatted',{},'plot_u_cur_formatted',{},'plot_v_cur_formatted',{});
            rois_f_ref_prelim = ncorr_class_roi.empty;
            rois_f_cur_prelim = ncorr_class_roi.empty;            
            if (strcmp(obj.data_dic.dispinfo.type,'regular'))
                % Gather reference and current images in the correct order
                imgs_ref = ncorr_class_img.empty;
                imgs_cur = ncorr_class_img.empty;
                for i = 0:length(obj.data_dic.dispinfo.imgcorr)-1
                    imgs_ref = horzcat(imgs_ref,repmat(imgs(obj.data_dic.dispinfo.imgcorr(i+1).idx_ref+1),1,obj.data_dic.dispinfo.imgcorr(i+1).idx_cur-obj.data_dic.dispinfo.imgcorr(i+1).idx_ref)); %#ok<AGROW>
                    imgs_cur = horzcat(imgs_cur,imgs(obj.data_dic.dispinfo.imgcorr(i+1).idx_ref+2:obj.data_dic.dispinfo.imgcorr(i+1).idx_cur+1)); %#ok<AGROW>
                end
                
                % Format step displacements ------------------------------%
                [plots_disp_f_step_prelim,rois_f_step_prelim,pixtounits_prelim,units_prelim,cutoff_corrcoef_prelim,lenscoef_prelim,outstate_formatdisp_step] = ncorr_gui_formatdisp(imgs_ref, ...
                                                                                                                                                                                    imgs_cur, ...
                                                                                                                                                                                    [obj.data_dic.displacements.roi_dic], ...                 
                                                                                                                                                                                    {obj.data_dic.displacements.plot_u_dic}, ...  
                                                                                                                                                                                    {obj.data_dic.displacements.plot_v_dic}, ...      
                                                                                                                                                                                    {obj.data_dic.displacements.plot_corrcoef_dic}, ...  
                                                                                                                                                                                    obj.data_dic.dispinfo.spacing, ...                                                                                                                                                                                
                                                                                                                                                                                    get(obj.handles_gui.figure,'OuterPosition'), ...
                                                                                                                                                                                    params_init);
                % Check if analysis was cancelled
                if (outstate_formatdisp_step ~= out.success)
                    return;
                end
                
                % Convert plots back to pixel displacements
                for i = 0:length(obj.current)-1
                    plots_disp_f_step_prelim(i+1).plot_u_formatted = plots_disp_f_step_prelim(i+1).plot_u_formatted./pixtounits_prelim;
                    plots_disp_f_step_prelim(i+1).plot_v_formatted = plots_disp_f_step_prelim(i+1).plot_v_formatted./pixtounits_prelim;  
                end

                % Add plots ----------------------------------------------%
                % Gather images that dont need to be added
                for i = obj.data_dic.dispinfo.imgcorr(1).idx_ref:obj.data_dic.dispinfo.imgcorr(1).idx_cur-1    
                    plots_disp_f_prelim(i+1).plot_u_ref_formatted = plots_disp_f_step_prelim(i+1).plot_u_formatted;  
                    plots_disp_f_prelim(i+1).plot_v_ref_formatted = plots_disp_f_step_prelim(i+1).plot_v_formatted;
                    rois_f_ref_prelim(i+1) = rois_f_step_prelim(i+1);
                end

                % Reduce reference ROI
                roi_ref_reduced = obj.reference.roi.reduce(obj.data_dic.dispinfo.spacing);
                
                % Add other plots
                plot_added_buffer = cell(0);
                % Cycle over imgcorr - skip first index in imgcorr
                for i = 1:length(obj.data_dic.dispinfo.imgcorr)-1    
                    for j = obj.data_dic.dispinfo.imgcorr(i+1).idx_ref:obj.data_dic.dispinfo.imgcorr(i+1).idx_cur-1                            
                        % Get idxs of displacement fields used to add
                        idx_dispadd = [[obj.data_dic.dispinfo.imgcorr(1:i).idx_cur]-1 j];

                        % Get added plots                             
                        [plot_added_buffer{1},outstate_add] = ncorr_alg_addanalysis({plots_disp_f_step_prelim(idx_dispadd+1).plot_u_formatted}, ...
                                                                                    {plots_disp_f_step_prelim(idx_dispadd+1).plot_v_formatted}, ...
                                                                                    rois_f_step_prelim(idx_dispadd+1), ...
                                                                                    obj.data_dic.dispinfo.spacing, ...
                                                                                    j+1, ...
                                                                                    length(obj.current));
                        % See if analysis was cancelled
                        if (outstate_add ~= out.success)
                            return;
                        end

                        % Store "added" plots
                        plots_disp_f_prelim(j+1).plot_u_ref_formatted = plot_added_buffer{1}.plot_u_added; 
                        plots_disp_f_prelim(j+1).plot_v_ref_formatted = plot_added_buffer{1}.plot_v_added; 
                        rois_f_ref_prelim(j+1) = roi_ref_reduced.get_union(plot_added_buffer{1}.plot_validpoints,0); 
                    end
                end 

                % Convert plots to Eulerian ------------------------------%
                % Get Eulerian displacements
                plot_eulerian_buffer = cell(0);
                roi_eulerian_buffer = cell(0);
                for i = 0:length(obj.current)-1                           
                    % Convert displacements from Lagrangian to Eulerian perspective
                    [plot_eulerian_buffer{1},roi_eulerian_buffer{1},outstate_convert] = ncorr_alg_convertanalysis(obj.current(i+1), ...
                                                                                                                  obj.reference, ...
                                                                                                                  {plots_disp_f_prelim(i+1).plot_u_ref_formatted}, ...
                                                                                                                  {plots_disp_f_prelim(i+1).plot_v_ref_formatted}, ...
                                                                                                                  rois_f_ref_prelim(i+1), ...
                                                                                                                  obj.data_dic.dispinfo.spacing, ...
                                                                                                                  i, ...
                                                                                                                  length(obj.current));

                    % See if analysis was cancelled or failed. Return for
                    % both.
                    if (outstate_convert ~= out.success)
                        return;
                    end

                    % Store
                    plots_disp_f_prelim(i+1).plot_u_cur_formatted = plot_eulerian_buffer{1}.plot_u_new;
                    plots_disp_f_prelim(i+1).plot_v_cur_formatted = plot_eulerian_buffer{1}.plot_v_new;
                    rois_f_cur_prelim(i+1) = roi_eulerian_buffer{1};
                end

                % Make sure to convert displacements back to real units from
                % pixels. Also multiply Eulerian displacements by -1 to 
                % make them truly Eulerian
                for i = 0:length(obj.current)-1
                    % Lagrangian
                    plots_disp_f_prelim(i+1).plot_u_ref_formatted = plots_disp_f_prelim(i+1).plot_u_ref_formatted.*pixtounits_prelim;
                    plots_disp_f_prelim(i+1).plot_v_ref_formatted = plots_disp_f_prelim(i+1).plot_v_ref_formatted.*pixtounits_prelim;  

                    % Eulerian
                    plots_disp_f_prelim(i+1).plot_u_cur_formatted = -plots_disp_f_prelim(i+1).plot_u_cur_formatted.*pixtounits_prelim;
                    plots_disp_f_prelim(i+1).plot_v_cur_formatted = -plots_disp_f_prelim(i+1).plot_v_cur_formatted.*pixtounits_prelim;  
                end
            else
                % Gather reference and current images in the correct order
                imgs_ref = imgs(2:end);
                imgs_cur = ncorr_class_img.empty;
                for i = length(obj.data_dic.dispinfo.imgcorr)-1:-1:0
                    imgs_cur = horzcat(imgs_cur,repmat(imgs(obj.data_dic.dispinfo.imgcorr(i+1).idx_cur+1),1,obj.data_dic.dispinfo.imgcorr(i+1).idx_ref-obj.data_dic.dispinfo.imgcorr(i+1).idx_cur)); %#ok<AGROW>
                end
                
                % Format step displacements ------------------------------%
                [plots_disp_f_step_prelim,rois_f_step_prelim,pixtounits_prelim,units_prelim,cutoff_corrcoef_prelim,lenscoef_prelim,outstate_formatdisp_step] = ncorr_gui_formatdisp(imgs_ref, ...
                                                                                                                                                                                    imgs_cur, ...
                                                                                                                                                                                    [obj.data_dic.displacements.roi_dic], ...                 
                                                                                                                                                                                    {obj.data_dic.displacements.plot_u_dic}, ...  
                                                                                                                                                                                    {obj.data_dic.displacements.plot_v_dic}, ...      
                                                                                                                                                                                    {obj.data_dic.displacements.plot_corrcoef_dic}, ...  
                                                                                                                                                                                    obj.data_dic.dispinfo.spacing, ...                                                                                                                                                                                 
                                                                                                                                                                                    get(obj.handles_gui.figure,'OuterPosition'), ...
                                                                                                                                                                                    params_init);
                % Check if analysis was cancelled
                if (outstate_formatdisp_step ~= out.success)
                    return;
                end
                
                % Convert plots to pixel displacements
                for i = 0:length(obj.current)-1
                    plots_disp_f_step_prelim(i+1).plot_u_formatted = plots_disp_f_step_prelim(i+1).plot_u_formatted./pixtounits_prelim;
                    plots_disp_f_step_prelim(i+1).plot_v_formatted = plots_disp_f_step_prelim(i+1).plot_v_formatted./pixtounits_prelim;  
                end                                      

                % Add plots ------------------------------------------%
                % Gather images that dont need to be added
                for i = obj.data_dic.dispinfo.imgcorr(end).idx_cur:obj.data_dic.dispinfo.imgcorr(end).idx_ref-1   
                    plots_disp_f_prelim(i+1).plot_u_cur_formatted = plots_disp_f_step_prelim(i+1).plot_u_formatted;  
                    plots_disp_f_prelim(i+1).plot_v_cur_formatted = plots_disp_f_step_prelim(i+1).plot_v_formatted;
                    rois_f_cur_prelim(i+1) = rois_f_step_prelim(i+1);
                end

                % Add other plots
                plot_added_buffer = cell(0);
                counter = obj.data_dic.dispinfo.imgcorr(end).idx_ref+1; 
                % Cycle over imgcorr - skip over last index
                for i = 0:length(obj.data_dic.dispinfo.imgcorr)-2    
                    for j = obj.data_dic.dispinfo.imgcorr(i+1).idx_ref-1:-1:obj.data_dic.dispinfo.imgcorr(i+1).idx_cur                            
                        % Get idxs of displacement fields used to add
                        idx_dispadd = [j [obj.data_dic.dispinfo.imgcorr(i+2:end).idx_ref]-1];

                        % Get added plots
                        [plot_added_buffer{1},outstate_add] = ncorr_alg_addanalysis({plots_disp_f_step_prelim(idx_dispadd+1).plot_u_formatted}, ...
                                                                                    {plots_disp_f_step_prelim(idx_dispadd+1).plot_v_formatted}, ...
                                                                                    rois_f_step_prelim(idx_dispadd+1), ...
                                                                                    obj.data_dic.dispinfo.spacing, ...
                                                                                    counter, ...
                                                                                    length(obj.current));
                        % See if analysis was cancelled
                        if (outstate_add ~= out.success)
                            return;
                        end

                        % Reduce current ROI
                        roi_cur_reduced = obj.current(j+1).roi.reduce(obj.data_dic.dispinfo.spacing);

                        % Store "added" plots
                        plots_disp_f_prelim(j+1).plot_u_cur_formatted = plot_added_buffer{1}.plot_u_added; 
                        plots_disp_f_prelim(j+1).plot_v_cur_formatted = plot_added_buffer{1}.plot_v_added; 
                        rois_f_cur_prelim(j+1) = roi_cur_reduced.get_union(plot_added_buffer{1}.plot_validpoints,0); 

                        % Update counter
                        counter = counter+1;
                    end
                end 

                % Convert plots to Lagrangian ----------------------------%
                % Get Lagrangian Displacements
                plot_lagrangian_buffer = cell(0);
                roi_lagrangian_buffer = cell(0);
                [plot_lagrangian_buffer{1},roi_lagrangian_buffer{1},outstate_convert] = ncorr_alg_convertanalysis(obj.reference, ...
                                                                                                                  obj.current, ...
                                                                                                                  {plots_disp_f_prelim.plot_u_cur_formatted}, ...
                                                                                                                  {plots_disp_f_prelim.plot_v_cur_formatted}, ...
                                                                                                                  rois_f_cur_prelim, ...
                                                                                                                  obj.data_dic.dispinfo.spacing, ...
                                                                                                                  0, ...
                                                                                                                  length(obj.current));
                % See if analysis was cancelled or failed. Return for
                % both.
                if (outstate_convert ~= out.success)
                    return;
                end

                % Make sure to convert displacements back to real units from 
                % pixels. Also multiply Eulerian displacements by -1 to 
                % make them truly
                % Eulerian
                for i = 0:length(obj.current)-1
                    % Lagrangian
                    plots_disp_f_prelim(i+1).plot_u_ref_formatted = plot_lagrangian_buffer{1}(i+1).plot_u_new.*pixtounits_prelim;
                    plots_disp_f_prelim(i+1).plot_v_ref_formatted = plot_lagrangian_buffer{1}(i+1).plot_v_new.*pixtounits_prelim;
                    rois_f_ref_prelim(i+1) = roi_lagrangian_buffer{1}(i+1);

                    % Eulerian
                    plots_disp_f_prelim(i+1).plot_u_cur_formatted = -plots_disp_f_prelim(i+1).plot_u_cur_formatted.*pixtounits_prelim;
                    plots_disp_f_prelim(i+1).plot_v_cur_formatted = -plots_disp_f_prelim(i+1).plot_v_cur_formatted.*pixtounits_prelim; 
                end
            end

            % Dont query overwrite, just do it since formatting
            % displacement data and calculating strain maps is quick
            obj.clear_downstream('formatdisp'); 

            % For regular analysis just store formatted rois
            % Package data:
            for i = 0:length(obj.current)-1
                displacements_template(i+1).plot_u_dic = obj.data_dic.displacements(i+1).plot_u_dic; %#ok<AGROW>
                displacements_template(i+1).plot_v_dic = obj.data_dic.displacements(i+1).plot_v_dic; %#ok<AGROW>
                displacements_template(i+1).plot_corrcoef_dic = obj.data_dic.displacements(i+1).plot_corrcoef_dic; %#ok<AGROW>
                displacements_template(i+1).roi_dic = obj.data_dic.displacements(i+1).roi_dic; %#ok<AGROW>
                displacements_template(i+1).plot_u_ref_formatted = plots_disp_f_prelim(i+1).plot_u_ref_formatted; %#ok<AGROW>
                displacements_template(i+1).plot_v_ref_formatted = plots_disp_f_prelim(i+1).plot_v_ref_formatted; %#ok<AGROW>
                displacements_template(i+1).roi_ref_formatted = rois_f_ref_prelim(i+1); %#ok<AGROW>
                displacements_template(i+1).plot_u_cur_formatted = plots_disp_f_prelim(i+1).plot_u_cur_formatted; %#ok<AGROW>
                displacements_template(i+1).plot_v_cur_formatted = plots_disp_f_prelim(i+1).plot_v_cur_formatted; %#ok<AGROW>
                displacements_template(i+1).roi_cur_formatted = rois_f_cur_prelim(i+1); %#ok<AGROW>
            end
            dispinfo_template.type = obj.data_dic.dispinfo.type;
            dispinfo_template.radius = obj.data_dic.dispinfo.radius; 
            dispinfo_template.spacing = obj.data_dic.dispinfo.spacing;
            dispinfo_template.cutoff_diffnorm = obj.data_dic.dispinfo.cutoff_diffnorm;
            dispinfo_template.cutoff_iteration = obj.data_dic.dispinfo.cutoff_iteration;
            dispinfo_template.total_threads = obj.data_dic.dispinfo.total_threads;
            dispinfo_template.stepanalysis = obj.data_dic.dispinfo.stepanalysis;
            dispinfo_template.subsettrunc = obj.data_dic.dispinfo.subsettrunc;
            dispinfo_template.imgcorr = obj.data_dic.dispinfo.imgcorr;
            dispinfo_template.pixtounits = pixtounits_prelim;
            dispinfo_template.units = units_prelim;
            dispinfo_template.cutoff_corrcoef = cutoff_corrcoef_prelim;
            dispinfo_template.lenscoef = lenscoef_prelim;

            % Store data:
            for i = 0:length(obj.current)-1
                obj.data_dic.displacements(i+1) = displacements_template(i+1);   
            end       
            obj.data_dic.dispinfo(1) = dispinfo_template;
        end
       
        function callback_topmenu_calcstrain(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for computing the strains after the
        % displacements have been formatted.
            % Get strain radius ------------------------------------------%
            % Check if strain has been calculated before
            params_init = [];
            if (~isempty(obj.data_dic.straininfo) && ~isempty(obj.data_dic.straininfo.radius) && ~isempty(obj.data_dic.straininfo.subsettrunc))
                params_init(1) = obj.data_dic.straininfo.radius;
                params_init(2) = obj.data_dic.straininfo.subsettrunc;
            end
            
            % Get strain parameters
            [radius_strain_prelim,subsettrunc_prelim,outstate_radius_strain] = ncorr_gui_setstrainradius(obj.reference.imginfo, ...
                                                                                                         [obj.current.imginfo], ... 
                                                                                                         [obj.data_dic.displacements.roi_ref_formatted], ...
                                                                                                         [obj.data_dic.displacements.roi_cur_formatted], ...
                                                                                                         {obj.data_dic.displacements.plot_u_ref_formatted}, ...
                                                                                                         {obj.data_dic.displacements.plot_v_ref_formatted}, ...
                                                                                                         {obj.data_dic.displacements.plot_u_cur_formatted}, ...
                                                                                                         {obj.data_dic.displacements.plot_v_cur_formatted}, ...
                                                                                                         obj.data_dic.dispinfo.spacing, ...
                                                                                                         obj.data_dic.dispinfo.pixtounits, ...
                                                                                                         obj.data_dic.dispinfo.units, ...
                                                                                                         get(obj.handles_gui.figure,'OuterPosition'), ...
                                                                                                         params_init);
                                                                                 
            % See if analysis was cancelled
            if (outstate_radius_strain ~= out.success)
                return;
            end
            
            % Calculate strains --------------------------------------%
            % Calculate displacement gradients through least squares
            % plane fit first 
            % Lagrangian
            plots_dispgrad_ref_prelim = struct('plot_dudx',{},'plot_dudy',{},'plot_dvdx',{},'plot_dvdy',{},'plot_validpoints',{});
            plots_dispgrad_cur_prelim = struct('plot_dudx',{},'plot_dudy',{},'plot_dvdx',{},'plot_dvdy',{},'plot_validpoints',{});

            % Lagrangian
            for i = 0:length(obj.current)-1
                [plots_dispgrad_ref_prelim(i+1),outstate_dispgrad_ref] = ncorr_alg_dispgrad(obj.data_dic.displacements(i+1).plot_u_ref_formatted, ...
                                                                                            obj.data_dic.displacements(i+1).plot_v_ref_formatted, ...
                                                                                            obj.data_dic.displacements(i+1).roi_ref_formatted.formatted(), ...
                                                                                            int32(radius_strain_prelim), ...
                                                                                            obj.data_dic.dispinfo.pixtounits, ...
                                                                                            int32(obj.data_dic.dispinfo.spacing), ...
                                                                                            logical(subsettrunc_prelim), ...
                                                                                            int32(i), ...
                                                                                            int32(length(obj.current)));
                % See if analysis was cancelled
                if (outstate_dispgrad_ref ~= out.success)
                    return;
                end 
            end  

            % Eulerian
            for i = 0:length(obj.current)-1
                % Eulerian
                [plots_dispgrad_cur_prelim(i+1),outstate_dispgrad_cur] = ncorr_alg_dispgrad(obj.data_dic.displacements(i+1).plot_u_cur_formatted, ...
                                                                                            obj.data_dic.displacements(i+1).plot_v_cur_formatted, ...
                                                                                            obj.data_dic.displacements(i+1).roi_cur_formatted.formatted(), ...
                                                                                            int32(radius_strain_prelim), ...
                                                                                            obj.data_dic.dispinfo.pixtounits, ...
                                                                                            int32(obj.data_dic.dispinfo.spacing), ...
                                                                                            logical(subsettrunc_prelim), ...
                                                                                            int32(i), ...
                                                                                            int32(length(obj.current)));
                % See if analysis was cancelled
                if (outstate_dispgrad_cur ~= out.success)
                    return;
                end 
            end     
                
            % Store
            % Dont query overwrite, just do it since calculating strain maps is quick      
            obj.clear_downstream('strains');

            % Package info:
            for i = 0:length(obj.current)-1
                % Green Lagrangian
                strains_template(i+1).plot_exx_ref_formatted = (1/2)*(2*plots_dispgrad_ref_prelim(i+1).plot_dudx+plots_dispgrad_ref_prelim(i+1).plot_dudx.^2+plots_dispgrad_ref_prelim(i+1).plot_dvdx.^2); %#ok<AGROW>
                strains_template(i+1).plot_exy_ref_formatted = (1/2)*(plots_dispgrad_ref_prelim(i+1).plot_dudy+plots_dispgrad_ref_prelim(i+1).plot_dvdx+plots_dispgrad_ref_prelim(i+1).plot_dudx.*plots_dispgrad_ref_prelim(i+1).plot_dudy+plots_dispgrad_ref_prelim(i+1).plot_dvdx.*plots_dispgrad_ref_prelim(i+1).plot_dvdy); %#ok<AGROW>
                strains_template(i+1).plot_eyy_ref_formatted = (1/2)*(2*plots_dispgrad_ref_prelim(i+1).plot_dvdy+plots_dispgrad_ref_prelim(i+1).plot_dudy.^2+plots_dispgrad_ref_prelim(i+1).plot_dvdy.^2); %#ok<AGROW>
                strains_template(i+1).roi_ref_formatted = obj.data_dic.displacements(i+1).roi_ref_formatted.get_union(plots_dispgrad_ref_prelim(i+1).plot_validpoints,0); %#ok<AGROW>

                % Eulerian-Almansi
                strains_template(i+1).plot_exx_cur_formatted = (1/2)*(2*plots_dispgrad_cur_prelim(i+1).plot_dudx-plots_dispgrad_cur_prelim(i+1).plot_dudx.^2-plots_dispgrad_cur_prelim(i+1).plot_dvdx.^2); %#ok<AGROW>
                strains_template(i+1).plot_exy_cur_formatted = (1/2)*(plots_dispgrad_cur_prelim(i+1).plot_dudy+plots_dispgrad_cur_prelim(i+1).plot_dvdx-plots_dispgrad_cur_prelim(i+1).plot_dudx.*plots_dispgrad_cur_prelim(i+1).plot_dudy-plots_dispgrad_cur_prelim(i+1).plot_dvdx.*plots_dispgrad_cur_prelim(i+1).plot_dvdy); %#ok<AGROW>
                strains_template(i+1).plot_eyy_cur_formatted = (1/2)*(2*plots_dispgrad_cur_prelim(i+1).plot_dvdy-plots_dispgrad_cur_prelim(i+1).plot_dudy.^2-plots_dispgrad_cur_prelim(i+1).plot_dvdy.^2); %#ok<AGROW>
                strains_template(i+1).roi_cur_formatted = obj.data_dic.displacements(i+1).roi_cur_formatted.get_union(plots_dispgrad_cur_prelim(i+1).plot_validpoints,0); %#ok<AGROW>
            end

            % Store info:
            for i = 0:length(obj.current)-1
                obj.data_dic.strains(i+1) = strains_template(i+1);
            end
            obj.data_dic.straininfo(1).radius = radius_strain_prelim;
            obj.data_dic.straininfo(1).subsettrunc = subsettrunc_prelim;

            % Tell user calculating strains was successful
            msgbox('Strains calculation complete. Press ok to finish.','WindowStyle','modal');
        end
        
        function callback_topmenu_viewplot(obj,hObject,eventdata,type_fig) %#ok<INUSL>
        % This is the callback for viewing the plots after the
        % displacements have been formatted and after the strains
        % have been computed.
            % Get Data
            handles_plot = getappdata(obj.handles_gui.figure,'handles_plot');   
            
            % Deactivate closerequestfcn - this prevents any closing of
            % figures when handles_plot is being updated
            for i = 0:length(handles_plot)-1
                 set(handles_plot(i+1),'CloseRequestFcn','');
            end      
            
            % Cycle over figures (it's possible type fig has more than one
            % type).
            for i = 0:length(type_fig)-1
                % Check if plot already exists
                plotexists = false;
                for j = 0:length(handles_plot)-1
                    if (strcmp(getappdata(handles_plot(j+1),'type_plot'),type_fig{i+1}))
                        % Plot already exists, bring it back to focus
                        figure(handles_plot(j+1));
                        plotexists = true;
                        break;
                    end
                end
                
                if (~plotexists)
                    % Get params_init;
                    params_init = cell(0); 
                    if (~isempty(handles_plot))
                        % All windows will have the same paramters so just
                        % grab the info from the first window
                        params_init{1} = getappdata(handles_plot(1),'num_cur');                                 % {1} = num_cur
                        % Get data
                        params_init{2} = getappdata(handles_plot(1),'scalebarlength');                          % {2} = scalebarlength
                        params_init{3} = getappdata(handles_plot(1),'val_checkbox_scalebar');                   % {3} = scalebar;
                        params_init{4} = getappdata(handles_plot(1),'val_checkbox_axes');                       % {4} = axes;
                        params_init{5} = getappdata(handles_plot(1),'val_popupmenu');                           % {5} = lagrangian or eulerian
                        params_init{6} = strcmp(get(getappdata(handles_plot(1),'handle_zoom'),'Enable'),'on');  % {6} = zoomed
                        params_init{7} = strcmp(get(getappdata(handles_plot(1),'handle_pan'),'Enable'),'on');   % {7} = panned
                    end
                    
                    % Go over friends list; find the right most figure and
                    % make this the figure whose coordinates are used for 
                    % pos_parent.
                    if (isempty(handles_plot))
                        % Use coordinates of ncorr
                        pos_parent = get(obj.handles_gui.figure,'OuterPosition');
                    else
                        % Plots already exist, fight right most plot
                        pos_parent = get(handles_plot(1),'OuterPosition');
                        for j = 1:length(handles_plot)-1
                            pos_buf =  get(handles_plot(j+1),'OuterPosition');
                            if (pos_buf(1) + pos_buf(3) > pos_parent(1) + pos_parent(3))
                                pos_parent = pos_buf;
                            elseif (pos_buf(1) + pos_buf(3) == pos_parent(1) + pos_parent(3))
                                % If they are equal, get the one that is
                                % lower
                                if (pos_buf(2) < pos_parent(2))
                                    pos_parent = pos_buf;
                                end
                            end
                        end
                    end
                    
                    % Create new plot - This function will return all
                    % gui handles. Store these as "friends" so plots can 
                    % remain consistent. The figure has already has its 
                    % closereqestfcn disabled
                    handles_gui_new = ncorr_gui_viewplots(obj.reference.imginfo, ...
                                                          [obj.current.imginfo], ...
                                                          obj.data_dic, ...
                                                          type_fig{i+1}, ...
                                                          pos_parent, ...
                                                          params_init);
                    
                    % Figure can be invalid if the ROI is empty
                    if (ishandle(handles_gui_new.figure))
                        % Set delete functions
                        set(handles_gui_new.figure,'DeleteFcn',ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata)callback_update_handles_plot(obj,hObject,eventdata),obj.handles_gui.figure));
                        
                        % Update friends list - this allows the things to
                        % sync
                        if (isempty(handles_plot))
                            % Check if a plot already exists - if it does then get
                            % the old friends list from the first plot
                            friends_old = [];
                        else
                            friends_old = getappdata(handles_plot(1),'friends');
                        end
                        friends_new = horzcat(friends_old,{handles_gui_new});
                        
                        % Link the axes
                        vec_friends = [friends_new{:}];
                        linkaxes([vec_friends.axes_formatplot],'xy');
                        
                        % Set Data
                        % First update handles_gui_new
                        setappdata(handles_gui_new.figure,'friends',friends_new);                        
                        % Now update the preexisting handles
                        for j = 0:length(handles_plot)-1
                            setappdata(handles_plot(j+1),'friends',friends_new);
                        end
                        
                        % Now store figure handle
                        handles_plot = horzcat(handles_plot,handles_gui_new.figure); %#ok<AGROW>
                    end
                end
            end      
            
            % Reactivate closerequestfcns
            for i = 0:length(handles_plot)-1
                 set(handles_plot(i+1),'CloseRequestFcn','closereq');
            end 
            
            % Store Data
            setappdata(obj.handles_gui.figure,'handles_plot',handles_plot);
        end        
        
        function callback_topmenu_closeplots(obj,hObject,eventdata) %#ok<INUSD>
        % This callback closes all the open plots.
            obj.close_dicplots('all'); 
        end      
                        
        %-----------------------------------------------------------------%
        % Other Callbacks ------------------------------------------------%
        %-----------------------------------------------------------------%               
               
        function callback_update_title(obj,hObject,eventdata) %#ok<INUSD>
        % This function will update the title of ncorr to display the
        % handle pointing to ncorr
            % Get Name
            name = util_get_name(obj);
            if (isempty(name))
                % No handle exists
                set(obj.handles_gui.figure,'Name','Ncorr'); 
            else
                % Handle exists
                set(obj.handles_gui.figure,'Name',['Ncorr - ' name]); 
            end
        end
               
        function callback_edit_imgnum(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for setting the image number directly
            % Get data
            num_cur = getappdata(obj.handles_gui.figure,'num_cur');
            
            % Get img num
            num_cur_prelim = str2double(get(obj.handles_gui.edit_imgnum,'string')); 
            if (ncorr_util_isintbb(num_cur_prelim,1,length(obj.current),'Current image number') == out.success)  
                % convert num_cur_prelim to zero based indexed
                num_cur = num_cur_prelim-1;
            end    
            
            % Set data
            setappdata(obj.handles_gui.figure,'num_cur',num_cur); 
        end
        
        function callback_button_left(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for the left button to display different
        % current images
            % Get data
            num_cur = getappdata(obj.handles_gui.figure,'num_cur');
            
            % Check for overshoot
            if (num_cur > 0)                
                num_cur = num_cur-1;
            end    
            
            % Set data
            setappdata(obj.handles_gui.figure,'num_cur',num_cur); 
        end
        
        function callback_button_right(obj,hObject,eventdata) %#ok<INUSD>
        % This is the callback for the right button to display different
        % current images
            % Get data
            num_cur = getappdata(obj.handles_gui.figure,'num_cur');  
            
            % Check for overshoot
            if (num_cur < length(obj.current)-1)         
                num_cur = num_cur+1;
            end      
            
            % Set data
            setappdata(obj.handles_gui.figure,'num_cur',num_cur);
        end
        
        function callback_update_handles_plot(obj,hObject,eventdata) %#ok<INUSD>
        % This function is called when a plot figure is closed; this
        % function removes the figure handle from handles_plot.
            % Get data
            handles_plot = getappdata(obj.handles_gui.figure,'handles_plot');
            
            % Disable closerequestfcns - ensure other plots aren't closed
            % while the handle is being updated.
            for i = 0:length(handles_plot)-1
                set(handles_plot(i+1),'CloseRequestFcn','');
            end
            
            % Determine index of closing figure
            closingfig = gcbf;              
            idx_closingfig = -1;
            for i = 0:length(handles_plot)-1
                if (handles_plot(i+1) == closingfig) % update handle array
                    idx_closingfig = i;
                    break;
                end
            end
            
            % Delete closing figure from handles_plot
            handles_plot(idx_closingfig+1) = [];     
            
            % Update Friends            
            for i = 0:length(handles_plot)-1
                % Get friends list for each plot handle. Delete the friend
                % that contains closing fig
                % Get data
                friends = getappdata(handles_plot(i+1),'friends');
                
                for j = 0:length(friends)-1                  
                    if (friends{j+1}.figure == closingfig)
                        friends(j+1) = [];
                        break;
                    end                    
                end
                
                % Set data
                setappdata(handles_plot(i+1),'friends',friends);
            end
            
            % Reenable CloseRequestfFcns
            for i = 0:length(handles_plot)-1
               set(handles_plot(i+1),'CloseRequestFcn','closereq');
            end
            
            % Set data
            setappdata(obj.handles_gui.figure,'handles_plot',handles_plot);   
        end 
        
        function callback_close_function(obj,hObject,eventdata,force) %#ok<INUSL>
            % Force will cause a force close if it's true - this situation
            % occurs when the user manually deletes the GUI handle
            closencorr = false;
            if (~force && ~isempty(obj.data_dic.dispinfo) && ~isempty(obj.data_dic.displacements))
%                 contbutton = questdlg('Prior DIC data has been detected. If you continue without saving, this data will be deleted. Do you want to proceed?','Continue Operation','Yes','No','Yes');
                contbutton = 'Yes'; % MultiDIC and DuoDIC modification to disable this warning
                if (strcmp(contbutton,'Yes'))
                    closencorr = true;
                end
            else
                % Force close if user calls "delete(handles_ncorr)" or if
                % dic analysis has not been run yet.
                closencorr = true;
            end   
            
            if (closencorr)
                % Close plot figures
                obj.close_dicplots('all'); 
                
                % Delete timer
                timer_ncorr = getappdata(obj.handles_gui.figure,'timer');
                if (isa(timer_ncorr,'timer') && isvalid(timer_ncorr))
                    stop(timer_ncorr);
                    delete(timer_ncorr);
                end
                
                % Close GUI
                delete(obj.handles_gui.figure);
                
                % Delete data
                delete(obj); 
            end
        end 
        
        %-----------------------------------------------------------------%
        % Other Functions ------------------------------------------------%
        %-----------------------------------------------------------------%
        
        function close_dicplots(obj,closetype)
        % This is the function for closing either the displacement plots, 
        % strain plots, or all the plots.                 
            
            % Displacement plots
            if (strcmp(closetype,'displacements') || strcmp(closetype,'all'))                  
                % Get data
                handles_plot = getappdata(obj.handles_gui.figure,'handles_plot');  
                
                handles_plot_close = [];
                for i = 0:length(handles_plot)-1
                    type_plot = getappdata(handles_plot(i+1),'type_plot');  
                    if (strcmp(type_plot,'u') || strcmp(type_plot,'v'))
                        handles_plot_close = horzcat(handles_plot_close,handles_plot(i+1)); %#ok<AGROW>
                    end                    
                end
                
                % Close figures
                close(handles_plot_close); % handles_plot updates after this
            end

            % Strain plots
            if (strcmp(closetype,'strains') || strcmp(closetype,'all'))  
                % Get Data
                handles_plot = getappdata(obj.handles_gui.figure,'handles_plot');  
                
                handles_plot_close = [];
                for i = 0:length(handles_plot)-1
                    type_plot = getappdata(handles_plot(i+1),'type_plot');  
                    if (strcmp(type_plot,'exx') || strcmp(type_plot,'exy') || strcmp(type_plot,'eyy'))
                        handles_plot_close = horzcat(handles_plot_close,handles_plot(i+1)); %#ok<AGROW>
                    end                    
                end
                
                % Close figures
                close(handles_plot_close); % handles_plot updates after this 
            end
        end       
        
        function overwrite = clear_downstream(obj,action) 
        % Middle of callback - alert user if data is going to be 
        % overwritten and then clear downstream data.
            % Initialize - set default to true.
            overwrite = true;             
            
            % See whether to query overwrite
            if (strcmp(action,'all') || strcmp(action,'set_ref') || strcmp(action,'set_cur') || strcmp(action,'set_roi') || strcmp(action,'set_dicparameters') || strcmp(action,'dicanalysis'))      
                % Only query if displacements have been set, since this is
                % usually the most time consuming process
                if (~isempty(obj.data_dic.displacements))
                    contbutton = questdlg('Prior DIC data has been detected. If you continue without saving, this data will be deleted. Do you want to proceed?','Continue Operation','Yes','No','Yes');
                    if (~strcmp(contbutton,'Yes'))
                        overwrite = false;
                    end                                
                end
            end
           
            % If overwrite is permitted, then clear downstream data.
            if (overwrite)
                % Get field names since they are used a lot
                fields_data_dic = fieldnames(obj.data_dic); 
                fields_dic_data_dispinfo = fieldnames(obj.data_dic.dispinfo); 
                
                % Handle each step accordingly
                if (strcmp(action,'all'))
                    % Close all DIC figs if they exist -------------------%
                    obj.close_dicplots('all'); 
                    % Clear downstream data ------------------------------%
                    % All info needs to be deleted
                    obj.reference(:) = [];
                    obj.current(:) = [];
                    setappdata(obj.handles_gui.figure,'num_cur',[]); 
                    for i = 0:length(fields_data_dic)-1
                        obj.data_dic.(fields_data_dic{i+1})(:) = [];
                    end  
                elseif (strcmp(action,'set_ref'))   
                    % Close all DIC figs if they exist -------------------%
                    obj.close_dicplots('all'); 
                    % Clear downstream data ------------------------------%
                    % Clear all reference info
                    obj.reference(:) = [];   
                    % For current images, clear all ROIs if analysis
                    % was regular. If it was backward, then clear all
                    % ROIs except for the last one
                    if (~isempty(obj.current) && ~isempty(obj.data_dic.dispinfo))                            
                        if (strcmp(obj.data_dic.dispinfo.type,'regular'))
                            for i = 0:length(obj.current)-1
                                obj.current(i+1).roi = ncorr_class_roi.empty;
                            end
                        else
                            for i = 0:length(obj.current)-2
                                obj.current(i+1).roi = ncorr_class_roi.empty;
                            end
                        end
                    end
                    % Clear everything in data_dic except for dispinfo
                    for i = 0:length(fields_data_dic)-1
                        if (~strcmp(fields_data_dic{i+1},'dispinfo'))
                            obj.data_dic.(fields_data_dic{i+1})(:) = [];
                        end                                    
                    end
                    % Delete everything except for the type in
                    % dispinfo. If the type was regular, then delete 
                    % that too                        
                    if (~isempty(obj.data_dic.dispinfo))
                        if (strcmp(obj.data_dic.dispinfo.type,'regular'))
                            obj.data_dic.dispinfo(:) = [];
                        else
                            for i = 0:length(fields_dic_data_dispinfo)-1
                                if (~strcmp(fields_dic_data_dispinfo{i+1},'type'))                                       
                                    obj.data_dic.dispinfo.(fields_dic_data_dispinfo{i+1})(:) = [];
                                end                                    
                            end
                        end
                    end
                elseif (strcmp(action,'set_cur'))
                    % Close all DIC figs if they exist -------------------%
                    obj.close_dicplots('all'); 
                    % Clear downstream data ------------------------------%
                    % Delete reference ROI if type is backward
                    if (~isempty(obj.reference) && ~isempty(obj.data_dic.dispinfo))
                        if (strcmp(obj.data_dic.dispinfo.type,'backward'))
                            obj.reference.roi = ncorr_class_roi.empty;
                        end
                    end
                    % Clear all current images
                    obj.current(:) = [];
                    setappdata(obj.handles_gui.figure,'num_cur',[]); 
                    % Clear everything in data_dic except for dispinfo
                    for i = 0:length(fields_data_dic)-1
                        if (~strcmp(fields_data_dic{i+1},'dispinfo'))
                            obj.data_dic.(fields_data_dic{i+1})(:) = [];
                        end                                    
                    end
                    % Delete everything except for the type in
                    % dispinfo. If the type was backward, then delete
                    % that too
                    if (~isempty(obj.data_dic.dispinfo))
                        if (strcmp(obj.data_dic.dispinfo.type,'backward'))
                            obj.data_dic.dispinfo(:) = [];
                        else
                            for i = 0:length(fields_dic_data_dispinfo)-1
                                if (~strcmp(fields_dic_data_dispinfo{i+1},'type'))                                       
                                    obj.data_dic.dispinfo.(fields_dic_data_dispinfo{i+1})(:) = [];
                                end                                    
                            end
                        end
                    end
                elseif (strcmp(action,'set_roi'))
                    % Close all DIC figs if they exist -------------------%
                    obj.close_dicplots('all'); 
                    % Clear downstream data ------------------------------%
                    % Clear all ROIs in reference image
                    if (~isempty(obj.reference))
                        obj.reference.roi = ncorr_class_roi.empty;
                    end
                    % Clear all ROIS in current image
                    if (~isempty(obj.current))
                        for i = 0:length(obj.current)-1
                            obj.current(i+1).roi = ncorr_class_roi.empty;
                        end
                    end
                    % Clear all data_dic
                    for i = 0:length(fields_data_dic)-1
                        obj.data_dic.(fields_data_dic{i+1})(:) = [];
                    end  
                 elseif (strcmp(action,'set_dicparameters'))
                    % Close DIC figs if they exist -----------------------%
                    obj.close_dicplots('all'); 
                    % Clear downstream data ------------------------------%
                    % If type is regular, delete all current ROIS. 
                    % If type is backward, delete all ROIS except for
                    % last current ROI
                    if (strcmp(obj.data_dic.dispinfo.type,'regular'))
                        for i = 0:length(obj.current)-1
                            obj.current(i+1).roi = ncorr_class_roi.empty;
                        end
                    else
                        obj.reference.roi = ncorr_class_roi.empty;
                        for i = 0:length(obj.current)-2
                            obj.current(i+1).roi = ncorr_class_roi.empty;
                        end
                    end
                    % Clear everything in data_dic except for dispinfo
                    for i = 0:length(fields_data_dic)-1
                        if (~strcmp(fields_data_dic{i+1},'dispinfo'))
                            obj.data_dic.(fields_data_dic{i+1})(:) = [];
                        end                                    
                    end
                    % Delete everything in dispinfo except for the type
                    fields_dic_data_dispinfo = fieldnames(obj.data_dic.dispinfo);
                    for i = 0:length(fields_dic_data_dispinfo)-1
                        if (~strcmp(fields_dic_data_dispinfo{i+1},'type'))                                       
                            obj.data_dic.dispinfo.(fields_dic_data_dispinfo{i+1})(:) = [];
                        end                                    
                    end
                 elseif (strcmp(action,'dicanalysis'))
                    % Close all DIC figs if they exist -------------------%
                    obj.close_dicplots('all'); 
                    % Clear downstream data ------------------------------%
                    % If type is regular, delete all current ROIs. 
                    % If type is backward, delete all ROIs except for
                    % last current ROI                        
                    if (strcmp(obj.data_dic.dispinfo.type,'regular'))
                        for i = 0:length(obj.current)-1
                            obj.current(i+1).roi = ncorr_class_roi.empty;
                        end
                    else
                        obj.reference.roi = ncorr_class_roi.empty;
                        for i = 0:length(obj.current)-2
                            obj.current(i+1).roi = ncorr_class_roi.empty;
                        end
                    end
                    % Delete all info except dispinfo
                    for i = 0:length(fields_data_dic)-1                            
                        if (~strcmp(fields_data_dic{i+1},'dispinfo'))
                            obj.data_dic.(fields_data_dic{i+1})(:) = [];
                        end
                    end
                    % Delete downstream data in data_dic.dispinfo:
                    fields_dic_data_dispinfo = fieldnames(obj.data_dic.dispinfo);
                    for i = 0:length(fields_dic_data_dispinfo)-1
                        if (~strcmp(fields_dic_data_dispinfo{i+1},'type') && ...
                            ~strcmp(fields_dic_data_dispinfo{i+1},'radius') && ...
                            ~strcmp(fields_dic_data_dispinfo{i+1},'spacing') && ...                               
                            ~strcmp(fields_dic_data_dispinfo{i+1},'cutoff_diffnorm') && ...
                            ~strcmp(fields_dic_data_dispinfo{i+1},'cutoff_iteration') && ...
                            ~strcmp(fields_dic_data_dispinfo{i+1},'total_threads') && ...     
                            ~strcmp(fields_dic_data_dispinfo{i+1},'stepanalysis') && ...
                            ~strcmp(fields_dic_data_dispinfo{i+1},'subsettrunc'))                       
                            obj.data_dic.dispinfo.(fields_dic_data_dispinfo{i+1})(:) = [];
                        end                                    
                    end
                elseif (strcmp(action,'formatdisp'))
                    % Close DIC figs if they exist -----------------------%
                    obj.close_dicplots('all'); 
                    % Clear downstream data ------------------------------%
                    % Delete all info except for displacements and
                    % dispinfo;
                    for i = 0:length(fields_data_dic)-1
                        if (~strcmp(fields_data_dic{i+1},'displacements') && ...
                            ~strcmp(fields_data_dic{i+1},'dispinfo'))
                            obj.data_dic.(fields_data_dic{i+1})(:) = [];
                        end                                    
                    end
                    % Delete downstream data in data_dic.dispinfo:
                    fields_dic_data_dispinfo = fieldnames(obj.data_dic.dispinfo);
                    for i = 0:length(fields_dic_data_dispinfo)-1
                        if (~strcmp(fields_dic_data_dispinfo{i+1},'type') && ...
                            ~strcmp(fields_dic_data_dispinfo{i+1},'radius') && ...
                            ~strcmp(fields_dic_data_dispinfo{i+1},'spacing') && ...                               
                            ~strcmp(fields_dic_data_dispinfo{i+1},'cutoff_diffnorm') && ...
                            ~strcmp(fields_dic_data_dispinfo{i+1},'cutoff_iteration') && ...
                            ~strcmp(fields_dic_data_dispinfo{i+1},'total_threads') && ...     
                            ~strcmp(fields_dic_data_dispinfo{i+1},'stepanalysis') && ...      
                            ~strcmp(fields_dic_data_dispinfo{i+1},'subsettrunc') && ...
                            ~strcmp(fields_dic_data_dispinfo{i+1},'imgcorr'))                       
                            obj.data_dic.dispinfo.(fields_dic_data_dispinfo{i+1})(:) = [];
                        end                                    
                    end
                    % Delete downstream data in data_dic.displacements:
                    dic_data_disp_fields = fieldnames(obj.data_dic.displacements);
                    for i = 0:length(dic_data_disp_fields)-1
                        if (~strcmp(dic_data_disp_fields{i+1},'plot_u_dic') && ...
                            ~strcmp(dic_data_disp_fields{i+1},'plot_v_dic') && ...
                            ~strcmp(dic_data_disp_fields{i+1},'plot_corrcoef_dic') && ...
                            ~strcmp(dic_data_disp_fields{i+1},'roi_dic'))                        
                            for j = 0:length(obj.data_dic.displacements)-1
                                obj.data_dic.displacements(j+1).(dic_data_disp_fields{i+1})(:) = [];
                            end
                        end                                    
                    end
                elseif (strcmp(action,'strains'))
                    % Close strain figs if they exist --------------------%
                    obj.close_dicplots('strains'); 
                    % Clear downstream data ------------------------------%                    
                    % Delete all info except dispinfo and displacements
                    for i = 0:length(fields_data_dic)-1                            
                        if (~strcmp(fields_data_dic{i+1},'dispinfo') && ...
                            ~strcmp(fields_data_dic{i+1},'displacements'))
                            obj.data_dic.(fields_data_dic{i+1})(:) = [];
                        end
                    end
                end
            end
        end        
    
        function freeze_menu(obj)
        % This function disables the menu when callbacks are executing
            set(obj.handles_gui.topmenu_file,'Enable','off'); 
            set(obj.handles_gui.topmenu_regionofinterest,'Enable','off'); 
            set(obj.handles_gui.topmenu_analysis,'Enable','off'); 
            set(obj.handles_gui.topmenu_plot,'Enable','off'); 
        end

        function unfreeze_menu(obj)
        % This function enables the menu when callbacks are done
            set(obj.handles_gui.topmenu_file,'Enable','on'); 
            set(obj.handles_gui.topmenu_regionofinterest,'Enable','on'); 
            set(obj.handles_gui.topmenu_analysis,'Enable','on'); 
            set(obj.handles_gui.topmenu_plot,'Enable','on'); 
        end

        function update_topmenu(obj)            
        % This function updates menu items depending on what is loaded
            % Reference image is loaded:
            if (~isempty(obj.reference))            
                set(obj.handles_gui.topmenu_setroi_ref,'Enable','on'); 
            else
                set(obj.handles_gui.topmenu_setroi_ref,'Enable','off'); 
            end
            
            % Current image is loaded:
            if (~isempty(obj.current))            
                set(obj.handles_gui.topmenu_setroi_cur,'Enable','on'); 
            else
                set(obj.handles_gui.topmenu_setroi_cur,'Enable','off'); 
            end
            
            % Reference image, current image, and ROI are loaded:
            % Must check if analysis is regular or backward because ROIs
            % can be set through the DIC analysis
            if (~isempty(obj.reference) && ~isempty(obj.current) && ~isempty(obj.data_dic.dispinfo) && ...
                ((strcmp(obj.data_dic.dispinfo.type,'regular') && ~isempty(obj.reference.roi)) || ...
                ((strcmp(obj.data_dic.dispinfo.type,'backward') && ~isempty(obj.current(end).roi)))))
                set(obj.handles_gui.topmenu_setdicparameters,'Enable','on');
            else
                set(obj.handles_gui.topmenu_setdicparameters,'Enable','off');
            end
            
            % Reference image, current image, and DIC parameters are loaded:
            if (~isempty(obj.reference) && ~isempty(obj.current) && ~isempty(obj.data_dic.dispinfo) && ...
                ~isempty(obj.data_dic.dispinfo.radius) && ...
                ~isempty(obj.data_dic.dispinfo.spacing) && ...
                ~isempty(obj.data_dic.dispinfo.cutoff_diffnorm) && ...
                ~isempty(obj.data_dic.dispinfo.cutoff_iteration) && ...
                ~isempty(obj.data_dic.dispinfo.total_threads) && ...
                ~isempty(obj.data_dic.dispinfo.stepanalysis) &&  ...
                ~isempty(obj.data_dic.dispinfo.subsettrunc))
                set(obj.handles_gui.topmenu_perfdic,'Enable','on');
            else
                set(obj.handles_gui.topmenu_perfdic,'Enable','off');
            end   
                                   
            % Reference image, current image, and DIC analysis are loaded:
            if (~isempty(obj.reference) && ~isempty(obj.current) && ~isempty(obj.data_dic.dispinfo) && ~isempty(obj.data_dic.displacements))
                set(obj.handles_gui.topmenu_savedata,'Enable','on');
                set(obj.handles_gui.topmenu_formatdisp,'Enable','on');  
            else
                set(obj.handles_gui.topmenu_savedata,'Enable','off');
                set(obj.handles_gui.topmenu_formatdisp,'Enable','off');  
            end   
            
            % Reference image, current image, DIC data, and formatted 
            % displacements are loaded:
            if (~isempty(obj.reference) && ~isempty(obj.current) && ~isempty(obj.data_dic.dispinfo) && ~isempty(obj.data_dic.displacements) && ...
                ~isempty(obj.data_dic.dispinfo.pixtounits) && ...
                ~isempty(obj.data_dic.dispinfo.units) && ...
                ~isempty(obj.data_dic.dispinfo.cutoff_corrcoef) && ...
                ~isempty(obj.data_dic.dispinfo.lenscoef))
                set(obj.handles_gui.topmenu_calcstrain,'Enable','on');
                set(obj.handles_gui.topmenu_viewdispplot_u,'Enable','on');  
                set(obj.handles_gui.topmenu_viewdispplot_v,'Enable','on');  
                set(obj.handles_gui.topmenu_viewdispplot_all,'Enable','on'); 
                set(obj.handles_gui.topmenu_closeplots,'Enable','on');   
            else
                set(obj.handles_gui.topmenu_calcstrain,'Enable','off');
                set(obj.handles_gui.topmenu_viewdispplot_u,'Enable','off'); 
                set(obj.handles_gui.topmenu_viewdispplot_v,'Enable','off');  
                set(obj.handles_gui.topmenu_viewdispplot_all,'Enable','off');    
                set(obj.handles_gui.topmenu_closeplots,'Enable','off');   
            end
                
            % Reference image, current image, DIC data, and strains are
            % loaded
            if (~isempty(obj.reference) && ~isempty(obj.current) && ~isempty(obj.data_dic.dispinfo)&& ~isempty(obj.data_dic.displacements) &&  ...
                ~isempty(obj.data_dic.straininfo) && ~isempty(obj.data_dic.strains))
                set(obj.handles_gui.topmenu_viewstrainplot_exx,'Enable','on');  
                set(obj.handles_gui.topmenu_viewstrainplot_exy,'Enable','on');  
                set(obj.handles_gui.topmenu_viewstrainplot_eyy,'Enable','on'); 
                set(obj.handles_gui.topmenu_viewstrainplot_all,'Enable','on');  
            else
                set(obj.handles_gui.topmenu_viewstrainplot_exx,'Enable','off');  
                set(obj.handles_gui.topmenu_viewstrainplot_exy,'Enable','off');  
                set(obj.handles_gui.topmenu_viewstrainplot_eyy,'Enable','off');  
                set(obj.handles_gui.topmenu_viewstrainplot_all,'Enable','off');  
            end
        end
        
        function update_dicstatetext(obj) 
        % This function updates state text depending on what is loaded
            % Reference image:
            if (~isempty(obj.reference))
                set(obj.handles_gui.text_refloaded_s,'String','SET','ForegroundColor','g','FontWeight','bold');            
            else
                set(obj.handles_gui.text_refloaded_s,'String','NOT SET','ForegroundColor','r','FontWeight','normal');
            end
            
            % Current image:
            if (~isempty(obj.current))
                set(obj.handles_gui.text_curloaded_s,'String','SET','ForegroundColor','g','FontWeight','bold');           
            else
                set(obj.handles_gui.text_curloaded_s,'String','NOT SET','ForegroundColor','r','FontWeight','normal');
            end
            
            % ROI: Must check if analysis is regular or backward because ROIs
            % can be set through the DIC analysis
            if (~isempty(obj.data_dic.dispinfo)&& ...
                ((strcmp(obj.data_dic.dispinfo.type,'regular') && ~isempty(obj.reference.roi)) || ...
                ((strcmp(obj.data_dic.dispinfo.type,'backward') && ~isempty(obj.current(end).roi)))))
                set(obj.handles_gui.text_roiloaded_s,'String','SET','ForegroundColor','g','FontWeight','bold');       
            else
                set(obj.handles_gui.text_roiloaded_s,'String','NOT SET','ForegroundColor','r','FontWeight','normal');
            end
            
            % DIC parameters:
            if (~isempty(obj.data_dic.dispinfo) && ...
                ~isempty(obj.data_dic.dispinfo.radius) && ...
                ~isempty(obj.data_dic.dispinfo.spacing) && ...
                ~isempty(obj.data_dic.dispinfo.cutoff_diffnorm) && ...
                ~isempty(obj.data_dic.dispinfo.cutoff_iteration) && ...
                ~isempty(obj.data_dic.dispinfo.total_threads) && ...
                ~isempty(obj.data_dic.dispinfo.stepanalysis) &&  ...
                ~isempty(obj.data_dic.dispinfo.subsettrunc))
                set(obj.handles_gui.text_dicparametersloaded_s,'String','SET','ForegroundColor','g','FontWeight','bold');    
            else
                set(obj.handles_gui.text_dicparametersloaded_s,'String','NOT SET','ForegroundColor','r','FontWeight','normal');    
            end      
            
            % DIC analysis:
            if (~isempty(obj.data_dic.dispinfo) && ~isempty(obj.data_dic.displacements))
                set(obj.handles_gui.text_dicloaded_s,'String','SET','ForegroundColor','g','FontWeight','bold');    
            else
                set(obj.handles_gui.text_dicloaded_s,'String','NOT SET','ForegroundColor','r','FontWeight','normal');    
            end       
            
            % Displacements formatted:
            if (~isempty(obj.data_dic.dispinfo) && ~isempty(obj.data_dic.displacements) && ...
                ~isempty(obj.data_dic.dispinfo.pixtounits) && ...
                ~isempty(obj.data_dic.dispinfo.units) && ...
                ~isempty(obj.data_dic.dispinfo.cutoff_corrcoef) && ...
                ~isempty(obj.data_dic.dispinfo.lenscoef))
                set(obj.handles_gui.text_disploaded_s,'String','SET','ForegroundColor','g','FontWeight','bold');   
            else
                set(obj.handles_gui.text_disploaded_s,'String','NOT SET','ForegroundColor','r','FontWeight','normal'); 
            end
            
            % Strains:
            if (~isempty(obj.data_dic.straininfo) && ~isempty(obj.data_dic.strains))
                set(obj.handles_gui.text_strainloaded_s,'String','SET','ForegroundColor','g','FontWeight','bold'); 
            else
                set(obj.handles_gui.text_strainloaded_s,'String','NOT SET','ForegroundColor','r','FontWeight','normal');    
            end    
        end
        
        function update_axes(obj,action)    
            % Get data
            num_cur = getappdata(obj.handles_gui.figure,'num_cur'); 
        
            % Updates axes, axes buttons, and axes texts below:  
            if (strcmp(action,'set'))
                % Reference axes
                if (isempty(obj.reference))
                    % Clear reference image
                    cla(obj.handles_gui.axes_ref);
                    set(obj.handles_gui.text_refname,'String','Name: ');
                    set(obj.handles_gui.text_refresolution,'String','Resolution: ');                               
                else
                    % Set reference image and format axes
                    imshow(obj.reference.imginfo.get_img(),[obj.reference.imginfo.min_gs obj.reference.imginfo.max_gs],'Parent',obj.handles_gui.axes_ref);
                    set(obj.handles_gui.axes_ref ,'Visible','off');
                    
                    % Show image text
                    set(obj.handles_gui.text_refname,'String',['Name: ' obj.reference.imginfo.name]);
                    set(obj.handles_gui.text_refresolution,'String',['Resolution: ' num2str(obj.reference.imginfo.width) ' x ' num2str(obj.reference.imginfo.height)]);
                end
                
                % Current axes
                if (isempty(obj.current))
                    % Disable buttons
                    set(obj.handles_gui.button_right,'Enable','off');
                    set(obj.handles_gui.button_left,'Enable','off');
                    set(obj.handles_gui.edit_imgnum,'Enable','off');
                    set(obj.handles_gui.edit_imgnum,'String',num2str(num_cur));
                    
                    % Reset current image:
                    cla(obj.handles_gui.axes_cur);
                    
                    % Update current image name and resolution:
                    set(obj.handles_gui.text_curname,'String','Name: ');
                    set(obj.handles_gui.text_curresolution,'String','Resolution: ');
                else                                
                    % Set current image and format axes
                    imshow(obj.current(num_cur+1).imginfo.get_img(),[obj.current(num_cur+1).imginfo.min_gs obj.current(num_cur+1).imginfo.max_gs],'Parent',obj.handles_gui.axes_cur);
                    set(obj.handles_gui.axes_cur ,'Visible','off');
                    
                    % Format/show text
                    set(obj.handles_gui.text_curname,'String',['Name: ' obj.current(num_cur+1).imginfo.name]);
                    set(obj.handles_gui.text_curresolution,'String',['Resolution: ' num2str(obj.current(num_cur+1).imginfo.width) ' x ' num2str(obj.current(num_cur+1).imginfo.height)]);
                      
                    % Set left/right buttons and editbox
                    set(obj.handles_gui.edit_imgnum,'String',num2str(num_cur+1));
                    if (length(obj.current) == 1)
                        set(obj.handles_gui.button_right,'Enable','off');
                        set(obj.handles_gui.button_left,'Enable','off');
                        set(obj.handles_gui.edit_imgnum,'Enable','off');
                    elseif (num_cur == 0)
                        set(obj.handles_gui.button_right,'Enable','on');
                        set(obj.handles_gui.button_left,'Enable','off');
                        set(obj.handles_gui.edit_imgnum,'Enable','on');
                    elseif (num_cur == length(obj.current)-1)
                        set(obj.handles_gui.button_right,'Enable','off');
                        set(obj.handles_gui.button_left,'Enable','on');
                        set(obj.handles_gui.edit_imgnum,'Enable','on');
                    else
                        set(obj.handles_gui.button_right,'Enable','on');
                        set(obj.handles_gui.button_left,'Enable','on');  
                        set(obj.handles_gui.edit_imgnum,'Enable','on');                                                                      
                    end
                end  

                % ROI axes
                if (~isempty(obj.data_dic.dispinfo) && strcmp(obj.data_dic.dispinfo.type,'regular'))
                    % Set reference ROI
                    imshow(obj.reference.roi.mask,'Parent',obj.handles_gui.axes_roi);
                    set(obj.handles_gui.axes_roi,'Visible','off');
                elseif (~isempty(obj.data_dic.dispinfo) && strcmp(obj.data_dic.dispinfo.type,'backward'))
                    % Set last current ROI
                    imshow(obj.current(end).roi.mask,'Parent',obj.handles_gui.axes_roi);
                    set(obj.handles_gui.axes_roi,'Visible','off');
                else
                    % Clear ROI
                    cla(obj.handles_gui.axes_roi);  
                end 
            end
        end
                        
        function outstate = check_matlabcompat(obj) %#ok<MANU>
        % This function checks the matlab compatibility. 
            % Initialize output
            outstate = out.failed;
        
            vm = datevec(version('-date'));
            if (vm(1) < 2009)
                % Pop error message but allow user to continue
                h_error = errordlg('Developed using MATLAB 2009; may not work for versions before this date. Ncorr will continue.','Error','modal');
                uiwait(h_error);
            end
            
            % Check to make sure correct toolboxes are installed ---------%
            % Requires the image processing toolbox and the statistics
            % toolbox            
            % Note that Matlab combined their statistics and machine
            % learning toolboxes so just search for "Statistics" to work
            % for both cases.
            v = ver;
            if (isempty(strfind([v.Name],'Image Processing Toolbox')) || isempty(strfind([v.Name],'Statistics')))  %#ok<STREMP>
                h_error = errordlg('Requires the image processing and statistics toolbox.','Error','modal');
                uiwait(h_error);                
                return;
            end
                        
            % Set output
            outstate = out.success;
        end
 
        function outstate = check_mexinstall(obj)
        % Check to see if mex files exist and loads ncorr_installinfo ----%
            % Initialize output
            outstate = out.cancelled;
            
            % Files to compile:
            % These are stand alone .cpp files compiled with '-c' flag:
            lib_ncorr_cpp = {'standard_datatypes','ncorr_datatypes','ncorr_lib'};
            % These are .cpp mex functions compiled with libraries
            func_ncorr_cpp = {'ncorr_alg_formmask', ...
                              'ncorr_alg_formregions', ...
                              'ncorr_alg_formboundary', ...
                              'ncorr_alg_formthreaddiagram', ...                              
                              'ncorr_alg_formunion', ...
                              'ncorr_alg_extrapdata', ...
                              'ncorr_alg_adddisp', ...
                              'ncorr_alg_convert', ...
                              'ncorr_alg_dispgrad'};
            % These are .cpp mex functions compiled with libraries and possibly openmp:
            func_ncorr_openmp_cpp = {'ncorr_alg_testopenmp','ncorr_alg_calcseeds','ncorr_alg_rgdic'};

            % Get object extension:                
            if (ispc) % pc
                objext = 'obj';
            elseif (isunix) % unix
                objext = 'o';
            end

            % See if these files have already been compiled:
            compile = false;
            
            % Check libraries first
            for i = 0:length(lib_ncorr_cpp)-1
                if (~exist([lib_ncorr_cpp{i+1} '.' objext],'file'))
                    compile = true;   
                    break;
                end
            end
            % Check ncorr functions next
            for i = 0:length(func_ncorr_cpp)-1
                if (~exist([func_ncorr_cpp{i+1} '.' mexext],'file'))
                    compile = true;    
                    break;                     
                end
            end                     
            % Check ncorr functions with openmp
            for i = 0:length(func_ncorr_openmp_cpp)-1
                if (~exist([func_ncorr_openmp_cpp{i+1} '.' mexext],'file'))
                    compile = true;  
                    break;                       
                end
            end
            
            if (compile)                
                % One of the compiled files is missing - tell user that ncorr
                % will recompile files.
                e = msgbox('One of the compiled files is missing. If using automatic installation, then mex files will compile; if using manual installation, please close Ncorr and make sure all the files compiled properly.','WindowStyle','modal'); 
                uiwait(e);
            end

            if (~compile)
                % Functions have already been compiled - check
                % "ncorr_installinfo.txt" file to verify openmp support. 
                % If "ncorr_installinfo.txt" has incorrect formatting, then 
                % files will be recompiled. If "ncorr_installinfo.txt" is 
                % empty, ncorr will take this as the standard way to 
                % reinstall mex files.
                try
                    % ncorr_installinfo should have the format of: 
                    % [support_openmp_prelim total_cores]
                    ncorr_installinfo = dlmread('ncorr_installinfo.txt');
                    if (size(ncorr_installinfo,1) == 1 && size(ncorr_installinfo,2) == 2)
                        % ncorr_installinfo.txt has the correct number of
                        % elements
                        if ((ncorr_installinfo(1) == 0 || ncorr_installinfo(1) == 1) && ...
                            ncorr_installinfo(2) >= 1)
                            % Assuming contents have not been modified
                            % by user, set them:
                            obj.support_openmp = logical(ncorr_installinfo(1));
                            obj.total_cores = ncorr_installinfo(2);
                        else
                            % ncorr_installinfo.txt has incorrect values
                            h_error = errordlg('For some reason "ncorr_installinfo.txt" has incorrect values. Program will recompile if using automatic installing. If manually installing Ncorr, please close it and make sure the "ncorr_installinfo.txt" was created correctly.','Error','modal');
                            uiwait(h_error);
                            
                            % recompile files
                            compile = true;
                        end
                    elseif (isempty(ncorr_installinfo))
                        % ncorr_installinfo is empty, this is the stanard
                        % way to reinstall Ncorr.
                        h_error = msgbox('"ncorr_installinfo.txt" file is empty. Program will recompile if using automatic installing. If manually installing Ncorr, please close it and make sure the "ncorr_installinfo.txt" was created correctly.','WindowStyle','modal'); 
                        uiwait(h_error);
                        
                        % recompile files
                        compile = true;
                    else
                        % ncorr_installinfo.txt has extra values which
                        % aren't supposed to be there
                        h_error = errordlg('For some reason the "ncorr_installinfo.txt" file has extra or missing values. Program will recompile if using automatic installing. If manually installing Ncorr, please close it and make sure the "ncorr_installinfo.txt" was created correctly.','Error','modal');
                        uiwait(h_error);
                        
                        % recompile files
                        compile = true;
                    end
                catch %#ok<CTCH>
                    % Returns error if file does not exist or format is
                    % wrong
                    h_error = msgbox('For some reason the "ncorr_installinfo.txt" file is missing or has improper format. Program will recompile if using automatic installing. If manually installing Ncorr, please close it and make sure the "ncorr_installinfo.txt" was created correctly.','WindowStyle','modal'); 
                    uiwait(h_error);
                    
                    % recompile files
                    compile = true;
                end    
            end
            
            if (~compile)
                % Check if previous installation was corrupted by OpenMP
                % support which does not actually exist
                if (obj.support_openmp)
                    if (~ncorr_alg_testopenmp())
                        % recompile files
                        compile = true;
                    end
                end
            end
            
            %-------------------------------------------------------------%
            % BEGIN COMPILE SECTION                                       %
            %-------------------------------------------------------------%
            % COMMENT THIS SECTION OUT AND MANUALLY COMPILE IF AUTOMATIC  % 
            % COMPILATION FAILS!!!                                        %
            % ------------------------------------------------------------%
            
            try       
                if (compile)
                    % Check if ncorr.m is in current directory  
                    listing = dir;            
                    if (any(strcmp('ncorr.m',{listing.name})))   
                        % Ask user about openmp support
                        [support_openmp_prelim,total_cores_prelim,outstate_install] = gui_install(get(obj.handles_gui.figure,'OuterPosition'));
                        
                        % See if openmp GUI was cancelled                        
                        if (outstate_install ~= out.success) 
                            return;
                        end
                        
                        % Get compiler flags for openmp support. Depends
                        % on both OS and compiler.
                        flags_f = {};  
                        if (support_openmp_prelim)
                            % Get compiler information
                            info_comp = mex.getCompilerConfigurations('cpp');
                            
                            if ((isfield(struct(info_comp),'Details') && isfield(struct(info_comp.Details),'CompilerExecutable')))
                                % Get compiler
                                compiler = info_comp.Details.CompilerExecutable;
                                
                                if (~isempty(strfind(compiler,'cl'))) %#ok<STREMP>
                                    % This is Microsoft visual studio
                                    % Openmp flag is /openmp
                                    flags_f = horzcat(flags_f{:},{'COMPFLAGS="$COMPFLAGS'},{'/openmp'},{'/DNCORR_OPENMP"'});
                                elseif (~isempty(strfind(compiler,'g++'))) %#ok<STREMP>
                                    % This is the GNU compiler; this is 
                                    % assumed to be Linux.
                                    % NOTE: GCC must manually link 
                                    % against -lgomp. -lgomp also has 
                                    % to be placed as the very last 
                                    % library or else it will not compile
                                    % correctly. 
                                    % Openmp flag is -fopenmp
                                    flags_f = horzcat(flags_f{:},{'CXXFLAGS="$CXXFLAGS'},{'-fopenmp'},{'-DNCORR_OPENMP"'},{'CXXLIBS="$CXXLIBS'},{'-lgomp"'});
                                else
                                    % C++ compiler is set up, but was not
                                    % recognized as either gcc or cl
                                    h_error = errordlg('C++ compiler is set up in mex, but was not recognized as G++ or Visual Studio. Either comment out the compile section in ncorr and manually compile the mex functions, or restart ncorr and disable openmp support.','Error','modal');
                                    uiwait(h_error);
                                    
                                    outstate = out.failed;
                                    return;
                                end                                
                            else
                                % C++ compiler not installed in mex. 
                                h_error = errordlg('C++ compiler was not found. If a C++ compiler was manually installed that supports openmp, then either comment out the compile section in ncorr and manually compile the functions or restart ncorr and disable openmp.','Error','modal');
                                uiwait(h_error);
                                
                                outstate = out.failed;
                                return;
                            end
                        end

                        % Compile libraries:
                        compile_lib_cpp_mex(lib_ncorr_cpp);   
                        % Compile functions:
                        compile_func_cpp_mex(func_ncorr_cpp,lib_ncorr_cpp,{});   
                        % Compile function with openmp:
                        compile_func_cpp_mex(func_ncorr_openmp_cpp,lib_ncorr_cpp,flags_f);  

                        % Write to ncorr_installinfo.txt file, this file will be
                        % used next time ncorr is started up to avoid
                        % having to reinput this information.
                        filename = 'ncorr_installinfo.txt';    
                        fid = fopen(filename,'w');            
                        if (fid ~= -1)    
                            fclose(fid);    

                            % Write info to file:
                            dlmwrite(filename,[support_openmp_prelim total_cores_prelim]);

                            % Set support variables:
                            obj.support_openmp = support_openmp_prelim;
                            obj.total_cores = total_cores_prelim;
                        else
                            % Error
                            h_error = errordlg('For some reason ncorr was not able to create "ncorr_installinfo.txt" file.','Error','modal');
                            uiwait(h_error);
                            
                            outstate = out.failed;
                            return;
                        end
                    else
                        % Error
                        h_error = errordlg('Please navigate to folder containing ncorr.m first before reinstalling.','Error','modal');                    
                        uiwait(h_error);                        
                        
                        outstate = out.failed;
                        return;
                    end    
                end    
            catch %#ok<CTCH>
                h_error = errordlg('Files did not compile properly. Make sure mex has a c++ compiler set and that the file names have not been altered or moved. If this is the case, then try reopening and installing ncorr without openmp. If problems persist, then try commenting out the compile section in ncorr and manually compiling the functions. Instructions for compilation are available in the manual at ncorr.com','Error','modal');
                uiwait(h_error);
                
                outstate = out.failed;
                return;
            end  
            
            %-------------------------------------------------------------%
            % END OF COMPILE SECTION                                      %
            %-------------------------------------------------------------%
            
            % Check if openmp is actually enabled by calling
            % ncorr_alg_testopenmp
            if (obj.support_openmp)
                if (~ncorr_alg_testopenmp())
                    % openmp was enabled and files were compiled, but openmp
                    % isn't actually supported.
                    h_error = errordlg('Files compiled, but it was determined that OpenMP is not actually supported. Please reinstall Ncorr in single threaded mode or use a compiler which supports OpenMP.','Error','modal');
                    uiwait(h_error);

                    outstate = out.failed;
                    return;
                end
            end
            
            % Set output
            outstate = out.success;
        end
        
        %-----------------------------------------------------------------%
        % Create UI Controls for GUI and Assign Callbacks ----------------%
        %-----------------------------------------------------------------% 
        
        function handles_gui = init_gui(obj)
            % Initialize GUI Figure --------------------------------------%
            % Start Ncorr in the middle of the screen
            set(0,'units','characters'); % Set this every time incase user changes units of root
            pos_screen = get(0,'screensize');
            
            height_ncorr = 25;
            width_ncorr= 167;
            left_ncorr = (pos_screen(3)-width_ncorr)/2;
            bottom_ncorr = (pos_screen(4)-height_ncorr)/2; 
            
            handles_gui.figure = figure( ...
                'Tag', 'figure', ...
                'Units', 'characters', ...
                'Position', [left_ncorr bottom_ncorr width_ncorr height_ncorr], ...
                'Name', 'Ncorr', ...
                'MenuBar', 'none', ...
                'NumberTitle', 'off', ...
                'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
                'Resize','off', ...
                'handlevisibility','off', ...
                'DockControls','off', ...    
                'Visible','off', ...
                'IntegerHandle','off', ...
                'Interruptible','off');      
            % Set closerequestfcn
            set(handles_gui.figure,'CloseRequestFcn',obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_close_function(obj,hObject,eventdata,false),handles_gui.figure)));
            
            % Menu Items -------------------------------------------------%
            % Under File
            handles_gui.topmenu_file = uimenu( ...
                'Parent', handles_gui.figure, ...
                'Tag', 'topmenu_file', ...
                'Label', 'File', ...
                'Checked', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_loadref = uimenu( ...
                'Parent', handles_gui.topmenu_file, ...
                'Tag', 'topmenu_loadref', ...
                'Label', 'Load Reference Image', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_loadref(obj,hObject,eventdata),handles_gui.figure)), ...
                'Interruptible','off');
            
            handles_gui.topmenu_loadcur = uimenu( ...
                'Parent', handles_gui.topmenu_file, ...
                'Tag', 'topmenu_loadcur', ...
                'Label', 'Load Current Image(s)', ...
                'Checked', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_loadcur_all = uimenu( ...
                'Parent', handles_gui.topmenu_loadcur, ...
                'Tag', 'topmenu_loadcur_all', ...
                'Label', 'Load All (memory heavy)', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_loadcur(obj,hObject,eventdata,false),handles_gui.figure)), ...
                'Interruptible','off');
            
            handles_gui.topmenu_loadcur_lazy = uimenu( ...
                'Parent', handles_gui.topmenu_loadcur, ...
                'Tag', 'topmenu_loadcur_lazy', ...
                'Label', 'Load Lazy (slower but less memory)', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_loadcur(obj,hObject,eventdata,true),handles_gui.figure)), ...
                'Interruptible','off');
            
            handles_gui.topmenu_loaddata = uimenu( ...
                'Parent', handles_gui.topmenu_file, ...
                'Tag', 'topmenu_loaddata', ...
                'Label', 'Load Data', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_loaddata(obj,hObject,eventdata),handles_gui.figure)), ...
                'Interruptible','off');
            
            handles_gui.topmenu_savedata = uimenu( ...
                'Parent', handles_gui.topmenu_file, ...
                'Tag', 'topmenu_savedata', ...
                'Label', 'Save Data', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_savedata(obj,hObject,eventdata),handles_gui.figure)), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_clear = uimenu( ...
                'Parent', handles_gui.topmenu_file, ...
                'Tag', 'topmenu_clear', ...
                'Label', 'Clear Data', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_clear(obj,hObject,eventdata),handles_gui.figure)), ...
                'Interruptible','off');
            
            handles_gui.topmenu_sethandle = uimenu( ...
                'Parent', handles_gui.topmenu_file, ...
                'Tag', 'topmenu_sethandle', ...
                'Label', 'Set Handle', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_sethandle(obj,hObject,eventdata),handles_gui.figure)), ...
                'Interruptible','off');
            
            handles_gui.topmenu_reinstall = uimenu( ...
                'Parent', handles_gui.topmenu_file, ...
                'Tag', 'topmenu_reinstall', ...
                'Label', 'Reinstall', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_reinstall(obj,hObject,eventdata),handles_gui.figure)), ...
                'Interruptible','off');
            
            handles_gui.topmenu_exit = uimenu( ...
                'Parent', handles_gui.topmenu_file, ...
                'Tag', 'topmenu_exit', ...
                'Label', 'Exit', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_exit_callback(obj,hObject,eventdata),handles_gui.figure)), ...
                'Interruptible','off');

            % Under region of interest        
            handles_gui.topmenu_regionofinterest = uimenu( ...
                'Parent', handles_gui.figure, ...
                'Tag', 'topmenu_regionofinterest', ...
                'Label', 'Region of Interest', ...
                'Checked', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_setroi_ref = uimenu( ...
                'Parent', handles_gui.topmenu_regionofinterest, ...
                'Tag', 'topmenu_setroi_ref', ...
                'Label', 'Set Reference ROI', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_set_roi_ref(obj,hObject,eventdata),handles_gui.figure)), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_setroi_cur = uimenu( ...
                'Parent', handles_gui.topmenu_regionofinterest, ...
                'Tag', 'topmenu_setroi_cur', ...
                'Label', 'Set Current ROI (For "Backward" Analysis)', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_set_roi_cur(obj,hObject,eventdata),handles_gui.figure)), ...
                'Enable', 'off', ...
                'Interruptible','off');

            % Under Analysis
            handles_gui.topmenu_analysis = uimenu( ...
                'Parent', handles_gui.figure, ...
                'Tag', 'topmenu_analysis', ...
                'Label', 'Analysis', ...
                'Checked', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_setdicparameters = uimenu( ...
                'Parent', handles_gui.topmenu_analysis, ...
                'Tag', 'topmenu_setdicparameters', ...
                'Label', 'Set DIC Parameters', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_setdicparameters(obj,hObject,eventdata),handles_gui.figure)), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_perfdic = uimenu( ...
                'Parent', handles_gui.topmenu_analysis, ...
                'Tag', 'topmenu_perfdic', ...
                'Label', 'Perform DIC Analysis', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_dic(obj,hObject,eventdata),handles_gui.figure)), ...
                'Enable', 'off', ...
                'Interruptible','off');
                       
            handles_gui.topmenu_formatdisp = uimenu( ...
                'Parent', handles_gui.topmenu_analysis, ...
                'Tag', 'topmenu_formatdisp', ...
                'Label', 'Format Displacements', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_formatdisp(obj,hObject,eventdata),handles_gui.figure)), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_calcstrain = uimenu( ...
                'Parent', handles_gui.topmenu_analysis, ...
                'Tag', 'topmenu_calcstrain', ...
                'Label', 'Calculate Strains', ...
                'Checked', 'off', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_calcstrain(obj,hObject,eventdata),handles_gui.figure)), ...
                'Enable', 'off', ...
                'Interruptible','off');

            % Under plot - Do not wrap with wrapcallbackall
            handles_gui.topmenu_plot = uimenu( ...
                'Parent', handles_gui.figure, ...
                'Tag', 'topmenu_plot', ...
                'Label', 'Plot', ...
                'Checked', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_viewdispplot = uimenu( ...
                'Parent', handles_gui.topmenu_plot, ...
                'Tag', 'topmenu_viewdispplot', ...
                'Label', 'View Displacement Plots', ...
                'Checked', 'off', ...
                'Enable', 'on', ...
                'Interruptible','off');
            
            handles_gui.topmenu_viewdispplot_u = uimenu( ...
                'Parent', handles_gui.topmenu_viewdispplot, ...
                'Tag', 'topmenu_viewdispplot_u', ...
                'Label', 'U', ...
                'Checked', 'off', ...
                'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_viewplot(obj,hObject,eventdata,{'u'}),handles_gui.figure), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_viewdispplot_v = uimenu( ...
                'Parent', handles_gui.topmenu_viewdispplot, ...
                'Tag', 'topmenu_viewdispplot_v', ...
                'Label', 'V', ...
                'Checked', 'off', ...
                'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_viewplot(obj,hObject,eventdata,{'v'}),handles_gui.figure), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_viewdispplot_all = uimenu( ...
                'Parent', handles_gui.topmenu_viewdispplot, ...
                'Tag', 'topmenu_viewdispplot_all', ...
                'Label', 'All', ...
                'Checked', 'off', ...
                'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_viewplot(obj,hObject,eventdata,{'u','v'}),handles_gui.figure), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_viewstrainplot = uimenu( ...
                'Parent', handles_gui.topmenu_plot, ...
                'Tag', 'topmenu_viewstrainpplot', ...
                'Label', 'View Strains Plots', ...
                'Checked', 'off', ...
                'Enable', 'on', ...
                'Interruptible','off');

            handles_gui.topmenu_viewstrainplot_exx = uimenu( ...
                'Parent', handles_gui.topmenu_viewstrainplot, ...
                'Tag', 'topmenu_viewstrainpplot_exx', ...
                'Label', 'Exx', ...
                'Checked', 'off', ...
                'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_viewplot(obj,hObject,eventdata,{'exx'}),handles_gui.figure), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_viewstrainplot_exy = uimenu( ...
                'Parent', handles_gui.topmenu_viewstrainplot, ...
                'Tag', 'topmenu_viewstrainpplot_exy', ...
                'Label', 'Exy', ...
                'Checked', 'off', ...
                'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_viewplot(obj,hObject,eventdata,{'exy'}),handles_gui.figure), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_viewstrainplot_eyy = uimenu( ...
                'Parent', handles_gui.topmenu_viewstrainplot, ...
                'Tag', 'topmenu_viewstrainpplot_eyy', ...
                'Label', 'Eyy', ...
                'Checked', 'off', ...
                'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_viewplot(obj,hObject,eventdata,{'eyy'}),handles_gui.figure), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_viewstrainplot_all = uimenu( ...
                'Parent', handles_gui.topmenu_viewstrainplot, ...
                'Tag', 'topmenu_viewstrainpplot_all', ...
                'Label', 'All', ...
                'Checked', 'off', ...
                'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_viewplot(obj,hObject,eventdata,{'exx','exy','eyy'}),handles_gui.figure), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            handles_gui.topmenu_closeplots = uimenu( ...
                'Parent', handles_gui.topmenu_plot, ...
                'Tag', 'topmenu_viewstrainpplot_all', ...
                'Label', 'Close All Plots', ...
                'Checked', 'off', ...
                'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_topmenu_closeplots(obj,hObject,eventdata),handles_gui.figure), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            % Panels -----------------------------------------------------%
            handles_gui.group_state = uibuttongroup( ...
                'Parent', handles_gui.figure, ...
                'Tag', 'group_state', ...
                'Units', 'characters', ...
                'Position', [2 13.8 35 10.7], ...
                'Title', 'Program State', ...
                'Interruptible','off');
            
            handles_gui.group_ref = uibuttongroup( ...
                'Parent', handles_gui.figure, ...
                'Tag', 'group_ref', ...
                'Units', 'characters', ...
                'Position', [39.0 0.75 62 23.75], ...
                'Title', 'Reference Image', ...
                'Interruptible','off');
            
            handles_gui.group_cur = uibuttongroup( ...
                'Parent', handles_gui.figure, ...
                'Tag', 'group_cur', ...
                'Units', 'characters', ...
                'Position', [103 0.75 62 23.75], ...
                'Title', 'Current Image(s)', ...
                'Interruptible','off');
            
            handles_gui.group_roi = uibuttongroup( ...
                'Parent', handles_gui.figure, ...
                'Tag', 'group_roi', ...
                'Units', 'characters', ...
                'Position', [2 0.75 35 12.3], ...
                'Title', 'Region of Interest', ...
                'Interruptible','off');

            % Axes -------------------------------------------------------%
            handles_gui.axes_ref = axes( ...
                'Parent', handles_gui.group_ref, ...
                'Tag', 'axes_ref', ...
                'Units', 'characters', ...
                'Position', [2 4 57 18], ...
                'Visible', 'off', ...
                'handlevisibility','off', ...
                'Interruptible','off');
            
            handles_gui.axes_cur = axes( ...
                'Parent', handles_gui.group_cur, ...
                'Tag', 'axes_cur', ...
                'Units', 'characters', ...
                'Position', [2 4 57 18], ...
                'Visible', 'off', ...
                'handlevisibility','off', ...
                'Interruptible','off');
            
            handles_gui.axes_roi = axes( ...
                'Parent', handles_gui.group_roi, ...
                'Tag', 'axes_roi', ...
                'Units', 'characters', ...
                'Position', [2 1 29.5 9.6],...
                'Visible', 'off', ...
                'handlevisibility','off', ...
                'Interruptible','off');

            % Image texts ------------------------------------------------%
            handles_gui.text_refname = uicontrol( ...
                'Parent', handles_gui.group_ref, ...
                'Tag', 'text_refname', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [2.1 2 56.9 1.3], ...
                'String', 'Name: ', ...
                'HorizontalAlignment', 'left', ...
                'Interruptible','off');
            
            handles_gui.text_refresolution = uicontrol( ...
                'Parent', handles_gui.group_ref, ...
                'Tag', 'text_refresolution', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [2.1 .5 56.9 1.3], ...
                'String', 'Resolution: ', ...
                'HorizontalAlignment', 'left', ...
                'Interruptible','off');
            
            handles_gui.text_curname = uicontrol( ...
                'Parent', handles_gui.group_cur, ...
                'Tag', 'text_curname', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [2.1 2 56.9 1.3], ...
                'String', 'Name: ', ...
                'HorizontalAlignment', 'left', ...
                'Interruptible','off');
            
            handles_gui.text_curresolution = uicontrol( ...
                'Parent', handles_gui.group_cur, ...
                'Tag', 'text_curresolution', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [2.1 .5 56.9 1.3], ...
                'String', 'Resolution: ', ...
                'HorizontalAlignment', 'left', ...
                'Interruptible','off');

            % State texts ------------------------------------------------%
            handles_gui.text_refloaded = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_refloaded', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [2.1 7.7 21 1.3], ...
                'String', 'Reference Image', ...
                'HorizontalAlignment', 'left', ...
                'Interruptible','off');
            
            handles_gui.text_curloaded = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_curloaded', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [2.1 6.5 21 1.3], ...
                'String', 'Current Image(s)', ...
                'HorizontalAlignment', 'left', ...
                'Interruptible','off');
            
            handles_gui.text_roiloaded = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_roiloaded', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [2.1 5.3 21 1.3], ...
                'String', 'Region of Interest', ...
                'HorizontalAlignment', 'left', ...
                'Interruptible','off');
            
            handles_gui.text_dicparamsloaded = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_dicparamsloaded', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [2.1 4.1 21 1.3], ...
                'String', 'DIC Parameters', ...
                'HorizontalAlignment', 'left', ...
                'Interruptible','off');
            
            handles_gui.text_dicloaded = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_dicloaded', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [2.1 2.9 21 1.3], ...
                'String', 'DIC Analysis', ...
                'HorizontalAlignment', 'left', ...
                'Interruptible','off');
            
            handles_gui.text_disploaded = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_disploaded', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [2.1 1.7 21 1.3], ...
                'String', 'Displacements', ...
                'HorizontalAlignment', 'left', ...
                'Interruptible','off');
            
            handles_gui.text_strainloaded = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_strainloaded', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [2.1 .5 21 1.3], ...
                'String', 'Strains', ...
                'HorizontalAlignment', 'left', ...
                'Interruptible','off');

            % More state texts -------------------------------------------%
            handles_gui.text_refloaded_s = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_refloaded_s', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [23 7.7 10 1.3], ...
                'String', 'NOT SET', ...
                'HorizontalAlignment', 'left', ...
                'ForegroundColor', 'r', ...
                'Interruptible','off');
            
            handles_gui.text_curloaded_s = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_curloaded_s', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [23 6.5 10 1.3], ...
                'String', 'NOT SET', ...
                'HorizontalAlignment', 'left', ...
                'ForegroundColor', 'r', ...
                'Interruptible','off');
            
            handles_gui.text_roiloaded_s = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_roiloaded_s', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [23 5.3 10 1.3], ...
                'String', 'NOT SET', ...
                'HorizontalAlignment', 'left', ...
                'ForegroundColor', 'r', ...
                'Interruptible','off');
            
            handles_gui.text_dicparametersloaded_s = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_dicparametersloaded_s', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [23 4.1 10 1.3], ...
                'String', 'NOT SET', ...
                'HorizontalAlignment', 'left', ...
                'ForegroundColor', 'r', ...
                'Interruptible','off');
            
            handles_gui.text_dicloaded_s = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_dicloaded_s', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [23 2.9 10 1.3], ...
                'String', 'NOT SET', ...
                'HorizontalAlignment', 'left', ...
                'ForegroundColor', 'r', ...
                'Interruptible','off');
            
            handles_gui.text_disploaded_s = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_disploaded_s', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [23 1.7 10 1.3], ...
                'String', 'NOT SET', ...
                'HorizontalAlignment', 'left', ...
                'ForegroundColor', 'r', ...
                'Interruptible','off');
            
            handles_gui.text_strainloaded_s = uicontrol( ...
                'Parent', handles_gui.group_state, ...
                'Tag', 'text_strainloaded_s', ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [23 0.5 10 1.3], ...
                'String', 'NOT SET', ...
                'HorizontalAlignment', 'left', ...
                'ForegroundColor', 'r', ...
                'Interruptible','off');   
            
            % Editbox ----------------------------------------------------%
            handles_gui.edit_imgnum = uicontrol( ...
                'Parent', handles_gui.group_cur, ...
                'Tag', 'edit_imgnum', ...
                'Style', 'edit', ...
                'Units', 'characters', ...
                'Position', [44 1.3 7 1.6], ...
                'String', '', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_edit_imgnum(obj,hObject,eventdata),handles_gui.figure)), ...
                'Enable', 'off', ...
                'Interruptible','off');
            
            % Pushbuttons ------------------------------------------------%
            handles_gui.button_left = uicontrol( ...
                'Parent', handles_gui.group_cur, ...
                'Tag', 'button_left', ...
                'Style', 'pushbutton', ...
                'Units', 'characters', ...
                'Position', [37 1.2 6 1.8], ...
                'String', '<', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_button_left(obj,hObject,eventdata),handles_gui.figure)), ...
                'Enable', 'off', ...
                'Interruptible','off');

            handles_gui.button_right = uicontrol( ...
                'Parent', handles_gui.group_cur, ...
                'Tag', 'button_right', ...
                'Style', 'pushbutton', ...
                'Units', 'characters', ...
                'Position', [52 1.2 6 1.8], ...
                'String', '>', ...
                'Callback', obj.util_wrapcallbackall(ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata) callback_button_right(obj,hObject,eventdata),handles_gui.figure)), ...
                'Enable', 'off', ...
                'Interruptible','off');
        end
    end
end

%-------------------------------------------------------------------------%
% Other functions --------------------------------------------------------%
%-------------------------------------------------------------------------%

function [handle_name,outstate] = gui_sethandle(pos_parent)  
% This is a GUI for setting the name of the variable that points to ncorr.
%
% Inputs -----------------------------------------------------------------%
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%
% Outputs ----------------------------------------------------------------%
%   handle_name - string; name of new handle
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.
  
    % Data ---------------------------------------------------------------%
    % Initialize outputs
    outstate = out.cancelled;
    handle_name = '';
    % Get GUI handles - Part of output
    handles_gui_sethandle = init_gui();  
    % Run c-tor
    feval(ncorr_util_wrapcallbacktrycatch(@constructor,handles_gui_sethandle.figure));

    % Callbacks and functions --------------------------------------------%
    function constructor()
        % Set Data
        setappdata(handles_gui_sethandle.figure,'handle_name_prelim','');

        % Set Visible
        set(handles_gui_sethandle.figure,'Visible','on'); 
    end

    function callback_edit_string_handle(hObject,eventdata) %#ok<INUSD>
        % Get data
        handle_name_prelim = getappdata(handles_gui_sethandle.figure,'handle_name_prelim');
        
        % Get Name
        handle_name_buffer = get(handles_gui_sethandle.edit_string_handle,'string'); 
        if (isvarname(handle_name_buffer))   
            handle_name_prelim = handle_name_buffer;
        else
            h_error = errordlg('Not a valid variable name.','Error','modal');
            uiwait(h_error);
        end

        % Set Data
        setappdata(handles_gui_sethandle.figure,'handle_name_prelim',handle_name_prelim);
            
        update_sidemenu();
    end

    function callback_button_finish(hObject,eventdata) %#ok<INUSD>
        % Get data 
        handle_name_prelim = getappdata(handles_gui_sethandle.figure,'handle_name_prelim');

        % Set output
        handle_name = handle_name_prelim;
        outstate = out.success;

        % Exit
        close(handles_gui_sethandle.figure);
    end

    function callback_button_cancel(hObject,eventdata) %#ok<INUSD>
        % Close
        close(handles_gui_sethandle.figure);
    end

    function update_sidemenu(hObject,eventdata) %#ok<INUSD>
        % Get data
        handle_name_prelim = getappdata(handles_gui_sethandle.figure,'handle_name_prelim');

        set(handles_gui_sethandle.edit_string_handle,'String',handle_name_prelim);
    end

    function handles_gui_sethandle = init_gui()
        % Form UIs
        handles_gui_sethandle.figure = figure( ...
            'Tag', 'figure', ...
            'Units', 'characters', ...
            'Position', ncorr_util_figpos(pos_parent,[6 40.5]), ...
            'Name', 'Set Image Size', ...
            'MenuBar', 'none', ...
            'NumberTitle', 'off', ...
            'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
            'handlevisibility','off', ...
            'DockControls','off', ...
            'Resize','off', ...
            'WindowStyle','modal', ...
            'Visible','off', ...
            'IntegerHandle','off', ...
            'Interruptible','off');    

        % Text
        handles_gui_sethandle.text_string_handle = uicontrol( ...
            'Parent', handles_gui_sethandle.figure, ...
            'Tag', 'text_string_handle', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.1 3.5 15 1.3], ...
            'String', 'Handle Name:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
        
        % Edit box
        handles_gui_sethandle.edit_string_handle = uicontrol( ...
            'Parent', handles_gui_sethandle.figure, ...
            'Tag', 'edit_string_handle', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [17.2 3.6 20.6 1.3], ...
            'HorizontalAlignment', 'left', ...            
            'String', '', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_string_handle,handles_gui_sethandle.figure), ...
            'Interruptible','off');  

        % Buttons
        handles_gui_sethandle.button_finish = uicontrol( ...
            'Parent', handles_gui_sethandle.figure, ...
            'Tag', 'button_finish', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 1 12 1.7], ...
            'String', 'Finish', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_finish,handles_gui_sethandle.figure), ...
            'Interruptible','off');

        handles_gui_sethandle.button_cancel = uicontrol( ...
            'Parent', handles_gui_sethandle.figure, ...
            'Tag', 'button_cancel', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [16 1 12 1.7], ...
            'String', 'Cancel', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_cancel,handles_gui_sethandle.figure), ...
            'Interruptible','off');
    end

    % Pause until figure is closed ---------------------------------------%
    waitfor(handles_gui_sethandle.figure);        
end
            
function [support_openmp,total_cores,outstate] = gui_install(pos_parent)
% This is a GUI to prompt user if openmp support exists.
%
% Inputs -----------------------------------------------------------------%
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%
% Outputs ----------------------------------------------------------------%
%   support_openmp - logical; tells whether openmp support exists
%   total_cores - integer; tells the number of cores that exist in the case
%   that openmp is supported
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.

    % Data ---------------------------------------------------------------%
    % Initialize outputs
    outstate = out.cancelled;
    support_openmp = false;
    total_cores = 1;
    % Get GUI handles
    handles_gui = init_gui();
    % Run c-tor
    feval(ncorr_util_wrapcallbacktrycatch(@constructor,handles_gui.figure));
    
    % Callbacks and functions --------------------------------------------%
    function constructor                       
        % Set data
        setappdata(handles_gui.figure,'total_cores_prelim',1); 
        setappdata(handles_gui.figure,'total_cores_min',1); 
        setappdata(handles_gui.figure,'total_cores_max',64); 
        setappdata(handles_gui.figure,'val_checkbox_openmp',false); 
        
        % Update        
        update_sidemenu();
        
        % Format window
        set(handles_gui.figure,'Visible','on');        
    end    
    
    % Callbacks and functions --------------------------------------------%
    function callback_checkbox_openmp(hObject,eventdata) %#ok<INUSD>
        % Get data
        total_cores_prelim = getappdata(handles_gui.figure,'total_cores_prelim');
        
        val_checkbox_openmp = get(handles_gui.checkbox_openmp,'Value');
        if (~val_checkbox_openmp)
            % Reset total cores to 1
            total_cores_prelim = 1;
        end
        
        % Set data
        setappdata(handles_gui.figure,'val_checkbox_openmp',val_checkbox_openmp);
        setappdata(handles_gui.figure,'total_cores_prelim',total_cores_prelim);
        
        update_sidemenu();
    end

    function callback_edit_total_cores(hObject,eventdata) %#ok<INUSD>
        % Get Data
        total_cores_prelim = getappdata(handles_gui.figure,'total_cores_prelim');
        total_cores_min = getappdata(handles_gui.figure,'total_cores_min'); 
        total_cores_max = getappdata(handles_gui.figure,'total_cores_max'); 
        
        % Get buffer
        total_cores_buffer = str2double(get(handles_gui.edit_total_cores,'string')); 
        if (ncorr_util_isintbb(total_cores_buffer,total_cores_min,total_cores_max,'Number of cores') == out.success)
            total_cores_prelim = total_cores_buffer;       
        end
        
        % Set data
        setappdata(handles_gui.figure,'total_cores_prelim',total_cores_prelim);
        
        update_sidemenu();
    end
        
    function callback_button_finish(hObject,eventdata) %#ok<INUSD>
        % Get data
        total_cores_prelim = getappdata(handles_gui.figure,'total_cores_prelim');
        support_openmp_prelim = getappdata(handles_gui.figure,'val_checkbox_openmp');       
        
        % Store
        total_cores = total_cores_prelim;
        support_openmp = support_openmp_prelim;
        outstate = out.success;       
                
        % Return
        close(handles_gui.figure);
    end

    function callback_button_cancel(hObject,eventdata) %#ok<INUSD>
        close(handles_gui.figure);
    end

    function update_sidemenu(hObject,eventdata) %#ok<INUSD>
        % Get data
        val_checkbox_openmp = getappdata(handles_gui.figure,'val_checkbox_openmp');
        total_cores_prelim = getappdata(handles_gui.figure,'total_cores_prelim');
        
        set(handles_gui.checkbox_openmp,'Value',val_checkbox_openmp);
        
        if (val_checkbox_openmp)
            set(handles_gui.edit_total_cores,'enable','on');
        else
            set(handles_gui.edit_total_cores,'enable','off');
        end
        
        set(handles_gui.edit_total_cores,'String',num2str(total_cores_prelim));
    end

    function handles_gui = init_gui()
        handles_gui.figure = figure( ...
            'Tag', 'figure', ...
            'Units', 'characters', ...
            'Position', ncorr_util_figpos(pos_parent,[9.5 70]), ...
            'Name', 'OpenMP support', ...
            'MenuBar', 'none', ...
            'NumberTitle', 'off', ...
            'Resize', 'off', ...
            'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
            'handlevisibility','off', ...
            'DockControls','off', ...
            'Visible','off', ...
            'IntegerHandle','off', ...
            'WindowStyle','modal', ...
            'Interruptible','off');      
        
        handles_gui.text_openmp = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'text_openmp', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [1.5 5.1 66.7 2.2], ...
            'HorizontalAlignment','left', ...
            'String', 'Requires multicore CPU and compiler which supports OpenMP installed through mex.', ...
            'Value', 0, ...
            'Interruptible','off');

        handles_gui.text_cores_total = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'text_core_total', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [10 3.2 18 1.2], ...
            'HorizontalAlignment','left', ...
            'String', 'Cores:', ...
            'Value', 0, ...
            'Interruptible','off');
        
        handles_gui.checkbox_openmp = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'checkbox_openmp', ...
            'Style', 'checkbox', ...
            'Units', 'characters', ...
            'Position', [1.5 7.5 66.7 1.3], ...
            'String', 'OpenMP Multithreading', ...
            'Value', 0, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_checkbox_openmp,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.edit_total_cores = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'edit_total_cores', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [22 3.3 8 1.2], ...
            'HorizontalAlignment','left', ...
            'String', '1', ...
            'Enable', 'off', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_total_cores,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_finish = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'button_finish', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [42.3 0.7 12 1.8], ...
            'String', 'Finish', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_finish,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_cancel = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'button_cancel', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [56 0.7 12 1.8], ...
            'String', 'Cancel', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_cancel,handles_gui.figure), ...
            'Interruptible','off');
    end

    % Pause until figure is closed ---------------------------------------%
    waitfor(handles_gui.figure);
end

function compile_lib_cpp_mex(lib_cpp)
% This function compiles cpp_function_lib as object files.
%
% Inputs -----------------------------------------------------------------%
%   lib_cpp - cell of strings; names of libraries to be compiled as object
%   files
%
% Returns error if files are not found.
    
    % Cycle over libraries and compile
    for i = 0:length(lib_cpp)-1           
        % First check if the cpp and header files exist
        if (exist([lib_cpp{i+1} '.cpp'],'file') && exist([lib_cpp{i+1} '.h'],'file'))
            disp(['Installing ' lib_cpp{i+1} '... Please wait']);

            % Generate compiler string
            string_compile = horzcat({'-c'},{[lib_cpp{i+1} '.cpp']});

            % Compile file
            mex(string_compile{:});
        else
            h_error = errordlg(['Files ' lib_cpp{i+1} '.cpp and ' lib_cpp{i+1} '.h were not found. Please find them and place them in the current directory before proceeding.'],'Error','modal');
            uiwait(h_error);
            
            % Spit error to invoke exception
            error('Compilation failed because file was not found.');
        end  
    end
end

function compile_func_cpp_mex(func_cpp,lib_cpp,flags_f)
% This function compiles func_cpp and linkes them with the lib_cpp object 
% files. If flags are provided then they will be appended when compiling.
%
% Inputs -----------------------------------------------------------------%
%   func_cpp - cell of strings; names of .cpp files to be
%   compiled
%   lib_cpp - cell of strings; names of object files to be linked
%   flags_f - cell of strings; formatted compiler flags.
%
% Returns error if library object files or source code files are not found.

    % Get the OS to find object extension
    if (ispc) % pc
        objext = 'obj';
    elseif (isunix) % unix
        objext = 'o';
    end

    % Check if lib files have been compiled first
    for i = 0:length(lib_cpp)-1           
        if (~exist([lib_cpp{i+1} '.' objext],'file'))
            % Object file doesnt exist, return error
            h_error = errordlg(['Library ' lib_cpp{i+1} ' has not yet been compiled. Please compile first and make sure the object file is located in the current directory before proceeding.'],'Error','modal');
            uiwait(h_error);
            
            % Spit error to invoke exception
            error('Library not compiled yet');
        end
    end

    % Cycle over mex files and compile
    for i = 0:length(func_cpp)-1
        % First check if the cpp file exists
        if (exist([func_cpp{i+1} '.cpp'],'file'))
            disp(['Installing ' func_cpp{i+1} '... Please wait']);

            % Generate compiler string
            string_compile = {[func_cpp{i+1} '.cpp']};   

            % Append libraries
            for j = 0:length(lib_cpp)-1         
                string_compile = horzcat(string_compile,{[lib_cpp{j+1} '.' objext]}); %#ok<AGROW>
            end               

            % Append flags                      
            string_compile = horzcat(string_compile,flags_f{:}); %#ok<AGROW>

            % Compile
            mex(string_compile{:});
        else
            h_error = errordlg(['File ' filename '.cpp was not found. Please find it and place it in the current directory before proceeding'],'Error','modal');
            uiwait(h_error);
            
            % Spit error to invoke exception
            error('Compilation failed because file was not found.');
        end
    end      
end

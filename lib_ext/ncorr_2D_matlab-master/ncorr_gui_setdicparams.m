function [radius,spacing,cutoff_diffnorm,cutoff_iteration,total_threads,stepanalysis,subsettrunc,outstate] = ncorr_gui_setdicparams(imgs,rois,support_openmp,total_cores,pos_parent,params_init)
% This is a GUI for setting the DIC parameters.
%
% Inputs -----------------------------------------------------------------%
%   imgs - ncorr_class_img; set of images. Generally, it's a single image;
%   usually the reference image or the last current image, but multiple
%   images are supported. Used to display the background image.
%   rois - ncorr_class_roi; ROI corresponding to each image, which is used
%   to display the region of interest over the image.
%   support_openmp - logical; indicates whether there is openmp support. If
%   not, then the number of threads is locked to one.
%   total_cores - integer; indicates the default number of cores specified
%   when ncorr was installed.
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%   params_init - cell; contains parameters {radius spacing cutoff_diffnorm 
%   cutoff_iteration total_threads stepanalysis subsettrunc} to intialize DIC
%   parameters if analysis has been done before.
%
% Outputs ----------------------------------------------------------------%
%   radius - integer; subset radius parameter
%   spacing - integer; spacing parameter
%   cutoff_diffnorm - double; cutoff for the norm of the difference vector 
%   cutoff_iteration - integer; number of iterations allowed before exiting
%   total_threads - integer; number of threads used in analysis
%   stepanalysis - struct; information about step analysis containing
%   struct('enabled',{},'type',{},'auto',{},'step',{})
%   subsettrunc - logical; true if subset truncation has been enabled for
%   discontinuous DIC analysis
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.
%
% Note that each region in rois must be non-empty. If called through Ncorr,
% then the regions will be guaranteed non-empty and contiguous at this point.

    % Data ---------------------------------------------------------------%
    % Initialize Outputs
    outstate = out.cancelled;
    radius = [];
    spacing = [];
    cutoff_diffnorm = [];
    cutoff_iteration = [];
    total_threads = [];
    stepanalysis = struct('enabled',{},'type',{},'auto',{},'step',{});
    subsettrunc = [];
    % Get GUI handles
    handles_gui = init_gui();
    % Run c-tor
    feval(ncorr_util_wrapcallbacktrycatch(@constructor,handles_gui.figure));
    
    % Callbacks and functions --------------------------------------------%
    function constructor()                
        % Mins and maxes - most are just heuristics
        min_radius = 10;
        min_spacing = 0;
        min_cutoff_diffnorm = 10^-8;
        min_cutoff_iteration = 5;
        min_total_threads = 1;   
        min_step = 1;
        max_radius = 200;
        max_spacing = 50;
        max_cutoff_diffnorm = 1e-2;
        max_cutoff_iteration = 100;
        max_total_threads = 64;
        max_step = 100;

        % Radius and spacing functions - these set the maximums for the
        % sliders, but higher values can be input through the edit box.
        get_radius = @(x)round(50*x+min_radius);
        get_spacing = @(x)round(10*x+min_spacing);
        
        % Get gs buffer - start off with the first image
        gs_buffer = imgs(1).get_gs();

        % Initialize 
        stepanalysis_prelim = struct('enabled',{},'type',{},'auto',{},'step',{});
        if (isempty(params_init))      
            % Initial radius and spacing based on heuristics and size of
            % the first image
            radius_prelim = min(max(round((sum(size(gs_buffer))/2)/35),min_radius),max_radius); 
            spacing_prelim = min(max(round(((sum(size(gs_buffer))-400)/2)/400),min_spacing),max_spacing); 
            cutoff_diffnorm_prelim = 10^-6;
            cutoff_iteration_prelim = 50;                
            if (support_openmp)
                % Initialize to the number of total cores if openmp is
                % supported. 
                total_threads_prelim = total_cores;
            else
                total_threads_prelim = 1;
            end
            
            % Set step analysis
            stepanalysis_prelim(1).enabled = false;
            stepanalysis_prelim.type = 'seed'; % Initialize to seed
            stepanalysis_prelim.auto = true; 
            stepanalysis_prelim.step = 5; 
            
            subsettrunc_prelim = false;
        else
            radius_prelim = params_init{1};
            spacing_prelim = params_init{2};            
            cutoff_diffnorm_prelim = params_init{3}; 
            cutoff_iteration_prelim = params_init{4};
            total_threads_prelim = params_init{5};
            stepanalysis_prelim(1) = params_init{6};
            subsettrunc_prelim = params_init{7};
        end

        % Set Data
        % min and max
        setappdata(handles_gui.figure,'min_radius',min_radius);
        setappdata(handles_gui.figure,'min_spacing',min_spacing);
        setappdata(handles_gui.figure,'min_cutoff_diffnorm',min_cutoff_diffnorm);
        setappdata(handles_gui.figure,'min_cutoff_iteration',min_cutoff_iteration);
        setappdata(handles_gui.figure,'min_total_threads',min_total_threads);
        setappdata(handles_gui.figure,'min_step',min_step);
        setappdata(handles_gui.figure,'max_radius',max_radius);
        setappdata(handles_gui.figure,'max_spacing',max_spacing);
        setappdata(handles_gui.figure,'max_cutoff_diffnorm',max_cutoff_diffnorm);
        setappdata(handles_gui.figure,'max_cutoff_iteration',max_cutoff_iteration);
        setappdata(handles_gui.figure,'max_total_threads',max_total_threads);
        setappdata(handles_gui.figure,'max_step',max_step);
        % Radius and spacing prelims
        setappdata(handles_gui.figure,'radius_prelim',radius_prelim);
        setappdata(handles_gui.figure,'spacing_prelim',spacing_prelim);
        setappdata(handles_gui.figure,'radius_slider_prelim',min(1,max(0, fzero(@(x) get_radius(x)-radius_prelim,0))));
        setappdata(handles_gui.figure,'spacing_slider_prelim',min(1,max(0, fzero(@(x) get_spacing(x)-spacing_prelim,0))));
        setappdata(handles_gui.figure,'cutoff_diffnorm_prelim',cutoff_diffnorm_prelim);
        setappdata(handles_gui.figure,'cutoff_iteration_prelim',cutoff_iteration_prelim);
        setappdata(handles_gui.figure,'total_threads_prelim',total_threads_prelim);
        setappdata(handles_gui.figure,'stepanalysis_prelim',stepanalysis_prelim);
        setappdata(handles_gui.figure,'subsettrunc_prelim',subsettrunc_prelim);
        % Image buffers
        setappdata(handles_gui.figure,'gs_buffer',gs_buffer);
        setappdata(handles_gui.figure,'preview_subsetloc',[]);
        setappdata(handles_gui.figure,'preview_subset',[]);
        % Images handles
        setappdata(handles_gui.figure,'handle_subsetloc',[]);
        setappdata(handles_gui.figure,'handle_subset',[]);
        % Impoint handles
        setappdata(handles_gui.figure,'handle_points',impoint.empty);
        % Info
        setappdata(handles_gui.figure,'num_img',0);
        setappdata(handles_gui.figure,'num_region',0);
        setappdata(handles_gui.figure,'pos_points',[]);
        % Radius and spacing functions
        setappdata(handles_gui.figure,'get_radius',get_radius);
        setappdata(handles_gui.figure,'get_spacing',get_spacing);

        % Update  
        update_axes('set');
        update_sidemenu();   
        
        % Set Visible
        set(handles_gui.figure,'Visible','on'); 
    end
    
    function callback_popupmenu(hObject,eventdata) %#ok<INUSD>
        % Get data
        stepanalysis_prelim = getappdata(handles_gui.figure,'stepanalysis_prelim');
        
        % Get popup menu value
        if (get(handles_gui.popupmenu,'Value') == 1)
            stepanalysis_prelim.type = 'seed';
        else
            stepanalysis_prelim.type = 'leapfrog';
        end
        
        % Set data
        setappdata(handles_gui.figure,'stepanalysis_prelim',stepanalysis_prelim);  
        
        % Update
        update_sidemenu();        
    end

    function callback_slider_radius(hObject,eventdata) %#ok<INUSD>
        % Get data
        get_radius = getappdata(handles_gui.figure,'get_radius');
        
        % Get slider value
        radius_slider_prelim = get(handles_gui.slider_radius,'value'); 
        
        % Get radius_prelim
        radius_prelim = get_radius(radius_slider_prelim);   
        
        % Set data
        setappdata(handles_gui.figure,'radius_prelim',radius_prelim);
        setappdata(handles_gui.figure,'radius_slider_prelim',radius_slider_prelim);    
        
        % Update
        update_axes('set');
        update_sidemenu();     
    end

    function callback_slider_spacing(hObject,eventdata) %#ok<INUSD>
        % Get data        
        get_spacing = getappdata(handles_gui.figure,'get_spacing');
        
        % Get slider value
        spacing_slider_prelim = get(handles_gui.slider_spacing,'value'); 
        
        % Get spacing_prelim
        spacing_prelim = get_spacing(spacing_slider_prelim);
        
        % Set data
        setappdata(handles_gui.figure,'spacing_prelim',spacing_prelim);
        setappdata(handles_gui.figure,'spacing_slider_prelim',spacing_slider_prelim);
        
        % Update
        update_axes('update');
        update_sidemenu();     
    end

    function callback_checkbox_stepanalysis(hObject,eventdata) %#ok<INUSD>
        % Get data
        stepanalysis_prelim = getappdata(handles_gui.figure,'stepanalysis_prelim');
        
        % Get checkbox value
        stepanalysis_prelim.enabled = get(handles_gui.checkbox_stepanalysis,'Value');
        
        % Set data
        setappdata(handles_gui.figure,'stepanalysis_prelim',stepanalysis_prelim);  
        
        % Update
        update_sidemenu();        
    end

    function callback_checkbox_seedauto(hObject,eventdata) %#ok<INUSD>
        % Get data
        stepanalysis_prelim = getappdata(handles_gui.figure,'stepanalysis_prelim');
        
        % Get checkbox value
        stepanalysis_prelim.auto = get(handles_gui.checkbox_seedauto,'Value');
        
        % Set data
        setappdata(handles_gui.figure,'stepanalysis_prelim',stepanalysis_prelim);  
        
        % Update
        update_sidemenu();        
    end

    function callback_checkbox_subsettrunc(hObject,eventdata) %#ok<INUSD>        
        % Get checkbox value
        subsettrunc_prelim = get(handles_gui.checkbox_subsettrunc,'Value');
        
        % Set data
        setappdata(handles_gui.figure,'subsettrunc_prelim',subsettrunc_prelim);  
        
        % Update
        update_axes('update');
        update_sidemenu();        
    end

    function callback_edit_radius(hObject,eventdata) %#ok<INUSD>
        % Get Data
        min_radius = getappdata(handles_gui.figure,'min_radius');
        max_radius = getappdata(handles_gui.figure,'max_radius');
        radius_prelim = getappdata(handles_gui.figure,'radius_prelim'); 
        radius_slider_prelim = getappdata(handles_gui.figure,'radius_slider_prelim');
        get_radius = getappdata(handles_gui.figure,'get_radius');
        
        % Get radius num
        radius_buffer = str2double(get(handles_gui.edit_radius,'string')); 
        if (ncorr_util_isintbb(radius_buffer,min_radius,max_radius,'Radius') == out.success)
            radius_prelim = radius_buffer;   
            radius_slider_prelim = min(1,max(0,fzero(@(x)get_radius(x)-radius_prelim,0)));
        end
            
        % Set data
        setappdata(handles_gui.figure,'radius_prelim',radius_prelim); 
        setappdata(handles_gui.figure,'radius_slider_prelim',radius_slider_prelim); 
                
        % Update
        update_axes('set');
        update_sidemenu();     
    end

    function callback_edit_spacing(hObject,eventdata) %#ok<INUSD>
        % Get Data
        min_spacing = getappdata(handles_gui.figure,'min_spacing');
        max_spacing = getappdata(handles_gui.figure,'max_spacing');
        spacing_prelim = getappdata(handles_gui.figure,'spacing_prelim'); 
        spacing_slider_prelim = getappdata(handles_gui.figure,'spacing_slider_prelim'); 
        get_spacing = getappdata(handles_gui.figure,'get_spacing');
        
        % Get spacing num
        spacing_buffer = str2double(get(handles_gui.edit_spacing,'string')); 
        if (ncorr_util_isintbb(spacing_buffer,min_spacing,max_spacing,'Spacing') == out.success)
            spacing_prelim = spacing_buffer;
            spacing_slider_prelim = min(1,max(0,fzero(@(x)get_spacing(x)-spacing_prelim,0)));
        end
        
        % Set data
        setappdata(handles_gui.figure,'spacing_prelim',spacing_prelim); 
        setappdata(handles_gui.figure,'spacing_slider_prelim',spacing_slider_prelim); 
                
        % Update
        update_axes('update');
        update_sidemenu();  
    end

    function callback_edit_cutoff_diffnorm(hObject,eventdata) %#ok<INUSD>
        % Get data
        min_cutoff_diffnorm = getappdata(handles_gui.figure,'min_cutoff_diffnorm');
        max_cutoff_diffnorm = getappdata(handles_gui.figure,'max_cutoff_diffnorm');
        cutoff_diffnorm_prelim = getappdata(handles_gui.figure,'cutoff_diffnorm_prelim');  
        
        % Get diffnorm
        cutoff_diffnorm_buffer = str2double(get(handles_gui.edit_cutoff_diffnorm,'string')); 
        if (ncorr_util_isrealbb(cutoff_diffnorm_buffer,min_cutoff_diffnorm,max_cutoff_diffnorm,'Difference norm cutoff') == out.success)
            cutoff_diffnorm_prelim = cutoff_diffnorm_buffer;
        end
        
        % Set data
        setappdata(handles_gui.figure,'cutoff_diffnorm_prelim',cutoff_diffnorm_prelim); 
        
        % Update
        update_sidemenu();  
    end

    function callback_edit_cutoff_iteration(hObject,eventdata) %#ok<INUSD>
        % Get data
        min_cutoff_iteration = getappdata(handles_gui.figure,'min_cutoff_iteration');
        max_cutoff_iteration = getappdata(handles_gui.figure,'max_cutoff_iteration');
        cutoff_iteration_prelim = getappdata(handles_gui.figure,'cutoff_iteration_prelim');
        
        % Get iterations
        cutoff_iteration_buffer = str2double(get(handles_gui.edit_cutoff_iteration,'string')); 
        if (ncorr_util_isintbb(cutoff_iteration_buffer,min_cutoff_iteration,max_cutoff_iteration,'Iteration number cutoff') == out.success)
            cutoff_iteration_prelim = cutoff_iteration_buffer;
        end
        
        % Set data
        setappdata(handles_gui.figure,'cutoff_iteration_prelim',cutoff_iteration_prelim);  
        
        % Update
        update_sidemenu();  
    end

    function callback_edit_total_threads(hObject,eventdata) %#ok<INUSD>
        % Get data
        min_total_threads = getappdata(handles_gui.figure,'min_total_threads');
        max_total_threads = getappdata(handles_gui.figure,'max_total_threads');
        total_threads_prelim = getappdata(handles_gui.figure,'total_threads_prelim');
        
        % Get total threads
        total_threads_buffer = str2double(get(handles_gui.edit_thread_total,'string')); 
        if (ncorr_util_isintbb(total_threads_buffer,min_total_threads,max_total_threads,'Total number of threads') == out.success)
            total_threads_prelim = total_threads_buffer;
        end
        
        % Set data
        setappdata(handles_gui.figure,'total_threads_prelim',total_threads_prelim);  
        
        % Update
        update_sidemenu();  
    end

    function callback_edit_leapfrog(hObject,eventdata) %#ok<INUSD>
        % Get data
        min_step = getappdata(handles_gui.figure,'min_step');
        max_step = getappdata(handles_gui.figure,'max_step');
        stepanalysis_prelim = getappdata(handles_gui.figure,'stepanalysis_prelim');
        
        % Get leapfrog
        step_buffer = str2double(get(handles_gui.edit_leapfrog,'string')); 
        if (ncorr_util_isintbb(step_buffer,min_step,max_step,'Step number') == out.success)
            stepanalysis_prelim.step = step_buffer;
        end
        
        % Set data
        setappdata(handles_gui.figure,'stepanalysis_prelim',stepanalysis_prelim);  
        
        % Update
        update_sidemenu();  
    end

    function callback_edit_imgnum(hObject,eventdata) %#ok<INUSD>
        % Get Data
        gs_buffer = getappdata(handles_gui.figure,'gs_buffer');
        num_img = getappdata(handles_gui.figure,'num_img');
        pos_points = getappdata(handles_gui.figure,'pos_points');

        % Get img num - this is one based indexing
        num_img_prelim = str2double(get(handles_gui.edit_imgnum,'string')); 
        if (ncorr_util_isintbb(num_img_prelim,1,length(imgs),'Image number') == out.success)  
            % Clear position of points and set new gs_buffer
            % convert num_img_prelim to zero based indexing
            num_img = num_img_prelim-1;
            gs_buffer = imgs(num_img+1).get_gs();
            pos_points = [];   
        end    
        
        % Set Data
        setappdata(handles_gui.figure,'gs_buffer',gs_buffer);
        setappdata(handles_gui.figure,'num_img',num_img);
        setappdata(handles_gui.figure,'pos_points',pos_points);
            
        % Update
        update_axes('set');
    end

    function callback_button_finish(hObject,eventdata) %#ok<INUSD>
        % Get data
        radius_prelim = getappdata(handles_gui.figure,'radius_prelim');
        spacing_prelim = getappdata(handles_gui.figure,'spacing_prelim');
        cutoff_diffnorm_prelim = getappdata(handles_gui.figure,'cutoff_diffnorm_prelim');
        cutoff_iteration_prelim = getappdata(handles_gui.figure,'cutoff_iteration_prelim');
        total_threads_prelim = getappdata(handles_gui.figure,'total_threads_prelim');
        stepanalysis_prelim = getappdata(handles_gui.figure,'stepanalysis_prelim');
        subsettrunc_prelim = getappdata(handles_gui.figure,'subsettrunc_prelim');
        
        % Check to make sure reduced ROI(s) are nonempty - can happen if
        % spacing is large; this really only happens if it is intended
        allfull = true;
        for i = 0:length(rois)-1
            roi_reduced_buffer = rois(i+1).reduce(spacing_prelim);
            if (roi_reduced_buffer.get_fullregions() == 0)
                allfull = false;
                break;
            end
        end        
        
        if (allfull)
            % Next, make sure the ROI(s) do not have points near the edges 
            % of the image; this can mess up the boundary updating algorithm
            % since updated boundary points near the image edges are
            % discarded.
            % "Near" is defined to be a distance of "radius_prelim" within the
            % boundary.        
            nearboundary = false;
            for i = 0:length(rois)-1
                if (any(any(rois(i+1).mask(1:min(radius_prelim,end),:))) || ...
                    any(any(rois(i+1).mask(max(end-radius_prelim+1,1):end,:))) || ...
                    any(any(rois(i+1).mask(:,1:min(radius_prelim,end)))) || ...
                    any(any(rois(i+1).mask(:,max(end-radius_prelim+1,1):end))))
                    nearboundary = true;
                    break;
                end
            end
            
            % Display error but let user continue
            if (nearboundary)
                h_error = errordlg('There are points in the ROI near the edges of the image; these points can cause problems when the ROI is updated. Please draw the ROI smaller, or possibly reduce the subset radius. You can proceed, but there may be problems when formatting the displacements.','Error','modal');
                uiwait(h_error);
            end
            
            % Ask user if inputs are correct
            if (stepanalysis_prelim.enabled)
                text_step = 'Enabled';
            else
                text_step = 'Disabled';
            end
            if (subsettrunc_prelim)
                text_st = 'Enabled';
            else
                text_st = 'Disabled';
            end
            contbutton = questdlg({['Radius is set to: ' num2str(radius_prelim) '.'], ...
                                   ['Spacing is set to: ' num2str(spacing_prelim) '.'], ... 
                                   ['Difference Norm Cutoff is set to: ' num2str(cutoff_diffnorm_prelim) '.'], ...
                                   ['Iteration Number Cutoff is set to: ' num2str(cutoff_iteration_prelim) '.'], ...
                                   ['Number of threads is set to: ' num2str(total_threads_prelim) '.'], ...
                                   ['Step analysis is: ' text_step '.'], ...
                                   ['Subset Truncation is: ' text_st '.'], ...
                                   'Is this correct?'},'Continue Operation','Yes','No','No');
            if (strcmp(contbutton,'Yes'))                
                % Set outputs
                radius = radius_prelim;
                spacing = spacing_prelim;
                cutoff_diffnorm = cutoff_diffnorm_prelim;
                cutoff_iteration = cutoff_iteration_prelim;
                total_threads = total_threads_prelim;
                stepanalysis = stepanalysis_prelim;
                subsettrunc = subsettrunc_prelim;
                outstate = out.success;

                % Exit                
                close(handles_gui.figure);
            end
        else
            h_error = errordlg('Spacing caused a region to become empty. Either decrease the spacing or redraw the ROI.','Error','modal');
            uiwait(h_error);
        end
    end

    function callback_button_cancel(hObject,eventdata) %#ok<INUSD>
        close(handles_gui.figure); 
    end

    function callback_button_left(hObject,eventdata) %#ok<INUSD>
        % Get Data
        gs_buffer = getappdata(handles_gui.figure,'gs_buffer');
        num_img = getappdata(handles_gui.figure,'num_img');
        pos_points = getappdata(handles_gui.figure,'pos_points');
        
        % Check for overshoot
        if (num_img > 0)
            % Clear position of points and set new gs_buffer
            num_img = num_img-1;
            gs_buffer = imgs(num_img+1).get_gs();
            pos_points = [];             
        end        
        
        % Set Data
        setappdata(handles_gui.figure,'gs_buffer',gs_buffer);
        setappdata(handles_gui.figure,'num_img',num_img);  
        setappdata(handles_gui.figure,'pos_points',pos_points);  
            
        % Update
        update_axes('set');
    end

    function callback_button_right(hObject,eventdata) %#ok<INUSD>
        % Get Data
        gs_buffer = getappdata(handles_gui.figure,'gs_buffer');
        num_img = getappdata(handles_gui.figure,'num_img');
        pos_points = getappdata(handles_gui.figure,'pos_points');
        
        % Check for overshoot
        if (num_img < length(imgs)-1)  
            % Clear position of points and set new gs_buffer
            num_img = num_img+1;
            gs_buffer = imgs(num_img+1).get_gs();
            pos_points = [];             
        end                 
        
        % Set Data
        setappdata(handles_gui.figure,'gs_buffer',gs_buffer);
        setappdata(handles_gui.figure,'num_img',num_img);
        setappdata(handles_gui.figure,'pos_points',pos_points);
            
        % Update
        update_axes('set');
    end

    function callback_impoint(pos,num_region)
        % Get data
        pos_points = getappdata(handles_gui.figure,'pos_points');
        
        % Store position point - convert to zero based indexing and round
        pos_points(num_region+1,:) = round(pos)-1;
        
        % Set Data 
        setappdata(handles_gui.figure,'pos_points',pos_points);
        setappdata(handles_gui.figure,'num_region',num_region);
        
        % Update
        update_axes('update');        
    end

    function update_sidemenu()
        % Get data        
        radius_prelim = getappdata(handles_gui.figure,'radius_prelim');
        spacing_prelim = getappdata(handles_gui.figure,'spacing_prelim');
        radius_slider_prelim = getappdata(handles_gui.figure,'radius_slider_prelim');
        spacing_slider_prelim = getappdata(handles_gui.figure,'spacing_slider_prelim');
        cutoff_diffnorm_prelim = getappdata(handles_gui.figure,'cutoff_diffnorm_prelim');
        cutoff_iteration_prelim = getappdata(handles_gui.figure,'cutoff_iteration_prelim');
        total_threads_prelim = getappdata(handles_gui.figure,'total_threads_prelim');
        stepanalysis_prelim = getappdata(handles_gui.figure,'stepanalysis_prelim');
        subsettrunc_prelim = getappdata(handles_gui.figure,'subsettrunc_prelim');
         
        % Disable total_threads if there is no open_mp support
        if (~support_openmp)
            set(handles_gui.edit_thread_total,'Enable','off');
        end
        
        % Enable seed analysis control panel  
        if (stepanalysis_prelim.enabled)            
            % Enable menu 
            set(handles_gui.popupmenu,'Enable','on');
            set(handles_gui.checkbox_seedauto,'Enable','on');
            if (strcmp(stepanalysis_prelim.type,'leapfrog'))
                set(handles_gui.edit_leapfrog,'Enable','on');    
            else
                set(handles_gui.edit_leapfrog,'Enable','off');  
            end
        else            
            % Disable menu
            set(handles_gui.popupmenu,'Enable','off');
            set(handles_gui.checkbox_seedauto,'Enable','off');
            set(handles_gui.edit_leapfrog,'Enable','off');          
        end
        
        % Set sliders
        set(handles_gui.slider_radius,'Value',radius_slider_prelim); 
        set(handles_gui.slider_spacing,'Value',spacing_slider_prelim); 

        % Set edit
        set(handles_gui.edit_radius,'String',num2str(radius_prelim));
        set(handles_gui.edit_spacing,'String',num2str(spacing_prelim));
        set(handles_gui.edit_cutoff_diffnorm,'String',num2str(cutoff_diffnorm_prelim));
        set(handles_gui.edit_cutoff_iteration,'String',num2str(cutoff_iteration_prelim));
        set(handles_gui.edit_thread_total,'String',num2str(total_threads_prelim));
        set(handles_gui.edit_leapfrog,'String',num2str(stepanalysis_prelim.step));  
                      
        % Set values
        set(handles_gui.checkbox_stepanalysis,'Value',stepanalysis_prelim.enabled);   
        set(handles_gui.checkbox_seedauto,'Value',stepanalysis_prelim.auto);  
        if (strcmp(stepanalysis_prelim.type,'seed'))
            set(handles_gui.popupmenu,'Value',1);
        else
            set(handles_gui.popupmenu,'Value',2);
        end
        set(handles_gui.checkbox_subsettrunc,'Value',subsettrunc_prelim);   
    end

    function update_axes(action)
        % Get data            
        radius_prelim = getappdata(handles_gui.figure,'radius_prelim');
        spacing_prelim = getappdata(handles_gui.figure,'spacing_prelim');
        subsettrunc_prelim = getappdata(handles_gui.figure,'subsettrunc_prelim');
        gs_buffer = getappdata(handles_gui.figure,'gs_buffer');
        preview_subsetloc = getappdata(handles_gui.figure,'preview_subsetloc');    
        preview_subset = getappdata(handles_gui.figure,'preview_subset');    
        handle_subsetloc = getappdata(handles_gui.figure,'handle_subsetloc');       
        handle_subset = getappdata(handles_gui.figure,'handle_subset');    
        handle_points = getappdata(handles_gui.figure,'handle_points');
        num_img = getappdata(handles_gui.figure,'num_img'); 
        num_region = getappdata(handles_gui.figure,'num_region'); 
        pos_points = getappdata(handles_gui.figure,'pos_points');
        
        if (strcmp(action,'set'))            
            % Form previews for subset location and the subset
            preview_subsetloc = gs_buffer;
            preview_subsetloc(rois(num_img+1).mask) = preview_subsetloc(rois(num_img+1).mask)+imgs(num_img+1).max_gs;     
            preview_subset = zeros(2*radius_prelim+1);        

            % Set initial subset location image
            handle_subsetloc = imshow(preview_subsetloc,[imgs(num_img+1).min_gs 3*imgs(num_img+1).max_gs],'Parent',handles_gui.axes_subsetloc);
            set(handles_gui.axes_subsetloc,'Visible','off');
            
            % Set initial subset image
            handle_subset = imshow(preview_subset,[imgs(num_img+1).min_gs 1.25*imgs(num_img+1).max_gs],'Parent',handles_gui.axes_subset);
            set(handles_gui.axes_subset,'Visible','off');
            
            % Form points (one for each region). Note that at this point 
            % each region is assumed to be full so each region will
            % have a gauranteed point. Also, if pos_points is empty, that
            % means the positions for the points needs to be reset; if it's 
            % not, then set the positions of the points to the positions 
            % stored in pos_points
            if (isempty(pos_points))
                setnewpos = true;
            else
                setnewpos = false;
            end
            
            handle_points = impoint.empty;
            for i = 0:length(rois(num_img+1).region)-1
                % Form point
                handle_points(i+1) = impoint(handles_gui.axes_subsetloc,round(size(gs_buffer,2)/2),round(size(gs_buffer,1)/2));            
                setColor(handle_points(i+1),'g');
                set(handle_points(i+1),'UIContextMenu','');
                
                % Set constraint function
                api_point = iptgetapi(handle_points(i+1));
                constrainfcn = ncorr_util_formregionconstraint(rois(num_img+1).region(i+1));
                api_point.setPositionConstraintFcn(constrainfcn);                
                
                % Set position of point
                if (setnewpos)
                    % Initialize point location - use 1 based indexing for input to constrainfcn  
                    x_center = round(rois(num_img+1).region(i+1).leftbound+(rois(num_img+1).region(i+1).rightbound-rois(num_img+1).region(i+1).leftbound)/2)+1;
                    y_center = round(rois(num_img+1).region(i+1).upperbound+(rois(num_img+1).region(i+1).lowerbound-rois(num_img+1).region(i+1).upperbound)/2)+1;
                
                    % Note constrainfcn returns a point using 1 based 
                    % indexing, so convert it to zero based to keep it
                    % consistent with the points stored in pos_points.
                    pos_points(i+1,1:2) = constrainfcn([x_center y_center])-1;      
                end        
                % Set position based on 1 based indexing
                setPosition(handle_points(i+1),pos_points(i+1,1:2)+1);
                
                % Assign Position callback
                addNewPositionCallback(handle_points(i+1),ncorr_util_wrapcallbacktrycatch(@(pos)callback_impoint(pos,i),handles_gui.figure));
            end
            
            % Set num_region to the last region if new positions were used.
            if (setnewpos)
                num_region = length(rois(num_img+1).region)-1;
            end       
                          
            % Format/show text
            set(handles_gui.text_name,'String',['Name: ' imgs(num_img+1).name(1:min(end,22))]);
            
            % Set buttons/editboxes
            set(handles_gui.edit_imgnum,'String',num2str(num_img+1));
            if (length(imgs) == 1)
                set(handles_gui.button_left,'Enable','off');
                set(handles_gui.button_right,'Enable','off');
                set(handles_gui.edit_imgnum,'Enable','off');
            elseif (num_img == 0)
                set(handles_gui.button_left,'Enable','off');
                set(handles_gui.button_right,'Enable','on');
                set(handles_gui.edit_imgnum,'Enable','on');
            elseif (num_img == length(imgs)-1)
                set(handles_gui.button_left,'Enable','on');
                set(handles_gui.button_right,'Enable','off');
                set(handles_gui.edit_imgnum,'Enable','on');
            else
                set(handles_gui.button_left,'Enable','on');
                set(handles_gui.button_right,'Enable','on');   
                set(handles_gui.edit_imgnum,'Enable','on');           
            end   
        end
        
        if (strcmp(action,'set') || strcmp(action,'update'))
            % Initialize CData for preview figures - Use buffer so original
            % preview_subsetloc doesn't get overwritten.
            preview_subsetloc_buffer = preview_subsetloc;
            preview_subset(:) = 0;          

            % Get cirroi - note that pos_points is zero based
            cirroi = rois(num_img+1).get_cirroi(pos_points(num_region+1,1),pos_points(num_region+1,2),radius_prelim,subsettrunc_prelim);    

            % Form subset and subset location images
            for i = 0:size(cirroi.region.noderange,1)-1
                x = i+(cirroi.x-cirroi.radius);
                for j = 0:2:cirroi.region.noderange(i+1)-1
                    y_vec = cirroi.region.nodelist(i+1,j+1):cirroi.region.nodelist(i+1,j+2);
                    preview_subsetloc_buffer(y_vec+1,x+1) = gs_buffer(y_vec+1,x+1) + 2*imgs(num_img+1).max_gs;
                    preview_subset(y_vec-(cirroi.y-cirroi.radius)+1,i+1) = gs_buffer(y_vec+1,x+1) + 0.25*imgs(num_img+1).max_gs;                            
                end     
            end    

            % Paint points over subset to illustrate subset spacing
            % Initialize column counter
            counter_column = mod(cirroi.x-cirroi.radius,spacing_prelim+1);
            if (counter_column == 0)
                counter_column = spacing_prelim; 
            else
                counter_column = counter_column-1;
            end
            
            for i = 0:size(cirroi.region.noderange,1)-1
                counter_column = counter_column+1;
                if (counter_column == spacing_prelim+1) % This a column contains a grid point - reset the counter and then use counter == 0 as a condition for a grid point
                    counter_column = 0; 
                end
                
                x = i+(cirroi.x-cirroi.radius);
                for j = 0:2:cirroi.region.noderange(i+1)-1
                    % Initialize row counter
                    if (counter_column == 0) % This is a column which contains a grid point
                        counter_row = mod(cirroi.region.nodelist(i+1,j+1),spacing_prelim+1);
                        if (counter_row == 0)
                            counter_row = spacing_prelim;
                        else
                            counter_row = counter_row-1;
                        end
                    end
                    
                    for k = cirroi.region.nodelist(i+1,j+1):cirroi.region.nodelist(i+1,j+2) 
                        y = k;
                        if (counter_column == 0) % This is a column which contacts a grid point
                            counter_row = counter_row + 1;
                            if (counter_row == spacing_prelim+1) % This is a row which contains a grid point   
                                preview_subset(y-(cirroi.y-cirroi.radius)+1,i+1) = gs_buffer(y+1,x+1) + 0.5*imgs(num_img+1).max_gs;
                                counter_row = 0;
                            end
                        end
                    end            
                end  
            end      
            
            % Set images
            set(handle_subsetloc,'CData',preview_subsetloc_buffer)
            set(handle_subset,'CData',preview_subset);
        end
        
        % Set data
        setappdata(handles_gui.figure,'preview_subsetloc',preview_subsetloc); 
        setappdata(handles_gui.figure,'preview_subset',preview_subset); 
        setappdata(handles_gui.figure,'handle_subsetloc',handle_subsetloc); 
        setappdata(handles_gui.figure,'handle_subset',handle_subset); 
        setappdata(handles_gui.figure,'handle_points',handle_points); 
        setappdata(handles_gui.figure,'num_region',num_region); 
        setappdata(handles_gui.figure,'pos_points',pos_points);
    end

    function handles_gui = init_gui()
    % GUI controls -------------------------------------------------------%
        % Figure
        handles_gui.figure = figure( ...
            'Tag', 'figure', ...
            'Units', 'characters', ...
            'Position', ncorr_util_figpos(pos_parent,[46.7 167]), ...
            'Name', 'Set DIC Parameters', ...
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

        % Panels
        handles_gui.group_subsetoptions = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_subsetoptions', ...
            'Units', 'characters', ...
            'Position', [2 36.7 35 9.4], ...
            'Title', 'Subset Options', ...
            'Interruptible','off');
        
        handles_gui.group_iterativeoptions = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_iterativeoptions', ...
            'Units', 'characters', ...
            'Position', [2 29.8 35 6.2], ...
            'Title', 'Iterative Solver Options', ...
            'Interruptible','off');
        
        handles_gui.group_multithreading = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_multithreading', ...
            'Units', 'characters', ...
            'Position', [2 24.8 35 4.3], ...
            'Title', 'Multithreading Options', ...
            'Interruptible','off');
        
        handles_gui.group_stepanalysis = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_stepanalysis', ...
            'Units', 'characters', ...
            'Position', [2 13.2 35 10.8], ...
            'Title', 'High Strain Analysis', ...
            'Interruptible','off');      
        
        handles_gui.group_subsettrunc = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_subsettrunc', ...
            'Units', 'characters', ...
            'Position', [2 8.2 35 4.3], ...
            'Title', 'Discontinuous Analysis', ...
            'Interruptible','off');  
        
        handles_gui.group_menu = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_menu', ...
            'Units', 'characters', ...
            'Position', [2 0.8 35 6.7], ...
            'Title', 'Menu', ...
            'Interruptible','off');

        handles_gui.group_preview = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_preview', ...
            'Units', 'characters', ...
            'Position', [39 0.75 126 45.3], ...
            'Title', 'Subset Radius and Spacing Preview', ...
            'Interruptible','off');

        % Axes
        handles_gui.axes_subsetloc = axes( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'axes_subsetloc', ...
            'Units', 'characters', ...
            'Position', [2 4 57.2 38.2], ...
            'Visible','off', ...
            'Interruptible','off');

        handles_gui.axes_subset = axes( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'axes_subset', ...
            'Units', 'characters', ...
            'Position', [65.4 4 57.2 38.2], ...
            'Visible','off', ...
            'Interruptible','off');
        
        % Drop-down Menu  
        handles_gui.popupmenu = uicontrol( ...
            'Parent', handles_gui.group_stepanalysis, ...
            'Tag', 'popupmenu', ...
            'Style', 'popupmenu', ...
            'Units', 'characters', ...
            'Position', [2.4 5.5 29.1 1.3], ...
            'BackgroundColor', [1 1 1], ...
            'String', {'Seed Propagation','Leapfrog'}, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_popupmenu,handles_gui.figure), ...
            'Interruptible','off');
        
        % Static Texts
        handles_gui.text_subsetloc = uicontrol( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'text_subsetloc', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.3 42.4 20 1.3], ...
            'String', 'Subset Location: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_subset = uicontrol( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'text_subset', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [65.4 42.4 20 1.3], ...
            'String', 'Subset: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_radius = uicontrol( ...
            'Parent', handles_gui.group_subsetoptions, ...
            'Tag', 'text_radius', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 6.0 16 1.3], ...
            'String', 'Subset Radius: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_spacing = uicontrol( ...
            'Parent', handles_gui.group_subsetoptions, ...
            'Tag', 'text_spacing', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 2.5 16 1.3], ...
            'String', 'Subset Spacing: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_cutoff_diffnorm = uicontrol( ...
            'Parent', handles_gui.group_iterativeoptions, ...
            'Tag', 'text_cutoff_diffnorm', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 2.8 16 1.3], ...
            'String', 'Diff Norm C/O:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_cutoff_iteration = uicontrol( ...
            'Parent', handles_gui.group_iterativeoptions, ...
            'Tag', 'text_cutoff_iteration', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 0.9 16 1.3], ...
            'String', 'Iteration # C/O:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_thread_total = uicontrol( ...
            'Parent', handles_gui.group_multithreading, ...
            'Tag', 'text_thread_total', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 0.9 16 1.3], ...
            'String', 'Num Threads:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');        
        
        handles_gui.text_name = uicontrol( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'text_name', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [4 1.3 56.9 1.3], ...
            'String', 'Name: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');   
        
        handles_gui.text_leapfrog = uicontrol( ...
            'Parent', handles_gui.group_stepanalysis, ...
            'Tag', 'text_leapfrog', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 0.9 20 1.3], ...
            'String', 'Step #: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');        
        
        % Sliders
        handles_gui.slider_radius = uicontrol( ...
            'Parent', handles_gui.group_subsetoptions, ...
            'Tag', 'slider_radius', ...
            'Style', 'slider', ...
            'Units', 'characters', ...
            'Position', [2.5 4.5 28.9 1.3], ...
            'BackgroundColor', [0.9 0.9 0.9], ...
            'String', {'Slider'}, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_slider_radius,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.slider_spacing = uicontrol( ...
            'Parent', handles_gui.group_subsetoptions, ...
            'Tag', 'slider_spacing', ...
            'Style', 'slider', ...
            'Units', 'characters', ...
            'Position', [2.5 1 28.9 1.3], ...
            'BackgroundColor', [0.9 0.9 0.9], ...
            'String', {'Slider'}, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_slider_spacing,handles_gui.figure), ...
            'Interruptible','off');
        
        % Checkbox
        handles_gui.checkbox_stepanalysis = uicontrol( ...
            'Parent', handles_gui.group_stepanalysis, ...
            'Tag', 'checkbox_stepanalysis', ...
            'Style', 'checkbox', ...
            'Units', 'characters', ...
            'Position', [2.3 7.5 29.7 1.3], ...
            'String', 'Enable Step Analysis', ...
            'Value', 0, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_checkbox_stepanalysis,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.checkbox_seedauto = uicontrol( ...
            'Parent', handles_gui.group_stepanalysis, ...
            'Tag', 'checkbox_seedauto', ...
            'Style', 'checkbox', ...
            'Units', 'characters', ...
            'Position', [2.3 3.0 26.7 1.3], ...
            'String', 'Auto Propagation', ...
            'Value', 0, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_checkbox_seedauto,handles_gui.figure), ...
            'Interruptible','off');     
        
        handles_gui.checkbox_subsettrunc = uicontrol( ...
            'Parent', handles_gui.group_subsettrunc, ...
            'Tag', 'checkbox_subsettrunc', ...
            'Style', 'checkbox', ...
            'Units', 'characters', ...
            'Position', [2.3 1.0 26.7 1.3], ...
            'String', 'Subset Truncation', ...
            'Value', 0, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_checkbox_subsettrunc,handles_gui.figure), ...
            'Interruptible','off');     
        
        % Edit
        handles_gui.edit_radius = uicontrol( ...
            'Parent', handles_gui.group_subsetoptions, ...
            'Tag', 'edit_radius', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 6.1 11.1 1.3], ...
            'String', '10', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_radius,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.edit_spacing = uicontrol( ...
            'Parent', handles_gui.group_subsetoptions, ...
            'Tag', 'edit_spacing', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 2.6 11.1 1.3], ...
            'String', '0', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_spacing,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.edit_cutoff_diffnorm = uicontrol( ...
            'Parent', handles_gui.group_iterativeoptions, ...
            'Tag', 'edit_cutoff_diffnorm', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 2.9 11.1 1.3], ...
            'String', '0.000001', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_cutoff_diffnorm,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.edit_cutoff_iteration = uicontrol( ...
            'Parent', handles_gui.group_iterativeoptions, ...
            'Tag', 'edit_cutoff_iteration', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 1.0 11.1 1.3], ...
            'String', '50', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_cutoff_iteration,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.edit_thread_total = uicontrol( ...
            'Parent', handles_gui.group_multithreading, ...
            'Tag', 'edit_thread_total', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 1.0 11.1 1.3], ...
            'String', '1', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_total_threads,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.edit_leapfrog = uicontrol( ...
            'Parent', handles_gui.group_stepanalysis, ...
            'Tag', 'edit_leapfrog', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 1.0 11.1 1.3], ...
            'String', '0', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_leapfrog,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.edit_imgnum = uicontrol( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'edit_imgnum', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [108.2 1.3 7 1.6], ...
            'String', '', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_imgnum,handles_gui.figure), ...
            'Enable', 'off', ...
            'Interruptible','off');
        
        % Pushbuttons
        handles_gui.button_finish = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'button_finish', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 3 29.7 1.7], ...
            'String', 'Finish', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_finish,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_cancel = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'button_cancel', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 1 29.7 1.7], ...
            'String', 'Cancel', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_cancel,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_left = uicontrol( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'button_left', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [101.2 1.2 6 1.8], ...
            'String', '<', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_left,handles_gui.figure), ...
            'Interruptible','off', ...
            'Enable', 'off');

        handles_gui.button_right = uicontrol( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'button_right', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [116.2 1.2 6 1.8], ...
            'String', '>', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_right,handles_gui.figure), ...
            'Interruptible','off', ...
            'Enable', 'off');
    end

    % Pause until figure is closed ---------------------------------------%
    waitfor(handles_gui.figure);
end

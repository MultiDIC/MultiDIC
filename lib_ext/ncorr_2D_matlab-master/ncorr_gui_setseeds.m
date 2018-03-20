function [seedinfo,threaddiagram,outstate] = ncorr_gui_setseeds(reference,current,roi,num_region,spacing,radius,cutoff_diffnorm,cutoff_iteration,total_threads,enabled_stepanalysis,subsettrunc,num_img,total_imgs,pos_parent)
% This is a GUI for setting the rgdic seeds. The number of seeds placed in
% each region depends on the number of threads. 
%
% Inputs -----------------------------------------------------------------%
%   reference - ncorr_class_img; used for displaying the background
%   image and calculations
%   current - ncorr_class_img(s); used for displaying the background image
%   and calculations
%   roi - ncorr_class_roi; roi corresponding to the reference image.
%   num_region - integer; number corresponding to the region; 
%   spacing - integer; spacing between subsets
%   radius - integer; radius of subsets
%   cutoff_diffnorm - double; cutoff of norm of the difference vector
%   cutoff_iteration - integer; cutoff for the number of iterations
%   total_threads - integer; total number of threads
%   enabled_stepanalysis - logical; if true, then process as many images as
%   possible. If false, process all images.
%   subsettrunc - logical; if true, then subset truncation is enabled
%   num_img - integer; current number of reference image
%   total_imgs - integer; total number of images being processed
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%
% Outputs ----------------------------------------------------------------%
%   seedinfo - struct; contains
%   struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints'
%   ,{}) 
%   threaddiagram - integer array; array the same size as the reduced
%   ROI. The elements indicate which points will be computed by which
%   threads
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.
%
% Note that if step analysis is enabled, it is gauranteed that at least one 
% image will be succesfully seeded if outstate returns true. If step 
% analysis is disabled, it is gauranteed that all images will be seeded or 
% if outstate returns true.

    % Data ---------------------------------------------------------------%     
    % Initialize Outputs 
    outstate = out.cancelled;
    threaddiagram = [];
    seedinfo = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{}); % paramvector = [x y u v du/dx du/dy dv/dx dv/dy corrcoef]
    % Get GUI handles    
    handles_gui = init_gui();
    % Run c-tor
    feval(ncorr_util_wrapcallbacktrycatch(@constructor,handles_gui.figure));
    
    % Callbacks and functions --------------------------------------------%
    function constructor()                        
        % Initialize buffers
        % Initialize reduced ref
        ref_reduced = reference.reduce(spacing);
        % Initialize roi_reduced
        roi_reduced = roi.reduce(spacing);
        % Get region mask
        regionmask = roi_reduced.get_regionmask(num_region);   
        % Initialize threaddiagram_prelim and threaddiagram preview - these
        % are modified in-place
        threaddiagram_prelim = -ones(size(roi_reduced.mask));
        preview_threaddiagram = zeros(size(roi_reduced.mask));
        
        % Store data   
        setappdata(handles_gui.figure,'ref_reduced',ref_reduced);  
        setappdata(handles_gui.figure,'roi_reduced',roi_reduced);     
        setappdata(handles_gui.figure,'handle_preview',[]);    
        setappdata(handles_gui.figure,'handle_points',[]);          
        setappdata(handles_gui.figure,'regionmask',regionmask);        
        setappdata(handles_gui.figure,'threaddiagram_prelim',threaddiagram_prelim);       
        setappdata(handles_gui.figure,'preview_threaddiagram',preview_threaddiagram);   
        
        % Update thread diagram and preview - must do this after setting the data
        get_threaddiagram([]);        
        
        % Update
        update_axes('set');
        update_sidemenu();
        
        % Set visible        
        set(handles_gui.figure,'Visible','on');         
    end    

    function callback_button_setseeds(hObject,eventdata) %#ok<INUSD>
        freeze_menu();
        
        % Get data
        roi_reduced = getappdata(handles_gui.figure,'roi_reduced');
        
        % Set seeds ------------------------------------------------------%        
        % Initialize point handle
        handle_points = impoint.empty;
        for i = 0:total_threads-1
            % Place seed and make sure it's within region
            seedfinished = false;                
            while (~seedfinished)   
                handle_point_buffer = impoint(handles_gui.axes_setseeds);  
                
                % impoint will be empty if esc is pressed
                if (~isempty(handle_point_buffer))   
                    % Get seed position - convert to zero based indexing
                    % and round
                    pos_seed = round(getPosition(handle_point_buffer))-1; 
                    
                    % Make sure point is within the correct region
                    if (roi_reduced.get_num_region(pos_seed(1),pos_seed(2),false(1,length(roi_reduced.region))) == num_region)       
                        % Store and format impoint
                        handle_points(i+1) = handle_point_buffer;
                        setColor(handle_points(i+1),'g');
                        set(handle_points(i+1),'UIContextMenu','');                    

                        % Set constraint function
                        api_point = iptgetapi(handle_points(i+1));
                        constrainfcn = ncorr_util_formregionconstraint(roi_reduced.region(num_region+1));
                        api_point.setPositionConstraintFcn(constrainfcn);    

                        % Assign Position callback
                        addNewPositionCallback(handle_points(i+1),ncorr_util_wrapcallbacktrycatch(@(pos)callback_impoint(pos),handles_gui.figure));
                        
                        % Update threaddiagram and threaddiagram preview           
                        get_threaddiagram(handle_points);

                        % Set data
                        setappdata(handles_gui.figure,'handle_points',handle_points);

                        % Update here to ensure preview is updated as the seeds
                        % are being placed.
                        update_sidemenu();
                        update_axes('update');                        
                        
                        % Exit
                        seedfinished = true;  
                    else
                        delete(handle_point_buffer);
                        
                        h_error = errordlg('Seed point not within region, try again.','Error','modal');     
                        uiwait(h_error);
                    end 
                else
                    % Escape was pressed - delete all points
                    delete(handle_points);
                    handle_points = [];    
                    
                    % Update threaddiagram and threaddiagram preview           
                    get_threaddiagram(handle_points);

                    % Set data
                    setappdata(handles_gui.figure,'handle_points',handle_points);
                    
                    % Update
                    update_sidemenu();
                    update_axes('update');
                    
                    % Return since analysis was cancelled
                    return;
                end
            end
        end
                                
        % Update
        update_sidemenu();
        
        unfreeze_menu();
    end    

    function callback_button_finish(hObject,eventdata) %#ok<INUSD>
        % Get data
        handle_points = getappdata(handles_gui.figure,'handle_points'); 
        threaddiagram_prelim = getappdata(handles_gui.figure,'threaddiagram_prelim');
        
        % Get position - and make sure no two seeds overlap
        pos_seed = [];
        for i = 0:length(handle_points)-1
            % Convert to zero based indexing and round
            pos_buffer = round(getPosition(handle_points(i+1)))-1; 
            pos_seed = vertcat(pos_seed,pos_buffer); %#ok<AGROW>
        end
        if (size(unique(pos_seed,'rows'),1) == size(pos_seed,1))
            % Convert pos_seed to regular coordinates since they are WRT to
            % the reduced ROI at this point
            pos_seed = pos_seed*(spacing+1);
            
            % Get seedinfo and preview of seeds. Note that if step analysis
            % is enabled, then ncorr_gui_seedpreview will return at least one
            % seed if it returns success. If step analysis is disabled, 
            % then ncorr_gui_seedpreview will provide seeds for 
            % all images if it returns true. 
            [seedinfo_prelim,outstate_preview] = ncorr_gui_seedpreview(reference, ...
                                                                       current, ...
                                                                       roi, ...
                                                                       num_region, ...
                                                                       pos_seed, ...
                                                                       radius, ...
                                                                       cutoff_diffnorm, ...
                                                                       cutoff_iteration, ...
                                                                       enabled_stepanalysis, ...
                                                                       subsettrunc, ...
                                                                       num_img, ...
                                                                       total_imgs, ...
                                                                       get(handles_gui.figure,'OuterPosition'));
            if (outstate_preview == out.success)
                % Set outputs
                seedinfo = seedinfo_prelim;
                % Determine number of compute points and append this to
                % seedinfo_buffer
                for i = 0:size(seedinfo,3)-1
                    for j = 0:size(seedinfo,1)-1
                        seedinfo(j+1,1,i+1).computepoints = length(find(threaddiagram_prelim == j));
                    end
                end
                % Assign thread diagram
                threaddiagram = threaddiagram_prelim;
                outstate = out.success;
                
                % Exit                
                close(handles_gui.figure);
            end
        else
            % Two seeds occupy same location, tell user to move
            % generators until they occupy unique locations
            h_error = errordlg('Two or more seeds occupy the same location, please move all seeds to unique locations.','Error','modal');
            uiwait(h_error);
        end
    end

    function callback_button_cancel(hObject,eventdata) %#ok<INUSD>
        close(handles_gui.figure);
    end

    function callback_impoint(pos) %#ok<INUSD>
        % Get data
        handle_points = getappdata(handles_gui.figure,'handle_points');
        
        % Update threaddiagram and threaddiagram preview           
        get_threaddiagram(handle_points);
        
        % Update
        update_axes('update');
    end

    function freeze_menu(hObject,eventdata) %#ok<INUSD>
        % Disable setseeds button
        set(handles_gui.button_setseeds,'Enable','off');  
    end

    function unfreeze_menu(hObject,eventdata) %#ok<INUSD>
    end

    function update_sidemenu()
        % Get data
        handle_points = getappdata(handles_gui.figure,'handle_points');
        
        % Update text
        set(handles_gui.text_seedcount,'String',['Seed(s) Set: ' num2str(length(handle_points)) ' of ' num2str(total_threads)]) 
        
        % Set buttons
        if (length(handle_points) == total_threads)
            % Seeds are finished
            set(handles_gui.button_setseeds,'Enable','off');
            set(handles_gui.button_finish,'Enable','on'); 
        elseif (~isempty(handle_points))
            % Seeds are being placed
            set(handles_gui.button_setseeds,'Enable','off');
            set(handles_gui.button_finish,'Enable','off'); 
        else
            % Seeds are not being placed
            set(handles_gui.button_setseeds,'Enable','on');
            set(handles_gui.button_finish,'Enable','off'); 
        end
    end

    function update_axes(action)
        % Get data    
        handle_preview = getappdata(handles_gui.figure,'handle_preview');  
        preview_threaddiagram = getappdata(handles_gui.figure,'preview_threaddiagram');    
        
        if (strcmp(action,'set'))
            handle_preview = imshow(preview_threaddiagram,[reference.min_gs 2*reference.max_gs],'Parent',handles_gui.axes_setseeds);
            set(handles_gui.axes_setseeds,'Visible','off');
            set(handles_gui.text_name,'String',['Name: ' reference.name(1:min(end,22))]);
        end
        
        if (strcmp(action,'set') || strcmp(action,'update'))
            % Update the overlay image
            set(handle_preview,'CData',preview_threaddiagram)
        end
        
        % Set data 
        setappdata(handles_gui.figure,'handle_preview',handle_preview);   
    end

    function get_threaddiagram(handle_points)
        % Get data
        ref_reduced = getappdata(handles_gui.figure,'ref_reduced');     
        regionmask = getappdata(handles_gui.figure,'regionmask');      
        threaddiagram_prelim = getappdata(handles_gui.figure,'threaddiagram_prelim');    
        preview_threaddiagram = getappdata(handles_gui.figure,'preview_threaddiagram');    
        
        % Draw thread diagram - get generators first
        generators = [];
        for i = 0:length(handle_points)-1
            % Convert to zero based indexing and round
            pos_buffer = round(getPosition(handle_points(i+1)))-1; 
            generators = vertcat(generators,pos_buffer); %#ok<AGROW>
        end
        
        % This function modifies the arrays in-place
        ncorr_alg_formthreaddiagram(threaddiagram_prelim,preview_threaddiagram,int32(generators),regionmask,ref_reduced.formatted());
    end

    function handles_gui = init_gui()
    % GUI controls -------------------------------------------------------%
        % Figure
        handles_gui.figure = figure( ...
            'Tag', 'figure', ...
            'Units', 'characters', ...
            'Position', ncorr_util_figpos(pos_parent,[35 133]), ...
            'Name', 'Set Seeds', ...
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
        handles_gui.group_seedoptions = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_seedoptions', ...
            'Units', 'characters', ...
            'Position', [2 28.5 35 6.0], ...
            'Title', 'Seed Options', ...
            'Interruptible','off');
        
        handles_gui.group_menu = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_menu', ...
            'Units', 'characters', ...
            'Position', [2 21.2 35 6.4], ...
            'Title', 'Menu', ...
            'Interruptible','off');
        
        handles_gui.group_setseeds = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_setseeds', ...
            'Units', 'characters', ...
            'Position', [39.0 0.75 92 33.75], ...
            'Title', 'Select Region', ...
            'Interruptible','off');

        % Axes
        handles_gui.axes_setseeds = axes( ...
            'Parent', handles_gui.group_setseeds, ...
            'Tag', 'axes_setseeds', ...
            'Units', 'characters', ...
            'Position', [2.4 3 86.2 28.8], ...
            'Visible', 'off', ...
            'Interruptible','off');

        % Static Texts
        handles_gui.text_seedcount = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'text_seedcount', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.2 2.8 29.7 1.5], ...
            'String', 'Seeds Set: 0 of ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
        
        handles_gui.text_name = uicontrol( ...
            'Parent', handles_gui.group_setseeds, ...
            'Tag', 'text_name', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [3 0.7 56.9 1.3], ...
            'String', 'Name: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');    

        % Pushbuttons
        handles_gui.button_setseeds = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'button_setseeds', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 1 29.7 1.7], ...
            'String', 'Set Seeds', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_setseeds,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.button_finish = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'button_finish', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 3 29.7 1.7], ...
            'String', 'Finish', ...
            'Enable', 'off', ...
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
    end
        
    % Pause until figure is closed ---------------------------------------%
    waitfor(handles_gui.figure);    
end

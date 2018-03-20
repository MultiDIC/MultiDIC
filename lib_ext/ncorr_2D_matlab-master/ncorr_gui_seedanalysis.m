function [seedinfo,threaddiagram,outstate] = ncorr_gui_seedanalysis(reference,current,roi,radius,spacing,cutoff_diffnorm,cutoff_iteration,total_threads,enabled_stepanalysis,subsettrunc,num_img,total_imgs,pos_parent)
% This is a GUI for selecting the region(s) to place seeds in. It provides 
% all the info needed to run the DIC analysis after it returns.
%
% Inputs -----------------------------------------------------------------%
%   reference - ncorr_class_img; used for displaying the background image 
%   and calculations
%   current - ncorr_class_img(s); used for displaying the background image
%   and calculations
%   roi - ncorr_class_roi; ROI corresponding to the reference image
%   radius - integer; radius of subsets 
%   spacing - integer; spacing used between subsets
%   cutoff_diffnorm - double; cutoff for norm of difference vector
%   cutoff_iteration - integer; cutoff for number of iterations
%   total_threads - integer; number of threads used in analysis
%   enabled_stepanalysis - logical; if true, then process as many seeds as
%   possible. If false, process all the seeds.
%   subsettrunc - logical; if true, then enable subset truncation
%   num_img - integer; reference image number
%   total_imgs - integer; total number of images being processed
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%
% Outputs ----------------------------------------------------------------%
%   seedinfo - struct; contains struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{}) 
%   which is information on the location of the seed, the deformation parameters
%   for that seed, the region the seed is located in, and the thread which
%   will process the seed.
%   threaddiagram - integer array; array the same size as the reduced
%   ROI. The array indicates which points will be computed by which
%   threads
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.
%
% Note that if step analysis is enabled, this function will return at least
% one seed if outstate is set to success. If step analysis is disabled then
% all images will be seeded if outstate is set to success. 

    % Data ---------------------------------------------------------------%     
    % Initialize outputs
    seedinfo = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{}); % paramvector = [x y u v du/dx du/dy dv/dx dv/dy corrcoef]
    threaddiagram = [];
    outstate = out.cancelled; 
    % Get GUI handles    
    handles_gui = init_gui();
    % Run c-tor
    feval(ncorr_util_wrapcallbacktrycatch(@constructor,handles_gui.figure));
    
    % Callbacks and functions --------------------------------------------%
    function constructor()       
        % Callback for escape key can be omitted since it is handled
        % properly within the callback
                
        % Initialize buffers
        % Initialize reduced ref
        ref_reduced = reference.reduce(spacing);
        % gs buffer
        gs_buffer_reduced = ref_reduced.get_gs();
        % Initialize reduced ROI
        roi_reduced = roi.reduce(spacing);
        % Initialize seedinfo buffer
        seedinfo_prelim = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{});
        % Initialize thread diagram buffer
        threaddiagram_prelim = -ones(size(roi_reduced.mask));
        % Initialize list for regions - this keeps track of regions that have
        % been seeded
        list_region = false(length(roi_reduced.region),1);    
        % Number of images successfully seeded
        num_imgs_success = inf;
        
        % Store data
        setappdata(handles_gui.figure,'ref_reduced',ref_reduced);  
        setappdata(handles_gui.figure,'gs_buffer_reduced',gs_buffer_reduced);
        setappdata(handles_gui.figure,'roi_reduced',roi_reduced);     
        setappdata(handles_gui.figure,'seedinfo_prelim',seedinfo_prelim); 
        setappdata(handles_gui.figure,'threaddiagram_prelim',threaddiagram_prelim);        
        setappdata(handles_gui.figure,'list_region',list_region);    
        setappdata(handles_gui.figure,'handle_preview',[]);
        setappdata(handles_gui.figure,'num_imgs_success',num_imgs_success); 

        % Update
        update_axes('set');
        update_sidemenu();

        % Set Visible
        set(handles_gui.figure,'Visible','on'); 
    end

    function callback_button_selectregion(hObject,eventdata) %#ok<INUSD>
        % Freeze menu
        freeze_menu();
        
        % Get Data       
        roi_reduced = getappdata(handles_gui.figure,'roi_reduced');  
        seedinfo_prelim = getappdata(handles_gui.figure,'seedinfo_prelim'); 
        threaddiagram_prelim = getappdata(handles_gui.figure,'threaddiagram_prelim');       
        list_region = getappdata(handles_gui.figure,'list_region');    
        num_imgs_success = getappdata(handles_gui.figure,'num_imgs_success');
        
        % Set loop condition - Exit when set to true
        regionfinished = false;
        while (~regionfinished)
            % Initialize point
            handle_point = impoint(handles_gui.axes_selectregion);            
            
            % handle_point can be empty if escape key is pressed
            if (~isempty(handle_point))          
                % Get position of point - convert to integer with zero based indexing
                pos_region = round(getPosition(handle_point))-1;

                % Make sure point is within reduced ROI 
                num_region_prelim = roi_reduced.get_num_region(pos_region(1),pos_region(2),list_region);  
                if (num_region_prelim ~= -1)         
                    % Format point
                    setColor(handle_point,'g');
                    set(handle_point,'UIContextMenu','');
                    set(handle_point,'ButtonDownFcn','');

                    % If step analysis is enabled, ncorr_gui_setseeds is
                    % gauranteed to return at least one succesfully seeded 
                    % image if it outstate is success. If step analysis is
                    % disabled, ncorr_gui_setseeds will return all seeds if
                    % outstate returns success
                    [seedinfo_buffer,threaddiagram_buffer,outstate_setseeds] = ncorr_gui_setseeds(reference, ...
                                                                                                  current, ...
                                                                                                  roi, ...
                                                                                                  num_region_prelim, ...
                                                                                                  spacing, ...
                                                                                                  radius, ...
                                                                                                  cutoff_diffnorm, ...
                                                                                                  cutoff_iteration, ...
                                                                                                  total_threads, ...
                                                                                                  enabled_stepanalysis, ...
                                                                                                  subsettrunc, ...
                                                                                                  num_img, ...
                                                                                                  total_imgs, ...
                                                                                                  get(handles_gui.figure,'OuterPosition')); 

                    if (outstate_setseeds ~= out.success) 
                        % Analysis for seedplacement was cancelled. Delete
                        % the point.
                        delete(handle_point);  
                        
                        % Exit loop
                        break;
                    end
                  
                    % Take minimum of num_imgs_success and
                    % seedinfo_buffer. Buffer can be more or less than
                    % num_imgs_success
                    num_imgs_success = min(num_imgs_success,size(seedinfo_buffer,3));

                    % Clear out other images in prelim and buffer
                    if (~isempty(seedinfo_prelim))
                        seedinfo_prelim = seedinfo_prelim(:,:,1:num_imgs_success);
                    end
                    seedinfo_buffer = seedinfo_buffer(:,:,1:num_imgs_success);

                    % Append seedinfo_buffer - append along 2nd dimension
                    seedinfo_prelim = horzcat(seedinfo_prelim,seedinfo_buffer); %#ok<AGROW>

                    % Merge threaddiagram from previous iteration
                    threaddiagram_prelim(threaddiagram_buffer ~= -1) = threaddiagram_buffer(threaddiagram_buffer ~= -1);

                    % Update region list
                    list_region(num_region_prelim+1) = true;

                    % Delete region point
                    delete(handle_point);      

                    % Exit loop
                    regionfinished = true;
                else   
                    % Point was not within ROI. Delete it and then let user
                    % place it again.
                    delete(handle_point);
                    
                    h_error = errordlg('Point not within ROI, try again.','Error','modal');    
                    uiwait(h_error);
                end  
            else
                % Escaped was pressed. 
                
                % Exit loop
                regionfinished = true;
            end
        end
        
        % Set data
        setappdata(handles_gui.figure,'seedinfo_prelim',seedinfo_prelim);
        setappdata(handles_gui.figure,'threaddiagram_prelim',threaddiagram_prelim);
        setappdata(handles_gui.figure,'list_region',list_region);    
        setappdata(handles_gui.figure,'num_imgs_success',num_imgs_success);
        
        % Update
        update_axes('set'); 
        update_sidemenu(); 
        
        unfreeze_menu();  
    end

    function callback_button_finish(hObject,eventdata) %#ok<INUSD>
        % Get Info      
        seedinfo_prelim = getappdata(handles_gui.figure,'seedinfo_prelim');    
        threaddiagram_prelim = getappdata(handles_gui.figure,'threaddiagram_prelim');   
        list_region = getappdata(handles_gui.figure,'list_region');     
        
        % See if there are still regions left; if there are, alert the
        % user; otherwise just accept 
        contbutton = 'Yes';
        if (sum(list_region) ~= roi.get_fullregions()) 
            contbutton = questdlg('There are still regions left. Do you want to finish?','Continue Operation','Yes','No','No');
        end
        
        if (strcmp(contbutton,'Yes'))
            % Set outputs
            for i = 0:size(seedinfo_prelim,1)-1
                for j = 0:size(seedinfo_prelim,2)-1
                    for k = 0:size(seedinfo_prelim,3)-1
                        seedinfo(i+1,j+1,k+1) = seedinfo_prelim(i+1,j+1,k+1);
                    end
                end
            end
            threaddiagram = threaddiagram_prelim;      
            outstate = out.success;
            
            % Exit            
            close(handles_gui.figure);
        end
    end

    function callback_button_cancel(hObject,eventdata) %#ok<INUSD>
        close(handles_gui.figure);
    end
        
    function freeze_menu()
        set(handles_gui.button_selectregion,'Enable','off'); 
    end

    function unfreeze_menu()    
    end

    function update_sidemenu()
        % Get Data     
        list_region = getappdata(handles_gui.figure,'list_region');  
                
        % Set number of regions analyzed - use get_fullregions() since some
        % regions might be empty.
        set(handles_gui.text_regioncount,'String',['Region(s) Set: ' num2str(sum(list_region)) ' of ' num2str(roi.get_fullregions())]);  
        
        % Update buttons - Make sure there is at least one seeded region
        % before enabling finish button, in the case that an empty ROI is
        % provided.
        if (any(list_region) && sum(list_region) == roi.get_fullregions())
            % All regions have been finished
            set(handles_gui.button_finish,'Enable','on'); 
            set(handles_gui.button_selectregion,'Enable','off'); 
        elseif (any(list_region))
            % Some regions finished
            set(handles_gui.button_finish,'Enable','on'); 
            set(handles_gui.button_selectregion,'Enable','on');
        else
            % No regions finished            
            set(handles_gui.button_finish,'Enable','off'); 
            set(handles_gui.button_selectregion,'Enable','on');
        end
    end

    function update_axes(action)
        % Get data  
        ref_reduced = getappdata(handles_gui.figure,'ref_reduced'); 
        gs_buffer_reduced = getappdata(handles_gui.figure,'gs_buffer_reduced');
        roi_reduced = getappdata(handles_gui.figure,'roi_reduced');       
        list_region = getappdata(handles_gui.figure,'list_region');    
        handle_preview = getappdata(handles_gui.figure,'handle_preview');       
        
        if (strcmp(action,'set'))
            % Draw figure with ROIs overlayed
            handle_preview = imshow(gs_buffer_reduced,[reference.min_gs 2*reference.max_gs],'Parent',handles_gui.axes_selectregion);
            set(handles_gui.axes_selectregion,'Visible','off');
            set(handles_gui.text_name,'String',['Name: ' reference.name(1:min(end,22))]);
        end
        
        if (strcmp(action,'set') || strcmp(action,'update'))
            preview_roi = gs_buffer_reduced;
            for i = 0:length(roi_reduced.region)-1
                if (list_region(i+1)) 
                    continue; % this region has been analyzed already
                else 
                    for j = 0:size(roi_reduced.region(i+1).noderange,1)-1
                        x = j+roi_reduced.region(i+1).leftbound;
                        for k = 0:2:roi_reduced.region(i+1).noderange(j+1)-1
                            vec_y = roi_reduced.region(i+1).nodelist(j+1,k+1):roi_reduced.region(i+1).nodelist(j+1,k+2);
                            preview_roi(vec_y+1,x+1) = preview_roi(vec_y+1,x+1)+ref_reduced.max_gs;            
                        end  
                    end
                end
            end        
            set(handle_preview,'CData',preview_roi)
        end
        
        % Set data
        setappdata(handles_gui.figure,'handle_preview',handle_preview);
    end

    function handles_gui = init_gui()
    % GUI controls -------------------------------------------------------%
        % Figure
        handles_gui.figure = figure( ...
            'Tag', 'figure', ...
            'Units', 'characters', ...
            'Position', ncorr_util_figpos(pos_parent,[25 103]), ...
            'Name', 'Select Region', ...
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
            'Position', [2 18.3 35 6.2], ...
            'Title', 'Region Options', ...
            'Interruptible','off');
        
        handles_gui.group_menu = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_menu', ...
            'Units', 'characters', ...
            'Position', [2 10.8 35 6.7], ...
            'Title', 'Menu', ...
            'Interruptible','off');
        
        handles_gui.group_selectregion = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_selectregion', ...
            'Units', 'characters', ...
            'Position', [39.0 0.75 62 23.75], ...
            'Title', 'Select Region', ...
            'Interruptible','off');

        % Axes
        handles_gui.axes_selectregion = axes( ...
            'Parent', handles_gui.group_selectregion, ...
            'Tag', 'axes_selectregion', ...
            'Units', 'characters', ...
            'Position', [2.4 3 56.2 18.8], ...
            'Visible', 'off', ...
            'Interruptible','off');

        % Static Texts
        handles_gui.text_regioncount = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'text_regioncount', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.2 2.8 29.7 1.5], ...
            'String', 'Regions Set: 0 of ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_name = uicontrol( ...
            'Parent', handles_gui.group_selectregion, ...
            'Tag', 'text_name', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [3 0.7 56.9 1.3], ...
            'String', 'Name: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');        
        
        % Pushbuttons
        handles_gui.button_selectregion = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'button_selectregion', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 1 29.7 1.7], ...
            'String', 'Select Region', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_selectregion,handles_gui.figure), ...
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

function [rois,outstate] = ncorr_gui_setrois(imgs,pos_parent,params_init)
% This is a GUI for setting the ROIs for either the reference or current
% image(s).
% 
% Inputs -----------------------------------------------------------------%
%   imgs - ncorr_class_img; the input images(s); Generally, it's a single
%   image, although multiple image are supported.
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%   params_init - ncorr_class_roi; contains ROI(s) if they have been set
%   before; otherwise it's empty. Number of input ROI(s) must be the same
%   length as imgs.
%
% Outputs ----------------------------------------------------------------%
%   rois - ncorr_class_roi; contains the ROI(s) for the analysis
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.

    % Data ---------------------------------------------------------------%
    % Initialize outputs
    outstate = out.cancelled;
    rois = ncorr_class_roi.empty;
    % Get GUI handles
    handles_gui = init_gui();       
    % Run c-tor
    feval(ncorr_util_wrapcallbacktrycatch(@constructor,handles_gui.figure));
    
    % Callbacks and functions --------------------------------------------%
    function constructor()                
        % Initialize buffers
        if (isempty(params_init))
            rois_prelim(1:length(imgs)) = ncorr_class_roi;
            list_rois = false(1,length(imgs)); % Used to keep track of which rois are set
        else
            rois_prelim = params_init;
            list_rois = true(1,length(imgs)); % Used to keep track of which rois are set
        end
        
        % Set Data
        setappdata(handles_gui.figure,'num_img',length(imgs)-1);
        setappdata(handles_gui.figure,'rois_prelim',rois_prelim);
        setappdata(handles_gui.figure,'list_rois',list_rois);
        
        % Update
        update_axes('set');
        update_sidemenu();
        
        % Format window
        set(handles_gui.figure,'Visible','on'); 
    end

    function callback_edit_imgnum(hObject,eventdata) %#ok<INUSD>
        freeze_menu();   
        
        % Get data
        num_img = getappdata(handles_gui.figure,'num_img');

        % Get img num
        num_img_prelim = str2double(get(handles_gui.edit_imgnum,'string')); 
        if (ncorr_util_isintbb(num_img_prelim,1,length(imgs),'Image number') == out.success)  
            % num_img_prelim will be zero based indexed
            num_img = num_img_prelim-1;
        end    

        % Set data
        setappdata(handles_gui.figure,'num_img',num_img); 
        
        % Update
        update_axes('set');
        update_sidemenu();
        
        unfreeze_menu();
    end

    function callback_button_loadroi(hObject,evendata) %#ok<INUSD>
        freeze_menu();
        
        % Get Data
        num_img = getappdata(handles_gui.figure,'num_img');
        rois_prelim = getappdata(handles_gui.figure,'rois_prelim');
        list_rois = getappdata(handles_gui.figure,'list_rois');
        
        % Get ROI
        [mask_img_prelim,outstate_roi] = ncorr_util_loadimgs(false);
        
        if (outstate_roi == out.success)
            % Convert to binary mask  
            mask_load_prelim = im2bw(mask_img_prelim.get_gs());         %#ok<IM2BW>

            % Make sure ROI and image are the same size and that the ROI 
            % is nonempty:
            if (isequal(size(mask_load_prelim),[imgs(num_img+1).height imgs(num_img+1).width]) && ~isempty(find(mask_load_prelim,1)))                    
                % Form roi_load_prelim
                roi_load_prelim = ncorr_class_roi;
                roi_load_prelim.set_roi('load',struct('mask',mask_load_prelim,'cutoff',2000)); 
                if (roi_load_prelim.get_fullregions() > 0)
                    % At this point, roi_load_prelim fits the critera for a
                    % ROI. Display it to the user.    
                    outstate_roi = ncorr_gui_loadroi(imgs(num_img+1), ...
                                                     roi_load_prelim, ...
                                                     get(handles_gui.figure,'OuterPosition'));

                    if (outstate_roi == out.success)
                        % Set ROI
                        rois_prelim(num_img+1) = roi_load_prelim;
                        list_rois(num_img+1) = true;
                    end
                else
                    h_error = errordlg('ROI must contain a large contiguous region. Small regions are filtered out by default.','Error','modal');
                    uiwait(h_error);
                end
            else
                h_error = errordlg('ROI was invalid - either not the same size as the reference images or empty.','Error','modal');
                uiwait(h_error);
            end   
        end
        
        % Set Data                        
        setappdata(handles_gui.figure,'rois_prelim',rois_prelim);
        setappdata(handles_gui.figure,'list_rois',list_rois);
              
        % Update
        update_axes('set');     
        update_sidemenu();
        
        unfreeze_menu();
    end

    function callback_button_drawroi(hObject,evendata) %#ok<INUSD>
        freeze_menu();
        
        % Get Data
        num_img = getappdata(handles_gui.figure,'num_img');
        rois_prelim = getappdata(handles_gui.figure,'rois_prelim');
        list_rois = getappdata(handles_gui.figure,'list_rois');
        
        % Check if ROI was drawn before
        params_init_sub = []; % must use sub postfix to prevent overwriting of params_init input
        if (list_rois(num_img+1) && strcmp(rois_prelim(num_img+1).type,'draw'))
            params_init_sub = rois_prelim(num_img+1).data.drawobjects;
        end
        
        [roi_draw_prelim,outstate_roi] = ncorr_gui_drawroi(imgs(num_img+1), ...
                                                           get(handles_gui.figure,'OuterPosition'), ...
                                                           params_init_sub);  
        
        if (outstate_roi == out.success)
            % Set ROI
            rois_prelim(num_img+1) = roi_draw_prelim;
            list_rois(num_img+1) = true;    
        end        
        
        % Set Data
        setappdata(handles_gui.figure,'rois_prelim',rois_prelim);
        setappdata(handles_gui.figure,'list_rois',list_rois); 
        
        % Update
        update_axes('set');     
        update_sidemenu();
        
        unfreeze_menu();
    end

    function callback_button_finish(hObject,evendata) %#ok<INUSD>
        % Get Data        
        rois_prelim = getappdata(handles_gui.figure,'rois_prelim');
        
        % Set outputs
        for i = 0:length(imgs)-1
            rois(i+1) = rois_prelim(i+1);
        end
        outstate = out.success;

        % Exit        
        close(handles_gui.figure);
    end

    function callback_button_cancel(hObject,evendata) %#ok<INUSD>
        close(handles_gui.figure);
    end

    function callback_button_left(hObject,evendata) %#ok<INUSD>
        freeze_menu();   
        
        % Get Data
        num_img = getappdata(handles_gui.figure,'num_img');
        
        % Check for overshoot
        if (num_img > 0)
            num_img = num_img-1;
        end    
        
        % Set Data
        setappdata(handles_gui.figure,'num_img',num_img);  
        
        % Update
        update_axes('set');   
        update_sidemenu();
        
        unfreeze_menu();
    end

    function callback_button_right(hObject,evendata) %#ok<INUSD>
        freeze_menu();   
        
        % Get Data
        num_img = getappdata(handles_gui.figure,'num_img'); 
        
        % Check for overshoot
        if (num_img < length(imgs)-1)
            num_img = num_img+1;
        end
        
        % Set Data
        setappdata(handles_gui.figure,'num_img',num_img);    
        
        % Update
        update_axes('set');   
        update_sidemenu();
        
        unfreeze_menu();
    end

    function freeze_menu()
        set(handles_gui.button_loadroi,'Enable','off'); 
        set(handles_gui.button_drawroi,'Enable','off');
        set(handles_gui.button_cancel,'Enable','off');
        set(handles_gui.button_finish,'Enable','off');
    end

    function unfreeze_menu()
        set(handles_gui.button_loadroi,'Enable','on'); 
        set(handles_gui.button_drawroi,'Enable','on');
        set(handles_gui.button_cancel,'Enable','on'); 
    end

    function update_sidemenu()
        % Get data
        list_rois = getappdata(handles_gui.figure,'list_rois');
        
        % Enable finish button if all ROIs have been set
        if (all(list_rois))
            set(handles_gui.button_finish,'Enable','on');
        end
        
        % Update text
        set(handles_gui.text_menu,'String',[num2str(sum(list_rois)) ' of ' num2str(length(list_rois)) ' ROIs set.']);   
    end

    function update_axes(action) 
        % Get Data
        num_img = getappdata(handles_gui.figure,'num_img');
        rois_prelim = getappdata(handles_gui.figure,'rois_prelim');      
        list_rois = getappdata(handles_gui.figure,'list_rois');  

        if (strcmp(action,'set'))
            % Show preview 
            if (list_rois(num_img+1))
                % Paint ROI over images
                preview_roi = imgs(num_img+1).get_gs();
                preview_roi(rois_prelim(num_img+1).mask) = preview_roi(rois_prelim(num_img+1).mask)+imgs(num_img+1).max_gs;
                imshow(preview_roi,[imgs(num_img+1).min_gs 2*imgs(num_img+1).max_gs],'Parent',handles_gui.axes_preview);
                set(handles_gui.axes_preview,'Visible','off');
            else
                % Paint only images
                imshow(imgs(num_img+1).get_gs(),[imgs(num_img+1).min_gs 2*imgs(num_img+1).max_gs],'Parent',handles_gui.axes_preview);
                set(handles_gui.axes_preview,'Visible','off');
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
    end

    function handles_gui = init_gui()
    % GUI controls -------------------------------------------------------%
        % Figure
        handles_gui.figure = figure( ...
            'Tag', 'figure', ...
            'Units', 'characters', ...
            'Position', ncorr_util_figpos(pos_parent,[35 126.7]), ...
            'Name', 'Draw ROI', ...
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
        handles_gui.group_roioptions = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_roioptions', ...
            'Units', 'characters', ...
            'Position', [2 26.1 35 8.3], ...
            'Title', 'ROI Options', ...
            'Interruptible','off');
        
        handles_gui.group_menu = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_menu', ...
            'Units', 'characters', ...
            'Position', [2 18.6 35 6.7], ...
            'Title', 'Menu', ...
            'Interruptible','off');

        handles_gui.group_preview = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_preview', ...
            'Units', 'characters', ...
            'Position', [38.7 0.75 85.8 33.6], ...
            'Title', 'ROI', ...
            'Interruptible','off');

        % Axes
        handles_gui.axes_preview = axes( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'axes_preview', ...
            'Units', 'characters', ...
            'Position', [2.4 4 79.7 27.6]', ...
            'Visible','off', ...
            'Interruptible','off');

        % Static Texts
        handles_gui.text_menu = uicontrol( ...
            'Parent', handles_gui.group_roioptions, ...
            'Tag', 'text_menu', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.3 5.0 23.5 1.3], ...
            'String', ['0 of ' num2str(length(imgs)) ' ROIs set.'], ...
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
            
        % Editbox 
        handles_gui.edit_imgnum = uicontrol( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'edit_imgnum', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [68.1 1.3 7 1.6], ...
            'String', '', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_imgnum,handles_gui.figure), ...
            'Enable', 'off', ...
            'Interruptible','off');

        % Push Buttons
        handles_gui.button_loadroi = uicontrol( ...
            'Parent', handles_gui.group_roioptions, ...
            'Tag', 'loadroi', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 3 29.7 1.7], ...
            'String', 'Load ROI', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_loadroi,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_drawroi = uicontrol( ...
            'Parent', handles_gui.group_roioptions, ...
            'Tag', 'draw', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 1 29.7 1.7], ...
            'String', 'Draw ROI', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_drawroi,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_finish = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'finish', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 3 29.7 1.7], ...
            'String', 'Finish', ...
            'Enable', 'off', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_finish,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_cancel = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'draw', ...
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
            'Position', [61.1 1.2 6 1.8], ...
            'String', '<', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_left,handles_gui.figure), ...
            'Interruptible','off', ...
            'Enable', 'off');

        handles_gui.button_right = uicontrol( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'button_right', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [76.1 1.2 6 1.8], ...
            'String', '>', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_right,handles_gui.figure), ...
            'Interruptible','off', ...
            'Enable', 'off');
    end

    % Pause until figure is closed ---------------------------------------%
    waitfor(handles_gui.figure);
end

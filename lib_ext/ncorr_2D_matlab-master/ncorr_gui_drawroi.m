function [roi,outstate] = ncorr_gui_drawroi(img,pos_parent,params_init)
% This is a GUI for drawing the ROI. 
%
% Inputs -----------------------------------------------------------------%
%   img - ncorr_class_img; used for displaying the background image
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%   params_init - struct; contains struct('pos_imroi',{},'type',{},'addorsub',{}) 
%   which intializes the drawing if it was done previously.
% 
% Outputs ----------------------------------------------------------------%
%   roi - ncorr_class_roi; contains the drawn ROI
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success. 

    % Data ---------------------------------------------------------------%
    % Initialize Outputs
    outstate = out.cancelled;
    roi = ncorr_class_roi.empty;
    % Get GUI handles
    handles_gui = init_gui();    
    % Run c-tor
    feval(ncorr_util_wrapcallbacktrycatch(@constructor,handles_gui.figure));
    
    % Callbacks and functions --------------------------------------------%
    function constructor()  
        % Set zoom and pan
        handle_zoom = zoom(handles_gui.figure);  
        handle_pan = pan(handles_gui.figure);   
        
        % Link axes for zoom
        linkaxes([handles_gui.axes_roi handles_gui.axes_preview],'xy');
        
        % Callback for escape key can be omitted since it is handled
        % properly within the callback
        
        % Get gs values
        gs_buffer = img.get_gs();
        
        % Set data
        setappdata(handles_gui.figure,'mask_prelim',false(size(gs_buffer)));
        setappdata(handles_gui.figure,'gs_buffer',gs_buffer);                                  
        setappdata(handles_gui.figure,'handle_preview',[]);
        setappdata(handles_gui.figure,'drawobjects',struct('imroi',{},'type',{},'addorsub',{}));  
        setappdata(handles_gui.figure,'imroi_active',[]);
        setappdata(handles_gui.figure,'handle_zoom',handle_zoom);
        setappdata(handles_gui.figure,'handle_pan',handle_pan);
                
        % Need to set axes before doing initialization because image must
        % be set before placing imrois
        update_axes('set');
        
        % Initialize imrois if drawing has already been done
        if (~isempty(params_init))
            for i = 0:length(params_init)-1
                % Convert position to 1 based indexing
                pos_init = params_init(i+1).pos_imroi;
                if (strcmp(params_init(i+1).type,'rect') || strcmp(params_init(i+1).type,'ellipse'))
                    pos_init(1:2) = pos_init(1:2)+1;
                elseif (strcmp(params_init(i+1).type,'poly'))
                    pos_init = pos_init+1;
                end            
                
                % Set imroi
                if (strcmp(params_init(i+1).type,'rect'))
                    imroi_init = imrect(handles_gui.axes_roi,pos_init);
                elseif (strcmp(params_init(i+1).type,'ellipse'))
                    imroi_init = imellipse(handles_gui.axes_roi,pos_init);
                elseif (strcmp(params_init(i+1).type,'poly'))
                    imroi_init = impoly(handles_gui.axes_roi,pos_init);
                end      
                
                % Format 
                format_imroi(imroi_init);            
                
                % Append
                append_imroi(imroi_init,params_init(i+1).type,params_init(i+1).addorsub);
            end            
        end   
        
        % Update mask - it is updated in-place
        get_mask();
        
        % Need to update the axes
        update_axes('update');
        update_sidemenu();
                
        % Set Visible
        set(handles_gui.figure,'Visible','on'); 
    end
    
    function callback_button_imroi(hObject,eventdata,type,addorsub) %#ok<INUSL>
    % This is the callback for the add/sub imroi buttons.
        freeze_menu();                
        
        % Placing imroi --------------------------------------------------%
        if (strcmp(type,'rect'))
            imroi = imrect(handles_gui.axes_roi);
        elseif (strcmp(type,'ellipse'))
            imroi = imellipse(handles_gui.axes_roi);
        else
            imroi = impoly(handles_gui.axes_roi);
        end
        
        % imroi can be empty if escape key is pressed
        if (~isempty(imroi)) 
            % Format
            format_imroi(imroi);    

            % Append
            append_imroi(imroi,type,addorsub);
            
            % Update mask
            get_mask();

            % Set this as the active ROI           
            setappdata(handles_gui.figure,'imroi_active',imroi);
        end          
                
        % Update
        update_axes('update');
        update_sidemenu();  
        
        % Unfreeze
        unfreeze_menu();     
    end

    function callback_button_clear(hObject,eventdata) %#ok<INUSD>
        % Get data
        drawobjects = getappdata(handles_gui.figure,'drawobjects');  
        
        % Delete imrois
        for i = 0:length(drawobjects)-1
            delete(drawobjects(i+1).imroi);      
        end                
        drawobjects = [];
        
        % Set data
        setappdata(handles_gui.figure,'drawobjects',drawobjects);  
        
        % Update mask
        get_mask();
        
        % Update
        update_axes('update');
        update_sidemenu();  
    end

    function callback_button_zoom(hObject,eventdata) %#ok<INUSD>
        % Get data
        handle_zoom = getappdata(handles_gui.figure,'handle_zoom');  
        
        if (strcmp(get(handle_zoom,'Enable'),'on'))
            % Zoom is already enabled; disable it
            set(handle_zoom,'Enable','off');   
        else
            % Zoom not enabled; enable it
            set(handle_zoom,'Enable','on');   
        end        
            
        % Set data
        setappdata(handles_gui.figure,'handle_zoom',handle_zoom);  
        
        % Update
        update_sidemenu()
    end

    function callback_button_pan(hObject,eventdata) %#ok<INUSD>
        % Get data
        handle_pan = getappdata(handles_gui.figure,'handle_pan');  
        
        if (strcmp(get(handle_pan,'Enable'),'on'))
            % Pan is already enabled; disable it
            set(handle_pan,'Enable','off');   
        else
            % Pan not enabled; enable it
            set(handle_pan,'Enable','on');   
        end        
            
        % Set data
        setappdata(handles_gui.figure,'handle_pan',handle_pan);  
        
        % Update
        update_sidemenu()
    end

    function callback_button_finish(hObject,eventdata) %#ok<INUSD>        
        % Get data
        drawobjects = getappdata(handles_gui.figure,'drawobjects');     
        mask_prelim = getappdata(handles_gui.figure,'mask_prelim');
        
        % Format data - store this version, as it stores the pos_imroi instead
        % of the imrois directly, which get deleted after the figure is closed
        drawobjects_formatted = format_drawobjects(drawobjects);         
      
        % Form preliminary ROI
        roi_prelim = ncorr_class_roi;
        roi_prelim.set_roi('draw',struct('mask',mask_prelim,'drawobjects',drawobjects_formatted,'cutoff',2000)); 
        if (roi_prelim.get_fullregions() > 0)
            % Set outputs
            roi(1) = roi_prelim;
            outstate = out.success;
            
            % Exit            
            close(handles_gui.figure);  
        else
            h_error = errordlg('ROI must contain a large contiguous region. Please make ROI larger.','Error','modal');
            uiwait(h_error);
        end
    end

    function callback_button_cancel(hObject,eventdata) %#ok<INUSD>
        close(handles_gui.figure);       
    end

    function callback_axes_click(hObject,eventdata) %#ok<INUSD>
    % This is the callback for when the axes is clicked. It essentially
    % deactivates the active roi.
        setappdata(handles_gui.figure,'imroi_active',[]);
        
        % Update
        update_axes('update');
    end 

    function callback_delete(imroi)
    % This is the callback for when a ROI is deleted
        % Get data
        drawobjects = getappdata(handles_gui.figure,'drawobjects');    
        
        % Scan drawobjects
        for i = 0:length(drawobjects)-1
            if (drawobjects(i+1).imroi == imroi) 
                delete(drawobjects(i+1).imroi);
                drawobjects(i+1) = [];
                break;
            end
        end         
        
        % Set data  
        setappdata(handles_gui.figure,'drawobjects',drawobjects);         
        
        % Update mask
        get_mask();
        
        % Update
        update_axes('update');
        update_sidemenu();  
    end

    function callback_imroi(pos,drawobject) %#ok<INUSL>
    % This callback receives the current drawobject and updates it. pos is
    % ignored. After the data is set the axes are updated
        setappdata(handles_gui.figure,'imroi_active',drawobject);
        
        % Update mask
        get_mask();
        
        % Update
        update_axes('update');
    end

    function freeze_menu()        
        set(handles_gui.button_clear,'Enable','off'); 
        set(handles_gui.button_finish,'Enable','off'); 
        set(handles_gui.button_addrect,'Enable','off');
        set(handles_gui.button_subrect,'Enable','off');
        set(handles_gui.button_addellipse,'Enable','off');
        set(handles_gui.button_subellipse,'Enable','off');
        set(handles_gui.button_addpoly,'Enable','off');
        set(handles_gui.button_subpoly,'Enable','off'); 
        set(handles_gui.button_zoom,'Enable','off'); 
        set(handles_gui.button_pan,'Enable','off');         
    end

    function unfreeze_menu()
        set(handles_gui.button_clear,'Enable','on'); 
        set(handles_gui.button_addrect,'Enable','on');
        set(handles_gui.button_subrect,'Enable','on');
        set(handles_gui.button_addellipse,'Enable','on');
        set(handles_gui.button_subellipse,'Enable','on');
        set(handles_gui.button_addpoly,'Enable','on');
        set(handles_gui.button_subpoly,'Enable','on'); 
        set(handles_gui.button_zoom,'Enable','on'); 
        set(handles_gui.button_pan,'Enable','on'); 
    end

    function update_sidemenu()
        % Get data
        mask_prelim = getappdata(handles_gui.figure,'mask_prelim');
        handle_zoom = getappdata(handles_gui.figure,'handle_zoom'); 
        handle_pan = getappdata(handles_gui.figure,'handle_pan');          
        
        % Check if mask_prelim is empty
        if (~isempty(find(mask_prelim,1)))
            set(handles_gui.button_finish,'Enable','on')
        else
            set(handles_gui.button_finish,'Enable','off')
        end
        
        % Zoom/pan
        if (strcmp(get(handle_pan,'Enable'),'on'))
            set(handles_gui.button_pan,'FontWeight','bold');   
        else
            set(handles_gui.button_pan,'FontWeight','normal');  
        end
        
        if (strcmp(get(handle_zoom,'Enable'),'on'))
            set(handles_gui.button_zoom,'FontWeight','bold');   
        else
            set(handles_gui.button_zoom,'FontWeight','normal');  
        end        
    end

    function update_axes(action)
        % Get data   
        mask_prelim = getappdata(handles_gui.figure,'mask_prelim');        
        gs_buffer = getappdata(handles_gui.figure,'gs_buffer');
        handle_preview = getappdata(handles_gui.figure,'handle_preview');  
        drawobjects = getappdata(handles_gui.figure,'drawobjects'); 
        imroi_active = getappdata(handles_gui.figure,'imroi_active');

        if (strcmp(action,'set'))
            % Display image
            set(imshow(gs_buffer,[img.min_gs img.max_gs],'Parent',handles_gui.axes_roi),'ButtonDownFcn',ncorr_util_wrapcallbacktrycatch(@callback_axes_click,handles_gui.figure),'Interruptible','off');
            set(handles_gui.axes_roi,'Visible','off');     
            set(handles_gui.text_name,'String',['Name: ' img.name(1:min(end,22))]);  

            % Display preview
            handle_preview = imshow(gs_buffer,[img.min_gs img.max_gs*2],'Parent',handles_gui.axes_preview);
            set(handles_gui.axes_preview,'visible','off');
        end
        
        if (strcmp(action,'set') || strcmp(action,'update'))   
            % Scan drawobjects
            for i = 0:length(drawobjects)-1
                if (drawobjects(i+1).imroi == imroi_active) 
                    % This is the active imroi
                    setColor(drawobjects(i+1).imroi,'b');
                else
                    % Inactivate other imrois
                    if (strcmp(drawobjects(i+1).addorsub,'add'))
                        setColor(drawobjects(i+1).imroi,'g');
                    else
                        setColor(drawobjects(i+1).imroi,'r');
                    end
                end
            end  
        
            % Form preview_roi 
            preview_roi = gs_buffer;
            preview_roi(mask_prelim) = preview_roi(mask_prelim) + img.max_gs;
            
            % Set preview
            set(handle_preview,'CData',preview_roi);
        end
        
        % Set data
        setappdata(handles_gui.figure,'handle_preview',handle_preview);
        
        % Menu must be updated
        update_sidemenu();  
    end

    function format_imroi(imroi)
    % This function adds callbacks for position updates and the right-click 
    % menu
        setColor(imroi,'b');
        addNewPositionCallback(imroi,ncorr_util_wrapcallbacktrycatch(@(pos)callback_imroi(pos,imroi),handles_gui.figure));
        
        % Append callback to patch - this activates imroi when it is
        % clicked
        patch_handle = findobj(get(imroi,'Children'),'Type','patch');
        patch_bdc = get(patch_handle,'ButtonDownFcn');
        set(patch_handle,'ButtonDownFcn',ncorr_util_wrapcallbacktrycatch(@(h,e)(cellfun(@(x)feval(x,h,e),{patch_bdc,@(h,e)callback_imroi([],imroi)})),handles_gui.figure),'Interruptible','off')
        
        % Disable UIControl Menus
        set(get(imroi,'Children'),'UIContextMenu','');
        
        % Set new context menu for delete option
        menu_context = uicontextmenu('parent',handles_gui.figure);
        menu_context_item1 = uimenu(menu_context,'Label','Delete','Callback',ncorr_util_wrapcallbacktrycatch(@(h,e)callback_delete(imroi),handles_gui.figure),'Interruptible','off'); %#ok<NASGU>
        set(patch_handle,'UIContextMenu',menu_context);
    end

    function append_imroi(imroi,type,addorsub)
    % This function appends this imroi to the drawobjects list
        % Get data
        drawobjects = getappdata(handles_gui.figure,'drawobjects');
              
        % Add drawobject
        drawobject_template.imroi = imroi;
        drawobject_template.type = type;
        drawobject_template.addorsub = addorsub;
        drawobjects = horzcat(drawobjects,drawobject_template); 
        
        % Set data
        setappdata(handles_gui.figure,'drawobjects',drawobjects);
    end

    function drawobjects_formatted = format_drawobjects(drawobjects)
    % This function formats the draw objects for the mex function and for
    % storing the info within an ncorr_class_roi.
        drawobjects_formatted = struct('pos_imroi',{},'type',{},'addorsub',{});
        for i = 0:length(drawobjects)-1
            % Get position and convert to 0 based indexing
            pos = getPosition(drawobjects(i+1).imroi);
            if (strcmp(drawobjects(i+1).type,'rect') || strcmp(drawobjects(i+1).type,'ellipse'))
                pos(1:2) = pos(1:2)-1;
            elseif (strcmp(drawobjects(i+1).type,'poly'))
                pos = pos-1;
            end            
            drawobjects_template(i+1).pos_imroi = pos; %#ok<AGROW>
            drawobjects_template(i+1).type = drawobjects(i+1).type; %#ok<AGROW>
            drawobjects_template(i+1).addorsub = drawobjects(i+1).addorsub; %#ok<AGROW>
        end   
        
        % Store data
        for i = 0:length(drawobjects)-1
            drawobjects_formatted(i+1) = drawobjects_template(i+1); 
        end
    end

    function get_mask() 
        % Get data
        mask_prelim = getappdata(handles_gui.figure,'mask_prelim'); 
        drawobjects = getappdata(handles_gui.figure,'drawobjects'); 

        % Format data 
        drawobjects_formatted = format_drawobjects(drawobjects);

        % Get mask; modified in-place.
        ncorr_alg_formmask(drawobjects_formatted,mask_prelim);    
    end

    function handles_gui = init_gui()
    % GUI controls -------------------------------------------------------%
        % Figure
        handles_gui.figure = figure( ...
            'Tag', 'figure', ...
            'Units', 'characters', ...
            'Position', ncorr_util_figpos(pos_parent,[35 214.4]), ...
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
        handles_gui.group_drawoptions = uibuttongroup( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_drawoptioncss', ...
            'Units', 'characters', ...
            'Position', [2 20.8 35 13.6], ...
            'Title', 'Drawing Options', ...
            'Interruptible','off');

        handles_gui.group_zoompan = uibuttongroup( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_zoompan', ...
            'Units', 'characters', ...
            'Position', [2 15.7 35 4.5], ...
            'Title', 'Zoom/Pan', ...
            'Interruptible','off');
        
        handles_gui.group_menu = uibuttongroup( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_menu', ...
            'Units', 'characters', ...
            'Position', [2 8.5 35 6.7], ...
            'Title', 'Menu', ...
            'Interruptible','off');
        
        handles_gui.group_roi = uibuttongroup( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_roi', ...
            'Units', 'characters', ...
            'Position', [38.7 0.75 85.8 33.6], ...
            'Title', 'Add/Remove Region', ...
            'Interruptible','off');
        
        handles_gui.group_preview = uibuttongroup( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_preview', ...
            'Units', 'characters', ...
            'Position', [126.4 0.75 85.8 33.6], ...
            'Title', 'ROI Preview', ...
            'Interruptible','off');

        % Axes
        handles_gui.axes_roi = axes( ...
            'Parent', handles_gui.group_roi, ...
            'Tag', 'axes_roi', ...
            'Units', 'characters', ...
            'Position', [2.4 3 80 28.6], ...
            'Visible','off', ...
            'Interruptible','off');

        handles_gui.axes_preview = axes( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'axes_preview', ...
            'Units', 'characters', ...
            'Position', [2.4 3 80 28.6], ...
            'Visible','off', ...
            'Interruptible','off');

        % Static texts    
        handles_gui.text_name = uicontrol( ...
            'Parent', handles_gui.group_roi, ...
            'Tag', 'text_name', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [3 0.7 56.9 1.3], ...
            'String', 'Name: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
                
        % Pushbuttons
        handles_gui.button_addrect = uicontrol( ...
            'Parent', handles_gui.group_drawoptions, ...
            'Tag', 'button_addrect', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 9.0 14.4 2.7], ...
            'String', '+ Rect', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata)callback_button_imroi(hObject,eventdata,'rect','add'),handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_subrect = uicontrol( ...
            'Parent', handles_gui.group_drawoptions, ...
            'Tag', 'button_subrect', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [17.4 9.0 14.4 2.7], ...
            'String', '- Rect', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata)callback_button_imroi(hObject,eventdata,'rect','sub'),handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_addellipse = uicontrol( ...
            'Parent', handles_gui.group_drawoptions, ...
            'Tag', 'button_addellipse', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 6.0 14.4 2.7], ...
            'String', '+ Ellipse', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata)callback_button_imroi(hObject,eventdata,'ellipse','add'),handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_subellipse = uicontrol( ...
            'Parent', handles_gui.group_drawoptions, ...
            'Tag', 'button_subellipse', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [17.4 6.0 14.4 2.7], ...
            'String', '- Ellipse', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata)callback_button_imroi(hObject,eventdata,'ellipse','sub'),handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_addpoly = uicontrol( ...
            'Parent', handles_gui.group_drawoptions, ...
            'Tag', 'button_addpoly', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 3.0 14.4 2.7], ...
            'String', '+ Poly', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata)callback_button_imroi(hObject,eventdata,'poly','add'),handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_subpoly = uicontrol( ...
            'Parent', handles_gui.group_drawoptions, ...
            'Tag', 'button_subpoly', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [17.4 3.0 14.4 2.7], ...
            'String', '- Poly', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata)callback_button_imroi(hObject,eventdata,'poly','sub'),handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_clear = uicontrol( ...
            'Parent', handles_gui.group_drawoptions, ...
            'Tag', 'button_clear', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 1 29.7 1.7], ...
            'String', 'Clear', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_clear,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.button_zoom = uicontrol( ...
            'Parent', handles_gui.group_zoompan, ...
            'Tag', 'button_zoom', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 1 14.4 1.7], ...
            'String', 'Zoom', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_zoom,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_pan = uicontrol( ...
            'Parent', handles_gui.group_zoompan, ...
            'Tag', 'button_pan', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [17.4 1 14.4 1.7], ...
            'String', 'Pan', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_pan,handles_gui.figure), ...
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

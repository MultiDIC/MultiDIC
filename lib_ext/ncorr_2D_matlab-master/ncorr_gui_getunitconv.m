function [pixtounits,units,outstate] = ncorr_gui_getunitconv(pos_parent)
% This is a GUI for getting the conversion factor from pixels to units.
% This assumes the image that is uploaded has square pixels, and that 
% the units per pixel is the same for all the images loaded. The units and 
% number of units are specified, then a line is draw to the number of
% units. Typically the image loaded will be an image of a ruler or the
% undeformed sample with known dimensions.
%
% Inputs -----------------------------------------------------------------%
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%
% Outputs ----------------------------------------------------------------%
%   pixtounits - double; units/pixels conversion factor.
%   units - string; units used in analysis.
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.

    % Data ---------------------------------------------------------------%     
    % Initialize Outputs 
    outstate = out.cancelled;
    pixtounits = [];
    units = '';
    % Get GUI handles    
    handles_gui = init_gui();
    % Run c-tor
    feval(ncorr_util_wrapcallbacktrycatch(@constructor,handles_gui.figure));
    
    % Callbacks and functions --------------------------------------------%
    function constructor()        
        % Set zoom and pan
        handle_zoom = zoom(handles_gui.figure);  
        handle_pan = pan(handles_gui.figure);   
        
        % Callback for escape key can be omitted since it is handled
        % properly within the callback
                
        % Initialize buffers
        img = ncorr_class_img.empty;
        units_prelim = 'mm';
        num_units_prelim = 10;
        pixtounits_prelim = [];
        
        % Store data   
        setappdata(handles_gui.figure,'img',img);  
        setappdata(handles_gui.figure,'units_prelim',units_prelim);  
        setappdata(handles_gui.figure,'num_units_prelim',num_units_prelim);  
        setappdata(handles_gui.figure,'pixtounits_prelim',pixtounits_prelim);  
        setappdata(handles_gui.figure,'handle_line',[]);  
        setappdata(handles_gui.figure,'handle_zoom',handle_zoom);
        setappdata(handles_gui.figure,'handle_pan',handle_pan);
        
        % Update
        update_sidemenu();
        
        % Set visible        
        set(handles_gui.figure,'Visible','on');         
    end
        
    function callback_edit_units(hObject,eventdata) %#ok<INUSD>        
        % Get units
        units_buffer = get(handles_gui.edit_units,'string');        
        
        % Set Data
        setappdata(handles_gui.figure,'units_prelim',units_buffer);
        
        % Update
        update_sidemenu();
    end

    function callback_edit_numunits(hObject,eventdata) %#ok<INUSD>
        % Get data
        num_units_prelim = getappdata(handles_gui.figure,'num_units_prelim');    
        pixtounits_prelim = getappdata(handles_gui.figure,'pixtounits_prelim');   
        handle_line = getappdata(handles_gui.figure,'handle_line');         
            
        num_units_buffer = str2double(get(handles_gui.edit_numunits,'string'));
        if (ncorr_util_isrealbb(num_units_buffer,1e-6,1e6,'# of units') == out.success)
            % Store in buffer
            num_units_prelim = num_units_buffer;
            % If line isnt set yet, get_pixtounits returns empty
            pixtounits_prelim = get_pixtounits(num_units_prelim,handle_line);
        end  
                
        % Store data
        setappdata(handles_gui.figure,'num_units_prelim',num_units_prelim);
        setappdata(handles_gui.figure,'pixtounits_prelim',pixtounits_prelim);
        
        % Update
        update_sidemenu();
    end

    function callback_button_setline(hObject,eventdata) %#ok<INUSD>
        freeze_menu();
          
        % Get data
        num_units_prelim = getappdata(handles_gui.figure,'num_units_prelim');    
        pixtounits_prelim = getappdata(handles_gui.figure,'pixtounits_prelim');
        handle_line = getappdata(handles_gui.figure,'handle_line');   
        
        % Place line
        handle_line_buffer = imline(handles_gui.axes_preview);              
            
        % Line could be empty if esc is pressed
        if (~isempty(handle_line_buffer)) 
            % Delete menu
            set(handle_line_buffer,'UIContextMenu','');

            % Assign position callback
            addNewPositionCallback(handle_line_buffer,ncorr_util_wrapcallbacktrycatch(@(pos)callback_line(pos),handles_gui.figure));
            
            % Store
            handle_line = handle_line_buffer;
            
            % Get new pixtounits parameter
            pixtounits_prelim = get_pixtounits(num_units_prelim,handle_line); 
        end
        
        % Set data
        setappdata(handles_gui.figure,'handle_line',handle_line);   
        setappdata(handles_gui.figure,'pixtounits_prelim',pixtounits_prelim);
            
        % Update
        update_sidemenu();  
        
        % Unfreeze
        unfreeze_menu();
    end

    function callback_button_loadimg(hObject,eventdata) %#ok<INUSD>
        freeze_menu();
        
        % Get data
        handle_line = getappdata(handles_gui.figure,'handle_line');   
        
        % false means lazy loading is disabled
        [img_prelim,outstate_loadimg] = ncorr_util_loadimgs(false);
        if (outstate_loadimg == out.success && length(img_prelim) == 1)     
            % Image sucessfully chosen; delete previous line if one was set
            if (~isempty(handle_line))
                delete(handle_line);
                handle_line = [];
            end
            
            % Store image
            img = img_prelim;
            
            % Set data
            setappdata(handles_gui.figure,'img',img);
            setappdata(handles_gui.figure,'pixtounits_prelim',[]); % Clear this
            setappdata(handles_gui.figure,'handle_line',handle_line); % Cleared
        
            % Update axes
            update_axes('set');
        elseif (outstate_loadimg == out.success && length(ref_prelim) > 1)
            h_error = errordlg('Please select only one calibration image.','Error','modal');  
            uiwait(h_error);
        end
        
        % Unfreeze menu
        unfreeze_menu();
        
        % Update
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
        % Get Data
        units_prelim = getappdata(handles_gui.figure,'units_prelim');  
        pixtounits_prelim = getappdata(handles_gui.figure,'pixtounits_prelim');  
        
        % Prompt user if inputs are correct
        string_disp = cell(0);
        string_disp{1} = ['Units/pixel is: ' num2str(pixtounits_prelim) '.'];
        string_disp{2} = ['Units are: ' units_prelim '.'];
        string_disp{3} = 'Is this correct?';
        contbutton = questdlg(string_disp,'Continue Operation','Yes','No','No');
        if (strcmp(contbutton,'Yes'))            
            units = units_prelim;
            pixtounits = pixtounits_prelim;
            outstate = out.success;
            
            % Exit
            close(handles_gui.figure);
        end
    end

    function callback_button_cancel(hObject,eventdata) %#ok<INUSD>
        close(handles_gui.figure);
    end

    function callback_line(pos) %#ok<INUSD>
        % Get data
        num_units_prelim = getappdata(handles_gui.figure,'num_units_prelim');     
        handle_line = getappdata(handles_gui.figure,'handle_line');         
            
        % Get update pixtounits
        pixtounits_prelim = get_pixtounits(num_units_prelim,handle_line);
                
        % Store data
        setappdata(handles_gui.figure,'pixtounits_prelim',pixtounits_prelim);
        
        % Update
        update_sidemenu();
    end

    function freeze_menu(hObject,eventdata) %#ok<INUSD>
        set(handles_gui.button_loadimg,'Enable','off');  
        set(handles_gui.button_setline,'Enable','off');  
        set(handles_gui.button_finish,'Enable','off');  
    end

    function unfreeze_menu(hObject,eventdata) %#ok<INUSD>
        set(handles_gui.button_loadimg,'Enable','on');  
    end

    function update_sidemenu()
        % Get data
        img = getappdata(handles_gui.figure,'img');
        units_prelim = getappdata(handles_gui.figure,'units_prelim');  
        num_units_prelim = getappdata(handles_gui.figure,'num_units_prelim');  
        pixtounits_prelim = getappdata(handles_gui.figure,'pixtounits_prelim');  
        handle_line = getappdata(handles_gui.figure,'handle_line'); 
        handle_zoom = getappdata(handles_gui.figure,'handle_zoom'); 
        handle_pan = getappdata(handles_gui.figure,'handle_pan');  
        
        % Set edit boxes
        set(handles_gui.edit_units,'String',units_prelim);
        set(handles_gui.edit_numunits,'String',num2str(num_units_prelim));
        
        % Set texts        
        set(handles_gui.text_pixtounits,'String',['Units/pixel : ' num2str(pixtounits_prelim)]);
        
        % Set buttons        
        if (~isempty(img) && isempty(handle_line))
            % Image has been loaded and line has not been set; enable the setting of a line            
            set(handles_gui.button_setline,'Enable','on');  
        else
            % Disable the line set button      
            set(handles_gui.button_setline,'Enable','off');        
        end
        
        if (~isempty(handle_line))
            % Line has been set, enable finish button     
            set(handles_gui.button_finish,'Enable','on');    
        else
            set(handles_gui.button_finish,'Enable','off');    
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
        img = getappdata(handles_gui.figure,'img');        
        
        if (strcmp(action,'set'))
            imshow(img.get_gs(),[img.min_gs img.max_gs],'Parent',handles_gui.axes_preview);
            set(handles_gui.axes_preview,'Visible','off');
            set(handles_gui.text_name,'String',['Name: ' img.name(1:min(end,22))]);
        end
    end

    function pixtounits_buf = get_pixtounits(num_units,handle_line)        
        % This calculates the pixtounits parameter. A postfix of buf is
        % used because pixtounits is an output argument with global scope.
        if (~isempty(handle_line))
            % Get Distance of line; units are in pixels
            pos_line = getPosition(handle_line);        
            dist_line = sqrt((pos_line(1,1)-pos_line(2,1))^2 + (pos_line(1,2)-pos_line(2,2))^2);

            % Divide number units by length of the line in pixels
            pixtounits_buf = num_units/dist_line;
        else
            pixtounits_buf = [];
        end
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
        handles_gui.group_units = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_units', ...
            'Units', 'characters', ...
            'Position', [2 22.5 35 12.0], ...
            'Title', 'Units options', ...
            'Interruptible','off');        
        
        handles_gui.group_zoompan = uibuttongroup( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_zoompan', ...
            'Units', 'characters', ...
            'Position', [2 17.4 35 4.5], ...
            'Title', 'Zoom/Pan', ...
            'Interruptible','off');
        
        handles_gui.group_menu = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_menu', ...
            'Units', 'characters', ...
            'Position', [2 10.4 35 6.4], ...
            'Title', 'Menu', ...
            'Interruptible','off');
        
        handles_gui.group_preview = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_preview', ...
            'Units', 'characters', ...
            'Position', [39.0 0.75 92 33.75], ...
            'Title', 'Line Preview', ...
            'Interruptible','off');

        % Axes
        handles_gui.axes_preview = axes( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'axes_setseeds', ...
            'Units', 'characters', ...
            'Position', [2.4 3 86.2 28.8], ...
            'Visible', 'off', ...
            'Interruptible','off');
        
        % Static Texts
        handles_gui.text_units = uicontrol( ...
            'Parent', handles_gui.group_units, ...
            'Tag', 'text_units', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.2 8.4 18 1.5], ...
            'String', 'Units: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
        
        handles_gui.text_numunits = uicontrol( ...
            'Parent', handles_gui.group_units, ...
            'Tag', 'text_numunits', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.2 6.7 18 1.5], ...
            'String', '# of Units: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');    
        
        handles_gui.text_pixtounits = uicontrol( ...
            'Parent', handles_gui.group_units, ...
            'Tag', 'text_pixtounits', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.2 5 29.7 1.5], ...
            'String', 'Units/pixel : ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');    

        handles_gui.text_name = uicontrol( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'text_ref', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [3 0.7 56.9 1.5], ...
            'String', 'Name: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
        
        % Edit
        handles_gui.edit_units = uicontrol( ...
            'Parent', handles_gui.group_units, ...
            'Tag', 'edit_units', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [21.2 8.7 10.1 1.3], ...
            'String', '0', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_units,handles_gui.figure), ...
            'Interruptible','off');  
        
        handles_gui.edit_numunits = uicontrol( ...
            'Parent', handles_gui.group_units, ...
            'Tag', 'edit_numunits', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [21.2 6.9 10.1 1.3], ...
            'String', '0', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_numunits,handles_gui.figure), ...
            'Interruptible','off');
        
        % Pushbuttons
        handles_gui.button_setline = uicontrol( ...
            'Parent', handles_gui.group_units, ...
            'Tag', 'button_setline', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 3 29.7 1.7], ...
            'String', 'Set Line', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_setline,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.button_loadimg = uicontrol( ...
            'Parent', handles_gui.group_units, ...
            'Tag', 'button_loadimg', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 1 29.7 1.7], ...
            'String', 'Load Calibration Image', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_loadimg,handles_gui.figure), ...
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

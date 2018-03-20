function handles_gui = ncorr_gui_viewplots(reference,current,data_dic,type_plot,pos_parent,params_init)
% This is a GUI for viewing the displacement/strain plots.
%
% Inputs -----------------------------------------------------------------%
%   reference - ncorr_class_img; used for displaying the background image.
%   current - ncorr_class_img; used for displaying the background image. 
%   data_dic - struct; contains displacement and strain info. This is the
%   main data structure used in Ncorr.
%   type_plot - string; specifies whether plot is u, v, exx, exy, or eyy.
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%   params_init - cell; {1} = num_cur, {2} = scalebarlength, 
%   {3} = if handle_scalebar is on, {4} = if handle_axes is on, {5} if
%   Lagrangian or Eulerian, {6} if zoomed; {7} if panned. Used to 
%   coordinate with any existing open plots.
%
% Outputs ----------------------------------------------------------------%
%   handles_gui - handles; these are the GUI handles which allow Ncorr to
%   manage the closing and synchronization of the open plots
%
% Note that this function checks both roi_ref_formatted and
% roi_cur_formatted to make sure all ROIs are non-empty. This is a very
% unlikely scenario, but if it happens this function will display and error
% dialogue and return.

    % Data ---------------------------------------------------------------%
    % Get GUI handles
    handles_gui = init_gui();  
    % Run c-tor
    feval(ncorr_util_wrapcallbacktrycatch(@constructor,handles_gui.figure));

    % Callbacks and functions --------------------------------------------%
    function constructor()                           
        % Check to make sure ROIs are full
        roifull = true;
        if (strcmp(type_plot,'u') || strcmp(type_plot,'v'))
            field = 'displacements';
        else
            field = 'strains';
        end
        for i = 0:length(current)-1 
            if (data_dic.(field)(i+1).roi_ref_formatted.get_fullregions == 0 || ...
                data_dic.(field)(i+1).roi_cur_formatted.get_fullregions == 0)
                roifull = false;
                break;
            end
        end  
               
        if (~roifull)            
            % One of the ROIs is empty, return; CloseRequestFcn is disabled,
            % in gui_init(), so reenable it before closing.
            set(handles_gui.figure,'CloseRequestFcn','closereq'); 
            h_error = errordlg('Some data plots are empty, please rerun analysis.','Error','modal');
            uiwait(h_error);
            close(handles_gui.figure);
            return;
        end      
        
        % Set zoom and pan
        handle_zoom = zoom(handles_gui.figure);  
        handle_pan = pan(handles_gui.figure);   

        % Initialize buffers -----------------------------------------%
        if (isempty(params_init))   
            % Set num_cur - by default show last image in set
            num_cur = length(current)-1;

            % Set scalebar length - initialize to 3/4 the width of the
            % reference image. This uses real units.
            scalebarlength_prelim = (3/4)*reference.width*data_dic.dispinfo.pixtounits;
            if (floor(scalebarlength_prelim) > 3)
                % Use integer
                scalebarlength = floor(scalebarlength_prelim);
            else
                % Use decimal 
                scalebarlength = scalebarlength_prelim;
            end                

            % handle_scalebar checkbox
            val_checkbox_scalebar = true;

            % handle_axes checkbox
            val_checkbox_axes = true;

            % Lagrangian or Eulerian; 1 == Lagrangian
            val_popupmenu = 1;
        else
            % Set num_cur
            num_cur = params_init{1};

            % Units and handle_scalebar buffer
            scalebarlength = params_init{2};

            % handle_scalebar checkbox
            val_checkbox_scalebar = params_init{3};

            % handle_axes checkbox
            val_checkbox_axes = params_init{4};

            % Lagrangian or Eulerian
            val_popupmenu = params_init{5};

            % Zoomed 
            if (params_init{6})
                set(handle_zoom,'Enable','on');
            end

            % Panned
            if (params_init{7})
                set(handle_pan,'Enable','on');
            end
        end

        % Transparency buffer
        transparency_prelim = 0.75;

        % Contour lines buffer
        contourlines_prelim = 20;
        max_contourlines = 100;            

        % Set slider buffer
        slider_buffer = struct('lagrangian',{},'eulerian',{});
        for i = 0:length(current)-1                 
            % Set parameters 
            % Get the data vector corresponding to this plot
            data_ref = get_data(type_plot,'lagrangian',i);                
            slider_buffer(i+1).lagrangian(1) = 0.5;  % Initialize to half
            slider_buffer(i+1).lagrangian(2) = 0.5;  % Initialize to half
            slider_buffer(i+1).lagrangian(3:7) = get_sliderparams(data_ref);

            % Get the data vector corresponding to this plot
            data_cur = get_data(type_plot,'eulerian',i);      
            slider_buffer(i+1).eulerian(1) = 0.5;  % Initialize to half
            slider_buffer(i+1).eulerian(2) = 0.5;  % Initialize to half
            slider_buffer(i+1).eulerian(3:7) = get_sliderparams(data_cur);
        end  
        max_upperbound = 1e10;
        min_lowerbound = -1e10;
        max_scalebarlength = 1e10;

        % Set data
        setappdata(handles_gui.figure,'num_cur',num_cur);  
        setappdata(handles_gui.figure,'transparency_prelim',transparency_prelim);  
        setappdata(handles_gui.figure,'contourlines_prelim',contourlines_prelim); 
        setappdata(handles_gui.figure,'max_contourlines',max_contourlines); 
        setappdata(handles_gui.figure,'scalebarlength',scalebarlength); 
        setappdata(handles_gui.figure,'slider_buffer',slider_buffer); 
        setappdata(handles_gui.figure,'max_upperbound',max_upperbound); 
        setappdata(handles_gui.figure,'min_lowerbound',min_lowerbound); 
        setappdata(handles_gui.figure,'max_scalebarlength',max_scalebarlength); 
        setappdata(handles_gui.figure,'type_plot',type_plot); % Store this explicitly so type and other parameters can be deduced from 'friends'
        setappdata(handles_gui.figure,'friends',[]); % This is used for modifying other viewplots; it stores their figure handles
        setappdata(handles_gui.figure,'text_info',[]); % This is the text box that displays when the cursor is over the data plot
        setappdata(handles_gui.figure,'val_checkbox_contour',false); 
        setappdata(handles_gui.figure,'val_checkbox_scalebar',val_checkbox_scalebar); 
        setappdata(handles_gui.figure,'val_checkbox_axes',val_checkbox_axes); 
        setappdata(handles_gui.figure,'val_checkbox_minmaxmarkers',false); 
        setappdata(handles_gui.figure,'val_popupmenu',val_popupmenu); 
        setappdata(handles_gui.figure,'handle_preview',[]);
        setappdata(handles_gui.figure,'handle_scalebar',[]);
        setappdata(handles_gui.figure,'handle_axes',[]);
        setappdata(handles_gui.figure,'handle_colorbar',[]);
        setappdata(handles_gui.figure,'handle_point_max',[]);
        setappdata(handles_gui.figure,'handle_point_min',[]);
        setappdata(handles_gui.figure,'handle_zoom',handle_zoom);
        setappdata(handles_gui.figure,'handle_pan',handle_pan);

        % Update - must send the GUI handles; this is because
        % update_axes can be used to update other figure handles stored
        % in 'friends'.
        update_axes('set',handles_gui);
        update_sidemenu(handles_gui);

        % Set resize and hover function callbacks; This GUI is
        % resizable and displays a data cursor when the cursor hovers
        % over the data plot.
        % Set resize function
        set(handles_gui.figure,'ResizeFcn',ncorr_util_wrapcallbacktrycatch(@callback_resizefunction,handles_gui.figure)); 
        % Set hover function
        set(handles_gui.figure,'WindowButtonMotionFcn',ncorr_util_wrapcallbacktrycatch(@callback_moveplot,handles_gui.figure));   

        
        % Set Visible
        set(handles_gui.figure,'Visible','on'); 
    end

    function callback_topmenu_save(hObject,eventdata,includeinfo) %#ok<INUSL>
        % Get info
        num_cur = getappdata(handles_gui.figure,'num_cur');

        % Form save figure -----------------------------------------------%
        handles_gui_savefig = form_savefig(num_cur,includeinfo,[]);

        % Disable close function so figure isnt inadvertently closed.
        set(handles_gui_savefig.figure,'CloseRequestFcn','');

        % Set size(s) ----------------------------------------------------%
        % gui_savesize is a local function with a GUI; it allows the user
        % to modify the size of the image
        [size_savefig,timedelay,outstate] = gui_savesize(handles_gui_savefig, ...
                                                         false, ...
                                                         get(handles_gui.figure,'OuterPosition')); %#ok<*ASGLU>
        
        if (outstate == out.success)
            % Save image -------------------------------------------------%   
            [filename,pathname] = uiputfile({'*.jpg';'*.png';'*.bmp';'*.tif'},'Save Image');
            
            if (~isequal(filename,0) && ~isequal(pathname,0))
                overwrite = true;
                if (exist(fullfile(pathname,filename),'file'))
                    contbutton = questdlg('File already exists. Do you want to overwrite?','Continue Operation','Yes','No','No');
                    if strcmp(contbutton,'No')
                        overwrite = false;
                    end
                end   
                if (overwrite)                
                    % Get image ------------------------------------------%
                    img_printscreen = getframe(handles_gui_savefig.figure); 
            
                    % Save the image
                    imwrite(img_printscreen.cdata,fullfile(pathname,filename));
                end 
            end
        end

        % Exit -----------------------------------------------------------%        
        set(handles_gui_savefig.figure,'CloseRequestFcn','closereq');
        close(handles_gui_savefig.figure);       
    end

    function callback_topmenu_savegif(hObject,eventdata) %#ok<INUSD>
        % Get data
        val_popupmenu = getappdata(handles_gui.figure,'val_popupmenu');
        
        % All current images must have the same size to use this feature
        % with the eulerian description
        size_cur = size(current(1).get_gs());
        samesize_cur = true;
        for i = 1:length(current)-1
            if (~isequal(size_cur,[current(i+1).height current(i+1).width]))
                samesize_cur = false;
                break;
            end
        end
        
        % Save image -----------------------------------------------------%
        % Note val_popupmenu equals 1 for the lagrangian perspective, and 2
        % for the eulerian perspective
        if (val_popupmenu == 1 || (val_popupmenu == 2 && samesize_cur))            
            % Form initial save figure -----------------------------------%
            handles_gui_savefig = form_savefig(0,0,[]);

            % Disable close function so figure isnt inadvertently closed
            set(handles_gui_savefig.figure,'CloseRequestFcn','');

            % Set Size(s) ------------------------------------------------%
            % gui_savesize is a local function with a GUI; it allows the user
            % to modify the size of the image
            [size_savefig,timedelay,outstate] = gui_savesize(handles_gui_savefig, ...
                                                             true, ...
                                                             get(handles_gui.figure,'OuterPosition'));
                                                        
            if (outstate == out.success)
                [filename,pathname] = uiputfile({'*.gif'},'Save Image');
                
                if (~isequal(filename,0) && ~isequal(pathname,0))
                    overwrite = true;
                    if (exist(fullfile(pathname,filename),'file'))
                        contbutton = questdlg('File already exists. Do you want to overwrite?','Continue Operation','Yes','No','No');
                        if strcmp(contbutton,'No')
                            overwrite = false;
                        end
                    end   
                    
                    if (overwrite)  
                        % Save initial image -----------------------------%
                        img_printscreen = getframe(handles_gui_savefig.figure); 
                        img_printscreen = frame2im(img_printscreen);
                        [imind,cm] = rgb2ind(img_printscreen,256);
                        
                        % Save
                        imwrite(imind,cm,[pathname filename],'gif','Loopcount',inf);
                        
                        % Cycle over other images to save
                        for i = 1:length(current)-1
                            % Close figure
                            set(handles_gui_savefig.figure,'CloseRequestFcn','closereq');
                            close(handles_gui_savefig.figure);  
                            
                            % Form updated save figure -------------------%
                            handles_gui_savefig = form_savefig(i,0,size_savefig);     

                            % Disable close function for now
                            set(handles_gui_savefig.figure,'CloseRequestFcn','');

                            % Get image ----------------------------------%
                            img_printscreen = getframe(handles_gui_savefig.figure); 
                            img_printscreen = frame2im(img_printscreen);
                            [imind,cm] = rgb2ind(img_printscreen,256);
                            
                            % Append
                            if (i ~= length(current)-1)
                                imwrite(imind,cm,[pathname filename],'gif','WriteMode','append','DelayTime',timedelay);   
                            else
                                % This is the last image, make the time
                                % delay longer so there's a pause at the
                                % end.
                                imwrite(imind,cm,[pathname filename],'gif','WriteMode','append','DelayTime',timedelay+0.5);  
                            end
                        end
                    end 
                end    
            end   
            
            % Close last figure
            set(handles_gui_savefig.figure,'CloseRequestFcn','closereq');
            close(handles_gui_savefig.figure);   
        else
            h_error = errordlg('All current images must have the same size to save gif with the Eulerian description.','Error','modal');
            uiwait(h_error);
        end
    end

    %---------------------------------------------------------------------%
    % These functions can potentially modify other viewplots -------------%
    %---------------------------------------------------------------------%   
    
    function callback_popupmenu(hObject,eventdata) %#ok<INUSD>
        % Get data
        friends = getappdata(handles_gui.figure,'friends');
        val_popupmenu = get(handles_gui.popupmenu,'Value');
        
        % Update other plots in friend list
        for i = 0:length(friends)-1
            if (~strcmp(getappdata(friends{i+1}.figure,'type_plot'),type_plot))
                % Must check to make sure figure has not been closed during
                % callback since protection is only gauranteed for this
                % figure.
                try                
                    % Set data
                    setappdata(friends{i+1}.figure,'val_popupmenu',val_popupmenu);

                    % Update
                    update_axes('set',friends{i+1});
                    update_sidemenu(friends{i+1});
                catch err
                    if (ishandle(friends{i+1}.figure))
                        rethrow(err);
                    end
                end   
            end
        end        
        
        % Set data
        setappdata(handles_gui.figure,'val_popupmenu',val_popupmenu); 

        % Update this plot
        update_axes('set',handles_gui);
        update_sidemenu(handles_gui);          
    end

    function callback_checkbox_scalebar(hObject,eventdata) %#ok<INUSD>
        % Get data
        friends = getappdata(handles_gui.figure,'friends');
        val_checkbox_scalebar = get(handles_gui.checkbox_scalebar,'Value');
        
        % Update other plots in friend list
        for i = 0:length(friends)-1
            if (~strcmp(getappdata(friends{i+1}.figure,'type_plot'),type_plot))
                % Must check to make sure figure has not been closed during
                % callback since protection is only gauranteed for this
                % figure.
                try                
                    % Set data
                    setappdata(friends{i+1}.figure,'val_checkbox_scalebar',val_checkbox_scalebar);

                    % Update
                    update_axes('update',friends{i+1});
                    update_sidemenu(friends{i+1});
                catch err
                    if (ishandle(friends{i+1}.figure))
                        rethrow(err);
                    end
                end     
            end
        end
        
        % Set data
        setappdata(handles_gui.figure,'val_checkbox_scalebar',val_checkbox_scalebar); 
        
        % Update this plot
        update_axes('update',handles_gui);
        update_sidemenu(handles_gui);   
    end

    function callback_checkbox_axes(hObject,eventdata) %#ok<INUSD>
        % Get data
        friends = getappdata(handles_gui.figure,'friends');
        val_checkbox_axes = get(handles_gui.checkbox_axes,'Value');
        
        % Update other plots in friend list
        for i = 0:length(friends)-1         
            if (~strcmp(getappdata(friends{i+1}.figure,'type_plot'),type_plot))  
                % Must check to make sure figure has not been closed during
                % callback since protection is only gauranteed for this
                % figure.
                try
                    % Set data
                    setappdata(friends{i+1}.figure,'val_checkbox_axes',val_checkbox_axes);

                    % Update
                    update_axes('update',friends{i+1});
                    update_sidemenu(friends{i+1});
                catch err
                    if (ishandle(friends{i+1}.figure))
                        rethrow(err);
                    end
                end  
            end
        end
        
        % Set data
        setappdata(handles_gui.figure,'val_checkbox_axes',val_checkbox_axes);

        % Update
        update_axes('update',handles_gui);
        update_sidemenu(handles_gui);
    end

    function callback_edit_scalebarlength(hObject,eventdata) %#ok<INUSD>
        % Get data
        scalebarlength = getappdata(handles_gui.figure,'scalebarlength');
        max_scalebarlength = getappdata(handles_gui.figure,'max_scalebarlength');
        friends = getappdata(handles_gui.figure,'friends');
        
        % Get Value
        scalebarlength_buffer = str2double(get(handles_gui.edit_scalebarlength,'string'));
        if (ncorr_util_isrealbb(scalebarlength_buffer,0,max_scalebarlength,'Scalebar Length') == out.success)
            % Update buffer
            scalebarlength = scalebarlength_buffer;
            
            % Update other plots in friend list
            for i = 0:length(friends)-1        
                if (~strcmp(getappdata(friends{i+1}.figure,'type_plot'),type_plot))   
                    % Must check to make sure figure has not been closed during
                    % callback since protection is only gauranteed for this
                    % figure.
                    try
                        % Set data
                        setappdata(friends{i+1}.figure,'scalebarlength',scalebarlength);

                        % Update
                        update_axes('update',friends{i+1});    
                        update_sidemenu(friends{i+1});
                    catch err
                        if (ishandle(friends{i+1}.figure))
                            rethrow(err);
                        end
                    end  
                end
            end
        end
        
        % Set data
        setappdata(handles_gui.figure,'scalebarlength',scalebarlength);

        % Update
        update_axes('update',handles_gui);    
        update_sidemenu(handles_gui);
    end

    function callback_edit_imgnum(hObject,eventdata) %#ok<INUSD>
        % Get data
        num_cur = getappdata(handles_gui.figure,'num_cur');   
        friends = getappdata(handles_gui.figure,'friends'); 
        
        % Get Value - uses one based indexing
        num_cur_prelim = str2double(get(handles_gui.edit_imgnum,'string'));
        if (ncorr_util_isintbb(num_cur_prelim,1,length(current),'Current Image number') == out.success)  
            % Store value - convert back to zero based indexing
            num_cur = num_cur_prelim-1;
            
            % Update other plots in friend list
            for i = 0:length(friends)-1        
                if (~strcmp(getappdata(friends{i+1}.figure,'type_plot'),type_plot))   
                    % Must check to make sure figure has not been closed during
                    % callback since protection is only gauranteed for this
                    % figure.
                    try
                        % Set data
                        setappdata(friends{i+1}.figure,'num_cur',num_cur);  
                        
                        % Update
                        update_axes('set',friends{i+1});
                        update_sidemenu(friends{i+1});                           
                    catch err
                        if (ishandle(friends{i+1}.figure))
                            rethrow(err);
                        end
                    end  
                end
            end        
        end    
        
        % Set data
        setappdata(handles_gui.figure,'num_cur',num_cur);  
        
        % Update this plot
        update_axes('set',handles_gui);
        update_sidemenu(handles_gui);  
    end

    function callback_button_zoom(hObject,eventdata) %#ok<INUSD>
        % Get data
        friends = getappdata(handles_gui.figure,'friends');
        handle_zoom = getappdata(handles_gui.figure,'handle_zoom');  
        
        if (strcmp(get(handle_zoom,'Enable'),'on'))
            % Zoom is already enabled; disable it
            val_zoom = false;
            set(handle_zoom,'Enable','off');  
        else
            % Zoom not enabled; enable it
            val_zoom = true;
            set(handle_zoom,'Enable','on');  
        end  
                
        % Update other plots in friend list
        for i = 0:length(friends)-1           
            if (~strcmp(getappdata(friends{i+1}.figure,'type_plot'),type_plot))   
                % Must check to make sure figure has not been closed during
                % callback since protection is only gauranteed for this
                % figure.
                try
                    % Get zoom handle from other plot
                    handle_zoom_sub = getappdata(friends{i+1}.figure,'handle_zoom');  

                    if (val_zoom)
                        set(handle_zoom_sub,'Enable','on');   
                    else
                        set(handle_zoom_sub,'Enable','off');   
                    end        

                    % Set data
                    setappdata(friends{i+1}.figure,'handle_zoom',handle_zoom_sub);  

                    % Update
                    update_sidemenu(friends{i+1});
                catch err
                    if (ishandle(friends{i+1}.figure))
                        rethrow(err);
                    end
                end  
            end
        end 

        % Set data
        setappdata(handles_gui.figure,'handle_zoom',handle_zoom);  

        % Update
        update_sidemenu(handles_gui);       
    end

    function callback_button_pan(hObject,eventdata) %#ok<INUSD>
        % Get data
        friends = getappdata(handles_gui.figure,'friends');
        handle_pan = getappdata(handles_gui.figure,'handle_pan');  
        
        if (strcmp(get(handle_pan,'Enable'),'on'))
            % Pan is already enabled; disable it
            val_pan = false;
            set(handle_pan,'Enable','off');  
        else
            % Pan not enabled; enable it
            val_pan = true;
            set(handle_pan,'Enable','on');   
        end  
                
        % Update other plots in friend list
        for i = 0:length(friends)-1        
            if (~strcmp(getappdata(friends{i+1}.figure,'type_plot'),type_plot))      
                % Must check to make sure figure has not been closed during
                % callback since protection is only gauranteed for this
                % figure.
                try
                    % Get Pan handle from other plot
                    handle_pan_sub = getappdata(friends{i+1}.figure,'handle_pan');  

                    if (val_pan)
                        set(handle_pan_sub,'Enable','on');   
                    else
                        set(handle_pan_sub,'Enable','off');   
                    end        

                    % Set data
                    setappdata(friends{i+1}.figure,'handle_pan',handle_pan_sub);  

                    % Update
                    update_sidemenu(friends{i+1});
                catch err
                    if (ishandle(friends{i+1}.figure))
                        rethrow(err);
                    end
                end  
            end
        end     

        % Set data
        setappdata(handles_gui.figure,'handle_pan',handle_pan);  

        % Update
        update_sidemenu(handles_gui);   
    end

    function callback_button_left(hObject,eventdata) %#ok<INUSD>
        % Get data
        num_cur = getappdata(handles_gui.figure,'num_cur');  
        friends = getappdata(handles_gui.figure,'friends'); 
        
        % Check for overshoot
        if (num_cur > 0)
            % Update other plots in friend list
            num_cur = num_cur-1;
            for i = 0:length(friends)-1  
                if (~strcmp(getappdata(friends{i+1}.figure,'type_plot'),type_plot))       
                    % Must check to make sure figure has not been closed during
                    % callback since protection is only gauranteed for this
                    % figure.
                    try
                        % Set data
                        setappdata(friends{i+1}.figure,'num_cur',num_cur);  
                        
                        % Update
                        update_axes('set',friends{i+1});    
                        update_sidemenu(friends{i+1});                
                    catch err
                        if (ishandle(friends{i+1}.figure))
                            rethrow(err);
                        end
                    end  
                end                
            end
        end    
        
        % Set data
        setappdata(handles_gui.figure,'num_cur',num_cur);  

        % Update this plot
        update_axes('set',handles_gui);     
        update_sidemenu(handles_gui);    
    end

    function callback_button_right(hObject,eventdata) %#ok<INUSD>
        % Get data
        num_cur = getappdata(handles_gui.figure,'num_cur');   
        friends = getappdata(handles_gui.figure,'friends'); 
        
        % Check for overshoot
        if (num_cur < length(current)-1)
            % Update other plots in friend list
            num_cur = num_cur+1;
            for i = 0:length(friends)-1    
                if (~strcmp(getappdata(friends{i+1}.figure,'type_plot'),type_plot))            
                    % Must check to make sure figure has not been closed during
                    % callback since protection is only gauranteed for this
                    % figure.
                    try
                        % Set data
                        setappdata(friends{i+1}.figure,'num_cur',num_cur);  
                        
                        % Update
                        update_axes('set',friends{i+1});
                        update_sidemenu(friends{i+1});                           
                    catch err
                        if (ishandle(friends{i+1}.figure))
                            rethrow(err);
                        end
                    end  
                end
            end   
        end    
        
        % Set data
        setappdata(handles_gui.figure,'num_cur',num_cur);  

        % Update this plot
        update_axes('set',handles_gui);
        update_sidemenu(handles_gui);   
    end

    function update_sidemenu(handles_gui_sub)
        % Get data - MUST USE HANDLES_GUI_SUB!!!
        num_cur = getappdata(handles_gui_sub.figure,'num_cur'); 
        scalebarlength = getappdata(handles_gui_sub.figure,'scalebarlength');  
        transparency_prelim = getappdata(handles_gui_sub.figure,'transparency_prelim');
        contourlines_prelim = getappdata(handles_gui_sub.figure,'contourlines_prelim');
        slider_buffer = getappdata(handles_gui_sub.figure,'slider_buffer'); 
        val_checkbox_scalebar = getappdata(handles_gui_sub.figure,'val_checkbox_scalebar');
        val_checkbox_axes = getappdata(handles_gui_sub.figure,'val_checkbox_axes');
        val_checkbox_contour = getappdata(handles_gui_sub.figure,'val_checkbox_contour');
        val_popupmenu = getappdata(handles_gui_sub.figure,'val_popupmenu');
        handle_zoom = getappdata(handles_gui_sub.figure,'handle_zoom'); 
        handle_pan = getappdata(handles_gui_sub.figure,'handle_pan');   
        
        if (val_popupmenu == 1)
            lore = 'lagrangian';
        else
            lore = 'eulerian';
        end
                
        % Sliders:
        set(handles_gui_sub.slider_transparency,'value',transparency_prelim);
        set(handles_gui_sub.slider_upperbound,'value',min(1,max(0,slider_buffer(num_cur+1).(lore)(1))));
        set(handles_gui_sub.slider_lowerbound,'value',min(1,max(0,slider_buffer(num_cur+1).(lore)(2))));
        
        % Edit:
        set(handles_gui_sub.edit_transparency,'String',num2str(transparency_prelim,'%6.4f'));
        set(handles_gui_sub.edit_upperbound,'String',num2str((slider_buffer(num_cur+1).(lore)(3)-slider_buffer(num_cur+1).(lore)(4))*slider_buffer(num_cur+1).(lore)(1)+slider_buffer(num_cur+1).(lore)(4),'%6.4f'));
        set(handles_gui_sub.edit_lowerbound,'String',num2str((slider_buffer(num_cur+1).(lore)(5)-slider_buffer(num_cur+1).(lore)(4))*slider_buffer(num_cur+1).(lore)(2)+slider_buffer(num_cur+1).(lore)(4),'%6.4f'));
        set(handles_gui_sub.edit_contour,'String',num2str(contourlines_prelim));            
        set(handles_gui_sub.edit_scalebarlength,'String',num2str(scalebarlength,'%6.2f'));
        
        % Popupmenu:
        set(handles_gui_sub.popupmenu,'Value',val_popupmenu);
        
        % Enable scalebar
        if (val_checkbox_scalebar)
            set(handles_gui_sub.edit_scalebarlength,'Enable','on','backgroundcolor',[0.9412 0.9412 0.9412]);
            set(handles_gui_sub.checkbox_scalebar,'Value',true);
        else
            set(handles_gui_sub.edit_scalebarlength,'Enable','off','backgroundcolor',[0.9412 0.9412 0.9412]);
            set(handles_gui_sub.checkbox_scalebar,'Value',false);
        end
        
        % Enable axes
        if (val_checkbox_axes)
            set(handles_gui_sub.checkbox_axes,'Value',true);
        else
            set(handles_gui_sub.checkbox_axes,'Value',false);
        end
        
        % Enable contour
        if (val_checkbox_contour)
            % Disable trans
            set(handles_gui_sub.edit_transparency,'Enable','off');
            set(handles_gui_sub.slider_transparency,'Enable','off');
            % Enable contours
            set(handles_gui_sub.edit_contour,'Enable','on','backgroundcolor',[0.9412 0.9412 0.9412]);
        else
            % Enable trans
            set(handles_gui_sub.edit_transparency,'Enable','on','backgroundcolor',[0.9412 0.9412 0.9412]);
            set(handles_gui_sub.slider_transparency,'Enable','on');
            % Disable contours
            set(handles_gui_sub.edit_contour,'Enable','off');
        end
           
        if (strcmp(get(handle_pan,'Enable'),'on'))
            set(handles_gui_sub.button_pan,'FontWeight','bold');   
        else
            set(handles_gui_sub.button_pan,'FontWeight','normal');  
        end
        
        if (strcmp(get(handle_zoom,'Enable'),'on'))
            set(handles_gui_sub.button_zoom,'FontWeight','bold');   
        else
            set(handles_gui_sub.button_zoom,'FontWeight','normal');  
        end
    end

    function update_axes(action,handles_gui_sub)
        % Get data - MUST USE HANDLES_GUI_SUB AND TYPE_PLOT_SUB!!!
        num_cur = getappdata(handles_gui_sub.figure,'num_cur');
        transparency_prelim = getappdata(handles_gui_sub.figure,'transparency_prelim');  
        contourlines_prelim = getappdata(handles_gui_sub.figure,'contourlines_prelim');
        scalebarlength = getappdata(handles_gui_sub.figure,'scalebarlength');
        slider_buffer = getappdata(handles_gui_sub.figure,'slider_buffer');
        type_plot_sub = getappdata(handles_gui_sub.figure,'type_plot');
        text_info = getappdata(handles_gui_sub.figure,'text_info');        
        val_checkbox_contour = getappdata(handles_gui_sub.figure,'val_checkbox_contour');
        val_checkbox_scalebar = getappdata(handles_gui_sub.figure,'val_checkbox_scalebar');
        val_checkbox_axes = getappdata(handles_gui_sub.figure,'val_checkbox_axes');
        val_checkbox_minmaxmarkers = getappdata(handles_gui_sub.figure,'val_checkbox_minmaxmarkers');
        val_popupmenu = getappdata(handles_gui_sub.figure,'val_popupmenu');
        handle_preview = getappdata(handles_gui_sub.figure,'handle_preview'); 
        handle_scalebar = getappdata(handles_gui_sub.figure,'handle_scalebar');
        handle_axes = getappdata(handles_gui_sub.figure,'handle_axes');
        handle_colorbar = getappdata(handles_gui_sub.figure,'handle_colorbar');
        handle_point_max = getappdata(handles_gui_sub.figure,'handle_point_max');
        handle_point_min = getappdata(handles_gui_sub.figure,'handle_point_min');
                
        if (val_popupmenu == 1)
            lore = 'lagrangian';
            img_bg = reference;
        else
            lore = 'eulerian';
            img_bg = current(num_cur+1);
        end
        
        % Get data
        data = get_data(type_plot_sub,lore,num_cur);
        plot_data = get_dataplot(type_plot_sub,lore,num_cur);
        roi_data = get_roi(type_plot_sub,lore,num_cur);
                
        if (strcmp(action,'set') || strcmp(action,'save'))
            % Get reduced img
            img_reduced = img_bg.reduce(data_dic.dispinfo.spacing);
            
            % Set Background Image
            imshow(img_reduced.get_img(),[img_reduced.min_gs img_reduced.max_gs],'Parent',handles_gui_sub.axes_formatplot);
            hold(handles_gui_sub.axes_formatplot,'on');
                        
            % Overlay plot - See if it's a contour plot
            if (val_checkbox_contour)
                % Format plot
                [grid_x,grid_y] = meshgrid(1:size(plot_data,2),1:size(plot_data,1));
                alphamap_nan = ~roi_data.mask;
                plot_data(alphamap_nan) = NaN;
                
                % Form contour plot
                [data_contour,handle_preview] = contourf(grid_x,grid_y,plot_data,contourlines_prelim,'Parent',handles_gui_sub.axes_formatplot);
                
                % Overlay image in NaN regions
                % This only works properly when opengl is disabled. When opengl
                % is enabled, stacking order is determined by projected Z
                % value. This is a problem because images and contourf functions
                % have z = 0 value, so sometimes the image will not overlay the
                % white regions.
                handle_trans = imshow(img_reduced.get_img(),[img_reduced.min_gs img_reduced.max_gs],'Parent',handles_gui_sub.axes_formatplot);
                set(handle_trans,'AlphaData',imdilate(alphamap_nan,true(3))); % This dilation covers the border
            else
                % Regular Plot
                handle_preview = imshow(plot_data,[],'Parent',handles_gui_sub.axes_formatplot);
            end
            % Set axes grid off
            set(handles_gui_sub.axes_formatplot,'Visible','off');
            
            % Place markers in min/max location        
            datamax = max(data);
            datamin = min(data);
            
            % Max marker
            [y_max,x_max] = find(plot_data == datamax,1);
            handle_point_max = impoint(handles_gui_sub.axes_formatplot,x_max,y_max);
            setColor(handle_point_max,'g');
            set(handle_point_max,'UIContextMenu','');
            set(handle_point_max,'ButtonDownFcn','');
            
            % Min Marker
            [y_min,x_min] = find(plot_data == datamin,1);
            handle_point_min = impoint(handles_gui_sub.axes_formatplot,x_min,y_min);
            setColor(handle_point_min,'g');
            set(handle_point_min,'UIContextMenu','');
            set(handle_point_min,'ButtonDownFcn','');        
            
            % Set invisible if markers are disabled
            if (~val_checkbox_minmaxmarkers)
                set(handle_point_max,'Visible','off');
                set(handle_point_min,'Visible','off');
            end
            
            % Turn hold off
            hold(handles_gui_sub.axes_formatplot,'off');
            
            % Set left/right buttons
            if (~strcmp(action,'save'))                
                set(handles_gui_sub.edit_imgnum,'String',num2str(num_cur+1));
                if (length(current) == 1)
                    set(handles_gui_sub.button_right,'Enable','off');
                    set(handles_gui_sub.button_left,'Enable','off');
                    set(handles_gui_sub.edit_imgnum,'Enable','off');
                elseif (num_cur == 0)
                    set(handles_gui_sub.button_right,'Enable','on');
                    set(handles_gui_sub.button_left,'Enable','off');
                    set(handles_gui_sub.edit_imgnum,'Enable','on');
                elseif (num_cur == length(current)-1)
                    set(handles_gui_sub.button_right,'Enable','off');
                    set(handles_gui_sub.button_left,'Enable','on');
                    set(handles_gui_sub.edit_imgnum,'Enable','on');
                else
                    set(handles_gui_sub.button_right,'Enable','on');
                    set(handles_gui_sub.button_left,'Enable','on');  
                    set(handles_gui_sub.edit_imgnum,'Enable','on');                                                                      
                end
            end
             
            % Static Texts:
            if (strcmp(action,'set') || (strcmp(action,'save') && getappdata(handles_gui_sub.figure,'includeinfo')))   
                set(handles_gui_sub.text_ref_name,'String',['Reference Name: ' reference.name(1:min(end,40))]);
                set(handles_gui_sub.text_cur_name,'String',['Current Name: ' current(num_cur+1).name(1:min(end,40))]);
                set(handles_gui_sub.text_type,'String',['Analysis type: ' data_dic.dispinfo.type]);
                
                if (strcmp(type_plot_sub,'u') || strcmp(type_plot_sub,'v'))
                    string_dicparams = ['RG-DIC Radius: ' num2str(data_dic.dispinfo.radius) ' | Subset Spacing: ' num2str(data_dic.dispinfo.spacing)];
                else    
                    string_dicparams = ['RG-DIC Radius: ' num2str(data_dic.dispinfo.radius) ' | Strain Radius: ' num2str(data_dic.straininfo.radius) ' | Subset Spacing: ' num2str(data_dic.dispinfo.spacing)];
                end
                set(handles_gui_sub.text_dicparams,'String', string_dicparams);
                
                set(handles_gui_sub.text_itparams,'String', ['Diffnorm Cutoff: ' num2str(data_dic.dispinfo.cutoff_diffnorm) ' | Iteration Cutoff: ' num2str(data_dic.dispinfo.cutoff_iteration) ' | Threads: ' num2str(data_dic.dispinfo.total_threads)]);
                
                if (data_dic.dispinfo.stepanalysis.enabled)
                    if (strcmp(data_dic.dispinfo.stepanalysis.type,'seed'))
                        string_stepanalysis = 'Step Analysis: Enabled | Type: Seed Propagation ';
                    else
                        string_stepanalysis = ['Step Analysis: Enabled | Type: Leap Frog | Step: ' num2str(data_dic.dispinfo.stepanalysis.step)];
                    end                    
                else    
                    string_stepanalysis = 'Step Analysis: Disabled';
                end
                set(handles_gui_sub.text_stepanalysis,'String',string_stepanalysis);
                
                if (data_dic.dispinfo.subsettrunc)
                    string_subsettrunc = 'RG-DIC Subset Truncation: Enabled';    
                else    
                    string_subsettrunc = 'RG-DIC Subset Truncation: Disabled';   
                end
                if (strcmp(type_plot_sub,'exx') || strcmp(type_plot_sub,'exy') || strcmp(type_plot_sub,'eyy'))
                    if (data_dic.straininfo.subsettrunc)
                        string_subsettrunc = horzcat(string_subsettrunc,' | Strain Subset Truncation: Enabled'); 
                    else
                        string_subsettrunc = horzcat(string_subsettrunc,' | Strain Subset Truncation: Disabled'); 
                    end
                end
                set(handles_gui_sub.text_subsettrunc,'String',string_subsettrunc);
                
                string_imgcorr = 'Image Correspondences: ';
                for i = 0:length(data_dic.dispinfo.imgcorr)-1 % Might be large
                    string_imgcorr = horzcat(string_imgcorr,['[' num2str(data_dic.dispinfo.imgcorr(i+1).idx_ref) ' ' num2str(data_dic.dispinfo.imgcorr(i+1).idx_cur) '] ']); %#ok<AGROW>
                end
                set(handles_gui_sub.text_imgcorr,'String',string_imgcorr);
                
                set(handles_gui_sub.text_pixtounits,'String', ['Units/pixels: ' num2str(data_dic.dispinfo.pixtounits) ' ' data_dic.dispinfo.units '/pixels']);
                set(handles_gui_sub.text_cutoff_corrcoef,'String', ['Correlation Coefficient Cutoff: ' num2str(data_dic.dispinfo.cutoff_corrcoef(num_cur+1),'%6.4f')]);
                set(handles_gui_sub.text_lenscoef,'String',['Radial Lens Distortion Coefficient: ' num2str(data_dic.dispinfo.lenscoef)]);
                
                if (strcmp(type_plot_sub,'u') || strcmp(type_plot_sub,'v'))
                    set(handles_gui_sub.text_minmedianmax,'String',['Max: ' num2str(slider_buffer(num_cur+1).(lore)(6),'%6.4f') ' ' data_dic.dispinfo.units ' | Median: ' num2str(slider_buffer(num_cur+1).(lore)(4),'%6.4f') ' ' data_dic.dispinfo.units ' | Min: ' num2str(slider_buffer(num_cur+1).(lore)(7),'%6.4f') ' ' data_dic.dispinfo.units]);
                elseif (strcmp(type_plot_sub,'exx') || strcmp(type_plot_sub,'exy') || strcmp(type_plot_sub,'eyy'))
                    set(handles_gui_sub.text_minmedianmax,'String',['Max: ' num2str(slider_buffer(num_cur+1).(lore)(6),'%6.4f') ' | Median: ' num2str(slider_buffer(num_cur+1).(lore)(4),'%6.4f') ' | Min: ' num2str(slider_buffer(num_cur+1).(lore)(7),'%6.4f')]);
                end 
            end
                        
            % Set empty text info
            text_info = text(0,0,2,'','Parent',handles_gui_sub.axes_formatplot,'HorizontalAlignment','left','VerticalAlignment','top','BackgroundColor',[1 1 1],'EdgeColor',[0.7 0.7 0.7],'tag','text_info','FontSize',8);
            
            % Update panel text if its strain, since it also displays the
            % strain tensor name
            if (~strcmp(action,'save') && (strcmp(type_plot_sub,'exx') || strcmp(type_plot_sub,'exy') || strcmp(type_plot_sub,'eyy')))
                if (strcmp(type_plot_sub,'exx'))
                    group_title = 'Exx ';
                elseif (strcmp(type_plot_sub,'exy'))
                    group_title = 'Exy ';  
                elseif (strcmp(type_plot_sub,'eyy'))
                    group_title = 'Eyy ';
                end     
                if (val_popupmenu == 1)
                    group_title = [group_title 'Green-Lagrangian'];
                else
                    group_title = [group_title 'Eulerian-Almansi'];  
                end
                % Update title
                set(handles_gui_sub.group_formataxes,'Title',group_title)
            end

            % Set colormap
            ncorr_util_colormap(handles_gui_sub.figure); 
        end        
        
        if (strcmp(action,'set') || strcmp(action,'update') || strcmp(action,'save'))
            % Get limits
            cmax = (slider_buffer(num_cur+1).(lore)(3)-slider_buffer(num_cur+1).(lore)(4))*slider_buffer(num_cur+1).(lore)(1)+slider_buffer(num_cur+1).(lore)(4);
            cmin = (slider_buffer(num_cur+1).(lore)(5)-slider_buffer(num_cur+1).(lore)(4))*slider_buffer(num_cur+1).(lore)(2)+slider_buffer(num_cur+1).(lore)(4);
            if (abs(cmin-cmax)<1e-4)
                cmax = cmax + 1e-5; % add a small increment to ensure they are not the same number
                cmin = cmin - 1e-5;
            end
            
            % Update Plot
            caxis(handles_gui_sub.axes_formatplot,[cmin cmax]);
            
            % Update transparency for regular plot
            if (~val_checkbox_contour)
                transparency = transparency_prelim;
                alphamap_plot = roi_data.mask.*transparency;
                set(handle_preview,'AlphaData',alphamap_plot);
            end
            
            % Set colorbar
            handle_colorbar = colorbar('peer',handles_gui_sub.axes_formatplot);        
            set(handle_colorbar,'UIContextMenu','');
            set(get(handle_colorbar,'child'),'YData',[cmin cmax]);
            set(handle_colorbar,'YLim',[cmin cmax]);
            set(handle_colorbar,'Units','Pixels');
            
            % Set max/min markers
            if (val_checkbox_minmaxmarkers)
                set(handle_point_max,'Visible','on');
                set(handle_point_min,'Visible','on');
            else
                set(handle_point_max,'Visible','off');
                set(handle_point_min,'Visible','off');
            end
        
            % Check for scale bar
            if (val_checkbox_scalebar)
                % Need to also test if handle is valid, since it is not
                % explicitly cleared when switching img_num
                if (isempty(handle_scalebar) || ~ishandle(handle_scalebar))
                    % Scalebar isnt present, so create it
                    handle_scalebar = form_scalebar(handles_gui_sub);
                else                                        
                    % Scalebar is already present - just update it
                    % Get display image dimensions                    
                    height_img = floor(img_bg.height/(data_dic.dispinfo.spacing+1));
                    width_img = floor(img_bg.width/(data_dic.dispinfo.spacing+1)); 
                    
                    line_scalebar = findobj(handle_scalebar,'tag','line');
                    width_sb = scalebarlength/(data_dic.dispinfo.pixtounits*(data_dic.dispinfo.spacing+1)); % Convert to pixels
                    height_sb = 0.015*height_img; % Height is 1.5% of img height        
                    offset_sb_left = 0.05*width_img; % Offset is 5% from left, and 95% from top
                    offset_sb_top = 0.95*height_img;
                    set(line_scalebar,'XData',[offset_sb_left offset_sb_left+width_sb offset_sb_left+width_sb offset_sb_left]);
                    set(line_scalebar,'YData',[offset_sb_top offset_sb_top offset_sb_top+height_sb offset_sb_top+height_sb]);

                    % Update bg
                    bg1 = findobj(handle_scalebar,'tag','bg1');
                    bg2 = findobj(handle_scalebar,'tag','bg2');
                    width_bg = (width_sb/width_img+0.06)*width_img;
                    height_bg = 0.12*height_img; % Height is a 12% of height        
                    offset_bg_left = 0.02*width_img; % Offset is 2% from left, and 86% from top
                    offset_bg_top = 0.86*height_img;
                    set(bg1,'XData',[offset_bg_left offset_bg_left+width_bg offset_bg_left+width_bg offset_sb_left+width_sb offset_sb_left+width_sb offset_sb_left offset_sb_left offset_bg_left]);
                    set(bg1,'YData',[offset_bg_top offset_bg_top offset_bg_top+height_bg offset_sb_top+height_sb offset_sb_top offset_sb_top offset_sb_top+height_sb offset_bg_top+height_bg]);
                    set(bg2,'XData',[offset_sb_left offset_sb_left+width_sb offset_bg_left+width_bg offset_bg_left]);
                    set(bg2,'YData',[offset_sb_top+height_sb offset_sb_top+height_sb offset_bg_top+height_bg offset_bg_top+height_bg]);

                    % Update text
                    text_scalebar = findobj(handle_scalebar,'tag','text_scalebar');
                    pos_x_text = offset_sb_left+width_sb/2;
                    pos_y_text = .905*height_img;
                    set(text_scalebar,'String',[num2str(scalebarlength) ' ' data_dic.dispinfo.units],'Position',[pos_x_text pos_y_text 1]);
                end
            else
                % See if scalebar exists
                if (~isempty(handle_scalebar) && ishandle(handle_scalebar))
                    delete(handle_scalebar);
                end
                handle_scalebar = [];
                
                if (~strcmp(action,'save'))
                    % Disable the edit box for editing scalebar length
                    set(handles_gui_sub.edit_scalebarlength,'Enable','off');
                end
            end

            % Check for the axes
            if (val_checkbox_axes)
                if (isempty(handle_axes) || ~ishandle(handle_axes))
                    % Create plot axes
                    handle_axes = form_plotaxes(handles_gui_sub);
                end
            else
                if (~isempty(handle_axes) && ishandle(handle_axes))
                    delete(handle_axes);
                end
                handle_axes = [];
            end
        end             
        
        % Set data
        setappdata(handles_gui_sub.figure,'handle_point_max', handle_point_max);
        setappdata(handles_gui_sub.figure,'handle_point_min', handle_point_min);
        setappdata(handles_gui_sub.figure,'handle_preview',handle_preview);
        setappdata(handles_gui_sub.figure,'handle_scalebar',handle_scalebar);
        setappdata(handles_gui_sub.figure,'handle_axes',handle_axes);
        setappdata(handles_gui_sub.figure,'handle_colorbar',handle_colorbar);
        setappdata(handles_gui_sub.figure,'text_info',text_info);
    end

    function handle_scalebar = form_scalebar(handles_gui_sub)
    % This function creates the scalebar
        % Get data
        num_cur = getappdata(handles_gui_sub.figure,'num_cur');
        scalebarlength = getappdata(handles_gui_sub.figure,'scalebarlength');        
        val_popupmenu = getappdata(handles_gui_sub.figure,'val_popupmenu');
                
        if (val_popupmenu == 1)
            img_bg = reference;
        else
            img_bg = current(num_cur+1);
        end

        % Get display image dimensions (these are reduced)
        height_img = floor(img_bg.height/(data_dic.dispinfo.spacing+1));
        width_img = floor(img_bg.width/(data_dic.dispinfo.spacing+1)); 
                
        % Form hggroup
        handle_scalebar = hggroup('parent',handles_gui_sub.axes_formatplot,'tag','handle_scalebar'); 
        
        % Form dimensions of scalebar - most importantly, width must be
        % converted back to pixels (in reduced coordinates). It must also be 
        % consistent among plots and based on the scalebarlength parameter. 
        % Any other parameter, like the scalebar height or offsets, can 
        % depend on this specific image size for display purposes.
        width_sb = scalebarlength/(data_dic.dispinfo.pixtounits*(data_dic.dispinfo.spacing+1));
        height_sb = 0.015*height_img; % Height is a 1.5% of im height        
        offset_sb_left = 0.05*width_img; % Offset is 5% from left, and 95% from top
        offset_sb_top = 0.95*height_img;    

        % BG - Form BG with two patches with a hole for the scale bar. 
        % For some reason alpha is applied to all overlapping patches 
        width_bg = (width_sb/width_img+0.06)*width_img;
        height_bg = 0.12*height_img; % Height is a 12% of height        
        offset_bg_left = 0.02*width_img; % Offset is 2% from left, and 86% from top
        offset_bg_top = 0.86*height_img;        
        patch([offset_bg_left offset_bg_left+width_bg offset_bg_left+width_bg offset_sb_left+width_sb offset_sb_left+width_sb offset_sb_left offset_sb_left offset_bg_left], ...
              [offset_bg_top offset_bg_top offset_bg_top+height_bg offset_sb_top+height_sb offset_sb_top offset_sb_top offset_sb_top+height_sb offset_bg_top+height_bg], ...
              [1 1 1 1 1 1 1 1], ...
              'k','parent',handle_scalebar,'linestyle','none','FaceAlpha',0.5,'tag','bg1'); 
        patch([offset_sb_left offset_sb_left+width_sb offset_bg_left+width_bg offset_bg_left], ...
              [offset_sb_top+height_sb offset_sb_top+height_sb offset_bg_top+height_bg offset_bg_top+height_bg], ...
              [1 1 1 1], ...
              'k','parent',handle_scalebar,'linestyle','none','FaceAlpha',0.5,'tag','bg2'); 

        % Line - The actual scale bar
        patch([offset_sb_left offset_sb_left+width_sb offset_sb_left+width_sb offset_sb_left], ...
              [offset_sb_top offset_sb_top offset_sb_top+height_sb offset_sb_top+height_sb], ...
              [1 1 1 1], ...
              'w','parent',handle_scalebar,'FaceAlpha',1,'tag','line');       

        % Text
        pos_x_text = offset_sb_left+width_sb/2;
        pos_y_text = .905*height_img;
        pos_img = get(handles_gui_sub.axes_formatplot,'Position');
        size_text = max(0.5*min(0.5*pos_img(4),0.14*pos_img(3)),8);
        text(pos_x_text,pos_y_text,1,[num2str(scalebarlength) ' ' data_dic.dispinfo.units],'parent',handle_scalebar,'color','w', ...
             'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','FontSize',size_text,'tag','text_scalebar');           
    end

    function handle_axes = form_plotaxes(handles_gui_sub)
    % This function forms the axes seen at the top left of the plot  
        % Get data
        num_cur = getappdata(handles_gui_sub.figure,'num_cur');       
        val_popupmenu = getappdata(handles_gui_sub.figure,'val_popupmenu');
                
        if (val_popupmenu == 1)
            img_bg = reference;
        else
            img_bg = current(num_cur+1);
        end
        
        % Do this WRT the display image 
        height_img = floor(img_bg.height/(data_dic.dispinfo.spacing+1));
        width_img = floor(img_bg.width/(data_dic.dispinfo.spacing+1));  
        
        % Form hggroup
        handle_axes = hggroup('parent',handles_gui_sub.axes_formatplot,'tag','handle_axes');        
        length_axes = 0.15*max(width_img,height_img);         
        width_axes = 0.012*max(width_img,height_img);    
        arrowwidth_axes = 1.0*width_axes;
        width_bg = width_axes*0.4;

        % BG
        patch([0 width_axes+width_bg width_axes+width_bg width_axes+arrowwidth_axes+2*width_bg 0],[0 0 length_axes-width_bg length_axes-width_bg length_axes+width_axes+arrowwidth_axes+width_bg],[1 1 1 1 1], ...
              'w','parent',handle_axes,'linestyle','none','tag','bg1');  
        patch([0 length_axes+width_axes+arrowwidth_axes+width_bg length_axes-width_bg length_axes-width_bg 0],[0 0 width_axes+arrowwidth_axes+2*width_bg width_axes+width_bg width_axes+width_bg],[1 1 1 1 1], ...
              'w','parent',handle_axes,'linestyle','none','tag','bg2'); 

        % Lines
        patch([0 width_axes width_axes width_axes+arrowwidth_axes 0],[0 0 length_axes length_axes length_axes+width_axes+arrowwidth_axes],[1 1 1 1 1], ...
              'k','parent',handle_axes,'linestyle','none','tag','arrow1');    
        patch([0 length_axes+width_axes+arrowwidth_axes length_axes length_axes 0],[0 0 width_axes+arrowwidth_axes width_axes width_axes],[1 1 1 1 1], ...
              'k','parent',handle_axes,'linestyle','none','tag','arrow2');   

        % Text BG
        spacing_textbg = 2.5*width_axes;
        offset_text = 3*width_axes;
        pos_text = length_axes+width_axes+arrowwidth_axes+width_bg+offset_text;
        pos_img = get(handles_gui_sub.axes_formatplot,'Position');
        size_text = max(0.5*min(0.5*pos_img(4),0.14*pos_img(3)),8);
        patch([offset_text-spacing_textbg offset_text+spacing_textbg offset_text+spacing_textbg offset_text-spacing_textbg],[pos_text-spacing_textbg pos_text-spacing_textbg pos_text+spacing_textbg pos_text+spacing_textbg],[1 1 1 1], ...
              'k','parent',handle_axes,'linestyle','none','FaceAlpha',0.5,'tag','text_bg1');   
        patch([pos_text-spacing_textbg pos_text+spacing_textbg pos_text+spacing_textbg pos_text-spacing_textbg],[offset_text-spacing_textbg offset_text-spacing_textbg offset_text+spacing_textbg offset_text+spacing_textbg],[1 1 1 1], ...
              'k','parent',handle_axes,'linestyle','none','FaceAlpha',0.5,'tag','text_bg2'); 

        % Text
        text(offset_text,pos_text,1,'Y','parent',handle_axes,'color','w', ...
             'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','FontSize',size_text,'tag','text_x');
        text(pos_text,offset_text,1,'X','parent',handle_axes,'color','w', ...
             'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold','FontSize',size_text,'tag','text_y');
    end

    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    %---------------------------------------------------------------------%
    
    function callback_slider_transparency(hObject,eventdata) %#ok<INUSD>
        % Get data
        transparency_prelim = get(handles_gui.slider_transparency,'value');
        
        % Set data
        setappdata(handles_gui.figure,'transparency_prelim',transparency_prelim);
        
        % Update
        update_axes('update',handles_gui);
        update_sidemenu(handles_gui);
    end

    function callback_slider_upperbound(hObject,eventdata) %#ok<INUSD>
        % Get data
        num_cur = getappdata(handles_gui.figure,'num_cur');
        slider_buffer = getappdata(handles_gui.figure,'slider_buffer');
        val_popupmenu = getappdata(handles_gui.figure,'val_popupmenu');
        
        if (val_popupmenu == 1)
            lore = 'lagrangian';
        else
            lore = 'eulerian';
        end
        
        % Update data
        slider_buffer(num_cur+1).(lore)(1) = get(handles_gui.slider_upperbound,'value');
        
        % Set data
        setappdata(handles_gui.figure,'slider_buffer',slider_buffer);
        
        % Update
        update_axes('update',handles_gui);
        update_sidemenu(handles_gui);
    end

    function callback_slider_lowerbound(hObject,eventdata) %#ok<INUSD>
        % Get data
        num_cur = getappdata(handles_gui.figure,'num_cur');
        slider_buffer = getappdata(handles_gui.figure,'slider_buffer');
        val_popupmenu = getappdata(handles_gui.figure,'val_popupmenu');
        
        if (val_popupmenu == 1)
            lore = 'lagrangian';
        else
            lore = 'eulerian';
        end
        
        % Update data
        slider_buffer(num_cur+1).(lore)(2) = get(handles_gui.slider_lowerbound,'value');
        
        % Set data
        setappdata(handles_gui.figure,'slider_buffer',slider_buffer);
        
        % Update
        update_axes('update',handles_gui);
        update_sidemenu(handles_gui);
    end
    
    function callback_checkbox_contour(hObject,eventdata) %#ok<INUSD>
        % Get data
        val_checkbox_contour = get(handles_gui.checkbox_contour,'Value');
        
        % Set data
        setappdata(handles_gui.figure,'val_checkbox_contour',val_checkbox_contour);
        
        % Update
        update_axes('set',handles_gui);  
        update_sidemenu(handles_gui);
    end
    
    function callback_checkbox_minmaxmarkers(hObject,eventdata) %#ok<INUSD>
        % Get data
        val_checkbox_minmaxmarkers = get(handles_gui.checkbox_minmaxmarkers,'Value');
        
        % Set data
        setappdata(handles_gui.figure,'val_checkbox_minmaxmarkers',val_checkbox_minmaxmarkers);
        
        % Update
        update_axes('update',handles_gui);
    end

    function callback_edit_transparency(hObject,eventdata) %#ok<INUSD>    
        % Get data
        transparency_prelim = getappdata(handles_gui.figure,'transparency_prelim');
        
        % Get Value
        transparency_buffer = str2double(get(handles_gui.edit_transparency,'string'));
        if (ncorr_util_isrealbb(transparency_buffer,0,1,'Transparency') == out.success)
            transparency_prelim = transparency_buffer;
        end
        
        % Set data
        setappdata(handles_gui.figure,'transparency_prelim',transparency_prelim);
        
        % Update
        update_axes('update',handles_gui);
        update_sidemenu(handles_gui);
    end
    
    function callback_edit_lowerbound(hObject,eventdata) %#ok<INUSD>
        % Get data
        num_cur = getappdata(handles_gui.figure,'num_cur');
        slider_buffer = getappdata(handles_gui.figure,'slider_buffer');
        min_lowerbound = getappdata(handles_gui.figure,'min_lowerbound');
        val_popupmenu = getappdata(handles_gui.figure,'val_popupmenu');
        
        if (val_popupmenu == 1)
            lore = 'lagrangian';
        else
            lore = 'eulerian';
        end
        
        % Get Value
        lowerbound_buffer = str2double(get(handles_gui.edit_lowerbound,'string'));
        if (ncorr_util_isrealbb(lowerbound_buffer,min_lowerbound,slider_buffer(num_cur+1).(lore)(4),'Lowerbound') == out.success)
            % Make sure denominator is not close to zero
            if (abs(slider_buffer(num_cur+1).(lore)(5)-slider_buffer(num_cur+1).(lore)(4)) <= 1e-10)
                slider_buffer(num_cur+1).(lore)(2) = 1;
            else
                slider_buffer(num_cur+1).(lore)(2) = (lowerbound_buffer-slider_buffer(num_cur+1).(lore)(4))/(slider_buffer(num_cur+1).(lore)(5)-slider_buffer(num_cur+1).(lore)(4));
            end
        end
            
        % Set data
        setappdata(handles_gui.figure,'slider_buffer',slider_buffer);
        
        % Update
        update_axes('update',handles_gui);
        update_sidemenu(handles_gui);
    end

    function callback_edit_upperbound(hObject,eventdata) %#ok<INUSD>
        % Get data
        num_cur = getappdata(handles_gui.figure,'num_cur');
        slider_buffer = getappdata(handles_gui.figure,'slider_buffer');
        max_upperbound = getappdata(handles_gui.figure,'max_upperbound');
        val_popupmenu = getappdata(handles_gui.figure,'val_popupmenu');
        
        if (val_popupmenu == 1)
            lore = 'lagrangian';
        else
            lore = 'eulerian';
        end
        
        % Get Value
        upperbound_buffer = str2double(get(handles_gui.edit_upperbound,'string'));
        if (ncorr_util_isrealbb(upperbound_buffer,slider_buffer(num_cur+1).(lore)(4),max_upperbound,'Upperbound') == out.success)
            % Make sure denominator is not close to zero
            if (abs(slider_buffer(num_cur+1).(lore)(3)-slider_buffer(num_cur+1).(lore)(4)) <= 1e-10)
                slider_buffer(num_cur+1).(lore)(1) = 1;
            else
                slider_buffer(num_cur+1).(lore)(1) = (upperbound_buffer-slider_buffer(num_cur+1).(lore)(4))/(slider_buffer(num_cur+1).(lore)(3)-slider_buffer(num_cur+1).(lore)(4));
            end                
        end
            
        % Set data
        setappdata(handles_gui.figure,'slider_buffer',slider_buffer);
        
        % Update
        update_axes('update',handles_gui);
        update_sidemenu(handles_gui);
    end
    
    function callback_edit_contour(hObject,eventdata) %#ok<INUSD>
        % Get data
        contourlines_prelim = getappdata(handles_gui.figure,'contourlines_prelim');
        max_contourlines = getappdata(handles_gui.figure,'max_contourlines');
        
        % Get Value
        contourlines_buffer = str2double(get(handles_gui.edit_contour,'string'));
        if (ncorr_util_isintbb(contourlines_buffer,1,max_contourlines,'Number of contour lines') == out.success)            
            contourlines_prelim = contourlines_buffer;
        end
        
        % Set data
        setappdata(handles_gui.figure,'contourlines_prelim',contourlines_prelim);
            
        % Update
        update_axes('set',handles_gui);
        update_sidemenu(handles_gui);
    end

    function callback_button_applytoall(hObject,eventdata) %#ok<INUSD>
        % Get data
        num_cur = getappdata(handles_gui.figure,'num_cur');
        slider_buffer = getappdata(handles_gui.figure,'slider_buffer');
        val_popupmenu = getappdata(handles_gui.figure,'val_popupmenu');
        
        if (val_popupmenu == 1)
            lore = 'lagrangian';
        else
            lore = 'eulerian';
        end
        
        % Get bounds
        upperbound_buffer = (slider_buffer(num_cur+1).(lore)(3)-slider_buffer(num_cur+1).(lore)(4))*slider_buffer(num_cur+1).(lore)(1)+slider_buffer(num_cur+1).(lore)(4);
        lowerbound_buffer = (slider_buffer(num_cur+1).(lore)(5)-slider_buffer(num_cur+1).(lore)(4))*slider_buffer(num_cur+1).(lore)(2)+slider_buffer(num_cur+1).(lore)(4);
        
        % Find cutoff for bounds.
        cutoff_buffer = vertcat(slider_buffer.(lore));
        cutoff_lowerbound = min(cutoff_buffer(:,4));
        cutoff_upperbound = max(cutoff_buffer(:,4));                        
        if (upperbound_buffer > cutoff_upperbound && lowerbound_buffer < cutoff_lowerbound)
            % Set bounds
            for i = 0:length(current)-1 
                % Set lowerbound
                if (abs(slider_buffer(i+1).(lore)(5)-slider_buffer(i+1).(lore)(4)) <= 1e-10)
                    slider_buffer(i+1).(lore)(2) = 1;
                else
                    slider_buffer(i+1).(lore)(2) = (lowerbound_buffer-slider_buffer(i+1).(lore)(4))/(slider_buffer(i+1).(lore)(5)-slider_buffer(i+1).(lore)(4));
                end
                
                % Set upperbound
                if (abs(slider_buffer(i+1).(lore)(3)-slider_buffer(i+1).(lore)(4)) <= 1e-10)
                    slider_buffer(i+1).(lore)(1) = 1;
                else
                    slider_buffer(i+1).(lore)(1) = (upperbound_buffer-slider_buffer(i+1).(lore)(4))/(slider_buffer(i+1).(lore)(3)-slider_buffer(i+1).(lore)(4));
                end  
            end              
        else
            h_error = errordlg(['Lowerbound must be a lower than ' num2str(cutoff_lowerbound) ' and upperbound must be greater than ' num2str(cutoff_upperbound) '.'],'Error','modal');
            uiwait(h_error);
        end
        
        % Set data
        setappdata(handles_gui.figure,'slider_buffer',slider_buffer);
        
        % Update
        update_axes('update',handles_gui);
        update_sidemenu(handles_gui);
    end

    function callback_button_resetdefaults(hObject,eventdata) %#ok<INUSD>
        % Get data
        slider_buffer = getappdata(handles_gui.figure,'slider_buffer');
        
        % Set slider buffer
        for i = 0:length(current)-1 
            % Set parameters
            data_ref = get_data(type_plot,'lagrangian',i);     
            slider_buffer(i+1).lagrangian(1) = 0.5;  % Initialize to half
            slider_buffer(i+1).lagrangian(2) = 0.5;  % Initialize to half
            slider_buffer(i+1).lagrangian(3:7) = get_sliderparams(data_ref);

            data_cur = get_data(type_plot,'eulerian',i);      
            slider_buffer(i+1).eulerian(1) = 0.5;  % Initialize to half
            slider_buffer(i+1).eulerian(2) = 0.5;  % Initialize to half
            slider_buffer(i+1).eulerian(3:7) = get_sliderparams(data_cur);
        end  
        
        % Set data
        setappdata(handles_gui.figure,'slider_buffer',slider_buffer);
        
        % Update
        update_axes('update',handles_gui);
        update_sidemenu(handles_gui);
    end

    function callback_resizefunction(hObject,eventdata) %#ok<INUSD>
    % This function constrains how the figure resizes its axes, buttons, etc
    % on resize. Be careful here because resize function will interrupt
    % callbacks even if 'interruptible' is set to 'off'
        % Get data
        handle_scalebar = getappdata(handles_gui.figure,'handle_scalebar');
        handle_axes = getappdata(handles_gui.figure,'handle_axes');
        pos_fig = get(handles_gui.figure,'Position');       
        val_checkbox_scalebar = getappdata(handles_gui.figure,'val_checkbox_scalebar');
        val_checkbox_axes = getappdata(handles_gui.figure,'val_checkbox_axes');
        
        % Set offsets
        offset_group_menu_y = 21.1;
        offset_group_view_y = 26.5;
        offset_group_scalebar_y = 33.8;
        offset_group_axes_y = 38.7;
        offset_group_zoompan_y = 43.7;
        offset_axes_x = 1.35;  
        offset_axes_y = 41.1;
        offset_img_x = 24.4;   
        offset_img_y = 80.8;
        offset_left_y = 66;
        offset_right_y = 51;  
        offset_edit_y = 59;  
        
        % Update menu panel
        pos_group_menu = get(handles_gui.group_menu,'Position');
        pos_group_menu(2) = pos_fig(4)-offset_group_menu_y;
        set(handles_gui.group_menu,'Position',pos_group_menu);
        
        % Update lore panel
        pos_group_view = get(handles_gui.group_viewoptions,'Position');
        pos_group_view(2) = pos_fig(4)-offset_group_view_y;
        set(handles_gui.group_viewoptions,'Position',pos_group_view);
        
        % Update scalebar panel
        pos_group_scalebar = get(handles_gui.group_scalebar,'Position');
        pos_group_scalebar(2) = pos_fig(4)-offset_group_scalebar_y;
        set(handles_gui.group_scalebar,'Position',pos_group_scalebar);
        
        % Update axes panel
        pos_group_axes = get(handles_gui.group_axes,'Position');
        pos_group_axes(2) = pos_fig(4)-offset_group_axes_y;
        set(handles_gui.group_axes,'Position',pos_group_axes); 
        
        % Update zoom/pan panel
        pos_group_zoompan = get(handles_gui.group_zoompan,'Position');
        pos_group_zoompan(2) = pos_fig(4)-offset_group_zoompan_y;
        set(handles_gui.group_zoompan,'Position',pos_group_zoompan); 

        % Update axes panel - make sure it has positive width and height
        if (pos_fig(3)-offset_axes_y > 0 && ...
            pos_fig(4)-offset_axes_x > 0)
            pos_axes = get(handles_gui.group_formataxes,'Position');
            pos_axes(3) = pos_fig(3)-offset_axes_y;
            pos_axes(4) = pos_fig(4)-offset_axes_x;
            set(handles_gui.group_formataxes,'Position',pos_axes);
        end
        
        % Update axes - make sure it has positive width and height
        if (pos_fig(3)-offset_img_y > 0 && ...
            pos_fig(4)-offset_img_x > 0)
            pos_img = get(handles_gui.axes_formatplot,'Position');
            pos_img(3) = pos_fig(3)-offset_img_y;
            pos_img(4) = pos_fig(4)-offset_img_x;
            set(handles_gui.axes_formatplot,'Position',pos_img);
        end
        
        % Update buttons
        if (pos_fig(3)-offset_left_y > 0)
            pos_left = get(handles_gui.button_left,'Position');
            pos_left(1) = pos_fig(3)-offset_left_y;
            set(handles_gui.button_left,'Position',pos_left);
        end            
        
        if (pos_fig(3)-offset_right_y > 0)
            pos_right = get(handles_gui.button_right,'Position');
            pos_right(1) = pos_fig(3)-offset_right_y;
            set(handles_gui.button_right,'Position',pos_right);
        end
        
        % Update imgnum edit        
        if (pos_fig(3)-offset_edit_y > 0)
            pos_right = get(handles_gui.edit_imgnum,'Position');
            pos_right(1) = pos_fig(3)-offset_edit_y;
            set(handles_gui.edit_imgnum,'Position',pos_right);
        end
        
        % Update handle_scalebar text if it it's enabled.
        % Only text is updated because text size is constant during resize;
        % patches will resize automatically with handle_axes. 
        % Note that resize function will interrupt callbacks, so check
        % explicitly that the scalebar exists first before resizing it (i.e.
        % its possible to interrupt the checkbox callback for the scalebar 
        % midway).
        if (val_checkbox_scalebar && ...
            ~isempty(handle_scalebar) && ishandle(handle_scalebar))
            text_scalebar = findobj(handle_scalebar,'tag','text_scalebar');
            pos_img = get(handles_gui.axes_formatplot,'Position');
            size_text = max(0.5*min(0.5*pos_img(4),0.14*pos_img(3)),8);
            set(text_scalebar,'FontSize',size_text);
        end
        
        % Update handle_axes text if it it's enabled.
        % Only text is updated because the size is constant during resize;
        % patches will resize automatically with handle_axes. 
        % Note that resize function will interrupt callbacks, so check
        % explicitly that the axes exists first before resizing it  (i.e.
        % its possible to interrupt the checkbox callback for the axes 
        % midway).
        if (val_checkbox_axes && ...
            ~isempty(handle_axes) && ishandle(handle_axes))
            text_x = findobj(handle_axes,'tag','text_x');
            text_y = findobj(handle_axes,'tag','text_y');
            pos_img = get(handles_gui.axes_formatplot,'Position');
            size_text = max(0.5*min(0.5*pos_img(4),0.14*pos_img(3)),8);
            set(text_x,'FontSize',size_text);
            set(text_y,'FontSize',size_text);
        end
    end 

    function callback_moveplot(hObject,eventdata) %#ok<INUSD>
    % This function tells the axes what to do if the mouse cursor hovers
    % over the plot axes
        % Get data
        num_cur = getappdata(handles_gui.figure,'num_cur');
        text_info = getappdata(handles_gui.figure,'text_info');
        handle_colorbar = getappdata(handles_gui.figure,'handle_colorbar'); 
        val_popupmenu = getappdata(handles_gui.figure,'val_popupmenu');
        handle_zoom = getappdata(handles_gui.figure,'handle_zoom'); 
        handle_pan = getappdata(handles_gui.figure,'handle_pan'); 
                
        if (val_popupmenu == 1)
            lore = 'lagrangian';
        else
            lore = 'eulerian';
        end
        
        % Get data
        plot_data = get_dataplot(type_plot,lore,num_cur);
        roi_data = get_roi(type_plot,lore,num_cur);
        
        % Test to see if mouse position is within ROI - uses 1 based
        % indexing. Also round it first.
        pos_mouse_axes = round(get(handles_gui.axes_formatplot,'CurrentPoint')); 
        if (pos_mouse_axes(1) >= 1 && pos_mouse_axes(3) >= 1 && pos_mouse_axes(1) <= size(roi_data.mask,2) && pos_mouse_axes(3) <= size(roi_data.mask,1))
            if (roi_data.mask(pos_mouse_axes(3),pos_mouse_axes(1)))
                % Get Plot info - note that strain is dimensionless
                if (strcmp(type_plot,'u'))
                    text_display = 'U-disp: ';
                    units = data_dic.dispinfo.units;
                elseif (strcmp(type_plot,'v'))
                    text_display = 'V-disp: ';
                    units = data_dic.dispinfo.units;
                elseif (strcmp(type_plot,'exx'))
                    text_display = 'Exx-strain: ';
                    units = '';
                elseif (strcmp(type_plot,'exy'))
                    text_display = 'Exy-strain: ';
                    units = '';
                elseif (strcmp(type_plot,'eyy'))
                    text_display = 'Eyy-strain: ';
                    units = '';
                end

                % Display text
                set(text_info,'String',{['x-pos: ' num2str(pos_mouse_axes(1))],['y-pos: ' num2str(pos_mouse_axes(3))],[text_display num2str(plot_data(pos_mouse_axes(3),pos_mouse_axes(1))) ' ' units]}, ...
                    'Position',[pos_mouse_axes(1) pos_mouse_axes(3) 2],'HorizontalAlignment','left','VerticalAlignment','top');

                % If zoom and pan arent being used, then use a cross for the
                % cursor
                if (strcmp(get(handle_zoom,'Enable'),'off') && strcmp(get(handle_pan,'Enable'),'off'))
                    set(handles_gui.figure,'Pointer','cross');
                end

                % Check if text overlaps colorbar or bottom of the axes; if so,
                % flip it so the text shows properly
                % Set all coordinates to pixels
                set(text_info,'units','pixels');
                set(handle_colorbar,'units','pixels');
                set(handles_gui.group_formataxes,'units','pixels');
                set(handles_gui.figure,'units','pixels');

                % Get positions
                pos_text = get(text_info,'extent');
                pos_colorbar = get(handle_colorbar,'Position');
                pos_group = get(handles_gui.group_formataxes,'position');
                pos_cursor = get(handles_gui.figure,'CurrentPoint');
                if (pos_cursor(1)+pos_text(3) > pos_group(1)+pos_colorbar(1) && pos_text(2) < 0)
                    set(text_info,'HorizontalAlignment','right','VerticalAlignment','bottom');
                elseif (pos_cursor(1)+pos_text(3) > pos_group(1)+pos_colorbar(1))
                    set(text_info,'HorizontalAlignment','right');
                elseif (pos_text(2) < 0) % Hangs below bottom
                    set(text_info,'VerticalAlignment','bottom');
                end

                % Convert back to original coordinates
                set(text_info,'units','data');
                set(handle_colorbar,'units','normalized');
                set(handles_gui.group_formataxes,'units','characters');
                set(handles_gui.figure,'units','characters');
            else
                % Clear text and set cursor back to arrow if zoom and pan
                % aren't enabled
                set(text_info,'String','');
                if (strcmp(get(handle_zoom,'Enable'),'off') && strcmp(get(handle_pan,'Enable'),'off'))
                    set(handles_gui.figure,'Pointer','arrow');
                end
            end
        else
            % Clear text since cursor is not on top of axes
            set(text_info,'String','');
            set(handles_gui.figure,'Pointer','arrow');
        end
    end

    function [data] = get_data(type_plot_sub,lore,num_cur)
    % This function returns displacement/strain data in vector from the data_dic structure    
        if (strcmp(lore,'lagrangian'))
            rorc = 'ref';
        else
            rorc = 'cur';
        end
    
        if (strcmp(type_plot_sub,'u'))
            % U:
            data = data_dic.displacements(num_cur+1).(['plot_u_' rorc '_formatted'])(data_dic.displacements(num_cur+1).(['roi_' rorc '_formatted']).mask);
        elseif (strcmp(type_plot_sub,'v'))
            % V:
            data = data_dic.displacements(num_cur+1).(['plot_v_' rorc '_formatted'])(data_dic.displacements(num_cur+1).(['roi_' rorc '_formatted']).mask);
        elseif (strcmp(type_plot_sub,'exx'))
            % Exx: 
            data = data_dic.strains(num_cur+1).(['plot_exx_' rorc '_formatted'])(data_dic.strains(num_cur+1).(['roi_' rorc '_formatted']).mask);
        elseif (strcmp(type_plot_sub,'exy'))
            % Exy: 
            data = data_dic.strains(num_cur+1).(['plot_exy_' rorc '_formatted'])(data_dic.strains(num_cur+1).(['roi_' rorc '_formatted']).mask);
        elseif (strcmp(type_plot_sub,'eyy'))
            % Eyy:
            data = data_dic.strains(num_cur+1).(['plot_eyy_' rorc '_formatted'])(data_dic.strains(num_cur+1).(['roi_' rorc '_formatted']).mask);
        end
    end

    function [plot_data] = get_dataplot(type_plot_sub,lore,num_cur)
    % This function returns plot_data from the data_dic structure    
        if (strcmp(lore,'lagrangian'))
            rorc = 'ref';
        else
            rorc = 'cur';
        end
    
        if (strcmp(type_plot_sub,'u'))
            % U:
            plot_data = data_dic.displacements(num_cur+1).(['plot_u_' rorc '_formatted']);
        elseif (strcmp(type_plot_sub,'v'))
            % V:
            plot_data = data_dic.displacements(num_cur+1).(['plot_v_' rorc '_formatted']);
        elseif (strcmp(type_plot_sub,'exx'))
            % Exx: 
            plot_data = data_dic.strains(num_cur+1).(['plot_exx_' rorc '_formatted']);
        elseif (strcmp(type_plot_sub,'exy'))
            % Exy: 
            plot_data = data_dic.strains(num_cur+1).(['plot_exy_' rorc '_formatted']);
        elseif (strcmp(type_plot_sub,'eyy'))
            % Eyy:
            plot_data = data_dic.strains(num_cur+1).(['plot_eyy_' rorc '_formatted']);
        end
    end

    function [roi_data] = get_roi(type_plot_sub,lore,num_cur)
    % This function returns the ROI from the data_dic structure      
        if (strcmp(lore,'lagrangian'))
            rorc = 'ref';
        else
            rorc = 'cur';
        end
        
        if (strcmp(type_plot_sub,'u') || strcmp(type_plot_sub,'v'))
            % Displacement
            roi_data = data_dic.displacements(num_cur+1).(['roi_' rorc '_formatted']);
        else
            % Strains
            roi_data = data_dic.strains(num_cur+1).(['roi_' rorc '_formatted']);
        end
    end  

    function sliderparams = get_sliderparams(data)
    % Returns slider params based on data input
        param_upperbound = prctile(data,99);
        param_median = prctile(data,50);
        param_lowerbound = prctile(data,1);
        
        sliderparams = [param_upperbound+(param_upperbound-param_median) ...
                        param_median ...
                        param_lowerbound-(param_median-param_lowerbound) ...
                        max(data) ...
                        min(data)];
    end

    function handles_gui_savefig = form_savefig(num_cur,includeinfo,params_init)
        % Get data - num_cur is specified
        transparency_prelim = getappdata(handles_gui.figure,'transparency_prelim');  
        contourlines_prelim = getappdata(handles_gui.figure,'contourlines_prelim');
        scalebarlength = getappdata(handles_gui.figure,'scalebarlength');
        slider_buffer = getappdata(handles_gui.figure,'slider_buffer'); 
        val_checkbox_contour = getappdata(handles_gui.figure,'val_checkbox_contour');
        val_checkbox_scalebar = getappdata(handles_gui.figure,'val_checkbox_scalebar');
        val_checkbox_axes = getappdata(handles_gui.figure,'val_checkbox_axes');
        val_checkbox_minmaxmarkers = getappdata(handles_gui.figure,'val_checkbox_minmaxmarkers');
        val_popupmenu = getappdata(handles_gui.figure,'val_popupmenu');
        handle_point_max = getappdata(handles_gui.figure,'handle_point_max');
        handle_point_min = getappdata(handles_gui.figure,'handle_point_min');
                             
        if (val_popupmenu == 1)
            img_bg = reference;
        else
            img_bg = current(num_cur+1);
        end
        
        % Form figure
        handles_gui_savefig.figure = figure( ...
            'Units','pixels', ...
            'Name', 'Save Preview', ...
            'MenuBar', 'none', ...
            'NumberTitle', 'off', ...
            'Color','white', ...
            'Resize','off', ...
            'handlevisibility','off', ...
            'DockControls','off', ...
            'Visible','off', ...
            'IntegerHandle','off', ...
            'Interruptible','off', ...
            'WindowStyle','modal');

        % Create axes
        handles_gui_savefig.axes_formatplot = axes('Parent',handles_gui_savefig.figure,'Units','pixels');  

        % If the user wants to include info
        if (includeinfo)
            % Set Info texts
            uicontrol( ...
                'Parent', handles_gui_savefig.figure, ...
                'Style', 'text', ...
                'Units', 'characters', ...
                'Position', [2.1 18.5 100 1.3], ...
                'String', ['Type: ' type_plot '-plot'], ...
                'HorizontalAlignment', 'left', ...
                'BackgroundColor','white', ...
                'Interruptible','off');

            % Now get info text which is below the plot
            handles_infotext = set_infotext(handles_gui_savefig.figure,'white');
            fields_infotext = fieldnames(handles_infotext);
            for i = 0:length(fields_infotext)-1
                handles_gui_savefig.(fields_infotext{i+1}) = handles_infotext.(fields_infotext{i+1});
            end

            % Set data height
            height_data = 255;
        else
            height_data = 0;
        end   

        % Transfer data, except for some handles, so that update_axes can set
        % the plot
        setappdata(handles_gui_savefig.figure,'num_cur',num_cur);
        setappdata(handles_gui_savefig.figure,'transparency_prelim',transparency_prelim);  
        setappdata(handles_gui_savefig.figure,'contourlines_prelim',contourlines_prelim);
        setappdata(handles_gui_savefig.figure,'scalebarlength',scalebarlength);
        setappdata(handles_gui_savefig.figure,'slider_buffer',slider_buffer); 
        setappdata(handles_gui_savefig.figure,'type_plot',type_plot); 
        setappdata(handles_gui_savefig.figure,'val_checkbox_contour',val_checkbox_contour);
        setappdata(handles_gui_savefig.figure,'val_checkbox_scalebar',val_checkbox_scalebar);
        setappdata(handles_gui_savefig.figure,'val_checkbox_axes',val_checkbox_axes);
        setappdata(handles_gui_savefig.figure,'val_checkbox_minmaxmarkers',val_checkbox_minmaxmarkers);
        setappdata(handles_gui_savefig.figure,'val_popupmenu',val_popupmenu);
        setappdata(handles_gui_savefig.figure,'handle_preview',[]); 
        setappdata(handles_gui_savefig.figure,'handle_scalebar',[]);
        setappdata(handles_gui_savefig.figure,'handle_axes',[]);
        setappdata(handles_gui_savefig.figure,'handle_point_max',handle_point_max);
        setappdata(handles_gui_savefig.figure,'handle_point_min',handle_point_min);   
        % Additional fields              
        setappdata(handles_gui_savefig.figure,'includeinfo',includeinfo);

        % Set plot
        update_axes('save',handles_gui_savefig);
        
        % Set default size
        if (isempty(params_init))
            width_img = floor(img_bg.width/(data_dic.dispinfo.spacing+1));
            height_img = floor(img_bg.height/(data_dic.dispinfo.spacing+1));
        else
            width_img = params_init(1);
            height_img = params_init(2);
        end
        border = 30;
        width_cb = 20;
        spacing_save_fig = 100;
        spacing_width_img = 4*border+width_cb;
        spacing_height_img = 2*border+height_data;
        
        % Set data in savefig - these are parameters which are modified by
        % gui_savesize
        setappdata(handles_gui_savefig.figure,'height_data',height_data); 
        setappdata(handles_gui_savefig.figure,'width_img',width_img);
        setappdata(handles_gui_savefig.figure,'height_img',height_img);
        setappdata(handles_gui_savefig.figure,'border',border);
        setappdata(handles_gui_savefig.figure,'width_cb',width_cb);
        setappdata(handles_gui_savefig.figure,'spacing_save_fig',spacing_save_fig);
        setappdata(handles_gui_savefig.figure,'spacing_width_img',spacing_width_img);
        setappdata(handles_gui_savefig.figure,'spacing_height_img',spacing_height_img);
        
        % Update savefig
        update_savefig(handles_gui_savefig);
        
        % Set plot visible 
        set(handles_gui_savefig.figure,'Visible','on');
    end
    
    function handles_gui = init_gui()        
        handles_gui.figure = figure( ...
            'Tag', 'figure', ...
            'Units', 'characters', ...
            'Position', ncorr_util_figpos(pos_parent,[44.5 150]), ...
            'Name', 'Data Plot', ...
            'MenuBar', 'none', ...
            'NumberTitle', 'off', ...
            'Color', get(0,'DefaultUicontrolBackgroundColor'), ... 
            'handlevisibility','off', ...
            'DockControls','off', ...
            'IntegerHandle','off', ...
            'Interruptible','off', ...
            'Visible','off', ...
            'CloseRequestFcn','');
        
        % Top Menu        
        handles_gui.topmenu_file = uimenu( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'topmenu_file', ...
            'Label', 'File', ...
            'Checked', 'off', ...
            'Interruptible','off');

        handles_gui.topmenu_save = uimenu( ...
            'Parent', handles_gui.topmenu_file, ...
            'Tag', 'topmenu_save', ...
            'Label', 'Save Image', ...
            'Checked', 'off', ...
            'Interruptible','off');        
        
        handles_gui.topmenu_save_withoutinfo = uimenu( ...
            'Parent', handles_gui.topmenu_save, ...
            'Tag', 'topmenu_save_withoutinfo', ...
            'Label', 'Save Image Without Info', ...
            'Checked', 'off', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata)callback_topmenu_save(hObject,eventdata,false),handles_gui.figure), ...
            'Interruptible','off');  
        
        handles_gui.topmenu_save_withinfo = uimenu( ...
            'Parent', handles_gui.topmenu_save, ...
            'Tag', 'topmenu_save_withinfo', ...
            'Label', 'Save Image With Info', ...
            'Checked', 'off', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata)callback_topmenu_save(hObject,eventdata,true),handles_gui.figure), ...
            'Interruptible','off');  
        
        handles_gui.topmenu_savegif = uimenu( ...
            'Parent', handles_gui.topmenu_file, ...
            'Tag', 'topmenu_savegif', ...
            'Label', 'Save GIF', ...
            'Checked', 'off', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@(hObject,eventdata)callback_topmenu_savegif (hObject,eventdata),handles_gui.figure), ...
            'Interruptible','off');  
        
        % Panels
        handles_gui.group_menu = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_menu', ...
            'Units', 'characters', ...
            'Position', [2 23.4 35.0 20.6], ...
            'Title', 'Local Plot Options', ...
            'Interruptible','off');
        
        handles_gui.group_viewoptions = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_viewoptions', ...
            'Units', 'characters', ...
            'Position', [2 18 35.0 4.8], ...
            'Title', 'View Options', ...
            'Interruptible','off');
        
        handles_gui.group_scalebar = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_scalebar', ...
            'Units', 'characters', ...
            'Position', [2 10.8 35.0 6.5], ...
            'Title', 'Scalebar Options', ...
            'Interruptible','off');
        
        handles_gui.group_axes = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_axes', ...
            'Units', 'characters', ...
            'Position', [2 5.8 35.0 4.4], ...
            'Title', 'Axes Options', ...
            'Interruptible','off');
        
        handles_gui.group_zoompan = uibuttongroup( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_zoompan', ...
            'Units', 'characters', ...
            'Position', [2 0.8 35 4.5], ...
            'Title', 'Zoom/Pan', ...
            'Interruptible','off');
        
        if (strcmp(type_plot,'u'))
            group_title = 'U-displacements';
        elseif (strcmp(type_plot,'v'))
            group_title = 'V-displacements';
        elseif (strcmp(type_plot,'exx'))
            group_title = 'Exx-strain';
        elseif (strcmp(type_plot,'exy'))
            group_title = 'Exy-strain';  
        elseif (strcmp(type_plot,'eyy'))
            group_title = 'Eyy-strain';
        end               
        handles_gui.group_formataxes = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_formataxes', ...
            'Units', 'characters', ...
            'Position', [39.0 0.8 109 43.2], ...
            'Title', group_title, ...
            'Interruptible','off');

        % Axes
        handles_gui.axes_formatplot = axes( ...
            'Parent', handles_gui.group_formataxes, ...
            'Tag', 'axes_formatplot', ...
            'Units', 'characters', ...
            'Position', [12.4 19.0 86.2 20], ...
            'Interruptible','off');

        % Drop-down Menu  
        handles_gui.popupmenu = uicontrol( ...
            'Parent', handles_gui.group_viewoptions, ...
            'Tag', 'popupmenu', ...
            'Style', 'popupmenu', ...
            'Units', 'characters', ...
            'Position', [2.4 1.4 29.1 1.3], ...
            'BackgroundColor', [1 1 1], ...
            'String', {'Lagrangian','Eulerian'}, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_popupmenu,handles_gui.figure), ...
            'Interruptible','off');
                
        % Static Texts
        handles_gui.text_transparency = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'text_transparency', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 13.5 16 1.3], ...
            'HorizontalAlignment', 'left', ...            
            'String', 'Transparency:', ...
            'Interruptible','off');

        handles_gui.text_upperbound = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'text_upperbound', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 9.9 16 1.3], ...
            'HorizontalAlignment', 'left', ...
            'String', 'Upperbound:', ...
            'Interruptible','off');

        handles_gui.text_lowerbound = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'text_lowerbound', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 6.6 16 1.3], ...
            'HorizontalAlignment', 'left', ...
            'String', 'Lowerbound:', ...
            'Interruptible','off');      
        
        handles_gui.text_scalebarlength = uicontrol( ...
            'Parent', handles_gui.group_scalebar, ...
            'Tag', 'text_scalebarlength', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 1 16 1.3], ...
            'HorizontalAlignment', 'left', ...            
            'String', 'Scale Bar Len:', ...
            'Interruptible','off');    
        
        % Now get info text which is below the plot
        handles_infotext = set_infotext(handles_gui.group_formataxes,'');
        fields_infotext = fieldnames(handles_infotext);
        for i = 0:length(fields_infotext)-1
            handles_gui.(fields_infotext{i+1}) = handles_infotext.(fields_infotext{i+1});
        end
        
        % Sliders
        handles_gui.slider_transparency = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'slider_transparency', ...
            'Style', 'slider', ...
            'Units', 'characters', ...
            'Position', [2.5 12.0 28.9 1.3], ...
            'BackgroundColor', [0.9 0.9 0.9], ...
            'String', 'Slider', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_slider_transparency,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.slider_upperbound = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'slider_upperbound', ...
            'Style', 'slider', ...
            'Units', 'characters', ...
            'Position', [2.5 8.5 28.9 1.3], ...
            'BackgroundColor', [0.9 0.9 0.9], ...
            'String', 'Slider', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_slider_upperbound,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.slider_lowerbound = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'slider_lowerbound', ...
            'Style', 'slider', ...
            'Units', 'characters', ...
            'Position', [2.5 5.0 28.9 1.3], ...
            'BackgroundColor', [0.9 0.9 0.9], ...
            'String', 'Slider', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_slider_lowerbound,handles_gui.figure), ...
            'Interruptible','off');
        
        % Check Box        
        handles_gui.checkbox_contour = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'checkbox_contour', ...
            'Style', 'checkbox', ...
            'Units', 'characters', ...
            'Position', [2.5 17.3 17 1.3], ...
            'String', 'Contour Plot: ', ...
            'Value', 0, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_checkbox_contour,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.checkbox_minmaxmarkers = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'checkbox_minmaxmarkers', ...
            'Style', 'checkbox', ...
            'Units', 'characters', ...
            'Position', [2.5 15.5 29.7 1.3], ...
            'String', 'Max/min markers', ...
            'Value', 0, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_checkbox_minmaxmarkers,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.checkbox_scalebar = uicontrol( ...
            'Parent', handles_gui.group_scalebar, ...
            'Tag', 'checkbox_scalebar', ...
            'Style', 'checkbox', ...
            'Units', 'characters', ...
            'Position', [2.5 3 17 1.3], ...
            'String', 'Scalebar ', ...
            'Value', 1, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_checkbox_scalebar,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.checkbox_axes = uicontrol( ...
            'Parent', handles_gui.group_axes, ...
            'Tag', 'checkbox_axes', ...
            'Style', 'checkbox', ...
            'Units', 'characters', ...
            'Position', [2.5 1.0 17 1.3], ...
            'String', 'Axes ', ...
            'Value', 1, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_checkbox_axes,handles_gui.figure), ...
            'Interruptible','off');
        
        % Edit
        handles_gui.edit_transparency = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'edit_transparency', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 13.6 11.1 1.3], ...
            'HorizontalAlignment', 'left', ...            
            'String', '', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_transparency,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.edit_upperbound = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'edit_upperbound', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 10.1 11.1 1.3], ...
            'HorizontalAlignment', 'left', ...            
            'String', '', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_upperbound,handles_gui.figure), ...
            'Interruptible','off');  
        
        handles_gui.edit_lowerbound = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'edit_lowerbound', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 6.6 11.1 1.3], ...
            'HorizontalAlignment', 'left', ...            
            'String', '', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_lowerbound,handles_gui.figure), ...
            'Interruptible','off'); 
        
        handles_gui.edit_contour = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'edit_contour', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 17.3 11.1 1.3], ...
            'HorizontalAlignment', 'left', ...            
            'String', '', ...
            'Enable', 'off', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_contour,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.edit_scalebarlength = uicontrol( ...
            'Parent', handles_gui.group_scalebar, ...
            'Tag', 'edit_scalebarlength', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 1.1 11.1 1.3], ...
            'HorizontalAlignment', 'left', ...            
            'String', '', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_scalebarlength,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.edit_imgnum = uicontrol( ...
            'Parent', handles_gui.group_formataxes, ...
            'Tag', 'edit_imgnum', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [91 1.1 7 1.6], ...
            'String', '', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_imgnum,handles_gui.figure), ...
            'Enable', 'off', ...
            'Interruptible','off');
        
        % Pushbuttons
        handles_gui.button_applytoall = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'button_applytoall', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.5 2.8 28.9 1.5], ...
            'String', 'Apply to All', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_applytoall,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.button_resetdefaults = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'button_resetdefaults', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.5 1 28.9 1.5], ...
            'String', 'Reset Defaults', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_resetdefaults,handles_gui.figure), ...
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
        
        handles_gui.button_left = uicontrol( ...
            'Parent', handles_gui.group_formataxes, ...
            'Tag', 'button_left', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [84 1.0 6 1.8], ...
            'String', '<', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_left,handles_gui.figure), ...
            'Interruptible','off', ...
            'Enable', 'off');

        handles_gui.button_right = uicontrol( ...
            'Parent', handles_gui.group_formataxes, ...
            'Tag', 'button_right', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [99 1.0 6 1.8], ...
            'String', '>', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_right,handles_gui.figure), ...
            'Interruptible','off', ...
            'Enable', 'off');
    end
end

%---------------------------------------------------------------------%
% SAVE GUI -----------------------------------------------------------%
%---------------------------------------------------------------------%

function [size_savefig,timedelay,outstate] = gui_savesize(handles_gui_savefig,isgif,pos_parent)  
% This is a GUI for inputting the height and width for the saved figure.
% This gui will update the figure in handles_gui_savefig.
% 
% Inputs -----------------------------------------------------------------%
%   handles_gui_savefig - handle; handle of the save figure. The position
%     of this figure is controlled by this GUI
%   isgif - logical; this tells whether or not an animated gif is being
%   saved. If a gif is being saved then the time delay can be adjusted
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%
% Outputs ----------------------------------------------------------------%
%   size_savefig - integer array; contains the size of the save fig
%   timedelay - double; this is the time delay between frames
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.

    % Data ---------------------------------------------------------------%
    % Initialize outputs
    outstate = out.cancelled;
    size_savefig = [];
    timedelay = [];
    % Get GUI handles - Part of output
    handles_gui = init_gui();  
    % Run c-tor
    feval(ncorr_util_wrapcallbacktrycatch(@constructor,handles_gui.figure));

    % Callbacks and functions --------------------------------------------%
    function constructor()
        % Get save size parameters; the screen size is needed to
        % prevent the user from resizing outside the screen, which will
        % result in cropping when using getframe()
        % Get data
        width_img = getappdata(handles_gui_savefig.figure,'width_img');
        height_img = getappdata(handles_gui_savefig.figure,'height_img');
        spacing_save_fig = getappdata(handles_gui_savefig.figure,'spacing_save_fig');
        spacing_width_img = getappdata(handles_gui_savefig.figure,'spacing_width_img');
        spacing_height_img = getappdata(handles_gui_savefig.figure,'spacing_height_img');

        set(0,'units','pixels'); % Set to pixels first 
        size_screen = get(0,'Screensize'); % width = size_screen(3); height = size_screen(4)        
        max_width_img = size_screen(3)-spacing_save_fig-spacing_width_img;
        max_height_img = size_screen(4)-spacing_save_fig-spacing_height_img;
        if (max_width_img/max_height_img > width_img/height_img)
            max_width_img = round(max_height_img*(width_img/height_img));
        else
            max_height_img = round(max_width_img/(width_img/height_img));
        end
        if (width_img > max_width_img || height_img > max_height_img)
            % Image is too large, resize it to fit within screen
            img_reduction = min(max_width_img/width_img,max_height_img/height_img);
            width_img = round(width_img*img_reduction);
            height_img = round(height_img*img_reduction);
        end

        timedelay_prelim = 0.04;
        
        % Set data
        setappdata(handles_gui.figure,'width_img',width_img);
        setappdata(handles_gui.figure,'height_img',height_img);
        setappdata(handles_gui.figure,'ratio_img',width_img/height_img);
        setappdata(handles_gui.figure,'max_width_img',max_width_img);
        setappdata(handles_gui.figure,'max_height_img',max_height_img);  
        setappdata(handles_gui.figure,'timedelay_prelim',timedelay_prelim);  

        % Update
        update_axes('set');
        update_sidemenu();

        % Set visible
        set(handles_gui.figure,'Visible','on'); 
    end

    function callback_edit_timedelay(hObject,eventdata) %#ok<INUSD>
        % Get data
        timedelay_prelim = getappdata(handles_gui.figure,'timedelay_prelim');

        % Max and min values
        max_timedelay = 1;
        min_timedelay = 1e-5;
        
        % Get width
        timedelay_buffer = str2double(get(handles_gui.edit_timedelay,'string')); 
        if (ncorr_util_isrealbb(timedelay_buffer,min_timedelay,max_timedelay,'Time delay') == out.success)
            timedelay_prelim = timedelay_buffer;
        end

        % Set data
        setappdata(handles_gui.figure,'timedelay_prelim',timedelay_prelim);        

        % Update preview image    
        update_axes('update');    
        update_sidemenu();
    end

    function callback_edit_width(hObject,eventdata) %#ok<INUSD>
        % Get data
        ratio_img = getappdata(handles_gui.figure,'ratio_img');
        width_img = getappdata(handles_gui.figure,'width_img');
        height_img = getappdata(handles_gui.figure,'height_img');
        max_width_img = getappdata(handles_gui.figure,'max_width_img');
        max_height_img = getappdata(handles_gui.figure,'max_height_img');

        % Get width
        width_img_buffer = str2double(get(handles_gui.edit_width,'string')); 
        if (ncorr_util_isintbb(width_img_buffer,1,max_width_img,'Width') == out.success)
            % Calculate new width and height - make sure they do
            % not go outside screen size
            height_img_buffer = round(width_img_buffer/ratio_img);
            if (ncorr_util_isintbb(height_img_buffer,1,max_height_img,'Height') == out.success)
                width_img = width_img_buffer;
                height_img = height_img_buffer;
            end
        end

        % Set data
        setappdata(handles_gui.figure,'width_img',width_img);
        setappdata(handles_gui.figure,'height_img',height_img);            

        % Update preview image    
        update_axes('update');    
        update_sidemenu();
    end

    function callback_edit_height(hObject,eventdata) %#ok<INUSD>
        % Get data
        ratio_img = getappdata(handles_gui.figure,'ratio_img');
        width_img = getappdata(handles_gui.figure,'width_img');
        height_img = getappdata(handles_gui.figure,'height_img');
        max_width_img = getappdata(handles_gui.figure,'max_width_img');
        max_height_img = getappdata(handles_gui.figure,'max_height_img');

        % Get height
        height_img_buffer = str2double(get(handles_gui.edit_height,'string'));     
        if (ncorr_util_isintbb(height_img_buffer,1,max_height_img,'Height') == out.success)
            % Calculate new width and height - make sure they do
            % not go outside screen size
            width_img_buffer = round(height_img_buffer*ratio_img);
            if (ncorr_util_isintbb(width_img_buffer,1,max_width_img,'Width') == out.success)     
                width_img = width_img_buffer;
                height_img = height_img_buffer;         
            end
        end

        % Set data
        setappdata(handles_gui.figure,'width_img',width_img);
        setappdata(handles_gui.figure,'height_img',height_img);            

        % Update
        update_axes('update');  
        update_sidemenu();
    end

    function callback_button_finish(hObject,eventdata) %#ok<INUSD>
        % Get data
        width_img = getappdata(handles_gui.figure,'width_img');
        height_img = getappdata(handles_gui.figure,'height_img');
        timedelay_prelim = getappdata(handles_gui.figure,'timedelay_prelim');
        
        % Set output
        size_savefig = [width_img height_img];
        timedelay = timedelay_prelim;
        outstate = out.success;

        % Exit
        close(handles_gui.figure);
    end

    function callback_button_cancel(hObject,eventdata) %#ok<INUSD>
        % Close
        close(handles_gui.figure);
    end

    function update_sidemenu()
        % Get data
        width_img = getappdata(handles_gui.figure,'width_img');
        height_img = getappdata(handles_gui.figure,'height_img');
        timedelay_prelim = getappdata(handles_gui.figure,'timedelay_prelim');

        set(handles_gui.edit_width,'String',num2str(width_img));
        set(handles_gui.edit_height,'String',num2str(height_img));
        if (isgif)
            set(handles_gui.edit_timedelay,'String',num2str(timedelay_prelim));
        else
            set(handles_gui.edit_timedelay,'enable','off');
        end
    end

    function update_axes(action)
        % Get data - savefig handles:
        width_img = getappdata(handles_gui.figure,'width_img');
        height_img = getappdata(handles_gui.figure,'height_img');

        if (strcmp(action,'set') || strcmp(action,'update'))    
            % Transfer data then update
            setappdata(handles_gui_savefig.figure,'width_img',width_img);
            setappdata(handles_gui_savefig.figure,'height_img',height_img);
            
            % Update
            update_savefig(handles_gui_savefig);
        end
    end

    function handles_gui = init_gui()        
        % Form UIs
        handles_gui.figure = figure( ...
            'Tag', 'figure', ...
            'Units', 'characters', ...
            'Position', ncorr_util_figpos(pos_parent,[10 30.5]), ...
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
        
        % Adjust position
        % Set this UI next to the save fig        
        pos_savefig_pixels = get(handles_gui_savefig.figure,'OuterPosition');
        
        % Get size of figure in pixels
        set(handles_gui.figure,'Units','pixels');
        pos_fig_pixels = get(handles_gui.figure,'Position');
        
        pos_fig_pixels(1) = pos_savefig_pixels(3)+1.5*pos_savefig_pixels(1);
        pos_fig_pixels(2) = pos_savefig_pixels(2)+pos_savefig_pixels(4)/2-pos_fig_pixels(4)/2;
        
        % Check to make sure position is not out of screen
        set(0,'units','pixels');
        pos_screen_pixels = get(0,'ScreenSize');
        
        % Check right - no need to check bottom
        screen_right = 20; % Use this as padding
        if (pos_fig_pixels(1)+pos_fig_pixels(3) > pos_screen_pixels(3)-screen_right)
            pos_fig_pixels(1) = pos_screen_pixels(3)-pos_fig_pixels(3)-screen_right;
        end
        
        % Set position
        set(handles_gui.figure,'Position',pos_fig_pixels);
        
        % Convert back to characters
        set(handles_gui.figure,'Units','characters');        

        % Text
        handles_gui.text_timedelay = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'text_timedelay', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.1 7.5 15 1.3], ...
            'String', 'Time Delay:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_width = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'text_width', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.1 5.5 15 1.3], ...
            'String', 'Width:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_height = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'text_height', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.1 3.5 15 1.3], ...
            'String', 'Height:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
        
        % Edit box
        handles_gui.edit_timedelay = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'edit_timedelay', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [18.2 7.5 9.6 1.3], ...
            'HorizontalAlignment', 'left', ...            
            'String', '', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_timedelay,handles_gui.figure), ...
            'Interruptible','off');  

        handles_gui.edit_width = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'edit_width', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [18.2 5.5 9.6 1.3], ...
            'HorizontalAlignment', 'left', ...            
            'String', '', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_width,handles_gui.figure), ...
            'Interruptible','off');  

        handles_gui.edit_height = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'edit_height', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [18.2 3.5 9.6 1.3], ...
            'HorizontalAlignment', 'left', ...            
            'String', '', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_height,handles_gui.figure), ...
            'Interruptible','off');

        % Buttons
        handles_gui.button_finish = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'button_finish', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 1 12 1.7], ...
            'String', 'Finish', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_finish,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_cancel = uicontrol( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'button_cancel', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [16 1 12 1.7], ...
            'String', 'Cancel', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_cancel,handles_gui.figure), ...
            'Interruptible','off');
    end

    % Pause until figure is closed -----------------------------------%
    waitfor(handles_gui.figure);        
end

function update_savefig(handles_gui_savefig)  
% This function is called whenever the savefig needs to be updated. It
% will resize the window and reposition the axes/colorbar.

    % Get data
    height_data = getappdata(handles_gui_savefig.figure,'height_data');
    width_img = getappdata(handles_gui_savefig.figure,'width_img');
    height_img = getappdata(handles_gui_savefig.figure,'height_img');
    border = getappdata(handles_gui_savefig.figure,'border');
    width_cb = getappdata(handles_gui_savefig.figure,'width_cb');
    spacing_save_fig = getappdata(handles_gui_savefig.figure,'spacing_save_fig');
    spacing_width_img = getappdata(handles_gui_savefig.figure,'spacing_width_img');
    spacing_height_img = getappdata(handles_gui_savefig.figure,'spacing_height_img');
    handle_colorbar = getappdata(handles_gui_savefig.figure,'handle_colorbar');
    handle_scalebar = getappdata(handles_gui_savefig.figure,'handle_scalebar');
    handle_axes = getappdata(handles_gui_savefig.figure,'handle_axes');    
    
    % Set image
    pos_savefig(1) = spacing_save_fig;
    pos_savefig(2) = spacing_save_fig;
    pos_savefig(3) = spacing_width_img+width_img; % Width
    pos_savefig(4) = spacing_height_img+height_img; % Height
    set(handles_gui_savefig.figure,'Position',pos_savefig);

    % Colorbar:
    set(handle_colorbar,'Position',[2*border+width_img border+height_data width_cb height_img]);

    % Plot:
    set(handles_gui_savefig.axes_formatplot,'Position',[border border+height_data width_img height_img]);                  

    % Update handle_scalebar text if it it's toggled
    % Only text is updated because the size is constant during resize;
    % patches will resize with handle_axes. Do not need to check with
    % 'ishandle' because the scalebar is either set once or not at all.
    if (~isempty(handle_scalebar))
        set(handles_gui_savefig.axes_formatplot,'units','characters');
        text_scalebar = findobj(handles_gui_savefig.axes_formatplot,'tag','text_scalebar');
        pos_img = get(handles_gui_savefig.axes_formatplot,'Position');
        size_text = max(0.5*max(0.5*pos_img(4),0.14*pos_img(3)),8);
        set(text_scalebar,'FontSize',size_text);
        set(handles_gui_savefig.axes_formatplot,'units','pixels');
    end

    % Update handle_axes text if it it's toggled
    % Only text is updated because the size is constant during resize;
    % patches will resize with handle_axes. Do not need to check with
    % 'ishandle' because the axes is either set once or not at all.
    if (~isempty(handle_axes))
        set(handles_gui_savefig.axes_formatplot,'units','characters');
        text_x = findobj(handles_gui_savefig.axes_formatplot,'tag','text_x');
        text_y = findobj(handles_gui_savefig.axes_formatplot,'tag','text_y');
        pos_img = get(handles_gui_savefig.axes_formatplot,'Position');
        size_text = max(0.5*max(0.5*pos_img(4),0.14*pos_img(3)),8);
        set(text_x,'FontSize',size_text);
        set(text_y,'FontSize',size_text);
        set(handles_gui_savefig.axes_formatplot,'units','pixels');
    end
    
    % Make sure figure is drawn to ensure getframe works properly.
    drawnow;
end

function handles_infotext = set_infotext(parent,bgcolor)
% This function sets information text
    handles_infotext.text_ref_name = uicontrol( ...
        'Parent', parent, ...
        'Tag', 'ref_name', ...
        'Style', 'text', ...
        'Units', 'characters', ...
        'Position', [2.1 17.0 100 1.3], ...
        'String', 'Reference Name: ', ...
        'HorizontalAlignment', 'left', ...
        'Interruptible','off');

    handles_infotext.text_cur_name = uicontrol( ...
        'Parent', parent, ...
        'Tag', 'cur_name', ...
        'Style', 'text', ...
        'Units', 'characters', ...
        'Position', [2.1 15.5 100 1.3], ...
        'String', 'Current Name: ', ...
        'HorizontalAlignment', 'left', ...
        'Interruptible','off');
    
    handles_infotext.text_type = uicontrol( ...
        'Parent', parent, ...
        'Tag', 'text_type', ...
        'Style', 'text', ...
        'Units', 'characters', ...
        'Position', [2.1 14.0 100 1.3], ...
        'String', 'Analysis Type: ', ...
        'HorizontalAlignment', 'left', ...
        'Interruptible','off');
    
    handles_infotext.text_dicparams = uicontrol( ...
        'Parent', parent, ...
        'Tag', 'text_dicparams', ...
        'Style', 'text', ...
        'Units', 'characters', ...
        'Position', [2.1 12.5 100 1.3], ...
        'String', 'RG-DIC Radius: | Subset Spacing:', ...
        'HorizontalAlignment', 'left', ...
        'Interruptible','off');
    
    handles_infotext.text_itparams = uicontrol( ...
        'Parent', parent, ...
        'Tag', 'text_itparams', ...
        'Style', 'text', ...
        'Units', 'characters', ...
        'Position', [2.1 11.0 100 1.3], ...
        'String', 'Diffnorm Cutoff: | Iteration Cutoff: | Total Threads: ', ...
        'HorizontalAlignment', 'left', ...
        'Interruptible','off');
    
    handles_infotext.text_stepanalysis = uicontrol( ...
        'Parent', parent, ...
        'Tag', 'text_stepanalysis', ...
        'Style', 'text', ...
        'Units', 'characters', ...
        'Position', [2.1 9.5 100 1.3], ...
        'String', 'Step Analysis: ', ...
        'HorizontalAlignment', 'left', ...
        'Interruptible','off');
   
    handles_infotext.text_subsettrunc = uicontrol( ...
        'Parent', parent, ...
        'Tag', 'text_subsettrunc', ...
        'Style', 'text', ...
        'Units', 'characters', ...
        'Position', [2.1 8.0 100 1.3], ...
        'String', 'Disp Subset Truncation: ', ...
        'HorizontalAlignment', 'left', ...
        'Interruptible','off');
    
    handles_infotext.text_imgcorr = uicontrol( ...
        'Parent', parent, ...
        'Tag', 'text_imgcorr', ...
        'Style', 'text', ...
        'Units', 'characters', ...
        'Position', [2.1 6.5 100 1.3], ...
        'String', 'Image Correspondences: ', ...
        'HorizontalAlignment', 'left', ...
        'Interruptible','off');
    
    handles_infotext.text_pixtounits = uicontrol( ...
        'Parent', parent, ...
        'Tag', 'text_pixtounits', ...
        'Style', 'text', ...
        'Units', 'characters', ...
        'Position', [2.1 5.0 100 1.3], ...
        'String', 'Reference Image FOV width: ', ...
        'HorizontalAlignment', 'left', ...
        'Interruptible','off');

    handles_infotext.text_cutoff_corrcoef = uicontrol( ...
        'Parent', parent, ...
        'Tag', 'text_cutoff_corrcoef', ...
        'Style', 'text', ...
        'Units', 'characters', ...
        'Position', [2.1 3.5 100 1.3], ...
        'String', 'Correlation Coefficient Cutoff: ', ...
        'HorizontalAlignment', 'left', ...
        'Interruptible','off');

    handles_infotext.text_lenscoef = uicontrol( ...
        'Parent', parent, ...
        'Tag', 'text_lenscoef', ...
        'Style', 'text', ...
        'Units', 'characters', ...
        'Position', [2.1 2.0 100 1.3], ...
        'String', 'Radial Lens Distortion Coefficient: ', ...
        'HorizontalAlignment', 'left', ...
        'Interruptible','off');

    handles_infotext.text_minmedianmax = uicontrol( ...
        'Parent', parent, ...
        'Tag', 'text_max', ...
        'Style', 'text', ...
        'Units', 'characters', ...
        'Position', [2.1 0.5 100 1.3], ...
        'String', 'Min: | Median: | Max: ', ...
        'HorizontalAlignment', 'left', ...
        'Interruptible','off'); 
    
    % Check if bgcolor was provided
    if (~isempty(bgcolor))
        fields_infotext = fieldnames(handles_infotext);
        for i = 0:length(fields_infotext)-1
            set(handles_infotext.(fields_infotext{i+1}),'BackgroundColor',bgcolor);
        end
    end
end

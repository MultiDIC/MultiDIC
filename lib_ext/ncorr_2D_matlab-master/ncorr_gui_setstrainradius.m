function [radius,subsettrunc,outstate] = ncorr_gui_setstrainradius(reference,current,rois_ref,rois_cur,plots_u_ref,plots_v_ref,plots_u_cur,plots_v_cur,spacing,pixtounits,units,pos_parent,params_init)
% This is a GUI for setting the strain radius.
%
% Inputs -----------------------------------------------------------------%
%   reference - ncorr_class_img; used for displaying the background image.
%   current - ncorr_class_img; used for displaying the background image.
%   rois_ref - ncorr_class_roi; contains ROIs WRT the reference image. Note 
%   that rois_ref has already been reduced.
%   rois_cur - ncorr_class_roi; contains ROIs WRT the current image. Note 
%   that rois_ref has already been reduced.
%   plots_u_ref - cell; contains displacement data WRT the reference image.
%   Uses real units. Note that displacement plots are reduced by default.
%   plots_v_ref - cell; contains displacement data WRT the reference image.
%   Uses real units. Note that displacement plots are reduced bym default.
%   plots_u_cur - cell; contains displacement data WRT the current image. 
%   Uses real units. Note that displacement plots are reduced by default.
%   plots_v_cur - cell; contains displacement data WRT the current image. 
%   Uses real units. Note that displacement plots are reduced by default.
%   spacing - double; spacing parameter
%   pixtounits - double; units/pixels conversion factor
%   units - string; this a string containing the units
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%   params_init - integer; contains [radius subsettrunc] if they have been set before.
%
% Outputs ----------------------------------------------------------------%
%   radius - integer; strain radius
%   subsettrunc - logical; true if subset truncation is set
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.
%
% Note that some regions for the input ROIs can be empty. This is handled
% properly.

    % Data ---------------------------------------------------------------%
    % Initialize outputs
    outstate = out.cancelled;
    radius = [];
    subsettrunc = [];
    % Get GUI handles
    handles_gui = init_gui();
    % Run c-tor
    feval(ncorr_util_wrapcallbacktrycatch(@constructor,handles_gui.figure));
    
    % Callbacks and functions --------------------------------------------%
    function constructor()
        % Set zoom and pan
        handle_zoom = zoom(handles_gui.figure);  
        handle_pan = pan(handles_gui.figure);   
        
        % Allow rotations - For some reason when rotations are enabled, the
        % impoint is no longer draggable, so only allow rotations when the
        % cursor is over the least squares plane fit axes.
        handle_rotate3d = rotate3d(handles_gui.axes_leastsquares);
        set(handles_gui.figure,'WindowButtonMotionFcn',ncorr_util_wrapcallbacktrycatch(@callback_moveplot,handles_gui.figure));
        
        % Set buffers        
        % Min radius and spacing
        min_radius = 1;
        max_radius = 100;
        
        % Create radius function this will take the input from the slider
        % and return a radius
        get_radius = @(x)round(x*50+min_radius);

        % Initialize strain radius and slider
        if (isempty(params_init))
            radius_prelim = 15;
            radius_slider_prelim = min(1,max(0,fzero(@(x) get_radius(x)-radius_prelim,0)));
            subsettrunc_prelim = false;
        else
            radius_prelim = params_init(1);
            radius_slider_prelim = min(1,max(0,fzero(@(x) get_radius(x)-radius_prelim,0)));
            subsettrunc_prelim = params_init(2);
        end
        
        % Set data
        setappdata(handles_gui.figure,'min_radius',min_radius);
        setappdata(handles_gui.figure,'max_radius',max_radius);
        setappdata(handles_gui.figure,'radius_prelim',radius_prelim);
        setappdata(handles_gui.figure,'radius_slider_prelim',radius_slider_prelim);
        setappdata(handles_gui.figure,'subsettrunc_prelim',subsettrunc_prelim);
        setappdata(handles_gui.figure,'get_radius',get_radius);
        setappdata(handles_gui.figure,'num_cur',length(current)-1);  
        setappdata(handles_gui.figure,'handle_points',[]);   
        setappdata(handles_gui.figure,'handle_preview',[]);
        setappdata(handles_gui.figure,'handle_rotate3d',handle_rotate3d);
        setappdata(handles_gui.figure,'val_popupmenu_disp',1);
        setappdata(handles_gui.figure,'val_popupmenu_lore',1);
        setappdata(handles_gui.figure,'pos_points',[]);
        setappdata(handles_gui.figure,'num_region',0);  
        setappdata(handles_gui.figure,'handle_zoom',handle_zoom);
        setappdata(handles_gui.figure,'handle_pan',handle_pan);
        
        % Update
        update_axes('set');  
        update_sidemenu();        
        
        % Set Visible
        set(handles_gui.figure,'Visible','on'); 
    end    

    function callback_popupmenu_lore(hObject,eventdata) %#ok<INUSD>
        % Get data
        val_popupmenu_lore = get(handles_gui.popupmenu_lore,'value');
        
        % Clear position of points
        pos_points = [];            
        
        % Set data
        setappdata(handles_gui.figure,'pos_points',pos_points);
        setappdata(handles_gui.figure,'val_popupmenu_lore',val_popupmenu_lore);
        
        % Update
        update_axes('set');
        update_sidemenu();
    end

    function callback_popupmenu_disp(hObject,eventdata) %#ok<INUSD>
        % Get data
        val_popupmenu_disp = get(handles_gui.popupmenu_disp,'value');
        
        % Set data
        setappdata(handles_gui.figure,'val_popupmenu_disp',val_popupmenu_disp);
        
        % Update
        update_axes('set');
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
        update_axes('update');  
        update_sidemenu();        
    end
    
    function callback_checkbox_subsettrunc(hObject,eventdata) %#ok<INUSD>        
        % Get value
        subsettrunc_prelim = get(handles_gui.checkbox_subsettrunc,'Value');
        
        % Set data
        setappdata(handles_gui.figure,'subsettrunc_prelim',subsettrunc_prelim);  
        
        % Update
        update_axes('update');
        update_sidemenu();        
    end

    function callback_edit_radius(hObject,eventdata) %#ok<INUSD>
        min_radius = getappdata(handles_gui.figure,'min_radius');
        max_radius = getappdata(handles_gui.figure,'max_radius');
        radius_prelim = getappdata(handles_gui.figure,'radius_prelim');
        radius_slider_prelim = getappdata(handles_gui.figure,'radius_slider_prelim');
        get_radius = getappdata(handles_gui.figure,'get_radius');
        
        % Get radius
        radius_buffer = str2double(get(handles_gui.edit_radius,'string')); 
        if (ncorr_util_isintbb(radius_buffer,min_radius,max_radius,'Radius') == out.success)
            % Set radius
            radius_prelim = radius_buffer;  
            radius_slider_prelim = min(1,max(0,fzero(@(x)get_radius(x)-radius_prelim,0)));
        end
        
        % Set data
        setappdata(handles_gui.figure,'radius_prelim',radius_prelim);
        setappdata(handles_gui.figure,'radius_slider_prelim',radius_slider_prelim); 
                        
        % Update
        update_axes('update');
        update_sidemenu(); 
    end

    function callback_edit_imgnum(hObject,eventdata) %#ok<INUSD>
        % Get Data
        num_cur = getappdata(handles_gui.figure,'num_cur');
        pos_points = getappdata(handles_gui.figure,'pos_points');

        % Get img num
        num_cur_prelim = str2double(get(handles_gui.edit_imgnum,'string')); 
        if (ncorr_util_isintbb(num_cur_prelim,1,length(current),'Current image number') == out.success)  
            % Clear position of points            
            pos_points = [];   
            
            % num_cur_prelim is zero based indexed
            num_cur = num_cur_prelim-1;
        end    
        
        % Set Data
        setappdata(handles_gui.figure,'num_cur',num_cur);
        setappdata(handles_gui.figure,'pos_points',pos_points);
            
        % Update
        update_axes('set');
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
        update_sidemenu();
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
        update_sidemenu();
    end

    function callback_button_finish(hObject,eventdata) %#ok<INUSD>
        % Get data
        radius_prelim = getappdata(handles_gui.figure,'radius_prelim');
        subsettrunc_prelim = getappdata(handles_gui.figure,'subsettrunc_prelim');
                
        if (subsettrunc_prelim)
            text_st = 'Enabled';
        else
            text_st = 'Disabled';
        end
        contbutton = questdlg({['Strain radius is set to: ' num2str(radius_prelim) '.'], ...
                               ['Subset Truncation is: ' text_st '.'], ...
                               'Is this correct?'}, ...
                               'Continue Operation','Yes','No','No');
        if (strcmp(contbutton,'Yes'))
            % Set Outputs
            radius = radius_prelim;
            subsettrunc = subsettrunc_prelim;
            outstate = out.success;
            
            % Exit                        
            close(handles_gui.figure);
        end
    end

    function callback_button_cancel(hObject,eventdata) %#ok<INUSD>
        close(handles_gui.figure); 
    end

    function callback_button_left(hObject,eventdata) %#ok<INUSD>
        % Get Data
        num_cur = getappdata(handles_gui.figure,'num_cur');
        pos_points = getappdata(handles_gui.figure,'pos_points');
        
        % Check for overshoot
        if (num_cur > 0)
            % Clear position of points
            pos_points = [];
            
            num_cur = num_cur-1;
        end        
        
        % Set Data
        setappdata(handles_gui.figure,'num_cur',num_cur);   
        setappdata(handles_gui.figure,'pos_points',pos_points);  
        
        % Update
        update_axes('set');
        update_sidemenu();
    end

    function callback_button_right(hObject,eventdata) %#ok<INUSD>
        % Get Data
        num_cur = getappdata(handles_gui.figure,'num_cur');
        pos_points = getappdata(handles_gui.figure,'pos_points');
        
        % Check for overshoot
        if (num_cur < length(current)-1)
            % Clear position of points
            pos_points = [];
            
            num_cur = num_cur+1;
        end                 
        
        % Set Data
        setappdata(handles_gui.figure,'num_cur',num_cur);  
        setappdata(handles_gui.figure,'pos_points',pos_points);  
            
        % Update
        update_axes('set');
        update_sidemenu();
    end   

    function callback_impoint(pos,num_region)
        % Get data
        pos_points = getappdata(handles_gui.figure,'pos_points');
        
        % Store position point
        pos_points(num_region+1,:) = round(pos)-1;
        
        % Set Data 
        setappdata(handles_gui.figure,'pos_points',pos_points);
        setappdata(handles_gui.figure,'num_region',num_region);
        
        % Update
        update_axes('update');        
    end

    function callback_moveplot(hObject,eventdata) %#ok<INUSD>
        % Get data
        handle_rotate3d = getappdata(handles_gui.figure,'handle_rotate3d');
        handle_zoom = getappdata(handles_gui.figure,'handle_zoom'); 
        handle_pan = getappdata(handles_gui.figure,'handle_pan'); 
                
        % Enable rotation if zoom and pan are disabled
        pos_mouse = get(handles_gui.figure,'CurrentPoint');
        pos_group = get(handles_gui.group_preview,'Position');
        pos_axes = get(handles_gui.axes_leastsquares,'Position');
        inleastsquares = (pos_mouse(1) > pos_group(1)+pos_axes(1) && pos_mouse(1) < pos_group(1)+pos_axes(1)+pos_axes(3) && pos_mouse(2) > pos_group(2)+pos_axes(2) && pos_mouse(2) < pos_group(2)+pos_axes(2)+pos_axes(4));
        if (inleastsquares)
            set(handle_rotate3d,'Enable','on');
        end
            
        % If pan and zoom are both disabled and cursor is outside of least
        % square axes, then disable rotation
        if (strcmp(get(handle_zoom,'Enable'),'off') && strcmp(get(handle_pan,'Enable'),'off') && ~inleastsquares)
            set(handle_rotate3d,'Enable','off');
        end
        
        % Update
        update_sidemenu();
    end
    
    function update_sidemenu()
        % Get data
        radius_prelim = getappdata(handles_gui.figure,'radius_prelim');
        radius_slider_prelim = getappdata(handles_gui.figure,'radius_slider_prelim');
        val_popupmenu_disp = getappdata(handles_gui.figure,'val_popupmenu_disp');
        val_popupmenu_lore = getappdata(handles_gui.figure,'val_popupmenu_lore');
        handle_zoom = getappdata(handles_gui.figure,'handle_zoom'); 
        handle_pan = getappdata(handles_gui.figure,'handle_pan'); 
        subsettrunc_prelim = getappdata(handles_gui.figure,'subsettrunc_prelim'); 
                
        % Update slider
        set(handles_gui.slider_radius,'value',radius_slider_prelim);
        
        % Update edit field
        set(handles_gui.edit_radius,'String',num2str(radius_prelim))
        
        % Update popupmenus
        set(handles_gui.popupmenu_disp,'value',val_popupmenu_disp);
        set(handles_gui.popupmenu_lore,'value',val_popupmenu_lore);
                       
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
        
        set(handles_gui.checkbox_subsettrunc,'Value',subsettrunc_prelim);   
    end
    
    function update_axes(action) 
        % Get data
        radius_prelim = getappdata(handles_gui.figure,'radius_prelim');
        num_cur = getappdata(handles_gui.figure,'num_cur'); 
        handle_points = getappdata(handles_gui.figure,'handle_points');   
        handle_preview = getappdata(handles_gui.figure,'handle_preview');
        val_popupmenu_disp = getappdata(handles_gui.figure,'val_popupmenu_disp');
        val_popupmenu_lore = getappdata(handles_gui.figure,'val_popupmenu_lore');
        pos_points = getappdata(handles_gui.figure,'pos_points');
        num_region = getappdata(handles_gui.figure,'num_region');
        subsettrunc_prelim = getappdata(handles_gui.figure,'subsettrunc_prelim');
        
        % Get img, data plot, and ROI
        if (val_popupmenu_disp == 1 && val_popupmenu_lore == 1)
            img_bg = reference;
            plot_data = plots_u_ref{num_cur+1};
            roi_data = rois_ref(num_cur+1);
        elseif (val_popupmenu_disp == 2 && val_popupmenu_lore == 1)
            img_bg = reference;
            plot_data = plots_v_ref{num_cur+1};
            roi_data = rois_ref(num_cur+1);
        elseif (val_popupmenu_disp == 1 && val_popupmenu_lore == 2)
            img_bg = current(num_cur+1);
            plot_data = plots_u_cur{num_cur+1};
            roi_data = rois_cur(num_cur+1);
        else
            img_bg = current(num_cur+1);
            plot_data = plots_v_cur{num_cur+1};
            roi_data = rois_cur(num_cur+1);
        end      
            
        if (strcmp(action,'set'))  
            % Get reduced img
            img_reduced = img_bg.reduce(spacing);
            
            % Initialize displacements
            imshow(img_reduced.get_img(),[img_reduced.min_gs img_reduced.max_gs],'Parent',handles_gui.axes_preview);
            hold(handles_gui.axes_preview,'on');
            
            % Get bounds
            cmin = prctile(plot_data(roi_data.mask),5)-1e-5;
            cmax = prctile(plot_data(roi_data.mask),95)+1e-5;
            if (abs(cmin-cmax)<1e-4)
                cmax = cmax + 1e-5; % add a small increment to insure they are not the same number
                cmin = cmin - 1e-5;
            end

            handle_preview = imshow(plot_data,'Parent',handles_gui.axes_preview); 
            set(handles_gui.axes_preview,'CLim',[cmin cmax]);            
            hold(handles_gui.axes_preview,'off');  
                                    
            % Set colorbar
            handle_colorbar = colorbar('peer',handles_gui.axes_preview);        
            set(handle_colorbar,'UIContextMenu','');
            set(get(handle_colorbar,'child'),'YData',[cmin cmax]);
            set(handle_colorbar,'YLim',[cmin cmax]);
            set(handle_colorbar,'Units','Pixels');
            
            % Set Impoints - only set impoints in regions which are full - use
            % a cell array here since it is not gauranteed that regions are
            % full.
            if (isempty(pos_points))
                setnewpos = true;
            else
                setnewpos = false;
            end
            
            handle_points = cell(1,length(roi_data.region));
            for i = 0:length(roi_data.region)-1
                if (roi_data.region(i+1).totalpoints > 0)
                    % Form point
                    handle_points{i+1} = impoint(handles_gui.axes_preview,round(size(roi_data.mask,2)/2),round(size(roi_data.mask,1)/2));
                    setColor(handle_points{i+1},'g');
                    set(handle_points{i+1},'UIContextMenu','');

                    % Set constraint function 
                    api_point = iptgetapi(handle_points{i+1});
                    constrainfcn = ncorr_util_formregionconstraint(roi_data.region(i+1));
                    api_point.setPositionConstraintFcn(constrainfcn);

                    % Initialize point location - use 1 based indexing for
                    % input to constrainfcn
                    x_center = round(roi_data.region(i+1).leftbound+(roi_data.region(i+1).rightbound-roi_data.region(i+1).leftbound)/2)+1;
                    y_center = round(roi_data.region(i+1).upperbound+(roi_data.region(i+1).lowerbound-roi_data.region(i+1).upperbound)/2)+1;
                    
                    % See if points already exist
                    if (setnewpos)
                        % constrainfcn returns points based on 1 based
                        % indexing, so convert it to zero based.
                        pos_points(i+1,1:2) = constrainfcn([x_center y_center])-1;
                    end
                    % Set position based on 1 based indexing
                    setPosition(handle_points{i+1},pos_points(i+1,1:2)+1);

                    % Assign Position callback
                    addNewPositionCallback(handle_points{i+1},ncorr_util_wrapcallbacktrycatch(@(pos)callback_impoint(pos,i),handles_gui.figure));
                    
                    % Set num_region
                    if (setnewpos)
                        num_region = i;
                    end  
                end
            end  
            
            % Set colormap
            ncorr_util_colormap(handles_gui.figure);
                        
            % Update buttons
            set(handles_gui.edit_imgnum,'String',num2str(num_cur+1));
            if (length(current) == 1)
                set(handles_gui.button_left,'Enable','off');
                set(handles_gui.button_right,'Enable','off');
                set(handles_gui.edit_imgnum,'Enable','off');
            elseif (num_cur == 0)
                set(handles_gui.button_left,'Enable','off');
                set(handles_gui.button_right,'Enable','on');
                set(handles_gui.edit_imgnum,'Enable','on');
            elseif (num_cur == length(current)-1)
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
            % Get cirroi
            cirroi = roi_data.get_cirroi(pos_points(num_region+1,1),pos_points(num_region+1,2),radius_prelim,subsettrunc_prelim);
            
            % Initialize matrix and vector to zeros
            mat_LS = zeros(3);
            z_vec_LS = zeros(3,1);

            % Initialize alphamask
            alpha_preview = 0.5*roi_data.mask;

            % Preallocate arrays
            plotpoints_x = zeros(cirroi.region.totalpoints,1);
            plotpoints_y = zeros(cirroi.region.totalpoints,1);
            plotpoints_z = zeros(cirroi.region.totalpoints,1);

            % Initialize counter
            counter = 0;      
            for i = 0:size(cirroi.region.noderange,1)-1
                x = i+(cirroi.x-cirroi.radius);
                for j = 0:2:cirroi.region.noderange(i+1)-1                                      
                    y_vec = cirroi.region.nodelist(i+1,j+1):cirroi.region.nodelist(i+1,j+2);

                    % Form cirroi.x,cirroi.y,z vectors                        
                    plotpoints_x(counter+1:counter+length(y_vec)) = x*ones(size(y_vec))*(pixtounits*(spacing+1));
                    plotpoints_y(counter+1:counter+length(y_vec)) = y_vec*(pixtounits*(spacing+1));
                    plotpoints_z(counter+1:counter+length(y_vec)) = plot_data(y_vec+1,x+1)';

                    % Do least squares analysis to find plane
                    x_LS = (x*ones(size(y_vec))-cirroi.x)*(pixtounits*(spacing+1));
                    y_LS = (y_vec-cirroi.y)*(pixtounits*(spacing+1));
                    z_LS = plot_data(y_vec+1,x+1)';

                    % Do the matrix first
                    mat_LS(1,1) = mat_LS(1,1)+sum(x_LS.^2);
                    mat_LS(2,1) = mat_LS(2,1)+sum(x_LS.*y_LS);
                    mat_LS(3,1) = mat_LS(3,1)+sum(x_LS);
                    mat_LS(1,2) = mat_LS(1,2)+sum(x_LS.*y_LS);
                    mat_LS(2,2) = mat_LS(2,2)+sum(y_LS.^2);
                    mat_LS(3,2) = mat_LS(3,2)+sum(y_LS);                                    
                    mat_LS(1,3) = mat_LS(1,3)+sum(x_LS);
                    mat_LS(2,3) = mat_LS(2,3)+sum(y_LS);
                    mat_LS(3,3) = mat_LS(3,3)+sum(ones(size(y_vec)));

                    % Do vec next
                    z_vec_LS(1) = z_vec_LS(1)+sum(x_LS.*z_LS);
                    z_vec_LS(2) = z_vec_LS(2)+sum(y_LS.*z_LS);
                    z_vec_LS(3) = z_vec_LS(3)+sum(z_LS);

                    % Update Preview Mask
                    alpha_preview(y_vec+1,x+1) = 1;

                    % Update counter
                    counter = counter + length(y_vec);
                end
            end

            % Set Alpha Mask
            set(handle_preview,'AlphaData',alpha_preview);
            
            if (abs(det(mat_LS)) > 10^-10) % Make sure matrix is not singular
                % Find dp/dx dp/cirroi.y 
                du_vec = linsolve(mat_LS,z_vec_LS);

                % Plot points first            
                plot3(plotpoints_x,plotpoints_y,plotpoints_z,'o','Parent',handles_gui.axes_leastsquares); 
                zlabel(handles_gui.axes_leastsquares,['Displacement - ' units]); 
                hold(handles_gui.axes_leastsquares,'on');
                
                % Plot the plane next
                [plane_x,plane_y] = meshgrid((-cirroi.radius:cirroi.radius)*(pixtounits*(spacing+1)),(-cirroi.radius:cirroi.radius)*(pixtounits*(spacing+1)));
                plane_z = du_vec(3)+du_vec(2)*plane_y+du_vec(1)*plane_x;
                mesh(plane_x+cirroi.x*(pixtounits*(spacing+1)),plane_y+cirroi.y*(pixtounits*(spacing+1)),plane_z,'Parent',handles_gui.axes_leastsquares);   
                
                hold(handles_gui.axes_leastsquares,'off');   
            else
                cla(handles_gui.axes_leastsquares);
            end
                        
            % Set Name
            set(handles_gui.text_ref,'String',['Reference Name: ' reference.name(1:min(end,22))]);
            set(handles_gui.text_cur,'String',['Current Name: ' current(num_cur+1).name(1:min(end,22))]);
        end
        
        % Set data
        setappdata(handles_gui.figure,'handle_points',handle_points);   
        setappdata(handles_gui.figure,'handle_preview',handle_preview);
        setappdata(handles_gui.figure,'pos_points',pos_points);
        setappdata(handles_gui.figure,'num_region',num_region); 
    end   

    function handles_gui = init_gui()
    % GUI controls -------------------------------------------------------%
        % Figure
        handles_gui.figure = figure( ...
            'Tag', 'figure', ...
            'Units', 'characters', ...
            'Position', ncorr_util_figpos(pos_parent,[32.5 167]), ...
            'Name', 'Set Strain Parameters', ...
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
        handles_gui.group_strainoptions = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_strainoptions', ...
            'Units', 'characters', ...
            'Position', [2 25.5 35 6.4], ...
            'Title', 'Strain Options', ...
            'Interruptible','off');
        
        handles_gui.group_viewoptions = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_viewoptions', ...
            'Units', 'characters', ...
            'Position', [2 17.8 35 7.2], ...
            'Title', 'View Options', ...
            'Interruptible','off');
        
        handles_gui.group_subsettrunc = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_subsettrunc', ...
            'Units', 'characters', ...
            'Position', [2 12.8 35 4.3], ...
            'Title', 'Discontinuous Analysis', ...
            'Interruptible','off');          

        handles_gui.group_zoompan = uibuttongroup( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_zoompan', ...
            'Units', 'characters', ...
            'Position', [2 8.0 35 4.5], ...
            'Title', 'Zoom/Pan', ...
            'Interruptible','off');   
        
        handles_gui.group_menu = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_menu', ...
            'Units', 'characters', ...
            'Position', [2 1 35 6.7], ...
            'Title', 'Menu', ...
            'Interruptible','off');

        handles_gui.group_preview = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_preview', ...
            'Units', 'characters', ...
            'Position', [39 0.75 126 31.2], ...
            'Title', 'Preview', ...
            'Interruptible','off');

        % Axes
        handles_gui.axes_preview = axes( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'axes_preview', ...
            'Units', 'characters', ...
            'Visible','off', ...
            'Position', [2 4 57 23.5], ...
            'Interruptible','off');

        handles_gui.axes_leastsquares = axes( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'axes_leastsquares', ...
            'Units', 'characters', ...
            'Visible', 'off', ...
            'Position', [74.8 6 38 19.5], ...
            'Interruptible','off');
        
        % Drop-down Menu  
        handles_gui.popupmenu_lore = uicontrol( ...
            'Parent', handles_gui.group_viewoptions, ...
            'Tag', 'popupmenu_lore', ...
            'Style', 'popupmenu', ...
            'Units', 'characters', ...
            'Position', [2.4 3.6 29.1 1.3], ...
            'BackgroundColor', [1 1 1], ...
            'String', {'Lagrangian','Eulerian'}, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_popupmenu_lore,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.popupmenu_disp = uicontrol( ...
            'Parent', handles_gui.group_viewoptions, ...
            'Tag', 'popupmenu_disp', ...
            'Style', 'popupmenu', ...
            'Units', 'characters', ...
            'Position', [2.4 1.5 29.1 1.3], ...
            'BackgroundColor', [1 1 1], ...
            'String', {'U-Displacement','V-Displacement'}, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_popupmenu_disp,handles_gui.figure), ...
            'Interruptible','off');
        
        % Static Texts
        handles_gui.text_radius = uicontrol( ...
            'Parent', handles_gui.group_strainoptions, ...
            'Tag', 'text_radius', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 2.6 16 1.3], ...
            'String', 'Strain Radius: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_ref = uicontrol( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'text_ref', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [3.0 1.8 56.9 1.3], ...
            'String', 'Reference Name: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
        
        handles_gui.text_cur = uicontrol( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'text_cur', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [3.0 0.6 56.9 1.3], ...
            'String', 'Current Name: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');        

        % Sliders
        handles_gui.slider_radius = uicontrol( ...
            'Parent', handles_gui.group_strainoptions, ...
            'Tag', 'slider_radius', ...
            'Style', 'slider', ...
            'Units', 'characters', ...
            'Position', [2.5 1.2 28.9 1.3], ...
            'BackgroundColor', [0.9 0.9 0.9], ...
            'String', {'Slider'}, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_slider_radius,handles_gui.figure), ...
            'Interruptible','off');
        
        % Checkbox
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

        % Edit texts
        handles_gui.edit_radius = uicontrol( ...
            'Parent', handles_gui.group_strainoptions, ...
            'Tag', 'edit_radius', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 2.8 11.1 1.3], ...
            'String', '', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_radius,handles_gui.figure), ...
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

function [plots_disp_f,rois_f,pixtounits,units,cutoff_corrcoef,lenscoef,outstate] = ncorr_gui_formatdisp(reference,current,roi,plots_u_ref,plots_v_ref,plots_corrcoef_ref,spacing,pos_parent,params_init)
% This is a GUI for formatting the displacements. This can filter out
% displacements with high correlation coefficient values, convert
% displacements from pixels to real units, and apply very basic corrections
% for lens distortion (although this is generally skipped). The number of
% reference images is equal to the number of current images.
%
% Inputs -----------------------------------------------------------------%
%   reference - ncorr_class_img; used for displaying the background image.
%   current - ncorr_class_img; used for displaying the background image.
%   roi - ncorr_class_roi; ROI WRT the reference image. Note that roi has 
%   already been reduced.
%   plots_u_ref - cell; cell array of u displacement plots WRT the reference
%   image. Units are pixels. Note that displacement plots are reduced by
%   default.
%   plots_v_ref - cell; cell array of v displacement plots WRT the reference
%   image. Units are pixels. Note that displacement plots are reduced by
%   default.
%   plots_corrcoef_ref - cell; cell array of correlation coefficient plots
%   WRT the reference image. Note that corrcoef plots are reduced by
%   default.
%   spacing - integer; spacing parameter
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%   params_init - cell; params_init{1} is the units/pixel factor, params_init{2} is
%   the units, params_init{3} is the correlation coefficient cutoff values,
%   and param_init{4} is the lens coefficient if these values have been set
%   already
%
% Outputs ----------------------------------------------------------------%
%   plots_disp_f - struct; contains the formatted displacement data
%   rois_f - ncorr_class_roi; updated ROI for displaying displacements
%   within the correlation coefficient range
%   pixtounits - double; ratio of units to pixels.
%   units - string; name of units
%   cutoff_corrcoef - double array; correlation coefficient cutoffs as determined
%   by the user
%   lenscoef - double; lens coefficient
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.

    % Data ---------------------------------------------------------------%
    % Initialize outputs
    outstate = out.cancelled;
    plots_disp_f = struct('plot_u_formatted',{},'plot_v_formatted',{});
    rois_f = ncorr_class_roi.empty;
    pixtounits = [];
    units = '';
    cutoff_corrcoef = [];
    lenscoef = [];
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
        linkaxes([handles_gui.axes_formatplot_u handles_gui.axes_formatplot_v],'xy');
                
        % Form buffers
        if (isempty(params_init))
            pixtounits_prelim = 1;    % Just keep pixels if nothing is specified
            units_prelim = 'pixels';
            lenscoef_prelim = 0;
        else
            pixtounits_prelim = params_init{1};
            units_prelim = params_init{2};
            lenscoef_prelim = params_init{4};
        end
        
        % Set slider buffer
        corrcoef_buffer = zeros(length(reference),3); % [corrcoef_slider_value corrcoef_max corrcoef_min]        
        for i = 0:length(reference)-1
            cc_max = prctile(plots_corrcoef_ref{i+1}(roi(i+1).mask),100);
            cc_min = prctile(plots_corrcoef_ref{i+1}(roi(i+1).mask),5); % Set this to 5 so user can't filter out all datapoints
            if (isempty(params_init) || abs(cc_max-cc_min) <= 1e-10)
                corrcoef_buffer(i+1,1) = 1;
            else
                corrcoef_buffer(i+1,1) = (params_init{3}(i+1)-cc_min)/(cc_max-cc_min);
            end
            corrcoef_buffer(i+1,2) = cc_max;
            corrcoef_buffer(i+1,3) = cc_min;
        end
        
        % Cutoffs
        min_lenscoef = -1;
        max_lenscoef = 1;    
        max_cutoff_corrcoef = 4;    
        
        % Initialize displacementss
        plots_disp_f_prelim = get_disp(pixtounits_prelim,lenscoef_prelim);
                               
        % Set data  
        setappdata(handles_gui.figure,'min_lenscoef',min_lenscoef);
        setappdata(handles_gui.figure,'max_lenscoef',max_lenscoef);
        setappdata(handles_gui.figure,'handle_point_umax',[]);
        setappdata(handles_gui.figure,'handle_point_umin',[]);
        setappdata(handles_gui.figure,'handle_point_vmax',[]);
        setappdata(handles_gui.figure,'handle_point_vmin',[]);
        setappdata(handles_gui.figure,'val_checkbox_minmaxmarkers',false);
        setappdata(handles_gui.figure,'num_img',length(reference)-1);
        setappdata(handles_gui.figure,'handle_u_preview',[]);
        setappdata(handles_gui.figure,'handle_v_preview',[]);
        setappdata(handles_gui.figure,'pixtounits_prelim',pixtounits_prelim);
        setappdata(handles_gui.figure,'units_prelim',units_prelim);
        setappdata(handles_gui.figure,'corrcoef_buffer',corrcoef_buffer);
        setappdata(handles_gui.figure,'max_cutoff_corrcoef',max_cutoff_corrcoef);
        setappdata(handles_gui.figure,'lenscoef_prelim',lenscoef_prelim);
        setappdata(handles_gui.figure,'plots_disp_f_prelim',plots_disp_f_prelim);
        setappdata(handles_gui.figure,'handle_zoom',handle_zoom);
        setappdata(handles_gui.figure,'handle_pan',handle_pan);
        
        % Update formatted displacement plots
        update_axes('set');
        update_sidemenu();
                
        % Set visible
        set(handles_gui.figure,'Visible','on'); 
    end    

    function callback_slider_corrcutoff(hObject,eventdata) %#ok<INUSD>
        % Get data
        corrcoef_buffer = getappdata(handles_gui.figure,'corrcoef_buffer');
        num_img = getappdata(handles_gui.figure,'num_img');
        
        % Store in buffer
        corrcoef_buffer(num_img+1,1) = get(handles_gui.slider_corrcutoff,'value');
        
        % Set data        
        setappdata(handles_gui.figure,'corrcoef_buffer',corrcoef_buffer);
        
        % Update 
        update_axes('update');  
        update_sidemenu();  
    end

    function callback_checkbox_maxminmarkers(hObject,eventdata) %#ok<INUSD>        
        % Get value
        val_checkbox_minmaxmarkers = get(handles_gui.checkbox_maxminmarkers,'Value');
        
        % Set data
        setappdata(handles_gui.figure,'val_checkbox_minmaxmarkers',val_checkbox_minmaxmarkers);   
        
        % Update
        update_axes('update');
    end

    function callback_edit_pixtounits(hObject,eventdata) %#ok<INUSD>       
        % Get data
        pixtounits_prelim = getappdata(handles_gui.figure,'pixtounits_prelim');
        lenscoef_prelim = getappdata(handles_gui.figure,'lenscoef_prelim');
        plots_disp_f_prelim = getappdata(handles_gui.figure,'plots_disp_f_prelim');
        
        % Get Value - set min value to right above zero.
        pixtounits_buffer = str2double(get(handles_gui.edit_pixtounits,'string'));        
        if (ncorr_util_isrealbb(pixtounits_buffer,1e-14,inf,'Pixel to units conversion factor') == out.success) 
            % Store in buffer
            pixtounits_prelim = pixtounits_buffer;
                        
            % Get displacements
            plots_disp_f_prelim = get_disp(pixtounits_prelim,lenscoef_prelim);
        end

        % Store data
        setappdata(handles_gui.figure,'pixtounits_prelim',pixtounits_prelim);
        setappdata(handles_gui.figure,'plots_disp_f_prelim',plots_disp_f_prelim);
            
        % Update
        update_axes('set'); 
        update_sidemenu();
    end
        
    function callback_edit_units(hObject,eventdata) %#ok<INUSD>
        % Get units
        units_buffer = get(handles_gui.edit_units,'string');        
        
        % Set Data
        setappdata(handles_gui.figure,'units_prelim',units_buffer);
        
        % Update
        update_sidemenu();
    end

    function callback_edit_corrcutoff(hObject,eventdata) %#ok<INUSD>
        % Get data
        corrcoef_buffer = getappdata(handles_gui.figure,'corrcoef_buffer');
        max_cutoff_corrcoef = getappdata(handles_gui.figure,'max_cutoff_corrcoef');
        num_img = getappdata(handles_gui.figure,'num_img');       
        
        cutoff_corrcoef_prelim = str2double(get(handles_gui.edit_corrcutoff,'string')); 
        if (ncorr_util_isrealbb(cutoff_corrcoef_prelim,corrcoef_buffer(num_img+1,3),max_cutoff_corrcoef,'Correlation Coefficient') == out.success)
            if (abs(corrcoef_buffer(num_img+1,2)-corrcoef_buffer(num_img+1,3)) <= 1e-10)
                corrcoef_prelim = 1;
            else
                corrcoef_prelim = (cutoff_corrcoef_prelim-corrcoef_buffer(num_img+1,3))/(corrcoef_buffer(num_img+1,2)-corrcoef_buffer(num_img+1,3));
            end
            
            % Store in buffer
            corrcoef_buffer(num_img+1,1) = corrcoef_prelim;
        end
        
        % Store data
        setappdata(handles_gui.figure,'corrcoef_buffer',corrcoef_buffer);

        % Update
        update_axes('update'); 
        update_sidemenu();   
    end

    function callback_edit_lenscoef(hObject,eventdata) %#ok<INUSD>
        % Get data
        min_lenscoef = getappdata(handles_gui.figure,'min_lenscoef');
        max_lenscoef = getappdata(handles_gui.figure,'max_lenscoef');  
        pixtounits_prelim = getappdata(handles_gui.figure,'pixtounits_prelim');  
        lenscoef_prelim = getappdata(handles_gui.figure,'lenscoef_prelim');
        plots_disp_f_prelim = getappdata(handles_gui.figure,'plots_disp_f_prelim');    
            
        % Get Value
        lenscoef_buffer = str2double(get(handles_gui.edit_lenscoef,'string')); 
        if (ncorr_util_isrealbb(lenscoef_buffer,min_lenscoef,max_lenscoef,'Lens distortion coefficient') == out.success)
            % Store in buffer
            lenscoef_prelim = lenscoef_buffer;

            % Get displacements
            plots_disp_f_prelim = get_disp(pixtounits_prelim,lenscoef_prelim);
        end           

        % Store data
        setappdata(handles_gui.figure,'lenscoef_prelim',lenscoef_prelim);
        setappdata(handles_gui.figure,'plots_disp_f_prelim',plots_disp_f_prelim);
                
        % Update
        update_axes('set'); 
        update_sidemenu();
    end

    function callback_edit_imgnum(hObject,eventdata) %#ok<INUSD>
        % Get Data
        num_img = getappdata(handles_gui.figure,'num_img');

        % Get img num
        num_img_prelim = str2double(get(handles_gui.edit_imgnum,'string')); 
        if (ncorr_util_isintbb(num_img_prelim,1,length(reference),'Image number') == out.success)  
            % num_img_prelim will be zero based indexed
            num_img = num_img_prelim-1;
        end    
        
        % Set Data
        setappdata(handles_gui.figure,'num_img',num_img);
            
        % Update
        update_axes('set');
        update_sidemenu();
    end

    function callback_button_getunitconv(hObject,eventdata) %#ok<INUSD>
        % Form the units conversion GUI
        [pixtounits_prelim,units_prelim,outstate_unitconv] = ncorr_gui_getunitconv(get(handles_gui.figure,'OuterPosition'));
        
        if (outstate_unitconv == out.success)       
            % Get data
            lenscoef_prelim = getappdata(handles_gui.figure,'lenscoef_prelim');

            % Get displacements
            plots_disp_f_prelim = get_disp(pixtounits_prelim,lenscoef_prelim);
            
            % Set Data
            setappdata(handles_gui.figure,'pixtounits_prelim',pixtounits_prelim);
            setappdata(handles_gui.figure,'units_prelim',units_prelim);
            setappdata(handles_gui.figure,'plots_disp_f_prelim',plots_disp_f_prelim);
            
            % Update
            update_axes('set'); 
            update_sidemenu();
        end
    end

    function callback_button_applytoall(hObject,eventdata) %#ok<INUSD>
        % Get data
        corrcoef_buffer = getappdata(handles_gui.figure,'corrcoef_buffer');   
        num_img = getappdata(handles_gui.figure,'num_img');
        
        % Get value        
        cutoff_corrcoef_prelim = (corrcoef_buffer(num_img+1,2)-corrcoef_buffer(num_img+1,3))*corrcoef_buffer(num_img+1,1)+corrcoef_buffer(num_img+1,3);
        
        % Make sure this value is above the minimum for every current
        % configuration - i.e. take the max of the minimums
        if (cutoff_corrcoef_prelim > max(corrcoef_buffer(:,3)))
            % Apply value to all, must calculate normalized value
            for i = 0:length(reference)-1
                if (abs(corrcoef_buffer(i+1,2)-corrcoef_buffer(i+1,3)) <= 1e-10)
                    corrcoef_prelim = 1;
                else
                    corrcoef_prelim = (cutoff_corrcoef_prelim-corrcoef_buffer(i+1,3))/(corrcoef_buffer(i+1,2)-corrcoef_buffer(i+1,3));
                end
                % Store in buffer
                corrcoef_buffer(i+1,1) = corrcoef_prelim;
            end             
        else
            h_error = errordlg(['Value is below the minimum allowed: ' num2str(max(corrcoef_buffer(:,3)))  '.'],'Error','modal');
            uiwait(h_error);
        end  
        
        % Set data; there's no need for an update since the other plots arent shown
        setappdata(handles_gui.figure,'corrcoef_buffer',corrcoef_buffer);  
    end

    function callback_button_getlenscoef(hObject,eventdata) %#ok<INUSD>
        % Might add callback to determine lens coefficient in the future -
        % this requires DIC for rigid body translated images; 
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
        pixtounits_prelim = getappdata(handles_gui.figure,'pixtounits_prelim');  
        units_prelim = getappdata(handles_gui.figure,'units_prelim');
        corrcoef_buffer = getappdata(handles_gui.figure,'corrcoef_buffer');  
        lenscoef_prelim = getappdata(handles_gui.figure,'lenscoef_prelim');
        plots_disp_f_prelim = getappdata(handles_gui.figure,'plots_disp_f_prelim');
        
        % Ask user if outputs are correct        
        cutoff_corrcoef_prelim = zeros(length(reference),1);
        for i = 0:length(reference)-1
            cutoff_corrcoef_prelim(i+1) = (corrcoef_buffer(i+1,2)-corrcoef_buffer(i+1,3))*corrcoef_buffer(i+1,1)+corrcoef_buffer(i+1,3);
        end
        
        % Prompt user if inputs are correct
        string_disp = cell(0);
        string_disp{1} = ['Units/pixel is: ' num2str(pixtounits_prelim) '.'];
        string_disp{2} = ['Units are: ' units_prelim '.'];
        % Show a max of max_imgs since this can cause the message to grow too
        % large
        max_imgs = 20;
        for i = 0:min(max_imgs-1,length(reference)-1)
            string_disp{end+1} = ['Correlation coefficent cutoff for image ' num2str(i+1) ' is : ' num2str(cutoff_corrcoef_prelim(i+1)) '.']; %#ok<AGROW>
        end
        if (length(reference) > 20)
            string_disp{end+1} = '...';
        end
        string_disp{end+1} = ['Lens distortion coefficient is : ' num2str(lenscoef_prelim) '.'];
        string_disp{end+1} = 'Is this correct?';
        contbutton = questdlg(string_disp,'Continue Operation','Yes','No','No');
        if (strcmp(contbutton,'Yes'))            
            % Set up waitbar
            h = waitbar(0,['Processing ROI 1 of ' num2str(length(reference)) '...'],'Name','Calculating...','WindowStyle','modal');
            set(h,'CloseRequestFcn','setappdata(gcbf,''canceling'',1)');
            setappdata(h,'canceling',0)
            
            % Assign outputs
            for i = 0:length(reference)-1                                
                % Get final mask
                mask_cc = roi(i+1).mask & plots_corrcoef_ref{i+1} <= cutoff_corrcoef_prelim(i+1);
                
                % Form the final plots
                plots_disp_f(i+1).plot_u_formatted = zeros(size(plots_disp_f_prelim(i+1).plot_u_formatted));
                plots_disp_f(i+1).plot_u_formatted(mask_cc) = plots_disp_f_prelim(i+1).plot_u_formatted(mask_cc);
                plots_disp_f(i+1).plot_v_formatted = zeros(size(plots_disp_f_prelim(i+1).plot_v_formatted));
                plots_disp_f(i+1).plot_v_formatted(mask_cc) = plots_disp_f_prelim(i+1).plot_v_formatted(mask_cc);
                
                % Get ROI - union will preserve region correspondences but
                % new regions may or may not be contiguous.
                rois_f(i+1) = roi(i+1).get_union(mask_cc,0);
                
                % See if analysis was cancelled by user
                if (getappdata(h,'canceling'))
                    delete(h);                    
                    % Exit
                    return;
                end
                % Update waitbar
                waitbar((i+1)/length(reference),h,['Processing ROI ' num2str(i+1) ' of ' num2str(length(reference)) '...']);
            end
            
            % Close wait bar     
            delete(h);
            
            % Set dispinfo parameters
            pixtounits = pixtounits_prelim;
            units = units_prelim;
            cutoff_corrcoef = cutoff_corrcoef_prelim;
            lenscoef = lenscoef_prelim;
                        
            % Set Success
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
        num_img = getappdata(handles_gui.figure,'num_img');
        
        % Check for overshoot
        if (num_img > 0)
            num_img = num_img-1;
        end        
        
        % Set Data
        setappdata(handles_gui.figure,'num_img', num_img);   
            
        % Update
        update_axes('set');
        update_sidemenu();
    end

    function callback_button_right(hObject,eventdata) %#ok<INUSD>
        % Get Data
        num_img = getappdata(handles_gui.figure,'num_img');
        
        % Check for overshoot
        if (num_img < length(reference)-1)
            num_img = num_img+1;
        end                 
        
        % Set Data
        setappdata(handles_gui.figure,'num_img', num_img);  
            
        % Update
        update_axes('set');
        update_sidemenu();
    end

    function update_sidemenu()
        % Get Data
        val_checkbox_minmaxmarkers = getappdata(handles_gui.figure,'val_checkbox_minmaxmarkers');
        num_img = getappdata(handles_gui.figure,'num_img');
        pixtounits_prelim = getappdata(handles_gui.figure,'pixtounits_prelim');
        units_prelim = getappdata(handles_gui.figure,'units_prelim');
        corrcoef_buffer = getappdata(handles_gui.figure,'corrcoef_buffer');
        lenscoef_prelim = getappdata(handles_gui.figure,'lenscoef_prelim');
        handle_zoom = getappdata(handles_gui.figure,'handle_zoom'); 
        handle_pan = getappdata(handles_gui.figure,'handle_pan');  
        
        % Set slider
        set(handles_gui.slider_corrcutoff,'value',min(1,max(0,corrcoef_buffer(num_img+1,1))));
        
        % Set checkbox
        set(handles_gui.checkbox_maxminmarkers,'Value',val_checkbox_minmaxmarkers);
        
        % Update edit texts
        set(handles_gui.edit_pixtounits,'String',num2str(pixtounits_prelim));
        set(handles_gui.edit_units,'String',units_prelim);
        set(handles_gui.edit_corrcutoff,'String',num2str((corrcoef_buffer(num_img+1,2)-corrcoef_buffer(num_img+1,3))*corrcoef_buffer(num_img+1,1)+corrcoef_buffer(num_img+1,3),'%6.4f'));
        set(handles_gui.edit_lenscoef,'String',num2str(lenscoef_prelim));
        
        % Update buttons
        if (length(reference) > 1)
            set(handles_gui.button_applytoall,'enable','on');
        end
        
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
        % Get Data
        handle_point_umax = getappdata(handles_gui.figure,'handle_point_umax');
        handle_point_umin = getappdata(handles_gui.figure,'handle_point_umin');
        handle_point_vmax = getappdata(handles_gui.figure,'handle_point_vmax');
        handle_point_vmin = getappdata(handles_gui.figure,'handle_point_vmin');
        val_checkbox_minmaxmarkers = getappdata(handles_gui.figure,'val_checkbox_minmaxmarkers');
        num_img = getappdata(handles_gui.figure,'num_img');
        handle_u_preview = getappdata(handles_gui.figure,'handle_u_preview');
        handle_v_preview = getappdata(handles_gui.figure,'handle_v_preview');
        corrcoef_buffer = getappdata(handles_gui.figure,'corrcoef_buffer');
        plots_disp_f_prelim = getappdata(handles_gui.figure,'plots_disp_f_prelim');
               
        if (strcmp(action,'set'))    
            % Get reduced reference image
            ref_reduced = reference(num_img+1).reduce(spacing);
            
            % Get displacement info     
            uplot = plots_disp_f_prelim(num_img+1).plot_u_formatted;
            vplot = plots_disp_f_prelim(num_img+1).plot_v_formatted;
          
            % U Plot:  
            max_u = max(uplot(roi(num_img+1).mask));
            min_u = min(uplot(roi(num_img+1).mask));
            imshow(ref_reduced.get_img(),[ref_reduced.min_gs ref_reduced.max_gs],'Parent',handles_gui.axes_formatplot_u);  
            hold(handles_gui.axes_formatplot_u,'on');
            handle_u_preview = imshow(uplot,[min_u max_u],'Parent',handles_gui.axes_formatplot_u);

            % Place markers on max and min locations - these positions get
            % overwritten in the update, so dont worry about making sure
            % the positions are within the ROI.
            [y_umax,x_umax] = find(uplot == max_u,1);
            [y_umin,x_umin] = find(uplot == min_u,1);

            handle_point_umax = impoint(handles_gui.axes_formatplot_u,x_umax,y_umax);
            setColor(handle_point_umax,'g');
            set(handle_point_umax,'UIContextMenu','');
            set(handle_point_umax,'ButtonDownFcn','');

            handle_point_umin = impoint(handles_gui.axes_formatplot_u,x_umin,y_umin);
            setColor(handle_point_umin,'g');
            set(handle_point_umin,'UIContextMenu','');
            set(handle_point_umin,'ButtonDownFcn','');

            % Turn hold off
            hold(handles_gui.axes_formatplot_u,'off');
            set(handles_gui.axes_formatplot_u,'Visible','off');
                             
            % V-Plot:
            max_v = max(vplot(roi(num_img+1).mask));
            min_v = min(vplot(roi(num_img+1).mask));
            imshow(ref_reduced.get_img(),[ref_reduced.min_gs ref_reduced.max_gs],'Parent',handles_gui.axes_formatplot_v);
            hold(handles_gui.axes_formatplot_v,'on');
            handle_v_preview = imshow(vplot,[min_v max_v],'Parent',handles_gui.axes_formatplot_v);

            % Place markers on max and min locations - these positions get
            % overwritten in the update, so dont worry about making sure
            % the positions are within the ROI.
            [y_vmax,x_vmax] = find(vplot == max_v,1);
            [y_vmin,x_vmin] = find(vplot == min_v,1); 

            handle_point_vmax = impoint(handles_gui.axes_formatplot_v,x_vmax,y_vmax);        
            setColor(handle_point_vmax,'g');
            set(handle_point_vmax,'UIContextMenu','');
            set(handle_point_vmax,'ButtonDownFcn','');

            handle_point_vmin = impoint(handles_gui.axes_formatplot_v,x_vmin,y_vmin);
            setColor(handle_point_vmin,'g');
            set(handle_point_vmin,'UIContextMenu','');
            set(handle_point_vmin,'ButtonDownFcn','');

            % Turn hold off
            hold(handles_gui.axes_formatplot_v,'off');
            set(handles_gui.axes_formatplot_v,'Visible','off');         

            % Set colormap
            ncorr_util_colormap(handles_gui.figure);
            
            % Update Texts
            set(handles_gui.text_cur,'String', ['Current Name: ' current(num_img+1).name(1:min(end,22))]);
            set(handles_gui.text_ref,'String', ['Reference Name: ' reference(num_img+1).name(1:min(end,22))]);

            % Update Buttons\Edit
            set(handles_gui.edit_imgnum,'String',num2str(num_img+1));
            if (length(reference) == 1)
                set(handles_gui.button_left,'Enable','off');
                set(handles_gui.button_right,'Enable','off');
                set(handles_gui.edit_imgnum,'Enable','off');
            elseif (num_img == 0)
                set(handles_gui.button_left,'Enable','off');
                set(handles_gui.button_right,'Enable','on');
                set(handles_gui.edit_imgnum,'Enable','on');
            elseif (num_img == length(reference)-1)
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
            % Get displacement info     
            uplot = plots_disp_f_prelim(num_img+1).plot_u_formatted;
            vplot = plots_disp_f_prelim(num_img+1).plot_v_formatted;
            
            % Calculate cutoff_corrcoef_prelim
            cutoff_corrcoef_prelim = (corrcoef_buffer(num_img+1,2)-corrcoef_buffer(num_img+1,3))*corrcoef_buffer(num_img+1,1)+corrcoef_buffer(num_img+1,3);

            % Get correlation coefficient unioned ROI
            mask_cc = roi(num_img+1).mask & plots_corrcoef_ref{num_img+1} <= cutoff_corrcoef_prelim;               

            % U-Plot:
            u_95 = prctile(uplot(mask_cc),95);
            u_5 = prctile(uplot(mask_cc),5);
            max_ulim = u_95+1e-5;
            min_ulim = u_5-1e-5;
            if (abs(min_ulim-max_ulim)<1e-4)
                max_ulim = max_ulim + 1e-5; % add a small increment to insure they are not the same number
                min_ulim = min_ulim - 1e-5;
            end
            % Update
            alphamap_u = mask_cc.*0.7;
            set(handle_u_preview,'AlphaData',alphamap_u);
            caxis(handles_gui.axes_formatplot_u,[min_ulim max_ulim]);
            
            % Set colorbar
            handle_colorbar_u = colorbar('peer',handles_gui.axes_formatplot_u);        
            set(handle_colorbar_u,'UIContextMenu','');
            set(get(handle_colorbar_u,'child'),'YData',[min_ulim max_ulim]);
            set(handle_colorbar_u,'YLim',[min_ulim max_ulim]);
            set(handle_colorbar_u,'Units','Pixels');

            % Update uplot points
            if (val_checkbox_minmaxmarkers)                
                % Get coordinates
                max_u = max(uplot(mask_cc));
                min_u = min(uplot(mask_cc));
                % Must use mask_cc to make sure coordinations are within
                % the ROI
                [y_umax,x_umax] = find(uplot == max_u & mask_cc,1);
                [y_umin,x_umin] = find(uplot == min_u & mask_cc,1);   
                
                % Set points visible
                set(handle_point_umax,'Visible','on');
                set(handle_point_umin,'Visible','on');
                
                % Set points position
                setPosition(handle_point_umax,[x_umax y_umax]);
                setPosition(handle_point_umin,[x_umin y_umin]);
            else
                set(handle_point_umax,'Visible','off');
                set(handle_point_umin,'Visible','off');  
            end
            
            % V-Plot:
            v_95 = prctile(vplot(mask_cc),95);
            v_5 = prctile(vplot(mask_cc),5);
            max_vlim = v_95+1e-5;
            min_vlim = v_5-1e-5;
            if (abs(min_vlim-max_vlim)<1e-4)
                max_vlim = max_vlim + 1e-5; % add a small increment to insure they are not the same number
                min_vlim = min_vlim - 1e-5;
            end
            
            % Update
            alphamap_v = mask_cc.*0.7;
            set(handle_v_preview,'AlphaData', alphamap_v);
            caxis(handles_gui.axes_formatplot_v,[min_vlim max_vlim]);

            % Set colorbar
            handle_colorbar_v = colorbar('peer',handles_gui.axes_formatplot_v);        
            set(handle_colorbar_v,'UIContextMenu','');
            set(get(handle_colorbar_v,'child'),'YData',[min_vlim max_vlim]);
            set(handle_colorbar_v,'YLim',[min_vlim max_vlim]);
            set(handle_colorbar_v,'Units','Pixels');
            
            % Update vplot points    
            if (val_checkbox_minmaxmarkers)
                % Get coordinates
                max_v = max(vplot(mask_cc));
                min_v = min(vplot(mask_cc));
                % Must use mask_cc to make sure coordinations are within
                % the ROI
                [y_vmax,x_vmax] = find(vplot == max_v & mask_cc,1);
                [y_vmin,x_vmin] = find(vplot == min_v & mask_cc,1);  
                
                % Set points visible
                set(handle_point_vmax,'Visible','on');
                set(handle_point_vmin,'Visible','on');
                
                % set points position
                setPosition(handle_point_vmax,[x_vmax y_vmax]);
                setPosition(handle_point_vmin,[x_vmin y_vmin]);
            else
                set(handle_point_vmax,'Visible','off');
                set(handle_point_vmin,'Visible','off');
            end                 
        end
        
        % Set data
        setappdata(handles_gui.figure,'handle_point_umax',handle_point_umax);
        setappdata(handles_gui.figure,'handle_point_umin',handle_point_umin);
        setappdata(handles_gui.figure,'handle_point_vmax',handle_point_vmax);
        setappdata(handles_gui.figure,'handle_point_vmin',handle_point_vmin);
        setappdata(handles_gui.figure,'handle_u_preview', handle_u_preview);
        setappdata(handles_gui.figure,'handle_v_preview', handle_v_preview);
    end

    function plots_disp_f_buf = get_disp(pixtounits_buf,lenscoef_buf)
        % This function performs corrections to all displacements given a lens
        % distortion coefficient and units/pixel conversion. This assumes
        % the distortion is radial and is centered WRT the center of the
        % reference image.
                
        % Copy input displacements first so they don't get overwritten
        plots_disp_f_buf = struct('plot_u_formatted',{},'plot_v_formatted',{});
        for i = 0:length(reference)-1
            plots_disp_f_buf(i+1).plot_u_formatted = plots_u_ref{i+1};
            plots_disp_f_buf(i+1).plot_v_formatted = plots_v_ref{i+1};
        end
        
        % Update formatted displacement plots based on lens distortion
        % coefficient. Do this WRT pixels.
        for i = 0:length(reference)-1
            % Correct plots for lens distortion. Assume it's radial WRT
            % the center of the reference image. Get coordinates of points
            % in pixels:
            [xcoords,ycoords] = meshgrid((0:spacing+1:reference(i+1).width-1)-reference(i+1).width/2, ...
                                         (0:spacing+1:reference(i+1).height-1)-reference(i+1).height/2);

            vec_x = xcoords(roi(i+1).mask);
            vec_y = ycoords(roi(i+1).mask);
            vec_x_tilda = vec_x+plots_disp_f_buf(i+1).plot_u_formatted(roi(i+1).mask); 
            vec_y_tilda = vec_y+plots_disp_f_buf(i+1).plot_v_formatted(roi(i+1).mask);

            % Eq.5 in Systematic errors in two-dimensional digital image 
            % correlation due to lens distortion. Also from fig.2.
            plots_disp_f_buf(i+1).plot_u_formatted(roi(i+1).mask) = plots_disp_f_buf(i+1).plot_u_formatted(roi(i+1).mask)+(lenscoef_buf*(vec_x_tilda.*(vec_x_tilda.^2+vec_y_tilda.^2)-vec_x.*(vec_x.^2+vec_y.^2)));
            plots_disp_f_buf(i+1).plot_v_formatted(roi(i+1).mask) = plots_disp_f_buf(i+1).plot_v_formatted(roi(i+1).mask)+(lenscoef_buf*(vec_y_tilda.*(vec_x_tilda.^2+vec_y_tilda.^2)-vec_y.*(vec_x.^2+vec_y.^2)));   
        end
            
        % Convert from pixels to real units
        for i = 0:length(reference)-1
            plots_disp_f_buf(i+1).plot_u_formatted = plots_disp_f_buf(i+1).plot_u_formatted.*pixtounits_buf;
            plots_disp_f_buf(i+1).plot_v_formatted = plots_disp_f_buf(i+1).plot_v_formatted.*pixtounits_buf;
        end
    end
    
    function handles_gui = init_gui()
    % GUI controls -------------------------------------------------------%
        % Figure
        handles_gui.figure = figure( ...
            'Tag', 'figure', ...
            'Units', 'characters', ...
            'Position', ncorr_util_figpos(pos_parent,[41 214.4]), ...
            'Name', 'Format Displacements', ...
            'MenuBar', 'none', ...
            'NumberTitle', 'off', ...
            'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
            'handlevisibility','off', ...
            'DockControls','off', ...
            'WindowStyle','modal', ...
            'Resize','off', ...
            'Visible','off', ...
            'IntegerHandle','off', ...
            'Interruptible','off');
            
        % Panels
        handles_gui.group_units = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_units', ...
            'Units', 'characters', ...
            'Position', [2 32.1 35 8.3], ...
            'Title', 'Units Options', ...
            'Interruptible','off');
        
        handles_gui.group_formatdispoptions = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_formatdispoptions', ...
            'Units', 'characters', ...
            'Position', [2 20.8 35 10.6], ...
            'Title', 'Formatting Options', ...
            'Interruptible','off');
        
        handles_gui.group_lenscoef = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_lenscoef', ...
            'Units', 'characters', ...
            'Position', [2 13.4 35 6.7], ...
            'Title', 'Lens Distortion Options', ...
            'Interruptible','off');

        handles_gui.group_zoompan = uibuttongroup( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_zoompan', ...
            'Units', 'characters', ...
            'Position', [2 8.1 35 4.5], ...
            'Title', 'Zoom/Pan', ...
            'Interruptible','off');
        
        handles_gui.group_menu = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_menu', ...
            'Units', 'characters', ...
            'Position', [2 0.75 35 6.7], ...
            'Title', 'Menu', ...
            'Interruptible','off');

        handles_gui.group_formatdisp = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_formatplot_u', ...
            'Units', 'characters', ...
            'Position', [38.7 0.75 173.8 39.6], ...
            'Title', 'Preview', ...
            'Interruptible','off');

        % Axes
        handles_gui.axes_formatplot_u = axes( ...
            'Parent', handles_gui.group_formatdisp, ...
            'Tag', 'axes_formatplot_u', ...
            'Units', 'characters', ...
            'Position', [2.4 4.1 75.8 31.5], ...
            'Visible','off', ...
            'Interruptible','off');

        handles_gui.axes_formatplot_v = axes( ...
            'Parent', handles_gui.group_formatdisp, ...
            'Tag', 'axes_formatplot_v', ...
            'Units', 'characters', ...
            'Position', [88.2 4.1 75.8 31.5], ...
            'Visible','off', ...
            'Interruptible','off');

        % Static Texts
        handles_gui.text_pixtounits = uicontrol( ...
            'Parent', handles_gui.group_units, ...
            'Tag', 'text_pixtounits', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 2.8 18 1.3], ...
            'String', 'Units/Pixel: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
        
        handles_gui.text_units = uicontrol( ...
            'Parent', handles_gui.group_units, ...
            'Tag', 'text_units', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 0.7 18 1.3], ...
            'String', 'Units:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
        
        handles_gui.text_corrcutoff = uicontrol( ...
            'Parent', handles_gui.group_formatdispoptions, ...
            'Tag', 'text_corrcutoff', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 4.9 18 1.3], ...
            'String', 'Corr-Coef Cutoff: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
        
        handles_gui.text_corrcutoff = uicontrol( ...
            'Parent', handles_gui.group_formatdispoptions, ...
            'Tag', 'text_corrcutoff', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 4.9 18 1.3], ...
            'String', 'Corr-Coef Cutoff: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_lenscoef = uicontrol( ...
            'Parent', handles_gui.group_lenscoef, ...
            'Tag', 'text_lenscoef', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.5 3.2 18 1.3], ...
            'String', 'Lens Coef: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
        
        handles_gui.text_ref = uicontrol( ...
            'Parent', handles_gui.group_formatdisp, ...
            'Tag', 'text_ref', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [3 1.8 56.9 1.5], ...
            'String', 'Reference Name: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
        
        handles_gui.text_cur = uicontrol( ...
            'Parent', handles_gui.group_formatdisp, ...
            'Tag', 'text_cur', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [3 0.6 56.9 1.5], ...
            'String', 'Current Name: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_plot_u = uicontrol( ...
            'Parent', handles_gui.group_formatdisp, ...
            'Tag', 'text_plot_u', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.4 36 29.4 1.5], ...
            'String', 'U-displacement:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_plot_v = uicontrol( ...
            'Parent', handles_gui.group_formatdisp, ...
            'Tag', 'text_plot_v', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [88.2 36 29.4 1.5], ...
            'String', 'V-displacement:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        % Sliders
        handles_gui.slider_corrcutoff = uicontrol( ...
            'Parent', handles_gui.group_formatdispoptions, ...
            'Tag', 'slider_corrcutoff', ...
            'Style', 'slider', ...
            'Units', 'characters', ...
            'Position', [2.5 3.5 28.9 1.3], ...
            'BackgroundColor', [0.9 0.9 0.9], ...
            'String', {'Slider'}, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_slider_corrcutoff,handles_gui.figure), ...
            'Interruptible','off');
                
        % Checkbox
        handles_gui.checkbox_maxminmarkers = uicontrol( ...
            'Parent', handles_gui.group_formatdispoptions, ...
            'Tag', 'checkbox_maxminmarkers', ...
            'Style', 'checkbox', ...
            'Units', 'characters', ...
            'Position', [2.5 7.2 29.7 1.3], ...
            'String', 'Max/min markers ', ...
            'Value', false, ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_checkbox_maxminmarkers,handles_gui.figure), ...
            'Interruptible','off');
        
        % Edit
        handles_gui.edit_pixtounits = uicontrol( ...
            'Parent', handles_gui.group_units, ...
            'Tag', 'edit_pixtounits', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [21.2 2.9 10.1 1.3], ...
            'String', '0', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_pixtounits,handles_gui.figure), ...
            'Interruptible','off');  
        
        handles_gui.edit_units = uicontrol( ...
            'Parent', handles_gui.group_units, ...
            'Tag', 'edit_units', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [21.2 0.9 10.1 1.3], ...
            'String', '0', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_units,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.edit_corrcutoff = uicontrol( ...
            'Parent', handles_gui.group_formatdispoptions, ...
            'Tag', 'edit_corrcutoff', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [21.2 5.1 10.1 1.3], ...
            'String', '0', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_corrcutoff,handles_gui.figure), ...
            'Interruptible','off');        
        
        handles_gui.edit_lenscoef = uicontrol( ...
            'Parent', handles_gui.group_lenscoef, ...
            'Tag', 'edit_lenscoef', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [21.2 3.4 10.1 1.3], ...
            'String', '0', ...
            'HorizontalAlignment', 'left', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_lenscoef,handles_gui.figure), ...
            'Interruptible','off');        
        
        handles_gui.edit_imgnum = uicontrol( ...
            'Parent', handles_gui.group_formatdisp, ...
            'Tag', 'edit_imgnum', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [156 1.1 7 1.6], ...
            'String', '', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_edit_imgnum,handles_gui.figure), ...
            'Enable', 'off', ...
            'Interruptible','off');
        
        % Pushbuttons
        handles_gui.button_getunitconv = uicontrol( ...
            'Parent', handles_gui.group_units, ...
            'Tag', 'button_getunitconv', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 4.8 29.7 1.7], ...
            'String', 'Get Unit Conversion', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_getunitconv,handles_gui.figure), ...
            'Interruptible','off', ...
            'Enable', 'on');
        
        handles_gui.button_applytoall = uicontrol( ...
            'Parent', handles_gui.group_formatdispoptions, ...
            'Tag', 'button_applytoall', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 1 29.7 1.7], ...
            'String', 'Apply To All', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_applytoall,handles_gui.figure), ...
            'Interruptible','off', ...
            'Enable', 'off');
        
        handles_gui.button_getlenscoef = uicontrol( ...
            'Parent', handles_gui.group_lenscoef, ...
            'Tag', 'button_lenscoef', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 1 29.7 1.7], ...
            'String', 'Get Lens Coef', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_getlenscoef,handles_gui.figure), ...
            'Interruptible','off', ...
            'Enable','off');        

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
            'Parent', handles_gui.group_formatdisp, ...
            'Tag', 'button_left', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [149 1 6 1.8], ...
            'String', '<', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_left,handles_gui.figure), ...
            'Interruptible','off', ...
            'Enable', 'off');

        handles_gui.button_right = uicontrol( ...
            'Parent', handles_gui.group_formatdisp, ...
            'Tag', 'button_right', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [164 1 6 1.8], ...
            'String', '>', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_right,handles_gui.figure), ...
            'Interruptible','off', ...
            'Enable', 'off');
    end

    % Pause until figure is closed ---------------------------------------%
    waitfor(handles_gui.figure);
end

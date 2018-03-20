function [plots_disp_new,rois_new,outstate] = ncorr_alg_convertanalysis(img_new,imgs_old,plots_u_old,plots_v_old,rois_old,spacing,num_img,total_imgs)
% This function is used to convert displacements from the old configuration
% to the new configuration. Can be Lagrangian to Eulerian or vice-versa, 
% which is why "old" and "new" are used instead of Lagrangian and Eulerian.
%
% Inputs -----------------------------------------------------------------%
%   img_new - struct; Contains struct('imginfo',{},'roi',{}) in the new 
%   configuration. 
%   imgs_old - struct; Contains struct('imginfo',{},'roi',{}) in the old 
%   configuration. 
%   plots_u_old - cell; contains "old" u displacement fields
%   plots_v_old - cell; contains "old" v displacement fields
%   rois_old - ncorr_class_roi; contains ROIs corresponding to the 
%   u and v plots. Note these are reduced by default.
%   spacing - integer; spacing parameter
%   num_img - integer; reference image number
%   total_imgs - integer; total number of images
%
% Outputs ----------------------------------------------------------------%
%   plots_disp_new - struct; contains struct('plot_u_new',{},'plot_v_new',{})
%   rois_new - struct; contains ROIs corresponding to the "new" u and v 
%   displacement plots
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.
%
% Returns failed if seed placement fails.

    % Initialize outputs
    outstate = out.cancelled;
    plots_disp_new = struct('plot_u_new',{},'plot_v_new',{});
    rois_new = ncorr_class_roi.empty;
    
    % Extrapolate displacement plots to improve interpolation near the
    % boundary points.
    % border_interp is the border added to the displacements when they are
    % extrapolated.
    border_interp = 20; % MUST BE GREATER THAN OR EQUAL TO 2
    plots_u_interp_old = cell(0);
    plots_v_interp_old = cell(0);
    % Note that ncorr_alg_extrapdata will return a separate extrapolated 
    % array for each region within a ROI. This is done to prevent 
    % displacements from adjacent regions from influencing each other.
    % Also note that the ordering returned forms a correspondence between
    % the regions stored in rois_old.
    for i = 0:length(imgs_old)-1
        plots_u_interp_old{i+1} = ncorr_alg_extrapdata(plots_u_old{i+1},rois_old(i+1).formatted(),int32(border_interp));
        plots_v_interp_old{i+1} = ncorr_alg_extrapdata(plots_v_old{i+1},rois_old(i+1).formatted(),int32(border_interp));
    end

    % Convert displacement plots to B-spline coefficients; these are used 
    % for interpolation. Make sure to do this for each region.
    for i = 0:length(imgs_old)-1
        for j = 0:length(rois_old(i+1).region)-1
            plots_u_interp_old{i+1}{j+1} = ncorr_class_img.form_bcoef(plots_u_interp_old{i+1}{j+1});
            plots_v_interp_old{i+1}{j+1} = ncorr_class_img.form_bcoef(plots_v_interp_old{i+1}{j+1});
        end
    end   

    % Set Seeds ----------------------------------------------------------%
    % This will automatically try to seed as many regions as possible. 
    roi_new_reduced = img_new.roi.reduce(spacing);
    convertseedinfo = cell(0);
    seedwindow = 1;
    for i = 0:length(imgs_old)-1
        % Will return either success or failed. Returns success if at
        % least one region is seeded.
        [convertseedinfo{i+1},outstate_convertseeds] = ncorr_alg_convertseeds(plots_u_old{i+1}, ...
                                                                              plots_v_old{i+1}, ...
                                                                              plots_u_interp_old{i+1}, ...
                                                                              plots_v_interp_old{i+1}, ...
                                                                              rois_old(i+1), ...
                                                                              roi_new_reduced, ...   
                                                                              seedwindow, ...                                                      
                                                                              spacing, ...
                                                                              border_interp);
                                                     
        % Must check to make sure seedinfo isnt empty
        if (outstate_convertseeds ~= out.success)
            % For now just fail the whole analysis since this will result 
            % in an empty plot.
            h_error = errordlg('Seeding failed for the conversion analysis. Please rerun DIC analysis and make sure regions have a large contiguous region.','Error','modal');
            uiwait(h_error);
            
            outstate = out.failed;
            break;
        end
    end
    
    if (outstate ~= out.failed) 
        % Find displacements in new configuration ------------------------%
        % Format seeds    
        convertseedinfo_f = convertseedinfo;
        for i = 0:length(convertseedinfo)-1
            for j = 0:length(convertseedinfo{i+1})-1
                convertseedinfo_f{i+1}(j+1).num_region_new = int32(convertseedinfo{i+1}(j+1).num_region_new);
                convertseedinfo_f{i+1}(j+1).num_region_old = int32(convertseedinfo{i+1}(j+1).num_region_old);
            end
        end  

        % Initialize structs for output 
        plots_disp_new_prelim = struct('plot_u_new',{},'plot_v_new',{});
        rois_new_prelim = ncorr_class_roi.empty;
        for i = 0:length(imgs_old)-1
            % ncorr_alg_convert returns either success or cancelled.
            [plot_disp_new_buffer,outstate_convert] = ncorr_alg_convert(plots_u_interp_old{i+1}, ...
                                                                        plots_v_interp_old{i+1}, ...
                                                                        rois_old(i+1).formatted(), ...
                                                                        roi_new_reduced.formatted(), ...
                                                                        convertseedinfo_f{i+1}, ...
                                                                        int32(spacing), ...
                                                                        int32(border_interp), ...
                                                                        int32(num_img+i), ...
                                                                        int32(total_imgs));
            % Check if analysis was cancelled
            if (outstate_convert ~= out.success)
                return;
            end

            % Store new displacement plots and then get the union of valid
            % points with roi_new_reduced.
            plots_disp_new_prelim(i+1).plot_u_new = plot_disp_new_buffer.plot_u_new;
            plots_disp_new_prelim(i+1).plot_v_new = plot_disp_new_buffer.plot_v_new;
            rois_new_prelim(i+1) = roi_new_reduced.get_union(plot_disp_new_buffer.plot_validpoints,0);
        end  

        % Set outputs
        for i = 0:length(plots_disp_new_prelim)-1
            plots_disp_new(i+1) = plots_disp_new_prelim(i+1);
            rois_new(i+1) = rois_new_prelim(i+1);
        end    
        outstate = out.success;
    end
end

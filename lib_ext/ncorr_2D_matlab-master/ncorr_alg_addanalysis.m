function [plot_added,outstate] = ncorr_alg_addanalysis(plots_u,plots_v,rois,spacing,num_img,total_imgs)
% This function performs the adding of displacement fields. It will cycle
% through every point in the mask of the first ROI and attempt to propogate 
% those points until the last displacement plot, so the added displacements 
% will be WRT to the first plot's perspective.
%
% Inputs -----------------------------------------------------------------%
%   plots_u - cell; u displacement plots
%   plots_v - cell; v displacement plots
%   rois - ncorr_class_roi; ROIs corresponding to the displacement plots. 
%   Note that rois will already be reduced by default.
%   spacing - integer; spacing parameter
%   num_img - integer; number of reference image being analyzed
%   total_imgs - integer; total number of images being analyzed
%
% Outputs ----------------------------------------------------------------%
%   plot_added - struct; contains
%   struct('plot_u_added',{},'plot_v_added',{},'plot_validpoints',{})
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.

    % Initialize outputs
    outstate = out.cancelled;    
    plot_added = struct('plot_u_added',{},'plot_v_added',{},'plot_validpoints',{});

    % Extrapolate displacement plots to improve interpolation near the
    % boundary points.
    % border_interp is the border added to the displacements when they are
    % extrapolated.
    border_interp = 20; % MUST BE GREATER THAN OR EQUAL TO 2
    plots_u_interp = cell(0);
    plots_v_interp = cell(0);
    % Note that ncorr_alg_extrapdata will return a separate extrapolated 
    % array for each region within a ROI. This is done to prevent 
    % displacements from adjacent regions from influencing each other.
    % Also note that the ordering returned forms a correspondence between
    % the regions stored in rois.
    for k = 0:length(plots_u)-1
        plots_u_interp{k+1} = ncorr_alg_extrapdata(plots_u{k+1},rois(k+1).formatted(),int32(border_interp));
        plots_v_interp{k+1} = ncorr_alg_extrapdata(plots_v{k+1},rois(k+1).formatted(),int32(border_interp));
    end

    % Convert displacement plots to B-spline coefficients; these are used 
    % for interpolation. Make sure to do this for each region.
    for k = 0:length(plots_u_interp)-1
        for l = 0:length(plots_u_interp{k+1})-1
            plots_u_interp{k+1}{l+1} = ncorr_class_img.form_bcoef(plots_u_interp{k+1}{l+1});
            plots_v_interp{k+1}{l+1} = ncorr_class_img.form_bcoef(plots_v_interp{k+1}{l+1});
        end
    end

    % Format ROIs                
    rois_formatted = ncorr_class_roi.empty;
    for k = 0:length(rois)-1
        rois_formatted(k+1) = rois(k+1).formatted();
    end

    % Add plots to get the overall displacement field. ncorr_alg_adddisp
    % returns either success or cancelled.
    plot_added_prelim = struct('plot_u_added',{},'plot_v_added',{},'plot_validpoints',{});
    [plot_added_prelim(1),outstate_add] = ncorr_alg_adddisp(plots_u_interp, ...
                                                            plots_v_interp, ...
                                                            rois_formatted, ...
                                                            int32(border_interp), ...
                                                            int32(spacing), ...
                                                            int32(num_img-1), ...
                                                            int32(total_imgs));
                                                        
    if (outstate_add == out.success)
        % Set outputs
        plot_added(1) = plot_added_prelim;
        outstate = out.success;
    end
end

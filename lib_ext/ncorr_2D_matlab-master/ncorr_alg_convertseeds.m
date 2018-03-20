function [convertseedinfo,outstate] = ncorr_alg_convertseeds(plot_u_old,plot_v_old,plot_u_interp_old,plot_v_interp_old,roi_old,roi_new,seedwindow,spacing,border_interp)
% This function obtains the "convert" seeds. It will try to seed as many
% "new" regions as possible (one seed per region). It is possible that no
% seeds are returned.
%
% Inputs -----------------------------------------------------------------%
%   plot_u_old - double array; u displacements WRT the "old" configuration. 
%   Units are pixels.
%   plot_v_old - double array; v displacements WRT the "old" configuration.
%   Units are pixels. 
%   plot_u_interp_old - cell; array of b-spline coefficients; one per
%   region.
%   plot_v_interp_old - cell; array of b-spline coefficients; one per
%   region;
%   roi_old - ncorr_class_roi; ROI corresponding to "old" displacements.
%   Note that ROI is already reduced by default.
%   roi_new - ncorr_class_roi; ROI corresponding to "new" displacements.
%   Note that ROI is already reduced by default.
%   seedwindow - integer; half width of window around seed that must contain
%   valid points before its processed. This prevents edge points from being
%   seeded.
%   spacing - integer; this is the spacing parameter.
%   border_interp - integer; the amount of padding used around the borders in
%   interpdata
%
% Outputs ----------------------------------------------------------------%
%   convertseedinfo - struct; contains
%   struct('paramvector',{},'num_region_new',{},'num_region_old',{}). The x
%   and y coords stored in convertseedinfo will be reduced while the
%   displacement values are in pixels.
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.
% 
% Note that the ordering of plot_u_interp_old and plot_v_interp_old must
% form a correspondence with the regions in roi_old. This means
% plot_u_interp_old{i} must correspond to roi_old.region(i). Regions in
% roi_old and roi_new do not need to form a direct correspondence, but they
% are assumed to be one-to-one. Returns failed if no seeds are found

    % Initialize outputs
    outstate = out.failed;
    convertseedinfo = struct('paramvector',{},'num_region_new',{},'num_region_old',{}); % paramvector = [x_new y_new x_old y_old u_old v_old distance]
        
    % Form convertseedinfo prelim
    convertseedinfo_prelim = struct('paramvector',{},'num_region_new',{},'num_region_old',{});

    % Keep track of the regions in roi_old that have already been analyzed.
    list_region_old = false(length(roi_old.region),1); 
    
    % Cycle over every "new" region and attempt to seed.
    for i = 0:length(roi_new.region)-1
        % Regions are not gauranteed contiguous at this point, so get the 
        % region corresponding to the largest contiguous area to make sure 
        % small areas arent seeded.
        regionmask_new_buffer = roi_new.get_regionmask(i);
        [region_new_buffer,removed] = ncorr_alg_formregions(regionmask_new_buffer,int32(0),false); %#ok<ASGLU>

        % Check if contiguous region(s) are empty or if there are more than one
        if (isempty(region_new_buffer))
            % Continue onto next region in roi_new if region_new_buffer is empty
            continue;
        elseif (length(region_new_buffer) > 1)
            % Select biggest region if there are more than one. This could 
            % possibly happen if a boundary is "pinched" or closes during 
            % deformation. Unlikely- but may happen.
            idx_max = find([region_new_buffer.totalpoints] == max([region_new_buffer.totalpoints]),1,'first');
            region_new_buffer = region_new_buffer(idx_max);
        end
        
        % Form convertseedinfo buffer
        convertseedinfo_buffer = struct('paramvector',{},'num_region_new',{},'num_region_old',{});        
        % Set num_region_new
        convertseedinfo_buffer(1).num_region_new = i;
        % Initialize
        successregion = false;
        for j = 0:size(region_new_buffer.noderange,1)-1
            x_new = j + region_new_buffer.leftbound;
            % Cycle over each point
            for k = 0:2:region_new_buffer.noderange(j+1)-1
                for l = region_new_buffer.nodelist(j+1,k+1):region_new_buffer.nodelist(j+1,k+2)
                    y_new = l;                    
                    % Make sure area of half-width seedwindow is valid
                    % around the seed before attempting to process.
                    if (all(all(regionmask_new_buffer(max(y_new-seedwindow+1,1):min(y_new+seedwindow+1,end), ...
                                                      max(x_new-seedwindow+1,1):min(x_new+seedwindow+1,end)))))
                        % Analyze point - outstate_calcpoint will either be
                        % success or failed
                        [convertseedinfo_buffer.paramvector,convertseedinfo_buffer.num_region_old,outstate_calcpoint] = calcpoint(x_new, ...
                                                                                                                                  y_new, ...
                                                                                                                                  plot_u_old, ...
                                                                                                                                  plot_v_old, ...
                                                                                                                                  plot_u_interp_old, ...
                                                                                                                                  plot_v_interp_old, ...
                                                                                                                                  roi_old, ...
                                                                                                                                  list_region_old, ...
                                                                                                                                  spacing, ...
                                                                                                                                  border_interp); 
                        % Check if nonlinear solver was successful
                        if (outstate_calcpoint == out.success)
                            % Check distance; make sure its less than a
                            % threshold and also make sure that x_old and
                            % y_old are in valid area
                            x_old = round(convertseedinfo_buffer.paramvector(3));
                            y_old = round(convertseedinfo_buffer.paramvector(4));                            
                            if (convertseedinfo_buffer.paramvector(7) < 0.005 && ...
                                x_old >= 0 && x_old < size(roi_old.mask,2) && ...
                                y_old >= 0 && y_old < size(roi_old.mask,1) && ...
                                roi_old.mask(y_old+1,x_old+1))  % Should really check regionmask, but for ease just check the whole mask
                                % Append
                                convertseedinfo_prelim = horzcat(convertseedinfo_prelim,convertseedinfo_buffer); %#ok<AGROW>
                                list_region_old(convertseedinfo_buffer.num_region_old+1) = true;
                                
                                % Break from this region
                                successregion = true;
                            end   
                        end 
                    end
                    % Break if seed was successful or all regions in the
                    % old configuration have been analyzed
                    if (successregion || all(list_region_old))
                        break;
                    end
                end
                % Break if seed was successful or all regions in the
                % old configuration have been analyzed
                if (successregion || all(list_region_old))
                    break;
                end
            end
            % Break if seed was successful or all regions in the
            % old configuration have been analyzed
            if (successregion || all(list_region_old))
                break;
            end
        end
    end
    
    % Assign Outputs
    if (~isempty(convertseedinfo_prelim))
        for i = 0:length(convertseedinfo_prelim)-1
            convertseedinfo(i+1) = convertseedinfo_prelim(i+1);
        end  
        outstate = out.success;
    end    
end

%-------------------------------------------------------------------------%
% Nonlinear solver equations ---------------------------------------------%
%-------------------------------------------------------------------------%

function [paramvector,num_region_old,outstate] = calcpoint(x_new,y_new,plot_u_old,plot_v_old,plot_u_interp_old,plot_v_interp_old,roi_old,list_region_old,spacing,border_interp)
% This function obtains the parameters for a "convert" seed.
%
% Inputs -----------------------------------------------------------------%
%   x_new - integer; x position WRT the "new" configuration. Note that
%   x_new is WRT reduced coordinates.
%   y_new - integer; y position WRT the "new" configuration. Note that
%   y_new is WRT reduced coordinates.
%   plot_u_old - double array; u displacements WRT the "old" configuration.
%   Units are pixels.
%   plot_v_old - double array; v displacements WRT the "old" configuration.
%   Units are pixels.
%   plot_u_interp_old - cell; array of b-spline coefficients; one per
%   region.
%   plot_v_interp_old - cell; array of b-spline coefficients; one per
%   region.
%   roi_old - ncorr_class_roi; ROI corresponding to "old" displacements.
%   Note that ROI is already reduced by default.
%   list_region_old - logical array; keeps track of which "old" regions have 
%   been analyzed
%   spacing - integer; this is the spacing parameter.
%   border_interp - integer; the amount of padding used around the borders in
%   interpdata
%
% Outputs ----------------------------------------------------------------%
%   paramvector - double array; [x_new y_new x_old y_old u_old v_old distance]
%   num_region_old - integer; The region number which seed resides in for 
%   the "old" configuration
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.

    % Initialize outputs
    paramvector = [];
    num_region_old = [];
    outstate = out.failed;
    
    % Perform global search - note: defvector = [x_old y_old] - reduced
    [defvector_init,num_region_old_prelim,outstate_initialguess] = initialguess(x_new, ...
                                                                                y_new, ...
                                                                                plot_u_old, ...
                                                                                plot_v_old, ...
                                                                                roi_old, ...
                                                                                list_region_old, ...
                                                                                spacing);
    
    if (outstate_initialguess == out.success)
        % Perform an iterative search
        [defvector,u_old,v_old,distance,outstate_iterative] = iterativesearch(x_new, ...
                                                                              y_new, ...
                                                                              defvector_init, ...
                                                                              plot_u_interp_old, ...
                                                                              plot_v_interp_old, ...
                                                                              roi_old, ... 
                                                                              num_region_old_prelim, ...
                                                                              spacing, ...
                                                                              border_interp);
                                                                             
        if (outstate_iterative == out.success)       
            % Set outputs
            paramvector = [x_new y_new defvector u_old v_old distance];
            num_region_old = num_region_old_prelim;
            outstate = out.success;
        end
    end
end

function [defvector_init,num_region_old,outstate] = initialguess(x_new,y_new,plot_u_old,plot_v_old,roi_old,list_region_old,spacing)
% This function finds the closest integer displacements as an initial
% guess.
%
% Inputs -----------------------------------------------------------------%
%   x_new - integer; x position WRT the "new" configuration. Note that
%   x_new is WRT reduced coordinates.
%   y_new - integer; y position WRT the "new" configuration. Note that
%   y_new is WRT reduced coordinates.
%   plot_u_old - double array; u displacements WRT the "old" configuration.
%   Units are pixels.
%   plot_v_old - double array; v displacements WRT the "old" configuration.
%   Units are pixels.
%   roi_old - ncorr_class_roi; ROI corresponding to "old" displacements.
%   Note that ROI is already reduced by default.
%   list_region_old - logical array; keeps track of which "old" regions have 
%   been analyzed
%   spacing - integer; this is the spacing parameter.
%
% Outputs ----------------------------------------------------------------%
%   defvector_init - integer array; [x_old y_old] - reduced
%   num_region_old - integer; number of region where x_old and y_old were found
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.

    % Initialize Outputs
    outstate = out.failed;
    defvector_init = [];
    num_region_old = [];       
        
    % Cycle through every point to insure a global minimum    
    x_old_prelim = -1;
    y_old_prelim = -1;
    num_region_old_prelim = -1;
    distance_prelim = inf; % arbitrarily large number            
    for i = 0:length(roi_old.region)-1 
        if (list_region_old(i+1)) % this ROI has been analyzed already
            continue;
        else 
            for j = 0:size(roi_old.region(i+1).noderange,1)-1
                x_old_buffer = j + roi_old.region(i+1).leftbound; 
                for k = 0:2:roi_old.region(i+1).noderange(j+1)-1
                    for l = roi_old.region(i+1).nodelist(j+1,k+1):roi_old.region(i+1).nodelist(j+1,k+2)
                        y_old_buffer = l;
                        % Must divide displacement by spacing so
                        % displacement is WRT reduced coordinates
                        u_old_buffer = plot_u_old(y_old_buffer+1,x_old_buffer+1)/(spacing+1); 
                        v_old_buffer = plot_v_old(y_old_buffer+1,x_old_buffer+1)/(spacing+1);
                        distance_buffer = sqrt((x_new-(x_old_buffer+u_old_buffer))^2+(y_new-(y_old_buffer+v_old_buffer))^2);
                        
                        % Check if this point is better
                        if (distance_buffer < distance_prelim)
                            x_old_prelim = x_old_buffer;
                            y_old_prelim = y_old_buffer;
                            num_region_old_prelim = i;
                            distance_prelim = distance_buffer;
                        end                                        
                    end            
                end 
            end 
        end
    end      
    
    if (x_old_prelim ~= -1 && y_old_prelim ~= -1)
        defvector_init = [x_old_prelim y_old_prelim];
        num_region_old = num_region_old_prelim;
        outstate = out.success;
    end
end

function [defvector,u_old,v_old,distance,outstate] = iterativesearch(x_new,y_new,defvector_init,plot_u_interp_old,plot_v_interp_old,roi_old,num_region_old,spacing,border_interp)
% This function uses Gauss-Newton iterations to find the subpixel x_old
% and y_old values.
%
% Inputs -----------------------------------------------------------------%
%   x_new - integer; x position WRT the "new" configuration. Note that
%   x_new is WRT reduced coordinates.
%   y_new - integer; y position WRT the "new" configuration. Note that
%   y_new is WRT reduced coordinates.
%   defvector_init - integer array; of form [x_old y_old]. It's integer
%   because these are the initial guesses from the global search.
%   plot_u_interp_old - cell; array of b-spline coefficients; one per
%   region.
%   plot_v_interp_old - cell; array of b-spline coefficients; one per
%   region.
%   roi_old - ncorr_class_roi; ROI corresponding to "old" displacements.
%   Note that ROI is already reduced by default.
%   num_region_old - integer; number of region being analyzed
%   spacing - integer; this is the spacing parameter.
%   border_interp - integer; the amount of padding used around the borders in
%   interpdata
%
% Outputs ----------------------------------------------------------------%
%   defvector - double array; [x_old y_old] - reduced
%   u_old - double; u displacement from x_old to x_new - pixels
%   v_old - double; v displacement from y_old to y_new - pixels
%   distance - double; distance between [x_new y_new] and [x_old y_old]
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.

    % Initialize outputs 
    outstate = out.failed;
    defvector = [];
    u_old = [];
    v_old = [];
    distance = []; 
    
    % Gauss Newton optimization - send only b-spline coefficients
    % corresponding to the num_region_old
    [defvector_prelim,u_old_prelim,v_old_prelim,distance_prelim,gradnorm,outstate_newton] = newton(x_new, ...
                                                                                                   y_new, ...
                                                                                                   defvector_init, ...
                                                                                                   plot_u_interp_old{num_region_old+1}, ...
                                                                                                   plot_v_interp_old{num_region_old+1}, ...
                                                                                                   roi_old.region(num_region_old+1), ...
                                                                                                   spacing, ...
                                                                                                   border_interp);
    counter = 1;    
    while (outstate_newton == out.success && gradnorm > 10^-5 && counter < 10)
        % Gauss Newton optimization - send only b-spline coefficients
        % corresponding to the num_region_old
        [defvector_prelim,u_old_prelim,v_old_prelim,distance_prelim,gradnorm,outstate_newton] = newton(x_new, ...
                                                                                                       y_new, ...
                                                                                                       defvector_prelim, ...
                                                                                                       plot_u_interp_old{num_region_old+1}, ...
                                                                                                       plot_v_interp_old{num_region_old+1}, ...
                                                                                                       roi_old.region(num_region_old+1), ...
                                                                                                       spacing, ...
                                                                                                       border_interp);
        counter = counter + 1;
    end

    if (outstate_newton == out.success)
        % Assign outputs
        defvector = defvector_prelim;
        u_old = u_old_prelim;
        v_old = v_old_prelim;
        distance = distance_prelim;
        outstate = out.success;
    end
end

function [defvector,u_old,v_old,distance,gradnorm,outstate] = newton(x_new,y_new,defvector_init,plot_u_interp_old,plot_v_interp_old,region_old,spacing,border_interp)
% This function actually performs the Gauss-Newton iteration
%
% Inputs -----------------------------------------------------------------%
%   x_new - integer; x position WRT the "new" configuration. Note that
%   x_new is WRT reduced coordinates.
%   y_new - integer; y position WRT the "new" configuration. Note that
%   y_new is WRT reduced coordinates.
%   defvector_init - double array; of form [x_old y_old]
%   plot_u_interp_old - double array; array of b-spline coefficients
%   corresponding to region_old
%   plot_v_interp_old - double array; array of b-spline coefficients
%   corresponding to region_old
%   region_old - struct; specific region being analyzed
%   spacing - double; this is the spacing parameter.
%   border_interp - integer; the amount of padding used around the borders in
%   interpdata
%
% Outputs ----------------------------------------------------------------%
%   defvector - double array; [x_old y_old] - reduced
%   u_old - double; u displacement from x_old to x_new - pixels
%   v_old - double; v displacement from y_old to y_new - pixels
%   distance - double; distance between [x_new y_new] and [x_old y_old]
%   gradnorm - double; norm of the gradient vector - should be close to 0
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.

    % Initialize inputs
    outstate = out.failed;
    defvector = [];
    u_old = [];
    v_old = [];
    gradnorm = [];
    distance = [];

    % Use interp function local to this m-file since it combines
    % interpolation of values and gradients.
    [interpvector,outstate_interp] = interpqbs_convert(defvector_init, ...
                                                       plot_u_interp_old, ...
                                                       plot_v_interp_old, ...
                                                       region_old, ...
                                                       border_interp);    

    if (outstate_interp == out.success)
        % Determine Gradient - note that interpolation found through
        % interpqbs_convert needs to be WRT reduced coordinates; However,
        % displacements use pixel units, so they need to be scaled in order
        % to be WRT reduced coordinates.
        gradient(1) = -2*((x_new-(defvector_init(1)+interpvector(1)/(spacing+1)))*(1+interpvector(3)/(spacing+1))+(y_new-(defvector_init(2)+interpvector(2)/(spacing+1)))*(interpvector(5)/(spacing+1)));
        gradient(2) = -2*((x_new-(defvector_init(1)+interpvector(1)/(spacing+1)))*(interpvector(4)/(spacing+1))+(y_new-(defvector_init(2)+interpvector(2)/(spacing+1)))*(1+interpvector(6)/(spacing+1)));
        
        % Determine Hessian
        hessian(1,1) = 2*((1+interpvector(3)/(spacing+1))^2+(interpvector(5)/(spacing+1))^2);
        hessian(2,1) = 2*((interpvector(4)/(spacing+1))*(1+interpvector(3)/(spacing+1))+(1+interpvector(6)/(spacing+1))*(interpvector(5)/(spacing+1)));
        hessian(1,2) = hessian(2,1); % symmetric
        hessian(2,2) = 2*((interpvector(4)/(spacing+1))^2+(1+interpvector(6)/(spacing+1))^2);
        
        det_hess = det(hessian);
        % Check to make sure hessian is positive definite
        % From :http://www.math.northwestern.edu/~clark/285/2006-07/handouts/pos-def.pdf
        % Make sure det(hess) > 0 and hess(1,1) > 0 
        if (det_hess > 0 && hessian(1) > 0)
            % Determine new coordinates
            defvector_prelim = (defvector_init'-hessian^-1*gradient')';

            % Calculate distance - we have to interpolate again with the new
            % coordinates
            [interpvector,outstate_interp] = interpqbs_convert(defvector_prelim, ...
                                                               plot_u_interp_old, ...
                                                               plot_v_interp_old, ...
                                                               region_old, ...
                                                               border_interp); 


            if (outstate_interp == out.success)            
                % Store outputs
                defvector = defvector_prelim;
                u_old = interpvector(1);
                v_old = interpvector(2);
                distance = sqrt((x_new-(defvector(1)+interpvector(1)/(spacing+1)))^2+(y_new-(defvector(2)+interpvector(2)/(spacing+1)))^2);                
                gradnorm = norm(gradient);
                outstate = out.success;
            end
        end
    end
end

function [interpvector,outstate] = interpqbs_convert(defvector,plot_u_interp_old,plot_v_interp_old,region_old,border_interp)
% This interpolation combines interpolation of values and gradients, so
% it's reimplemented here because it can save some time by combining the
% interpolation.
%
% Inputs -----------------------------------------------------------------%
%   defvector - double array; [x_old y_old] - reduced
%   plot_u_interp_old - double array; array of b-spline coefficients. Units
%   are pixels.
%   plot_v_interp_old - double array; array of b-spline coefficients. Units
%   are pixels.
%   region_old - struct; specific region being analyzed
%   border_interp - integer; the amount of padding used around the edges in
%   interpdata
%
% Outputs ----------------------------------------------------------------%
%   interpvector - double array; [u v du/dx du/dy dv/dx dv/dy]
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.

    % Initialize output 
    outstate = out.failed;
    interpvector = [];

    % Biquintic Kernel Matrix
    QK = [1/120  13/60  11/20 13/60  1/120 0;
         -1/24   -5/12    0    5/12  1/24  0;
          1/12    1/6   -1/2   1/6   1/12  0;
         -1/12    1/6     0   -1/6   1/12  0;
          1/24   -1/6    1/4  -1/6   1/24  0;
         -1/120   1/24  -1/12  1/12 -1/24 1/120];  
     
    % Interpolate if in bounds
    x_tilda = defvector(1);
    y_tilda = defvector(2);
    
    x_tilda_floor = floor(x_tilda);
    y_tilda_floor = floor(y_tilda);
        
    % Make sure top, left, bottom, and right are within the b-spline 
    % coefficient array. top, left, bottom and right are the bounding 
    % box of the b-spline coefficients used for interpolation of this
    % point;
    top = y_tilda_floor-region_old.upperbound+border_interp-2;
    left = x_tilda_floor-region_old.leftbound+border_interp-2;
    bottom = y_tilda_floor-region_old.upperbound+border_interp+3;
    right = x_tilda_floor-region_old.leftbound+border_interp+3;
        
    if (top >= 0 && ...
        left >= 0 && ...
        bottom < size(plot_u_interp_old,1) && ...
        right < size(plot_u_interp_old,2))
        % Set coords
        x_tilda_delta = x_tilda-x_tilda_floor;
        y_tilda_delta = y_tilda-y_tilda_floor;

        x_vec(1) = 1.0;
        x_vec(2) = x_tilda_delta;
        x_vec(3) = x_vec(2)*x_tilda_delta;
        x_vec(4) = x_vec(3)*x_tilda_delta;
        x_vec(5) = x_vec(4)*x_tilda_delta;
        x_vec(6) = x_vec(5)*x_tilda_delta;

        y_vec(1) = 1.0;
        y_vec(2) = y_tilda_delta;
        y_vec(3) = y_vec(2)*y_tilda_delta;
        y_vec(4) = y_vec(3)*y_tilda_delta;
        y_vec(5) = y_vec(4)*y_tilda_delta;
        y_vec(6) = y_vec(5)*y_tilda_delta;
        
        x_vec_dx(1) = 0.0;
        x_vec_dx(2) = 1.0;
        x_vec_dx(3) = 2.0*x_vec(2);
        x_vec_dx(4) = 3.0*x_vec(3);
        x_vec_dx(5) = 4.0*x_vec(4);
        x_vec_dx(6) = 5.0*x_vec(5);

        y_vec_dy(1) = 0.0;
        y_vec_dy(2) = 1.0;
        y_vec_dy(3) = 2.0*y_vec(2);
        y_vec_dy(4) = 3.0*y_vec(3);
        y_vec_dy(5) = 4.0*y_vec(4);
        y_vec_dy(6) = 5.0*y_vec(5);
                
        % Precompute
        QKMAT_u_plot = QK*plot_u_interp_old(top+1:bottom+1,left+1:right+1)*QK';
        QKMAT_v_plot = QK*plot_v_interp_old(top+1:bottom+1,left+1:right+1)*QK';

        % Get interpolated value
        interpvector(1) = y_vec*QKMAT_u_plot*x_vec';        % u
        interpvector(2) = y_vec*QKMAT_v_plot*x_vec';        % v
        interpvector(3) = y_vec*QKMAT_u_plot*x_vec_dx';     % du/dx
        interpvector(4) = y_vec_dy*QKMAT_u_plot*x_vec';     % du/dy
        interpvector(5) = y_vec*QKMAT_v_plot*x_vec_dx';     % dv/dx
        interpvector(6) = y_vec_dy*QKMAT_v_plot*x_vec';     % dv/dy
        outstate = out.success;
    end
end

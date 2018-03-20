classdef ncorr_class_roi < handle
% This is the class definition for the region of interest.
    
    % Properties ---------------------------------------------------------%
    properties(SetAccess = private)
        type;               % string
        mask;               % logical array
        region;             % struct('nodelist',{},'noderange',{},'leftbound',{},'rightbound',{},'upperbound',{},'lowerbound',{},'totalpoints',{})
        boundary;           % struct('add',{},'sub',{})
        data;               % varying struct
    end
        
    % Methods ------------------------------------------------------------%
    methods(Access = public)
        % Constructor
        function obj = ncorr_class_roi
            obj.type = '';
            obj.mask = false(0);
            obj.region = struct('nodelist',{},'noderange',{},'leftbound',{},'rightbound',{},'upperbound',{},'lowerbound',{},'totalpoints',{});
            obj.boundary = struct('add',{},'sub',{});
            obj.data = struct();
        end
        
        function set_roi(obj,type_i,data_i)       
        % This function sets the region of interest.
        %
        % Inputs ---------------------------------------------------------%
        %   type_i - string; describes how ROI was made. Supported types
        %   shown below:
        %       'load' - mask was loaded from an image file. data_i is 
        %       struct('mask',{},'cutoff',{}).
        %       'draw' - mask was drawn. data_i is
        %       struct('mask',{}','drawobjects',{},'cutoff',{}).
        %       'region' - region and size of the mask are provided. 
        %       data_i is struct('region',{},'size_mask',{}).
        %       'boundary' - boundary and size of the mask are provided. 
        %       data_i is struct('boundary',{},'size_mask',{}).
        %   data_i - struct; contains input data for the type of ROI.
        %
        % Outputs --------------------------------------------------------%
        %   None
        %
        % Returns error if type is wrong.
        
            if (~strcmp(type_i,'load') && ~strcmp(type_i,'draw') && ~strcmp(type_i,'region') && ~strcmp(type_i,'boundary'))  
                error('Incorrect type provided.');
            end
            
            if (strcmp(type_i,'load') || strcmp(type_i,'draw'))
                % This is the "load" and "draw" types, which mean a
                % mask for the ROI is provided directly.

                % Get mask - mask is directly provided for 'load' and 'draw'
                mask_prelim = data_i.mask;

                % Create region(s) - these are 4-way connected/contiguous
                [region_prelim,removed] = ncorr_alg_formregions(mask_prelim,int32(data_i.cutoff),false);

                % Get 20 largest contiguous regions. This is to prevent
                % boundary routine from running very slowly if there are a
                % bunch of small regions
                cutoff_regions = 20;
                [val_sorted,idx_sorted] = sort([region_prelim.totalpoints],'descend'); %#ok<ASGLU>
                idx_sorted = idx_sorted-1; % Convert to zero based
                if (length(idx_sorted) > cutoff_regions)
                    region_prelim = region_prelim(idx_sorted(1:cutoff_regions)+1);
                    removed = true;
                end

                % Update mask if regions were removed
                if (removed)
                    mask_prelim(:) = false;
                    for i = 0:length(region_prelim)-1
                        for j = 0:size(region_prelim(i+1).noderange,1)-1
                            x = j + region_prelim(i+1).leftbound;
                            for k = 0:2:region_prelim(i+1).noderange(j+1)-1
                                vec_y = region_prelim(i+1).nodelist(j+1,k+1):region_prelim(i+1).nodelist(j+1,k+2);
                                mask_prelim(vec_y+1,x+1) = true;
                            end
                        end
                    end
                end 

                % Now get the boundaries - Do the "add boundaries"
                % first. One per region.
                boundary_prelim = struct('add',cell(length(region_prelim),1),'sub',cell(length(region_prelim),1));
                for i = 0:length(region_prelim)-1
                    % Get mask corresponding to region - must do this
                    % because boundary tracing algorithm is 8-way
                    % connected, but the region is 4-way connected. I
                    % used the 8-way connected boundary because it is
                    % smoother. Note that if two 4-way connected
                    % regions are just touching at corners they are
                    % technically considered 8-way connected, which is
                    % why each 4-way region must be analyzed separately.
                    regionmask_buffer = false(size(mask_prelim));               
                    for j = 0:size(region_prelim(i+1).noderange,1)-1
                        x = j + region_prelim(i+1).leftbound;
                        for k = 0:2:region_prelim(i+1).noderange(j+1)-1
                            vec_y = region_prelim(i+1).nodelist(j+1,k+1):region_prelim(i+1).nodelist(j+1,k+2);
                            regionmask_buffer(vec_y+1,x+1) = true;
                        end
                    end

                    % Get add boundary:
                    % Use the top of the left most point in the region.
                    % Send this coordinate and the regionmask_buffer 
                    % corresponding to this region to determine the 
                    % boundary
                    boundary_buffer = ncorr_alg_formboundary(int32([region_prelim(i+1).leftbound region_prelim(i+1).nodelist(1,1)]),int32(0),regionmask_buffer);
                    % Append boundary
                    boundary_prelim(i+1).add = boundary_buffer;

                    % Get inverse of regionmask_buffer, and then do logical
                    % "and" with the filled in outer boundary. This
                    % leaves only internal regions 
                    regionmask_inv_buffer = false(size(mask_prelim));
                    % Get filled-in outter boundary - note that this
                    % does not fill the boundary exactly since the
                    % algorithm uses double precision points to
                    % estimate the boundary. 
                    ncorr_alg_formmask(struct('pos_imroi',boundary_prelim(i+1).add,'type','poly','addorsub','add'),regionmask_inv_buffer);                            
                    regionmask_inv_buffer = regionmask_inv_buffer & ~regionmask_buffer;

                    % Form mask_boundary - this keeps track of the sub
                    % boundaries already analyzed so we dont count one
                    % twice.
                    mask_boundary = false(size(mask_prelim)); % Make a copy to keep track of which points have been used to analyze the boundary

                    % Get sub boundaries:   
                    boundary_prelim(i+1).sub = cell(0);
                    for j = 0:size(region_prelim(i+1).noderange,1)-1
                        x = j + region_prelim(i+1).leftbound;
                        % Make sure noderange is greater than 2, or
                        % else this portion doesnt have a hole
                        if (region_prelim(i+1).noderange(j+1) > 2)
                            for k = 0:2:region_prelim(i+1).noderange(j+1)-3 % Dont test last node pair
                                % Only use bottom nodes - only test
                                % bottom nodes since these are the
                                % first points adjacent to a hole
                                y_bottom = region_prelim(i+1).nodelist(j+1,k+2);

                                % Test one pixel below bottom node. Make 
                                % sure this pixel hasnt already been analyzed, 
                                % and is also within the inverse regionmask
                                % buffer
                                if (regionmask_inv_buffer(y_bottom+2,x+1) && ~mask_boundary(y_bottom+2,x+1))
                                    % This is a boundary point which
                                    % hasn't been analyzed yet
                                    boundary_prelim(i+1).sub{end+1} = ncorr_alg_formboundary(int32([x y_bottom+1]),int32(0),regionmask_inv_buffer);
                                    % update mask_boundary
                                    mask_boundary(sub2ind(size(mask_boundary),boundary_prelim(i+1).sub{end}(:,2)+1,boundary_prelim(i+1).sub{end}(:,1)+1)) = true;
                                end
                            end
                        end
                    end   
                end
            elseif (strcmp(type_i,'region'))
                % This is the 'region' option; regions are supplied
                % directly. It's possible that analysis yielded regions
                % directly (like unioning).
                
                % Form mask
                mask_prelim = false(data_i.size_mask);

                % Form region
                region_prelim = data_i.region; 

                % Update mask from region:
                for i = 0:length(region_prelim)-1
                    for j = 0:size(region_prelim(i+1).noderange,1)-1
                        x = j + region_prelim(i+1).leftbound;
                        for k = 0:2:region_prelim(i+1).noderange(j+1)-1
                            vec_y = region_prelim(i+1).nodelist(j+1,k+1):region_prelim(i+1).nodelist(j+1,k+2);
                            mask_prelim(vec_y+1,x+1) = true;   
                        end
                    end
                end

                % For the custom case, do not form a real boundary, since
                % regions may not longer be contigous. Do form one boundary
                % per region though.
                for i = 0:length(region_prelim)-1
                    boundary_prelim(i+1).add = [-1 -1]; %#ok<AGROW>
                    boundary_prelim(i+1).sub = {}; %#ok<AGROW>
                end
            elseif (strcmp(type_i,'boundary'))
                % This is the 'boundary' option; a boundary is supplied
                % directly. This allows Ncorr to update the mask based 
                % on displacement values at the boundary point. Note 
                % there will be one region per "add" boundary, this 
                % preserves the correspondences between a region and 
                % boundary.

                % Form mask
                mask_prelim = false(data_i.size_mask);

                % Form boundary
                boundary_prelim = data_i.boundary;

                % Fill a region for each boundary -> Get the region 
                % corresponding to that mask -> append all
                % regions. 
                % Set drawobjects which are sent to ncorr_alg_formmask
                region_prelim = struct('nodelist',{},'noderange',{},'leftbound',{},'rightbound',{},'upperbound',{},'lowerbound',{},'totalpoints',{});   
                for i = 0:length(boundary_prelim)-1
                    % Initialize counter and draw objects
                    drawobjects = struct('pos_imroi',{},'type',{},'addorsub',{});
                    counter = 0;

                    drawobjects(counter+1).pos_imroi = boundary_prelim(i+1).add;
                    drawobjects(counter+1).type = 'poly';
                    drawobjects(counter+1).addorsub = 'add';      
                    counter = counter+1;
                    for j = 0:length(boundary_prelim(i+1).sub)-1
                        drawobjects(counter+1).pos_imroi = boundary_prelim(i+1).sub{j+1};
                        drawobjects(counter+1).type = 'poly';
                        drawobjects(counter+1).addorsub = 'sub';      
                        counter = counter+1;
                    end

                    % Get mask corresponding to this boundary
                    ncorr_alg_formmask(drawobjects,mask_prelim);

                    % Get region
                    [region_buffer,removed] = ncorr_alg_formregions(mask_prelim,int32(0),false); %#ok<ASGLU>

                    % There must be one region per boundary to preserve 
                    % their correspondence. It's  possible for the region 
                    % for a boundary to "pinch" or form more than one. 
                    if (isempty(region_buffer))
                        % Form an empty ROI list if region_buffer is
                        % empty - this keeps the correspondence between
                        % the boundary and region
                        region_buffer = struct('nodelist',[-1 -1],'noderange',0,'leftbound',0,'rightbound',0,'upperbound',0,'lowerbound',0,'totalpoints',0);
                    elseif (length(region_buffer) > 1)
                        % Select biggest ROI if there are more than
                        % one. This could possibly happen if a boundary
                        % is "pinched" or closes. Unlikely- but may
                        % happen.
                        idx_max = find([region_buffer.totalpoints] == max([region_buffer.totalpoints]),1,'first');
                        region_buffer = region_buffer(idx_max);
                    end

                    % Merge           
                    region_prelim(i+1) = region_buffer;
                end

                % Update mask from region:
                for i = 0:length(region_prelim)-1
                    for j = 0:size(region_prelim(i+1).noderange,1)-1
                        x = j + region_prelim(i+1).leftbound;
                        for k = 0:2:region_prelim(i+1).noderange(j+1)-1
                            vec_y = region_prelim(i+1).nodelist(j+1,k+1):region_prelim(i+1).nodelist(j+1,k+2);
                            mask_prelim(vec_y+1,x+1) = true;   
                        end
                    end
                end
            end

            % Set properties
            obj.type = type_i;
            obj.mask = mask_prelim;
            obj.region = region_prelim;
            obj.boundary = boundary_prelim;
            obj.data = data_i;  
        end
        
        function roi_reduced = reduce(obj,spacing) 
        % This function reduces the ROI. Its possible some regions will 
        % be deleted after the reduction; a placeholder is still used to 
        % preserve the number of regions.
        %
        % Inputs ---------------------------------------------------------%
        %   spacing - integer; spacing parameter
        %
        % Outputs --------------------------------------------------------%
        %   roi_reduced - ncorr_class_roi; reduced ROI        
        % 
        % Returns error if ROI has not been set yet.
        
            if (isempty(obj.type))  
                error('ROI has not been set yet');
            end
            
            % Initialize output
            roi_reduced = ncorr_class_roi.empty;

            if (spacing == 0)
                % Don't do anything, the ROI is not actually being reduced. 
                % Just make a deep copy of this ROI and return;
                roi_reduced(1) = ncorr_class_roi;
                p = properties(obj);
                for i = 0:length(p)-1
                    roi_reduced.(p{i+1}) = obj.(p{i+1});
                end                
            else 
                % Form reduced region and its template
                region_reduced = struct('nodelist',{},'noderange',{},'leftbound',{},'rightbound',{},'upperbound',{},'lowerbound',{},'totalpoints',{});  
                region_reduced_template = struct('nodelist',{},'noderange',{},'leftbound',{},'rightbound',{},'upperbound',{},'lowerbound',{},'totalpoints',{});  
                for i = 0:length(obj.region)-1
                    % Find new bounds - left and right bounds are the
                    % upper limit.
                    region_reduced_template(1).leftbound = ceil(obj.region(i+1).leftbound/(spacing+1));
                    region_reduced_template.rightbound = floor(obj.region(i+1).rightbound/(spacing+1));
                    region_reduced_template.upperbound = inf; % Initialize to infinitely high so it can updated
                    region_reduced_template.lowerbound = -inf; % Initialize to infinitely low so it can updated

                    % Initialize nodelist, noderange, and totalpoints
                    region_reduced_template.nodelist = -ones(max(1,region_reduced_template.rightbound-region_reduced_template.leftbound+1),size(obj.region(i+1).nodelist,2));
                    region_reduced_template.noderange = zeros(max(1,region_reduced_template.rightbound-region_reduced_template.leftbound+1),1);
                    region_reduced_template.totalpoints = 0;

                    % Scan columns - calculate nodelist, noderange,
                    % and totalpoints
                    for j = region_reduced_template.leftbound*(spacing+1)-obj.region(i+1).leftbound:spacing+1:region_reduced_template.rightbound*(spacing+1)-obj.region(i+1).leftbound
                        % This is the column which lies on the grid;
                        % Form the reduced nodelist and noderange for
                        % this column
                        x_reduced = (j-(region_reduced_template.leftbound*(spacing+1)-obj.region(i+1).leftbound))/(spacing+1);
                        for k = 0:2:obj.region(i+1).noderange(j+1)-1
                            % Find upper and lower bounds for each node pair -
                            % it is possible that entire node pair is skipped
                            % over by the spacing
                            node_top_reduced = ceil(obj.region(i+1).nodelist(j+1,k+1)/(spacing+1));
                            node_bottom_reduced = floor(obj.region(i+1).nodelist(j+1,k+2)/(spacing+1));                        
                            if (node_bottom_reduced >= node_top_reduced)
                                % Add this node pair                                        
                                region_reduced_template.nodelist(x_reduced+1,region_reduced_template.noderange(x_reduced+1)+1) = node_top_reduced;
                                region_reduced_template.nodelist(x_reduced+1,region_reduced_template.noderange(x_reduced+1)+2) = node_bottom_reduced;

                                % Update node range
                                region_reduced_template.noderange(x_reduced+1) = region_reduced_template.noderange((j-(region_reduced_template.leftbound*(spacing+1)-obj.region(i+1).leftbound))/(spacing+1)+1)+2;

                                % Update totalpoints
                                region_reduced_template.totalpoints = region_reduced_template.totalpoints + (node_bottom_reduced-node_top_reduced+1);

                                % Update upper and lower bounds
                                if (node_top_reduced < region_reduced_template.upperbound)
                                    region_reduced_template.upperbound = node_top_reduced;
                                end
                                if (node_bottom_reduced > region_reduced_template.lowerbound)
                                    region_reduced_template.lowerbound = node_bottom_reduced;
                                end
                            end      
                        end
                    end                         

                    % See if region is empty - if so use a placeholder.
                    if (region_reduced_template.totalpoints == 0)    
                        region_reduced_template.nodelist = [-1 -1];
                        region_reduced_template.noderange = 0;
                        region_reduced_template.leftbound = 0;     
                        region_reduced_template.rightbound = 0;   
                        region_reduced_template.upperbound = 0;    
                        region_reduced_template.lowerbound = 0;   
                    end

                    % Insert template into reduced region
                    region_reduced(i+1) = region_reduced_template;
                end

                % Create new ROI
                roi_reduced(1) = ncorr_class_roi;
                roi_reduced.set_roi('region',struct('region',region_reduced,'size_mask',size(obj.mask(1:spacing+1:end,1:spacing+1:end))));
            end
        end
        
        function roi_union = get_union(obj,mask_i,spacing)  
        % This function returns the union of the current ROI and the 
        % inputted mask (which may be reduced).
        %
        % Inputs ---------------------------------------------------------%
        %   mask_i - logical array; mask used to union
        %   spacing - integer; indicates how much mask_i has been reduced.
        % 
        % Outputs --------------------------------------------------------%
        %   roi_union - ncorr_class_roi; union of this ROI with mask_i.
        % 
        % Returns error if ROI has not been set yet or if the sizes dont
        % match.
        
            if (isempty(obj.type))  
                error('ROI has not been set yet');
            elseif (~isequal(size(obj.mask(1:spacing+1:end,1:spacing+1:end)),size(mask_i)))
                error('Reduced mask and input mask are not the same size');
            end
            
            % Initialize output
            roi_union = ncorr_class_roi;

            % Get reduced ROI
            roi_reduced = obj.reduce(spacing);

            % Get unioned mask
            mask_union = obj.mask(1:spacing+1:end,1:spacing+1:end) & mask_i;

            % Get unioned region, this will not necessarily be contiguous anymore                    
            region_union = ncorr_alg_formunion(roi_reduced.formatted(),mask_union);

            % Create new ROI
            roi_union.set_roi('region',struct('region',region_union,'size_mask',size(mask_union)));
        end
        
        function roi_update = update_roi(obj,plot_u,plot_v,roi_plot,size_mask_update,spacing,radius)
        % This function returns an updated ROI, based on the inputted
        % displacement fields. This works by updating the boundary and
        % redrawing the mask. This assumes boundary can be subpixel, so it
        % interpolates a new position based on the displacement field.
        % Furthermore, displacement fields are reduced by a spacing factor,
        % so it's required as an input.
        %
        % Inputs ---------------------------------------------------------%
        %   plot_u - double array; u displacement plot
        %   plot_v - double array; v displacement plot
        %   roi_plot - ncorr_class_roi; ROI corresponding to displacement
        %   plots
        %   size_mask_update - integer array; size of the updated mask. Can
        %     be a different size than obj.mask if different sized image are 
        %     used
        %   spacing - integer; spacing parameter 
        %   radius - integer; subset radius
        % 
        % Outputs --------------------------------------------------------%
        %   roi_update - ncorr_class_roi; updated ROI.
        % 
        % Returns error if ROI has not been set yet.
        
            if (isempty(obj.type))  
                error('ROI has not been set yet');
            end
            
            % Extrapolate displacement plots to improve interpolation near the
            % boundary points.
            % border_interp is the border added to the displacements when they are
            % extrapolated.
            border_interp = 20; % MUST BE GREATER THAN OR EQUAL TO 2!!!!
            plot_u_interp = ncorr_alg_extrapdata(plot_u,roi_plot.formatted(),int32(border_interp));
            plot_v_interp = ncorr_alg_extrapdata(plot_v,roi_plot.formatted(),int32(border_interp));

            % Convert displacement plots to B-spline coefficients; these are used 
            % for interpolation. Make sure to do this for each region.
            for i = 0:length(roi_plot.region)-1
                plot_u_interp{i+1} = ncorr_class_img.form_bcoef(plot_u_interp{i+1});
                plot_v_interp{i+1} = ncorr_class_img.form_bcoef(plot_v_interp{i+1});
            end

            % Initialize new boundary - preserve the number and
            % correspondences of the "add" boundaries, even if a
            % boundary is empty
            boundary_update = struct('add',{},'sub',{});
            for i = 0:length(obj.boundary)-1
                boundary_update(i+1).add = interp_border(obj.boundary(i+1).add, ...
                                                         plot_u_interp{i+1}, ...
                                                         plot_v_interp{i+1}, ...
                                                         roi_plot, ...
                                                         i, ...
                                                         size_mask_update, ...
                                                         border_interp, ...
                                                         spacing, ...
                                                         radius);     
                for j = 0:length(obj.boundary(i+1).sub)-1
                    boundary_update(i+1).sub{j+1} = interp_border(obj.boundary(i+1).sub{j+1}, ...
                                                                  plot_u_interp{i+1}, ...
                                                                  plot_v_interp{i+1}, ...
                                                                  roi_plot, ...
                                                                  i, ...
                                                                  size_mask_update, ...
                                                                  border_interp, ...
                                                                  spacing, ...
                                                                  radius);
                end
            end

            % Set new ROI  - note that size_mask_update must be used
            % instead of the size of the current mask. This is because
            % different sized current images can be used.
            roi_update = ncorr_class_roi;
            roi_update.set_roi('boundary',struct('boundary',boundary_update,'size_mask',size_mask_update));   
        end
        
        function roi_f = formatted(obj)
        % This function returns the formatted ROI, which can be inputted to a
        % mex function. Mex functions can receive ncorr_class_roi as either
        % a class or structure.
        %
        % Inputs ---------------------------------------------------------%
        %   none;
        % 
        % Outputs --------------------------------------------------------%
        %   roi_f - ncorr_class_roi; formatted ROI.
        % 
        % Returns error if ROI has not been set yet.
        
            if (isempty(obj.type))  
                error('ROI has not been set yet');
            end
            
            % Must make deep copy first
            roi_f = ncorr_class_roi;
            p = properties(obj);
            for i = 0:length(p)-1
                roi_f.(p{i+1}) = obj.(p{i+1});
            end     

            % Now format properties
            for i = 0:length(obj.region)-1
               roi_f.region(i+1).nodelist = int32(roi_f.region(i+1).nodelist);
               roi_f.region(i+1).noderange = int32(roi_f.region(i+1).noderange);
               roi_f.region(i+1).leftbound = int32(roi_f.region(i+1).leftbound);
               roi_f.region(i+1).rightbound = int32(roi_f.region(i+1).rightbound);
               roi_f.region(i+1).upperbound = int32(roi_f.region(i+1).upperbound);
               roi_f.region(i+1).lowerbound = int32(roi_f.region(i+1).lowerbound);
               roi_f.region(i+1).totalpoints = int32(roi_f.region(i+1).totalpoints);        
            end
        end        
        
        % ----------------------------------------------------------------%
        % Get functions --------------------------------------------------%
        % ----------------------------------------------------------------%
                
        function [num_region,idx_nodelist] = get_num_region(obj,x,y,list_region)
        % This function determines which region x and y lie in, excluding
        % the regions indicated by list_region
        %
        % Inputs ---------------------------------------------------------%
        %   x - integer; x_coordinate
        %   y - integer; y_coordinate
        %   list_region - logical array; indices of regions to search
        %   if x and y are contained within it. If list_region(i) == true,
        %   then skip over region(i).
        %
        % Output ---------------------------------------------------------%
        %   num_region - integer; index of region that contains x
        %   and y. Returns -1 if x and y are not within the regions
        %   specified by region.
        %   idx_nodelist - integer; index of node pair that contains y
        % 
        % Returns error if ROI has not been set yet.
            
            if (isempty(obj.type))  
                error('ROI has not been set yet');
            end
            
            % Initialize outputs
            num_region = -1;
            idx_nodelist = 0;

            % Initialize withinroi to false
            withinroi = false;
            for i = 0:length(obj.region)-1
                if (list_region(i+1)) % Skip this ROI
                    continue;
                else                    
                    if (x >= obj.region(i+1).leftbound && x <= obj.region(i+1).rightbound) % Check to make sure x coordinate is within bounds
                        idx_region = x-obj.region(i+1).leftbound;
                        for j =  0:2:obj.region(i+1).noderange(idx_region+1)-1 % Check to see if point is within node pairs
                            if (y >= obj.region(i+1).nodelist(idx_region+1,j+1) && y <= obj.region(i+1).nodelist(idx_region+1,j+2))
                                withinroi = true;
                                idx_nodelist = j;
                                num_region = i;
                                break;
                            end
                        end
                    end
                    if (withinroi)
                        break;
                    end
                end
            end 
        end
        
        function regionmask = get_regionmask(obj,num_region)
        % This function returns a mask corresponding only to the region
        % indicated by num_region.
        %
        % Inputs ---------------------------------------------------------%
        %   num_region - integer; number of region to use to form mask.
        %
        % Outputs --------------------------------------------------------%
        %   regionmask - logical array; mask corresponding to num_region
        % 
        % Returns error if ROI has not been set yet.
            
            if (isempty(obj.type))  
                error('ROI has not been set yet');
            end
            
            regionmask = false(size(obj.mask));                
            for i = 0:size(obj.region(num_region+1).noderange,1)-1
                x = i + obj.region(num_region+1).leftbound;
                for j = 0:2:obj.region(num_region+1).noderange(i+1)-1
                    vec_y = obj.region(num_region+1).nodelist(i+1,j+1):obj.region(num_region+1).nodelist(i+1,j+2);
                    regionmask(vec_y+1,x+1) = true;
                end
            end       
        end
                    
        function total_fullregions = get_fullregions(obj)
        % This function returns the number of regions which are
        % nonempty. Generally used to see if ROI is empty or not.
        %
        % Inputs ---------------------------------------------------------%
        %   none;
        %
        % Outputs --------------------------------------------------------%
        %   total_fullregions - integer; number of nonemtpy regions
        % 
        % Returns error if ROI has not been set yet.
            
            if (isempty(obj.type))  
                error('ROI has not been set yet');
            end
            
            total_fullregions = 0;         
            for i = 0:length(obj.region)-1
                if (obj.region(i+1).totalpoints > 0)
                    total_fullregions = total_fullregions+1;
                end
            end        
        end
        
        function cirroi = get_cirroi(obj,x,y,radius,subsettrunc)
        % This function returns a contiguous circular region for the subset
        % located at x and y.
        %
        % Inputs ---------------------------------------------------------%
        %   x - integer; x_coordinate of subset
        %   y - integer; y_coordinate of subset
        %   radius - integer; radius of subset
        %   subsettrunc - logical; indicates whether to enable subset
        %   truncation or not
        % 
        % Outputs --------------------------------------------------------%
        %   cirroi - circular subset; this is the circular subset
        %   corresponding to the points specified at x and y. Contains 
        %   struct('mask',{},'region',{},'boundary',{},'x',{},'y',{},'radius',{})
        % 
        % Returns error if ROI has not been set yet or if x and y are not
        % within the ROI. Note that the boundary field in cirroi is 
        % currently not being used.
        
            % Find idx_nodelist and num_region
            [num_region,idx_nodelist] = obj.get_num_region(x,y,false(size(obj.region)));  
            if (isempty(obj.type))  
                error('ROI has not been set yet');
            elseif (num_region == -1)
                error('x and y coordinates are not within the ROI');
            end
                
            % Initialize cirroi structure
            cirroi = struct('mask',{},'region',{},'boundary',{},'x',{},'y',{},'radius',{});
            cirroi_template = struct('mask',{},'region',{},'boundary',{},'x',{},'y',{},'radius',{});

            % Store x,y, and radius
            cirroi_template(1).mask = false(2*radius+1);  
            cirroi_template.region = struct('nodelist',{},'noderange',{},'leftbound',{},'rightbound',{},'upperbound',{},'lowerbound',{},'totalpoints',{});
            cirroi_template.boundary = [];
            cirroi_template.x = x;
            cirroi_template.y = y;
            cirroi_template.radius = radius;  

            % Find max_nodewidth_buffer - this is the maximum width
            % of the nodelist.
            max_nodewidth_buffer = 0;
            for i = 0:length(obj.region)-1  
                if (size(obj.region(i+1).nodelist,2) > max_nodewidth_buffer)
                    max_nodewidth_buffer = size(obj.region(i+1).nodelist,2);
                end
            end

            % Initialize cirroi nodelist and noderange
            cirroi_template.region(1).nodelist = -ones(cirroi_template.radius*2+1,max_nodewidth_buffer);     
            cirroi_template.region.noderange = zeros(cirroi_template.radius*2+1,1); 

            % Initialize totalpoints
            cirroi_template.region.totalpoints = 0;

            % Now find circle nodes
            % Fill must be contiguous and begins at the centerpoint
            queue_nodelist = -ones(size(cirroi_template.region.nodelist,1)*size(cirroi_template.region.nodelist,2),1); % Initialize to max size to prevent resizing
            queue_nodeindex = zeros(size(cirroi_template.region.nodelist,1)*size(cirroi_template.region.nodelist,2)/2,1); % Index of nodepair with respect to cirroi;
            length_queue = 0; %#ok<NASGU> % length of the queue_nodelist

            queue_nodelist_buffer = -ones(1,2); 
            queue_nodeindex_buffer = 0; %#ok<NASGU>

            % Initialize parameter which indicates whether the subset
            % has been truncated
            circ_untrunc = true;

            % Node pair which contains the center point is added to queue first; this
            % guarantees the circular region is contiguous WRT the centerpoint
            if (obj.region(num_region+1).nodelist(cirroi_template.x-obj.region(num_region+1).leftbound+1,idx_nodelist+1) < cirroi_template.y-cirroi_template.radius) % make sure top node is within the circle
                node_top = cirroi_template.y-cirroi_template.radius; 
            else
                node_top = obj.region(num_region+1).nodelist(cirroi_template.x-obj.region(num_region+1).leftbound+1,idx_nodelist+1);
                circ_untrunc = false;
            end
            if (obj.region(num_region+1).nodelist(cirroi_template.x-obj.region(num_region+1).leftbound+1,idx_nodelist+2) > cirroi_template.y+cirroi_template.radius) % make sure bottom node is within the circle
                node_bottom = cirroi_template.y+cirroi_template.radius; 
            else
                node_bottom = obj.region(num_region+1).nodelist(cirroi_template.x-obj.region(num_region+1).leftbound+1,idx_nodelist+2);
                circ_untrunc = false;
            end

            % Update mask
            cirroi_template.mask(node_top-(cirroi_template.y-cirroi_template.radius)+1:node_bottom-(cirroi_template.y-cirroi_template.radius)+1,cirroi_template.radius+1) = true;

            % Update queue
            queue_nodelist(1:2) = [node_top node_bottom];
            queue_nodeindex(1) = cirroi_template.radius;
            length_queue = 2;

            % Form activelines and then inactivate nodepair from nodelist
            % Activelines keeps track of which node pairs have been
            % analyzed
            activelines = true(size(obj.region(num_region+1).nodelist,1),size(obj.region(num_region+1).nodelist,2)/2);   
            activelines((cirroi_template.x-obj.region(num_region+1).leftbound)+1,idx_nodelist/2+1) = false;

            % Enter while loop. Exit when queue is empty        
            while (length_queue ~= 0)
                % STEPS:
                % 1) Load nodepair from queue
                % 2) Update queue
                % 3) Add nodepair to cirroi, sort, and then update
                % 4) Compare nodepair to nodepairs left and right and add nodes to queue
                % 5) Update totalpoints

                % Step 1) Pop queue
                queue_nodelist_buffer(1:2) = queue_nodelist(length_queue-1:length_queue); 
                queue_nodeindex_buffer = queue_nodeindex(length_queue/2);

                % Step 2)
                length_queue = length_queue-2;

                % Step 3)   
                if (cirroi_template.region.noderange(queue_nodeindex_buffer+1) == 0)
                    % Just place directly
                    cirroi_template.region.nodelist(queue_nodeindex_buffer+1,1:2) = queue_nodelist_buffer(1:2);
                else
                    % Merge nodes
                    inserted = false;
                    for i = 0:2:cirroi_template.region.noderange(queue_nodeindex_buffer+1)-1
                        if (queue_nodelist_buffer(2) < cirroi_template.region.nodelist(queue_nodeindex_buffer+1,i+1))
                            nodelist_buffer = cirroi_template.region.nodelist(queue_nodeindex_buffer+1,i+1:cirroi_template.region.noderange(queue_nodeindex_buffer+1));
                            cirroi_template.region.nodelist(queue_nodeindex_buffer+1,i+1:i+2) = queue_nodelist_buffer;
                            cirroi_template.region.nodelist(queue_nodeindex_buffer+1,i+3:cirroi_template.region.noderange(queue_nodeindex_buffer+1)+2) = nodelist_buffer;
                            inserted = true;
                            break
                        end
                    end
                    if (~inserted)
                        % Append at the end
                        cirroi_template.region.nodelist(queue_nodeindex_buffer+1,cirroi_template.region.noderange(queue_nodeindex_buffer+1)+1:cirroi_template.region.noderange(queue_nodeindex_buffer+1)+2) = queue_nodelist_buffer;
                    end
                end
                cirroi_template.region.noderange(queue_nodeindex_buffer+1) = cirroi_template.region.noderange(queue_nodeindex_buffer+1)+2;

                % Step 4)
                idx_froi = queue_nodeindex_buffer+(cirroi_template.x-cirroi_template.radius)-obj.region(num_region+1).leftbound; 

                % Compare to node pairs LEFT
                if (queue_nodeindex_buffer > 0 && idx_froi > 0) % makes sure previous node is within cirroi
                    for i = 0:2:obj.region(num_region+1).noderange(idx_froi)-1 
                        if (activelines(idx_froi,i/2+1) == 0) % Nodes are inactive
                            continue
                        end

                        upperlim = ceil(-sqrt(cirroi_template.radius^2-(((queue_nodeindex_buffer-1)+(cirroi_template.x-cirroi_template.radius))-cirroi_template.x)^2)+cirroi_template.y);
                        lowerlim = floor(sqrt(cirroi_template.radius^2-(((queue_nodeindex_buffer-1)+(cirroi_template.x-cirroi_template.radius))-cirroi_template.x)^2)+cirroi_template.y); 

                        if (obj.region(num_region+1).nodelist(idx_froi,i+2) < queue_nodelist_buffer(1)) % lower node comes before upper node of buffer 
                            continue;
                        elseif (obj.region(num_region+1).nodelist(idx_froi,i+1) <= queue_nodelist_buffer(2) && obj.region(num_region+1).nodelist(idx_froi,i+2) >= queue_nodelist_buffer(1)) % nodes interact
                            if (obj.region(num_region+1).nodelist(idx_froi,i+1) < upperlim) % make sure upper node is within the circle
                                node_top = upperlim; 
                            else
                                node_top = obj.region(num_region+1).nodelist(idx_froi,i+1);
                                circ_untrunc = false;
                            end

                            if (obj.region(num_region+1).nodelist(idx_froi,i+2) > lowerlim) % make sure lower node is within the circle
                                node_bottom = lowerlim; 
                            else
                                node_bottom = obj.region(num_region+1).nodelist(idx_froi,i+2);
                                circ_untrunc = false;
                            end

                            if (node_top > node_bottom || node_top > queue_nodelist_buffer(2) || node_bottom < queue_nodelist_buffer(1)) % This can occur along diagonal boundaries by virtue of how the four way direction contiguous algorithm works. If this happens, break, because it is no longer contiugous
                                continue;
                            end

                            % Add to queue
                            queue_nodelist(length_queue+1:length_queue+2) = [node_top node_bottom]; % add to queue
                            queue_nodeindex(length_queue/2+1) = queue_nodeindex_buffer-1; 
                            length_queue = length_queue + 2;

                            % Make node pair inactive
                            activelines((idx_froi-1)+1,i/2+1) = false;

                            % Update mask
                            cirroi_template.mask(node_top-(cirroi_template.y-cirroi_template.radius)+1:node_bottom-(cirroi_template.y-cirroi_template.radius)+1,queue_nodeindex_buffer) = true;
                        else
                            break;
                        end        
                    end
                end

                % Compare to node pairs RIGHT
                if (queue_nodeindex_buffer < 2*cirroi_template.radius && idx_froi < obj.region(num_region+1).rightbound-obj.region(num_region+1).leftbound) % makes sure next node is within cirroiint
                    for i = 0:2:obj.region(num_region+1).noderange(idx_froi+2)-1 
                        if (activelines(idx_froi+2,i/2+1) == 0)
                            continue
                        end

                        upperlim = ceil(-sqrt(cirroi_template.radius^2-(((queue_nodeindex_buffer+1)+(cirroi_template.x-cirroi_template.radius))-cirroi_template.x)^2)+cirroi_template.y); 
                        lowerlim = floor(sqrt(cirroi_template.radius^2-(((queue_nodeindex_buffer+1)+(cirroi_template.x-cirroi_template.radius))-cirroi_template.x)^2)+cirroi_template.y);

                        if (obj.region(num_region+1).nodelist(idx_froi+2,i+2) < queue_nodelist_buffer(1)) % lower node comes before upper node of buffer
                            continue;
                        elseif (obj.region(num_region+1).nodelist(idx_froi+2,i+1) <= queue_nodelist_buffer(2) && obj.region(num_region+1).nodelist(idx_froi+2,i+2) >= queue_nodelist_buffer(1)) % nodes interact
                            if (obj.region(num_region+1).nodelist(idx_froi+2,i+1) < upperlim) % make sure upper node is within the circle
                                node_top = upperlim; 
                            else
                                node_top = obj.region(num_region+1).nodelist(idx_froi+2,i+1);
                                circ_untrunc = false;
                            end

                            if (obj.region(num_region+1).nodelist(idx_froi+2,i+2) > lowerlim) % make sure lower node is within the circle
                                node_bottom = lowerlim; 
                            else
                                node_bottom = obj.region(num_region+1).nodelist(idx_froi+2,i+2);
                                circ_untrunc = false;
                            end

                            if (node_top > node_bottom || node_top > queue_nodelist_buffer(2) || node_bottom < queue_nodelist_buffer(1)) % This can occur along diagonal boundaries by virtue of how the four way direction contiguous algorithm works. If this happens, break, because it is no longer contiugous
                                continue;
                            end

                            % Add to queue
                            queue_nodelist(length_queue+1:length_queue+2) = [node_top node_bottom]; % add to queue
                            queue_nodeindex(length_queue/2+1) = queue_nodeindex_buffer+1; % 
                            length_queue = length_queue + 2;

                            % Make node pair inactive
                            activelines((idx_froi+1)+1,i/2+1) = false;

                            % Update mask
                            cirroi_template.mask(node_top-(cirroi_template.y-cirroi_template.radius)+1:node_bottom-(cirroi_template.y-cirroi_template.radius)+1,queue_nodeindex_buffer+2) = true;
                        else
                            break;
                        end        
                    end
                end

                % Update totalpoints
                cirroi_template.region.totalpoints = cirroi_template.region.totalpoints + (queue_nodelist_buffer(2)-queue_nodelist_buffer(1)+1);   
            end                                       

            % Set Bounds - since noderange length is always
            % 2*radius+1, set left and right bounds to:
            cirroi_template.region.leftbound = cirroi_template.x - cirroi_template.radius;
            cirroi_template.region.rightbound = cirroi_template.x + cirroi_template.radius;

            % Find upper and lower bounds:
            cirroi_template.region.upperbound = size(obj.mask,1);
            cirroi_template.region.lowerbound = 0;
            for i = 0:size(cirroi_template.region.noderange,1)-1 % Find max and min x-values of cirroi
                if (cirroi_template.region.noderange(i+1) > 0)
                    if (cirroi_template.region.upperbound > cirroi_template.region.nodelist(i+1,1))
                       cirroi_template.region.upperbound = cirroi_template.region.nodelist(i+1,1);
                    end
                    if (cirroi_template.region.lowerbound < cirroi_template.region.nodelist(i+1,cirroi_template.region.noderange(i+1)))
                       cirroi_template.region.lowerbound = cirroi_template.region.nodelist(i+1,cirroi_template.region.noderange(i+1));
                    end
                end
            end

            % If subset is truncated, do some special analysis to construct 
            % the cirroi. At this point, it's only contiguous. I've added 
            % this portion to truncate a part of the subset if this option 
            % is enabled. This is helpful for analyzing cracks so the subset 
            % does not wrap around the crack tip. Also make sure to check 
            % if any of the noderanges are zero, as this will not trigger a 
            % truncation.
            if (subsettrunc && (~circ_untrunc || any(cirroi_template.region.noderange == 0)))   
                % Get boundary
                point_topleft = [-1 -1]; % Initialize
                for i = 0:size(cirroi_template.mask,2)-1
                    for j = 0:size(cirroi_template.mask,1)-1
                        if (cirroi_template.mask(j+1,i+1))
                            point_topleft = [i j]; % This is the first point.
                            break;
                        end
                    end
                    if (cirroi_template.mask(j+1,i+1))
                        break;
                    end
                end                         

                boundary_cirroi = ncorr_alg_formboundary(int32(point_topleft),int32(7),cirroi_template.mask);                        

                % Find closest point -> Find line equation -> Find 
                % number of points to the left and right of the 
                % line; the side with the least amount of points 
                % will be discarded
                dists = (boundary_cirroi(:,1)-cirroi_template.radius).^2 + (boundary_cirroi(:,2)-cirroi_template.radius).^2;
                [val_min,idx_min] = min(dists);
                idx_min = idx_min-1; % Zero based indexing

                % Make sure closest point is not on the boundary,
                % which can happen since the boundary isnt
                % perfectly circular
                if (ceil(sqrt(val_min)) < cirroi_template.radius)
                    % Nonlinear solver
                    idx_space = 3;
                    idx_plus = mod(idx_min+idx_space,size(boundary_cirroi,1));
                    idx_minus = mod(idx_min-idx_space,size(boundary_cirroi,1));

                    % Initialize points to calculate derivatives with
                    x_plus_f = 0;
                    x_min_f = 0; 
                    x_minus_f = 0;
                    y_plus_f = 0;
                    y_min_f = 0; 
                    y_minus_f = 0;
                    % Set filt length
                    length_filt = 2;
                    for i = -length_filt:length_filt 
                        x_plus_f = x_plus_f + boundary_cirroi(mod(idx_plus+i,size(boundary_cirroi,1))+1,1);
                        x_min_f = x_min_f + boundary_cirroi(mod(idx_min+i,size(boundary_cirroi,1))+1,1);
                        x_minus_f = x_minus_f + boundary_cirroi(mod(idx_minus+i,size(boundary_cirroi,1))+1,1);

                        y_plus_f = y_plus_f + boundary_cirroi(mod(idx_plus+i,size(boundary_cirroi,1))+1,2);
                        y_min_f = y_min_f + boundary_cirroi(mod(idx_min+i,size(boundary_cirroi,1))+1,2);
                        y_minus_f = y_minus_f + boundary_cirroi(mod(idx_minus+i,size(boundary_cirroi,1))+1,2);
                    end
                    % Divide by length
                    x_plus_f = x_plus_f/(2*length_filt+1);
                    x_min_f = x_min_f/(2*length_filt+1);
                    x_minus_f = x_minus_f/(2*length_filt+1);
                    y_plus_f = y_plus_f/(2*length_filt+1);
                    y_min_f = y_min_f/(2*length_filt+1); 
                    y_minus_f = y_minus_f/(2*length_filt+1);                        

                    % Get derivatives
                    dx_di_f = (x_plus_f - x_minus_f)/(2*idx_space);
                    d2x_di2_f = (x_plus_f-2*x_min_f+x_minus_f)/(idx_space^2);
                    dy_di_f = (y_plus_f - y_minus_f)/(2*idx_space);
                    d2y_di2_f = (y_plus_f-2*y_min_f+y_minus_f)/(idx_space^2);

                    % Do one iteration to find change in index
                    deltai = ((-x_min_f*dx_di_f+cirroi_template.radius*dx_di_f) + (-y_min_f*dy_di_f+cirroi_template.radius*dy_di_f))/((dx_di_f^2+x_min_f*d2x_di2_f-cirroi_template.radius*d2x_di2_f)+(dy_di_f^2+y_min_f*d2y_di2_f-cirroi_template.radius*d2y_di2_f));

                    % Test deltai, make sure its between subsequent spacings; If
                    % not then set the points minus and plus points      
                    if (abs(deltai) < idx_space) 
                        % Calculate two points based on dx_di_i and dy_di_i
                        x_i = x_min_f + dx_di_f*deltai + (1/2)*d2x_di2_f*deltai^2; 
                        y_i = y_min_f + dy_di_f*deltai + (1/2)*d2y_di2_f*deltai^2; 
                        dx_di_i = dx_di_f + d2x_di2_f*deltai;
                        dy_di_i = dy_di_f + d2y_di2_f*deltai;

                        % Set stepsize - the magnitude doesnt matter
                        stepsize = 0.5;

                        % Get points
                        p0 = [x_i-dx_di_i*stepsize y_i-dy_di_i*stepsize];
                        p1 = [x_i+dx_di_i*stepsize y_i+dy_di_i*stepsize];
                    else
                        % Set points to averaged points on either side
                        p0 = [x_minus_f y_minus_f];
                        p1 = [x_plus_f y_plus_f];                            
                    end

                    % Find which side to clear: p_subset is a point defined 
                    % to be part of the subset, so the side opposite of 
                    % p_subset needs to be cleared. If the center of the 
                    % subset does not lie on boundary, then find displacement 
                    % from closest point to p0. Add this displacement to 
                    % the center, and then determine which side the center 
                    % is on.
                    p_subset = [0 0]; % Initialize
                    if (isequal(boundary_cirroi(idx_min+1,:),[cirroi_template.radius cirroi_template.radius]))
                        % Center is the closest point. Find valid points
                        % around the boundary and get the centroid of 
                        % these points.
                        width_win = 1; % This determines window of points collected
                        counter = 0; % Counts number of points added so a proper average can be taken
                        for i = -width_win:width_win
                            x_mask = boundary_cirroi(idx_min+1,1)+i;
                            for j = -width_win:width_win
                                y_mask = boundary_cirroi(idx_min+1,2)+j;

                                % Make sure points are within mask and
                                % that the mask is valid here.
                                if (x_mask >= 0 && x_mask <= 2*cirroi_template.radius && y_mask >= 0 && y_mask <= 2*cirroi_template.radius && ...
                                    cirroi_template.mask(y_mask+1,x_mask+1))
                                    p_subset = p_subset + [x_mask y_mask];
                                    counter = counter+1;
                                end
                            end
                        end                           
                        % Divide by counter to get the average
                        p_subset = p_subset/counter;

                        % Debug ------%
                        %{
                        figure(1); imshow(cirroi_template.mask); hold on;
                        plot(p0(1,1),p0(1,2),'rx');
                        plot(p1(1,1),p1(1,2),'rx');
                        plot(p_subset(1,1),p_subset(1,2),'gx');
                        hold off;
                        %}
                        % ------------%                                                                            
                    else     
                        % Set the center point
                        p_subset = [cirroi_template.radius cirroi_template.radius];

                        % Debug ------%
                        %{
                        figure(1); imshow(cirroi_template.mask); hold on;
                        plot(p0(1,1),p0(1,2),'rx');
                        plot(p1(1,1),p1(1,2),'rx');
                        plot(p_subset(1,1),p_subset(1,2),'gx');
                        hold off;
                        %}
                        % ------------%   
                    end           

                    % Normalize point against line since it can shift
                    disp_p0 = p0 - boundary_cirroi(idx_min+1,:); 
                    p_subset = p_subset + disp_p0;

                    % Now determine which side p_subset lies and
                    % clear the other side
                    sign_clear = -sign((p1(1)-p0(1))*(p_subset(2)-p0(2))-(p_subset(1)-p0(1))*(p1(2)-p0(2)));  
                    mask_cirroi_prelim = cirroi_template.mask; % Make a copy
                    for i = 0:size(cirroi_template.region.noderange,1)-1
                        for j = 0:2:cirroi_template.region.noderange(i+1)-1
                            for k = cirroi_template.region.nodelist(i+1,j+1):cirroi_template.region.nodelist(i+1,j+2)
                                p2 = [i k-(cirroi_template.y-cirroi_template.radius)];
                                if (sign((p1(1)-p0(1))*(p2(2)-p0(2))-(p2(1)-p0(1))*(p1(2)-p0(2))) == sign_clear)
                                    % Clear points on this side
                                    mask_cirroi_prelim(p2(2)+1,p2(1)+1) = false;                                        
                                end
                            end
                        end
                    end

                    % Get region again from the mask, then select the
                    % biggest region. This is required because subset 
                    % truncation does not preserve the cirroi's contiguous 
                    % nature. Setting the last parameter to true for 
                    % formregion will make the length of the noderange the 
                    % same as the width of the mask
                    [region_cirroi_prelim,removed] = ncorr_alg_formregions(mask_cirroi_prelim,int32(0),true); %#ok<ASGLU>

                    % Its possible for region to be empty. Must check it. If 
                    % its empty, then abandon subset truncation.
                    if (~isempty(region_cirroi_prelim))
                        % Only keep the region with the max number of
                        % totalpoints
                        region_cirroi_prelim = region_cirroi_prelim([region_cirroi_prelim.totalpoints] == max([region_cirroi_prelim.totalpoints]));

                        % Update mask then store it
                        mask_cirroi_prelim(:) = false;
                        for i = 0:size(cirroi_template.region.noderange,1)-1
                            for j = 0:2:region_cirroi_prelim.noderange(i+1)-1
                                for k = region_cirroi_prelim.nodelist(i+1,j+1):region_cirroi_prelim.nodelist(i+1,j+2)
                                    mask_cirroi_prelim(k+1,i+1) = true;   
                                end
                            end
                        end

                        % Store mask
                        cirroi_template.mask = mask_cirroi_prelim;

                        % Store region with greatest number of nodes
                        cirroi_template.region = region_cirroi_prelim;

                        % Bring back to regular coords
                        for i = 0:size(cirroi_template.region.noderange,1)-1
                            if (cirroi_template.region.noderange(i+1) > 0)
                                cirroi_template.region.nodelist(i+1,1:cirroi_template.region.noderange(i+1)) = cirroi_template.region.nodelist(i+1,1:cirroi_template.region.noderange(i+1)) + (cirroi_template.y-cirroi_template.radius);
                            end
                        end
                        cirroi_template.region.leftbound = cirroi_template.region.leftbound + (cirroi_template.x-cirroi_template.radius);
                        cirroi_template.region.rightbound = cirroi_template.region.rightbound + (cirroi_template.x-cirroi_template.radius);
                        cirroi_template.region.upperbound = cirroi_template.region.upperbound + (cirroi_template.y-cirroi_template.radius);
                        cirroi_template.region.lowerbound = cirroi_template.region.lowerbound + (cirroi_template.y-cirroi_template.radius);                        
                    end


                    % Debug
                    %{
                    figure(2), imshow(cirroi_template.mask);
                    hold on;
                    plot(boundary_cirroi(idx_min+1,1)+1,boundary_cirroi(idx_min+1,2)+1,'bx');
                    plot(p0(1)+1,p0(2)+1,'gs');
                    plot(p1(1)+1,p1(2)+1,'gs');
                    hold off;
                    %}
                end
            end         
            %-----------------------------------------------------%
            %-----------------------------------------------------%
            %-----------------------------------------------------%

            % Assign outputs
            cirroi(1) = cirroi_template;
        end   
                
        %-----------------------------------------------------------------%
        % For debugging --------------------------------------------------%
        %-----------------------------------------------------------------%
        function debug_region(obj)
            if (isempty(obj.type))  
                error('ROI has not been set yet');
            end
            
            mask_debug = false(size(obj.mask));
            for i = 0:length(obj.region)-1
                regionmask_debug = false(size(obj.mask));
                for j = 0:size(obj.region(i+1).noderange,1)-1
                    x = j + obj.region(i+1).leftbound;
                    for k = 0:2:obj.region(i+1).noderange(j+1)-1
                        vec_y = obj.region(i+1).nodelist(j+1,k+1):obj.region(i+1).nodelist(j+1,k+2);
                        mask_debug(vec_y+1,x+1) = true;   
                        regionmask_debug(vec_y+1,x+1) = true; 
                    end
                end
                % figure, imshow(regionmask_debug); - This can potential be a lot of figures
            end
            figure, imshow(mask_debug);
            if (isequal(mask_debug,obj.mask))
                disp('SUCCESS - mask matches region');
            else
                error('FAILURE - for some reason mask and region are different');
            end
        end 
        
        function debug_boundary(obj)
            if (isempty(obj.type))  
                error('ROI has not been set yet');
            end
            
            figure, imshow(obj.mask); 
            hold on;
            for i = 0:length(obj.boundary)-1
                plot(obj.boundary(i+1).add(:,1)+1,obj.boundary(i+1).add(:,2)+1,'go');
                for j = 0:length(obj.boundary(i+1).sub)-1
                    plot(obj.boundary(i+1).sub{j+1}(:,1)+1,obj.boundary(i+1).sub{j+1}(:,2)+1,'rx');
                end
            end
            hold off;
        end
    end
end

function boundary_update = interp_border(boundary,plot_u_interp,plot_v_interp,roi,num_region,size_mask_update,border_interp,spacing,radius)
% Updates boundary based on old boundary and displacement plots. 
% Does not preserve length of the boundary. If points are found to travel 
% outside a radius length around the size_mask_update, they are discarded.
%
% Inputs -----------------------------------------------------------------%
%   boundary - double array; boundary we're trying to update - uses pixel
%   units.
%   plot_u_interp - double array; b-spline coefficients of u displacement
%   field
%   plot_v_interp - double array; b-spline coefficients of v displacement
%   field
%   roi - ncorr_class_roi; ROI corresponding to displacement fields. Note
%   that roi is reduced.
%   num_region - integer; number corresponding to region
%   size_mask_update - integer array; size of the new mask (can be different 
%   from the mask in ROI if different sized current images are used).
%   border_interp - integer; border around b-spline coefficients
%   spacing - integer; spacing parameter
%   radius - integer; subset radius
%
% Outputs ----------------------------------------------------------------%
%   boundary_update - double array; updated boundary coordinates
% 
% Note, if boundary_update is found to be empty after analysis, it will 
% return [-1 -1] so it is not empty. This is a temporary workaround because
% the ncorr_alg_formmask function does not support empty polygons.
        
    % Get scaled coordinates since displacement plots are reduced.
    boundary_scaled = boundary/(spacing+1);
    
    % Get interpolated u displacements - can return NaN if out of range
    interp_u = ncorr_alg_interpqbs(boundary_scaled,plot_u_interp,roi.region(num_region+1).leftbound,roi.region(num_region+1).upperbound,border_interp);
    
    % Get interpolated v displacements - can return NaN if out of range
    interp_v = ncorr_alg_interpqbs(boundary_scaled,plot_v_interp,roi.region(num_region+1).leftbound,roi.region(num_region+1).upperbound,border_interp);
    
    % Update boundary
    boundary_update = boundary;
    boundary_update(:,1) = boundary_update(:,1) + interp_u;
    boundary_update(:,2) = boundary_update(:,2) + interp_v;
        
    % Set points near boundary to NaN;
    for i = 0:size(boundary_update,1)-1
        if (boundary_update(i+1,1) < radius || boundary_update(i+1,1) > size_mask_update(2)-radius || ...
            boundary_update(i+1,2) < radius || boundary_update(i+1,2) > size_mask_update(1)-radius)
            boundary_update(i+1,:) = NaN;
        end
    end
        
    % Remove all coordinates that have NaN
    counter = 0;
    while (counter < size(boundary_update,1))
        if (any(isnan(boundary_update(counter+1,:))))
            boundary_update(counter+1,:) = [];
            counter = counter-1;
        end
        counter = counter+1;
    end
    
    % Check to see if boundary is empty
    if (isempty(boundary_update))
        boundary_update = [-1 -1];
    end
end

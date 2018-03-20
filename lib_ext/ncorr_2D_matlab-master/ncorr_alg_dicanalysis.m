function [displacements,rois_dic,seedinfo,outstate] = ncorr_alg_dicanalysis(imgs,radius,spacing,cutoff_diffnorm,cutoff_iteration,total_threads,enabled_stepanalysis,subsettrunc,num_img,total_imgs,pos_parent,params_init)
% This function performs RG-DIC.
%
% Inputs -----------------------------------------------------------------%
%   imgs - struct; Contains struct('imginfo',{},'roi',{}). Contains the 
%   reference and current images packaged together. The imginfo field
%   contains an ncorr_class_img and the ROI field contains an ncorr_class_roi.
%   radius - integer; subset radius
%   spacing - integer; subset spacing
%   cutoff_diffnorm - double; cutoff for the norm of difference vector
%   cutoff_iteration - integer; cutoff for the number of IC-GN iterations
%   total_threads - integer; total number of threads
%   enabled_stepanalysis - logical; if true, then process as many seeds as
%   possible. If false, process all the seeds.
%   subsettrunc - logical; if true, then subset truncation is enabled
%   num_img - integer; reference image number
%   total_imgs - integer; total number of images
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%   params_init - struct; contains struct('paramvector',{},'num_region',{},
%   'num_thread',{},'computepoints',{}). If its not empty, it contains last 
%   set of seeds used for the previous iteration
%
% Outputs ----------------------------------------------------------------%
%   displacements - struct; contains
%   struct('plot_u',{},'plot_v',{},'plot_corrcoef',{},'plot_validpoints',{}) 
%   rois_dic - ncorr_class_img; ROIs after being unioned with reference roi 
%   and validpoints from DIC analysis.
%   seedinfo - struct; contains struct('paramvector',{},'num_region',{},
%   'num_thread',{},'computepoints',{})
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.
%
% outstate will only return failed if ncorr_alg_rgdic throws an exception.

    % Initialize outputs
    outstate = out.cancelled;
    displacements = struct('plot_u',{},'plot_v',{},'plot_corrcoef',{},'plot_validpoints',{});    
    rois_dic = ncorr_class_roi.empty;
    seedinfo = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{});
    
    % --------------------------------------------------------------------%
    % Get Seeds ----------------------------------------------------------%
    % --------------------------------------------------------------------%
    
    if (isempty(params_init))
        % Get new seeds manually through GUI. outstate_seeds will either be
        % success or cancelled.
        [seedinfo_prelim,threaddiagram,outstate_seeds] = ncorr_gui_seedanalysis(imgs(1).imginfo, ...
                                                                                [imgs(2:end).imginfo], ...
                                                                                imgs(1).roi, ...
                                                                                radius, ...
                                                                                spacing, ...
                                                                                cutoff_diffnorm, ...
                                                                                cutoff_iteration, ...
                                                                                total_threads, ...
                                                                                enabled_stepanalysis, ...
                                                                                subsettrunc, ...
                                                                                num_img, ...
                                                                                total_imgs, ...
                                                                                pos_parent);
        % See if analysis was cancelled
        if (outstate_seeds ~= out.success)
            return;
        end            
    else
        % Use params init to place seeds directly               
        % Initialize seedinfo_prelim
        seedinfo_prelim = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{});
            
        % Initialize Thread diagram buffers - some buffers are unused
        ref_reduced = imgs(1).imginfo.reduce(spacing); % Unused
        roi_reduced = imgs(1).roi.reduce(spacing); % Unused
        preview_threaddiagram = zeros(size(roi_reduced.mask)); % Unused
        threaddiagram_buffer = -ones(size(roi_reduced.mask));
        threaddiagram = -ones(size(roi_reduced.mask)); % This must be negative 1
        
        % Calculate seeds - cycle over regions
        num_imgs_success = inf; % Initialize to impossibly high value
        manualseed = false; % Set this parameter true if there is trouble. This causes the user to reselect the seeds manually through a GUI.
        % Cycle over regions
        for i = 0:size(params_init,2)-1 
            % Get updated seed positions based on seed locations of the
            % last successfully seeded image for this region
            % Initialize
            pos_seed = zeros(size(params_init,1),2);
            % Cycle over threads
            for j = 0:size(params_init,1)-1
                % Must round displacements so that the updated position
                % lies properly on "spaced grid" - THIS IS IMPORTANT!!!
                pos_seed(j+1,1) = params_init(j+1,i+1,end).paramvector(1) + round(params_init(j+1,i+1,end).paramvector(3)/(spacing+1))*(spacing+1);
                pos_seed(j+1,2) = params_init(j+1,i+1,end).paramvector(2) + round(params_init(j+1,i+1,end).paramvector(4)/(spacing+1))*(spacing+1);
            end
            
            % Get num_region - just use seed from the first thread
            num_region = params_init(1,i+1,end).num_region;
            
            % Make sure seed positions are within the ROI and are unique.
            % If they arent then prompt the user and have them replace the
            % seeds manually.
            regionmask = roi_reduced.get_regionmask(num_region); 
            if (size(pos_seed,1) == size(unique(pos_seed,'rows'),1) && ...
                all(pos_seed(:,1)./(spacing+1) >= 0) && all(pos_seed(:,1)./(spacing+1) < size(regionmask,2)) && ...
                all(pos_seed(:,2)./(spacing+1) >= 0) && all(pos_seed(:,2)./(spacing+1) < size(regionmask,1)) && ...
                all(regionmask(sub2ind(size(regionmask),pos_seed(:,2)/(spacing+1)+1,pos_seed(:,1)/(spacing+1)+1))))
                % Get seeds - outstate will be successful if at least one
                % image is seeded
                [seedinfo_buffer,convergence_buffer,outstate_seeds] = ncorr_alg_seedanalysis(imgs(1).imginfo, ...
                                                                                             [imgs(2:end).imginfo], ...
                                                                                             imgs(1).roi, ...
                                                                                             num_region, ...
                                                                                             pos_seed, ...
                                                                                             radius, ...
                                                                                             cutoff_diffnorm, ...
                                                                                             cutoff_iteration, ...
                                                                                             enabled_stepanalysis, ...
                                                                                             subsettrunc, ...
                                                                                             num_img, ...
                                                                                             total_imgs); %#ok<ASGLU>
                % See if analysis was cancelled
                if (outstate_seeds == out.cancelled)
                    return;
                end
                
                % See if analysis was successful or failed
                if (outstate_seeds == out.success)     
                    % Form thread diagram for these seeds
                    ncorr_alg_formthreaddiagram(threaddiagram_buffer,preview_threaddiagram,int32(pos_seed/(spacing+1)),regionmask,ref_reduced.formatted());

                    % Get compute points
                    for j = 0:size(seedinfo_buffer,3)-1
                        for k = 0:size(seedinfo_buffer,1)-1
                            seedinfo_buffer(k+1,1,j+1).computepoints = length(find(threaddiagram_buffer == k));
                        end
                    end    
                else
                    % Seed analysis failed. Alert user and then ask him/her
                    % to manually place seeds.
                    h_error = errordlg('Not a single image was seeded correctly; please replace seeds manually.','Error','modal');
                    uiwait(h_error);
                    
                    manualseed = true;
                end
            else
                % Alert user
                h_error = errordlg('One or more seeds went outside of the ROI **OR** converged on top of each other, please replace manually. If this happens often then dont place seeds near each other or near the boundary.','Error','modal');
                uiwait(h_error);
                manualseed = true;
            end            
            
            if (~manualseed)            
                % Take minimum of num_imgs_success and
                % seedinfo. Buffer can be more or less than
                % num_imgs_success
                num_imgs_success = min(num_imgs_success, size(seedinfo_buffer,3));

                % Clear out other images in buffer and prelim
                if (~isempty(seedinfo_prelim))
                    seedinfo_prelim = seedinfo_prelim(:,:,1:num_imgs_success);
                end
                seedinfo_buffer = seedinfo_buffer(:,:,1:num_imgs_success);

                % Append seedinfo_buffer - append along 2nd dimension 
                seedinfo_prelim = horzcat(seedinfo_prelim,seedinfo_buffer); %#ok<AGROW>  

                % Merge threaddiagram from previous iteration
                threaddiagram(threaddiagram_buffer ~= -1) = threaddiagram_buffer(threaddiagram_buffer ~= -1);
            else
                % If there's a problem, just ask user to manually reset
                % seeds. Outstate will either be success or cancelled.
                [seedinfo_prelim,threaddiagram,outstate_seeds] = ncorr_gui_seedanalysis(imgs(1).imginfo, ...
                                                                                        [imgs(2:end).imginfo], ...
                                                                                        imgs(1).roi, ...
                                                                                        radius, ...
                                                                                        spacing, ...
                                                                                        cutoff_diffnorm, ...
                                                                                        cutoff_iteration, ...
                                                                                        total_threads, ...
                                                                                        enabled_stepanalysis, ...
                                                                                        subsettrunc, ...
                                                                                        num_img, ...
                                                                                        total_imgs, ...
                                                                                        pos_parent);
                                                                                    
                % Check if analysis was cancelled
                if (outstate_seeds ~= out.success)
                    return;
                end
                
                % Break since seeds for all regions are placed in GUI.
                break;
            end
        end
        
        % Debug ----------------------------------------------------------%
        %{
        figure, imshow(threaddiagram,[]); hold on;
        for i = 0:size(seedinfo_prelim,2)-1
            for j = 0:size(seedinfo_prelim,1)-1
                plot(seedinfo_prelim(j+1,i+1,1).paramvector(1)/(spacing+1),seedinfo_prelim(j+1,i+1,1).paramvector(2)/(spacing+1),'ro');
            end
        end
        hold off;
        %}
        % ----------------------------------------------------------------%
    end                                                                     
                                        
    % --------------------------------------------------------------------%
    % Perform DIC Analysis -----------------------------------------------%
    % --------------------------------------------------------------------%
    
    % Format Seeds ---------------------------------------------------%
    seedinfo_prelim_f = seedinfo_prelim;
    for i = 0:size(seedinfo_prelim,1)-1
        for j = 0:size(seedinfo_prelim,2)-1
            for k = 0:size(seedinfo_prelim,3)-1
                seedinfo_prelim_f(i+1,j+1,k+1).num_region = int32(seedinfo_prelim(i+1,j+1,k+1).num_region);
                seedinfo_prelim_f(i+1,j+1,k+1).num_thread = int32(seedinfo_prelim(i+1,j+1,k+1).num_thread);
                seedinfo_prelim_f(i+1,j+1,k+1).computepoints = int32(seedinfo_prelim(i+1,j+1,k+1).computepoints);
            end
        end
    end    

    % Begin DIC analysis ---------------------------------------------%
    displacements_prelim = struct('plot_u',{},'plot_v',{},'plot_corrcoef',{},'plot_validpoints',{});    
    rois_dic_prelim = ncorr_class_roi.empty;
    for i = 0:size(seedinfo_prelim,3)-1
        % Put in try block because DIC can run
        % out of memory                                
        try
            tic
            % No need to send in the number of regions or total number
            % of threads because this is encoded in the size of the
            % seeds. 
            % FYI: seeinfo_prelim(i,j,k) where i refers to the thread
            % number, j refers to the region, and k refers to the image.
            % ncorr_alg_rgdic will either return success or cancelled or
            % return an exception.
            [displacements_prelim(i+1),outstate_dic] = ncorr_alg_rgdic(imgs(1).imginfo.formatted(), ...
                                                                       imgs(i+2).imginfo.formatted(), ...
                                                                       imgs(1).roi.formatted(), ...
                                                                       seedinfo_prelim_f(:,:,i+1), ...
                                                                       int32(threaddiagram), ...
                                                                       int32(radius), ...
                                                                       int32(spacing), ...
                                                                       cutoff_diffnorm, ...
                                                                       int32(cutoff_iteration), ...
                                                                       logical(subsettrunc), ...
                                                                       int32(num_img+i), ...
                                                                       int32(total_imgs));       
            toc
        catch %#ok<CTCH>
            % Only time an exception is returned is either if rgdic has
            % incorrect inputs, ran out of memory, or there is a bug in 
            % the code. The first option never happens if rgdic
            % is called through the program. The second option can 
            % happen and is thus the only real handleable exception 
            % which will get thrown by rgdic. The last case was not
            % intended, so it is not handled here.
            h_error = errordlg('Ncorr most likely ran out of memory while performing DIC. Please clear memory in workspace, restart Ncorr, or crop/use smaller images before doing analysis.','Error','modal');
            uiwait(h_error);
            return;
        end  

        % See if analysis was cancelled
        if (outstate_dic ~= out.success)
            return;
        end

        % Take union of reference ROI with validpoints to get the new
        % ROI.
        rois_dic_prelim(i+1) = imgs(1).roi.get_union(displacements_prelim(i+1).plot_validpoints,spacing);
    end

    % Set outputs
    for i = 0:length(displacements_prelim)-1
        displacements(i+1) = displacements_prelim(i+1);
        rois_dic(i+1) = rois_dic_prelim(i+1);
    end
    for i = 0:size(seedinfo_prelim,1)-1
        for j = 0:size(seedinfo_prelim,2)-1
            for k = 0:size(seedinfo_prelim,3)-1
                seedinfo(i+1,j+1,k+1) = seedinfo_prelim(i+1,j+1,k+1);
            end
        end
    end
    outstate = out.success;
end

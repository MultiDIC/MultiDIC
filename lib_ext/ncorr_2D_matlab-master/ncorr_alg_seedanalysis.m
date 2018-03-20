function [seedinfo,convergence,outstate] = ncorr_alg_seedanalysis(reference,current,roi,num_region,pos_seed,radius,cutoff_diffnorm,cutoff_iteration,enabled_stepanalysis,subsettrunc,num_img,total_imgs)
% This function attempts to calculate seeds for the region specified. 
%
% Inputs -----------------------------------------------------------------%
%   reference - ncorr_class_img; used for calculations.
%   current - ncorr_class_img; used for calculations.
%   roi - ncorr_class_roi; ROI corresponding to the reference image.
%   num_region - integer; number corresponding to the region.
%   pos_seed - integer array; 1st column is the x-position of the seed; 2nd
%   is the y-position. Each row corresponds to a thread in ascending order
%   (i.e. the first row is the first thread).
%   radius - integer; radius of subset
%   cutoff_diffnorm - double; cutoff of norm of the difference vector
%   cutoff_iteration - integer; cutoff for number of iterations
%   enabled_stepanalysis - logical; if true, then process as many seeds as
%   possible. If false, process all the seeds.
%   subsettrunc - logical; if true, then enabled subset truncation
%   num_img - integer; number of reference image being analyzed
%   total_imgs - integer; total number of images being analyzed
%
% Outputs ----------------------------------------------------------------%
%   seedinfo - struct; contains
%   struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{})
%   convergence - struct; contains
%   struct('num_iterations',{},'diffnorm',{})
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.
% 
% Note that if step analysis is enabled, this function will return at least
% one seed if outstate is set to success. If step analysis is disabled then
% all images will be seeded if outstate is set to success.

    % Initialze outputs      
    outstate = out.cancelled;
    seedinfo = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{}); % Contains seed information for the RG-DIC analysis
    convergence = struct('num_iterations',{},'diffnorm',{});                                % Contains convergence info  

    % Prepare buffers
    seedinfo_prelim = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{}); % Contains seed information for the RG-DIC analysis
    convergence_prelim = struct('num_iterations',{},'diffnorm',{});                                % Contains convergence info
    
    % Set up waitbar
    h = waitbar((num_img)/total_imgs,['Processing seed(s) for image ' num2str(num_img+1) ' of ' num2str(total_imgs) '...'],'Name','Calculating...','WindowStyle','modal');
    set(h,'CloseRequestFcn','setappdata(gcbf,''canceling'',1)');
    setappdata(h,'canceling',0);
        
    % Iterate over current images and place seeds
    tic
    for i = 0:length(current)-1
        % Calculate seeds for this image and region. Note that for
        % ncorr_alg_calcseeds, outstate will either return success if all
        % seeds are successfully calculated for this image/region, or it 
        % will return failed.
        [seedinfo_buffer,convergence_buffer,outstate_seeds] = ncorr_alg_calcseeds(reference.formatted(), ...
                                                                                  current(i+1).formatted(), ...
                                                                                  roi.formatted(), ...
                                                                                  int32(num_region), ...
                                                                                  int32(pos_seed), ...
                                                                                  int32(radius), ...
                                                                                  cutoff_diffnorm, ...
                                                                                  int32(cutoff_iteration), ...
                                                                                  logical(enabled_stepanalysis), ...
                                                                                  logical(subsettrunc));      
                
        if (enabled_stepanalysis)
            % If seed placement fails or if the diffnorm or corrcoef for any 
            % of the seeds in this image and region are greater than the 
            % absolute cutoff or if this is the 2nd+ image and the number 
            % of iterations is saturated for any of the seeds, then break. 
            % Choose the 2nd image because if the number of iterations is
            % saturated, that doesnt necessarily mean its analyzed
            % incorrectly, but it means that the match is probably poor.        
            
            % Set absolute cutoff for diffnorm and corrcoef
            cutoff_max_diffnorm = 0.1;
            cutoff_max_corrcoef = 0.5;
            
            corrcoef_buffer = vertcat(seedinfo_buffer.paramvector);
            corrcoef_buffer = corrcoef_buffer(:,9);
            if (outstate_seeds ~= out.success || ...
                any([convergence_buffer.diffnorm] > cutoff_max_diffnorm) || ...
                any(corrcoef_buffer > cutoff_max_corrcoef) || ...
                (i > 0 && any([convergence_buffer.num_iterations] == cutoff_iteration)))
                % Check if the number of seeds is zero for step analysis; 
                % if so, then fail the analysis
                if (isempty(seedinfo_prelim))
                    outstate = out.failed;
                end
                
                break;
            end
        else
            % Step analysis isn't enabled. If any seed placements fail then
            % fail the whole analysis
            if (outstate_seeds ~= out.success)
                h_error = errordlg('Some seeds could not be analyzed. If high strain is anticipated, then try enabling the high strain step analysis.','Error','modal');
                uiwait(h_error);
                
                outstate = out.failed;
                break;
            end
        end
        
        % Append seeds
        seedinfo_prelim(:,1,i+1) = seedinfo_buffer;
        convergence_prelim(:,1,i+1) = convergence_buffer;
        
        % See if analysis was cancelled by user
        if (getappdata(h,'canceling'))
            delete(h);         
            % Exit
            return;
        end

        % Update waitbar
        waitbar((num_img+i+1)/total_imgs,h,['Processing seed(s) for image ' num2str(num_img+i+1) ' of ' num2str(total_imgs) '...']);    
    end        
    toc
    
    % Close wait bar
    delete(h);
    
    if (outstate ~= out.failed)
        % Assign outputs
        seedinfo = seedinfo_prelim;
        convergence = convergence_prelim;
        outstate = out.success;
    end
end
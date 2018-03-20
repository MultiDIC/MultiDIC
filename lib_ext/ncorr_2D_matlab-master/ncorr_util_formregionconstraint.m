function handle_constraintfcn = ncorr_util_formregionconstraint(region)
% This function returns handle_constraintfcn which, given a point, will
% return the closest point inside the region. This is mainly used for the 
% constrainfcn for impoints which are restricted to be inside a region.
%
% Inputs -----------------------------------------------------------------%
%   region - struct; region info used to determine the position of returned
%   point
%
% Outputs ----------------------------------------------------------------%
%   handle_constraintfcn - function handle; 
%
% Returns error if region is empty.

    if (region.totalpoints == 0)   
        error('Region cannot be empty when forming constraint function.');  
    end
        
    handle_constraintfcn = @constraintfcn;

    function pos_nearest = constraintfcn(pos)     
    % This function will generally be called within a callback
    % automatically which will use 1 based indexing, so pos_nearest is
    % returned using 1 based indexing. However, for processing, 0 based
    % indexing is used, so convert back and forth as necessary.
    
        % Convert pos to 0 based indexing and round it:
        pos = round(pos)-1;        
                       
        % Find bounding box - must find this since left and right bounds 
        % are based on the length of the noderange/nodelist which may not
        % necessarily be the bounding box.
        left = region.rightbound;
        right = 0;        
        firstpoint = false;
        for i = 0:size(region.noderange,1)-1
            if (region.noderange(i+1) > 0 && ~firstpoint)
                firstpoint = true;
                left = i + region.leftbound;
            end

            if (region.noderange(i+1) > 0 && i+region.leftbound > right)
                right = i + region.leftbound;
            end
        end

        % Check if above, between, or beneath bounds
        if (pos(1) < left) 
            % Left of bounds - set x to left and then find the nearest node
            % for y
            pos_nearest(1) = left+1;
            pos_nearest(2) = nearestnode(region.nodelist(left-region.leftbound+1,:),region.noderange(left-region.leftbound+1),pos(2))+1;
        elseif (pos(1) >= left && pos(1) <= right) 
            % Between bounds - must check if nodes exist at this x position
            if (region.noderange(pos(1)-region.leftbound+1) > 0)
                pos_nearest(1) = pos(1)+1;                    
            else
                % Initialize
                leftfound = false;
                leftdist = 0;
                rightfound = false;
                rightdist = 0;

                % Cycle left until finding a positive noderange
                for i = pos(1)-1:-1:left
                    if (region.noderange(i-region.leftbound+1) > 0)
                        leftdist = pos(1)-i;
                        leftfound = true;
                        break;
                    end
                end

                % Cycle right until finding a positive noderange
                for i = pos(1)+1:right
                    if (region.noderange(i-region.leftbound+1) > 0)
                        rightdist = i-pos(1);
                        rightfound = true;
                        break;
                    end
                end

                % Now set x position
                if (leftfound && rightfound)
                    if (leftdist < rightdist)
                        pos_nearest(1) = pos(1)-leftdist+1;
                    else
                        pos_nearest(1) = pos(1)+rightdist+1;
                    end
                elseif (leftfound)
                    pos_nearest(1) = pos(1)-leftdist+1;
                else
                    pos_nearest(1) = pos(1)+rightdist+1;
                end                   
            end           

            pos_nearest(2) = nearestnode(region.nodelist((pos_nearest(1)-1)-region.leftbound+1,:),region.noderange((pos_nearest(1)-1)-region.leftbound+1),pos(2))+1;  
        else
            % Right of bounds - set x to left and then find the nearest node
            % for y
            pos_nearest(1) = right+1;
            pos_nearest(2) = nearestnode(region.nodelist(right-region.leftbound+1,:),region.noderange(right-region.leftbound+1),pos(2))+1;
        end
    end
end

function y_nearest = nearestnode(nodelist,noderange,y_coord)    
% Finds the nearest y value given the nodelist, noderange and y_coord. Note
% that nodelist and noderange are a single column for a specific x position. 
% Uses zero based indexing.

    % Initialize    
    y_nearest = -1;
    for i = 0:2:noderange-1
        % Test if node is before, within, or after a node pair
        if (y_coord < nodelist(i+1)) 
            if (i == 0) 
                % Above first node
                y_nearest = nodelist(i+1);
                break;
            else 
                % Above inner node pair
                % Find out if y_coord is closer to bottom node of previous
                % pair or top node of current pair
                if ((y_coord-nodelist(i)) < (nodelist(i+1)-y_coord)) 
                    % y_coord is closer to bottom node of previous pair
                    y_nearest = nodelist(i);
                    break;
                else 
                    % y_coord is closer to top node of current pair
                    y_nearest = nodelist(i+1);
                    break;
                end
            end            
        elseif (y_coord >= nodelist(i+1) && y_coord <= nodelist(i+2)) 
            % y_coord is between node pairs, so just return it
            y_nearest = y_coord;
            break;
        else
            % y_coord is after noderange - test to see if this is the last
            % node pair.
            if (i == noderange-2) 
                % This is last node pair, so y_coord is the below node of
                % last node pair
                y_nearest = nodelist(i+2);
                break;
            else 
                % Skip to next node pair
                continue;
            end
        end       
    end
end

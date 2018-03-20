function pos_child = ncorr_util_figpos(pos_parent,size_child)
% This function receives the position of the parent, size of the child, and
% then returns the position of the child. pos_parent should be obtained through 
% the outerposition property. Units for inputs and calculations are characters.
%
% Inputs -----------------------------------------------------------------%
%   pos_parent - integer array; position of parent figure - 
%   [left bottom width height].
%   size_child - integer array; size of the child figure - [height width];
%
% Outputs ----------------------------------------------------------------%
%   pos_child - integer array; position of the child figure

    % Get size of screen in character units
    set(0,'units','characters'); % Set this every time incase user changes units of root
    pos_screen = get(0,'screensize');

    % Get child figure position - Make it an offset from the parent
    offset_x = 5; 
    offset_y = 6;
    
    pos_child(1) = pos_parent(1)+offset_x;
    pos_child(2) = pos_parent(2)+pos_parent(4)-size_child(1)-offset_y;
    pos_child(3) = size_child(2);
    pos_child(4) = size_child(1);

    % Check if right side extends beyond the screen
    % Only check right and bottom since figures are shifted to the bottom
    % right of their parent (assuming parent is on screen)
    right_screen = 10; % Cushion from right of the screen
    bottom_screen = 6; % Cushion from bottom of the screen
    
    % Right
    if (pos_child(1)+pos_child(3) > pos_screen(3)-right_screen)
        % Set so right side touches side of the screen minus the offset
        pos_child(1) = pos_screen(3)-size_child(2)-right_screen;
    end

    % Bottom
    if (pos_child(2) < bottom_screen)
        % Set so bottom side touches bottom of the screen plus the offset
        pos_child(2) = bottom_screen;
    end        
end

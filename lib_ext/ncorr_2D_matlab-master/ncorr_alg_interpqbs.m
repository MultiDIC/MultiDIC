function vec_interp = ncorr_alg_interpqbs(coords,plot_bcoef,offset_x,offset_y,border_bcoef)
% This function performs biquintic interpolation.
% 
% Inputs -----------------------------------------------------------------%
%   coords - double array; nx2 array of coordinates. In the form of [x y]
%   plot_bcoef - double array; plot of bspline coefficients
%   offset_x - integer; x offset of bspline coef plot from the origin
%   offset_y - integer; y offset of bspline coef plot from the origin
%   border_bcoef - integer; border used when interpolating b-spline
%   coefficients. A border is usually applied to help mitigate the effects
%   of ringing near the edges since the bspline coefficient array is formed
%   using an FFT.
%
% Outputs ----------------------------------------------------------------%
%   vec_interp - vector of interpolated values the same length as the input 
%   coords. Values which lie outside the b-spline coefficient plot are 
%   returned as NaNs.
%
% Note that this function interpolates all coordinates that can be
% interpolated within the b-spline array (i.e. no ROI or mask is used).

    % Initialize vec_interp - points that can't be interpolated will stay
    % NaNs
    vec_interp = nan(size(coords,1),1);

    % Biquintic Kernel Matrix
    QK = [1/120  13/60  11/20 13/60  1/120 0;
         -1/24   -5/12    0    5/12  1/24  0;
          1/12    1/6   -1/2   1/6   1/12  0;
         -1/12    1/6     0   -1/6   1/12  0;
          1/24   -1/6    1/4  -1/6   1/24  0;
         -1/120   1/24  -1/12  1/12 -1/24 1/120];  
     
    % Cycle over coordinates
    for i = 0:size(coords,1)-1
        x_tilda = coords(i+1,1);
        y_tilda = coords(i+1,2);

        y_tilda_floor = floor(y_tilda);
        x_tilda_floor = floor(x_tilda);

        % Make sure top, left, bottom, and right are within the b-spline 
        % coefficient array. top, left, bottom and right are the bounding 
        % box of the b-spline coefficients used for interpolation of this
        % point;
        top = y_tilda_floor-offset_y+border_bcoef-2;
        left = x_tilda_floor-offset_x+border_bcoef-2;
        bottom = y_tilda_floor-offset_y+border_bcoef+3;
        right = x_tilda_floor-offset_x+border_bcoef+3;
        if (top >= 0 && ...
            left >= 0 && ...
            bottom < size(plot_bcoef,1) && ...
            right < size(plot_bcoef,2))
            % Set coords
            y_tilda_delta = y_tilda-y_tilda_floor;
            x_tilda_delta = x_tilda-x_tilda_floor;

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

            % Get interpolated value
            vec_interp(i+1) = y_vec*QK*plot_bcoef(top+1:bottom+1,left+1:right+1)*QK'*x_vec'; 
        end
    end
end
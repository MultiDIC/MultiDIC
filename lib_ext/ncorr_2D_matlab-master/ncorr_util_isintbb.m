function outstate = ncorr_util_isintbb(x,low,high,name)
% Tests if x is an integer between and including bounds low and high. 
%
% Inputs -----------------------------------------------------------------%
%   x - double; can be NaN, inf, a special character, complex, or a regular
%   real number.
%   low - integer; lower bound
%   high - integer; upper bound
%   name - string; used to display errordlg if number is not between bounds
%
% Outputs ----------------------------------------------------------------%
%   outstate - integer; returns either out.cancelled, out.failed, or out.success. 
%
% Note that x needs to be the output from the str2double function.

    % Initialize output
    outstate = out.failed;
    
    if (isfinite(x) && isreal(x) && mod(x,1) == 0 && x >= low && x <= high)        
        outstate = out.success;
    else
        h_error = errordlg([name ' must be an integer greater than or equal to ' num2str(low) ' and less than or equal to ' num2str(high) '.'],'Error','modal');
        uiwait(h_error);
    end
end


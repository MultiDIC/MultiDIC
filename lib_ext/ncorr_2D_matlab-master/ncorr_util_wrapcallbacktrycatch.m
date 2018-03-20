function handle_wrapcallbacktrycatch = ncorr_util_wrapcallbacktrycatch(handle_callback,handle_figure)
% This function is a wrapper that wraps the callback in a try-catch
% statement. If an error is thrown after the figure handle is closed, then
% disregard the error; if not, then rethrow it. The reason for this function
% is that closing a figure always interrupts a function, which can cause it
% to throw an error. Make sure to only throw the error if the figure which
% the function is attached to is still open.
%
% Inputs -----------------------------------------------------------------%
%   handle_callback - function handle;
%   handle_figure - figure handle;
%
% Outputs ----------------------------------------------------------------%
%   handle_wrapcallbacktrycatch - function handle;
    
    handle_wrapcallbacktrycatch = @wrapcallbacktrycatch;

    function wrapcallbacktrycatch(varargin)
        try
            handle_callback(varargin{:});
        catch err
            if (ishandle(handle_figure))
                rethrow(err);
            end
        end
    end
end

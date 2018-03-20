classdef out
% This is basically acting as an enumeration and should be called like: 
% (outstate == out.success) to see if output of a function was either 
% successful, failed, or was cancelled.

    properties(Constant,Access = private)
        s = 1;  % success
        f = 0;  % failed
        c = -1; % cancelled
    end

    methods(Static)
        function i = success()
            i = out.s;
        end       
        
        function i = failed()
            i = out.f;
        end
        
        function i = cancelled()
            i = out.c;
        end       
    end
end

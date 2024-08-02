classdef gtime_t
    % class to define the time
    
    properties
        time = 0;
        sec = 0.0;
    end
    
    methods
        function obj = gtime_t(time, sec)
            if nargin > 0
                obj.time = time;
            end
            if nargin > 1
                obj.sec = sec;
            end
        end
    end
end


classdef Timer < handle
    properties
        time_elapsed = 0;
    end
    methods 
        function reset(self)
            self.time_elapsed = 0;
        end

        function start(self)
            tic;
        end

        function pause(self)
            t = toc;
            self.time_elapsed = self.time_elapsed + t;
        end

        function time_elapsed = getTime(self)
            time_elapsed  = self.time_elapsed;
        end
    end
end
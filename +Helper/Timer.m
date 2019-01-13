classdef Timer < handle
    properties
        time_elapsed = 0;
    end
    methods 
        function reset(self)
            self.time_elapsed = 0;
        end

        function startTimer(self)
            tic;
        end

        function endTimer(self)
            t = toc;
            self.time_elapsed = self.time_elapsed + t;
        end

        function time_elapsed = getTime(self)
            time_elapsed  = self.time_elapsed;
        end
    end
end
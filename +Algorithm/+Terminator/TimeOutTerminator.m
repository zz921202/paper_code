classdef MaxIterTerminator < Algorithm.Terminator.TerminatorInterface
    properties
        MAXTIME = 1200;
    end

    methods
        function flag = terminate(self)
            cur_time = self.fo_algorithm.time_elapsed;
            flag = false;
            if cur_time >= self.MAXTIME
                str = sprintf('TIMEOUT, time elapsed:%s', num2str(cur_time));
                disp(str);
                flag = true;
            end
        end
    end
end
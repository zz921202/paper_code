classdef EstGapTerminator < Algorithm.Terminator.TerminatorInterface
    properties
        GAP = 1e-1;
    end

    methods
        function flag = terminate(self)
            cur_gap = self.fo_algorithm.cur_est_gap;
            flag = false;
            if cur_gap <= self.GAP
                str = sprintf('Using  Estimated Gap Criterion for Termination, est_gap: %s', num2str(cur_gap));
                disp(str);
                flag = true;
            end
        end
    end
end
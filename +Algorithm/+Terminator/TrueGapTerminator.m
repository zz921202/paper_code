classdef TrueGapTerminator < Algorithm.Terminator.TerminatorInterface
    properties
        GAP = 1e-2;
    end

    methods
        function flag = terminate(self)
            cur_gap = self.fo_algorithm.cur_true_gap;
            flag = false;
            if cur_gap <= self.GAP
                str = sprintf('Ref Gap Reached, true_gap: %s', num2str(cur_gap));
                disp(str);
                flag = true;
            end
        end
    end
end
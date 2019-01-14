classdef MaxIterTerminator < Algorithm.Terminator.TerminatorInterface
    properties
        MAXITERATION = 500;
    end

    methods
        function flag = terminate(self)
            cur_iter_count = self.num_iters
            flag = false;
            if cur_iter_count >= self.MAXITERATION
                str = sprintf('MAXITER reached, iter:%s', num2str(cur_iter_count));
                disp(str);
                flag = true;
            end
        end
    end
end
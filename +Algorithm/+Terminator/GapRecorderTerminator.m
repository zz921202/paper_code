classdef GapRecorderTerminator < Algorithm.Terminator.TerminatorInterface
    properties
        gap_ratios = [1e-1, 1e-2, 5e-3, 1e-3];
        next_gap_idx = 1;
        time_elapsed = [];
        iters = [];
        target_val;
        MAXITER = 5000; %TODO
    end

    methods

        function self = GapRecorderTerminator( objval)
            self.target_val = (1 + self.gap_ratios) * objval;
        end


        function flag = terminate(self)
            cur_obj = self.fo_algorithm.cur_obj;
            flag = false;
            cur_iter = self.fo_algorithm.num_iters;
            if cur_iter > self.MAXITER
                flag = true 
                return
            end

            if self.next_gap_idx > length(self.gap_ratios)
                flag = true;
                return
            end
            cur_requirement = self.target_val(self.next_gap_idx);
            
            if cur_obj <= cur_requirement
                
                cur_time = self.fo_algorithm.time_elapsed;
                
                fprintf('Reaching gap %s in %d iters and %s sec \n', num2str(self.gap_ratios(self.next_gap_idx)), cur_iter, cur_time);
                self.time_elapsed = [self.time_elapsed, cur_time];
                self.iters = [self.iters, cur_iter];

                if self.next_gap_idx == length(self.gap_ratios)
                    flag = true;
                end
                self.next_gap_idx = self.next_gap_idx + 1;
                
            end

        end

        function reset(self)
            self.next_gap_idx = 1;
            self.time_elapsed = [];
            self.iters = [];
        end

        function str = getResult(self)
            target_len = length(self.gap_ratios);
            real_len = length(self.time_elapsed);
            iter_str = ''; time_str = '';
            for idx = 1: target_len
                if idx <= real_len
                    iter_str = [iter_str, num2str(self.iters(idx)) 'iter '];
                    time_str = [time_str, num2str(self.time_elapsed(idx)) 'sec '];
                else
                    iter_str = [iter_str, 'NA '];
                    time_str = [time_str, 'NA '];
                end

            end
            str = [iter_str, ',', time_str];
        end
    end
end
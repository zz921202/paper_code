classdef SingleAlgorithm 
    properties
        alg_handle ;
        raw_writer;
        best_writer;
    end
    methods
        function self =SingleAlgorithm(alg_handle, raw_writer, best_writer)
            self.alg_handle = alg_handle;
            self.raw_writer = raw_writer;
            self.best_writer = best_writer;
        end

        function runAndRecord(self, single_problem)
            ref_problem  = single_problem.generateData();
            terminator = single_problem.getTerminator();
            data_str = single_problem.getProblemStr();
            alg = self.alg_handle(ref_problem, terminator);
            best_val = Inf;
            best_gap = Inf;
            best_iter = Inf;
            ind = 0;
            best_ind = 1;
            while alg.nextGridParam()
                ind  = ind + 1;
                param_str = alg.showGridParam(ind);

                [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run();
                res_str = sprintf('%d, %d, %s, %s, %s, %s sec, obj: %s', ind, num_iters, alg.getName(), num2str(true_gap), num2str(est_gap), num2str(time_elapsed), num2str(cur_obj_val) );
                self.raw_writer.write(sprintf( '%s,%s, %s', data_str, param_str, res_str));

                if cur_obj_val < best_val && num_iters == best_iter
                    best_val = cur_obj_val;
                    
                    best_ind = ind;
                end
                if num_iters < best_iter -1
                    best_iter = num_iters;
                    best_ind = ind;
                end
            end

            terminator = single_problem.getTerminator();
            terminator.all_terminators{1}.MAXITERATION = 100;
            alg = self.alg_handle(ref_problem, terminator);
            alg.setGridParam(best_ind);
            param_str = alg.showGridParam(best_ind);
            [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run();
            res_str = sprintf('%d, %d, %s, %s, %s, %s sec, obj: %s', best_ind, num_iters, alg.getName(), num2str(true_gap), num2str(est_gap), num2str(time_elapsed), num2str(cur_obj_val) );
            self.best_writer.write( sprintf('%s,%s, %s', data_str, param_str, res_str));
            self.raw_writer.write( sprintf('%s,%s, %s', data_str, param_str, res_str))
            fprintf('%s, gap: %s best_val: %s', alg.getName(), num2str(true_gap), num2str(cur_obj_val));
        end


    end
end
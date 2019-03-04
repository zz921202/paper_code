classdef SingleTuningAlgorithm 
    properties
        alg_handle ;
        raw_writer;
        res_writer;
        tunning_iters = 200;

    end
    methods
        function self =SingleTuningAlgorithm(alg_handle, res_writer, tuning_iters)
            self.alg_handle = alg_handle;            
            self.res_writer = res_writer;
            self.tunning_iters = tuning_iters;
        end

        function runAndRecord(self, single_problem)
            [eu_problem, en_problem]  = single_problem.generateData();
            
            data_str = single_problem.getProblemStr();
            
            rec_terminator = Algorithm.Terminator.GapRecorderTerminator(eu_problem.opt_val);

            eu_alg = self.alg_handle(eu_problem, rec_terminator);
            eu_alg.setTuningIter(self.tunning_iters);
            eu_alg.startTuning();
            [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = eu_alg.run();
            res_str1 = rec_terminator.getResult();
            res_str2 = sprintf('%d, %d, %s, %s, %s, %s sec, obj: %s', eu_alg.best_ind, num_iters, eu_alg.getName(), num2str(true_gap), num2str(est_gap), num2str(time_elapsed), num2str(cur_obj_val) );
            self.res_writer.write( sprintf('%s,%s, %s, %s', data_str, 'Euclidean', res_str1, res_str2));
            
            fprintf('%s, gap: %s best_val: %s, Euclidean', eu_alg.getName(), num2str(true_gap), num2str(cur_obj_val));
            
            
            rec_terminator = Algorithm.Terminator.GapRecorderTerminator(en_problem.opt_val);
            en_alg = self.alg_handle(en_problem, rec_terminator);
            if en_alg.solveEntropy()
                en_alg.setTuningIter(self.tunning_iters);
                en_alg.startTuning();
                [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = en_alg.run();
                res_str1 = rec_terminator.getResult();
                res_str2 = sprintf('%d, %d, %s, %s, %s, %s sec, obj: %s', en_alg.best_ind, num_iters, en_alg.getName(), num2str(true_gap), num2str(est_gap), num2str(time_elapsed), num2str(cur_obj_val) );
                self.res_writer.write( sprintf('%s,%s, %s, %s', data_str, 'Entropy', res_str1, res_str2));

                fprintf('%s, gap: %s best_val: %s, Entropy', en_alg.getName(), num2str(true_gap), num2str(cur_obj_val));
            end
        end


    end
end
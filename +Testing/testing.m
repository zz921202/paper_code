ms = {10, 20,  100};
ks = ms;
rng_seeds = {100};
max_iters = {10, 50, 100};
algs = {@Algorithm.Subgradient, @Algorithm.FixedNestrov, @Algorithm.DynamicNestrov, @Algorithm.PDDAlgorithm, ...
        @Algorithm.DSLAlgorithm, @Algorithm.ABPAlgorithm};

algs_names = {'subGrad', 'fixNes', 'dynNes', 'pdd', 'dsl', 'abp'};



for m_ind = 1 : 1%length(ms)
    for k_ind = 1: 1!%length(ks)
        for seed_ind = 1: length(rng_seeds)

            cur_seed = rng_seeds{seed_ind}; cur_m = ms{m_ind}; cur_k = ks{k_ind};
            
            mat_gen = DataGenerator.SimpleCompleteRecourseRandomData(cur_m,cur_m*2);
            ref_problem_eu = Problem.LinearRatioUncertainty(mat_gen,  'Euclidean');
            ref_problem_en = Problem.LinearRatioUncertainty(mat_gen,  'Entropy');
            rng(cur_seed);
            ref_problem_en.generateData(cur_k);
            rng(cur_seed);
            ref_problem_eu.generateData(cur_k);


            best_f = fopen('best_data.txt', 'a');
            raw_f = fopen('raw_data.txt', 'a');
            
            for alg_ind = 2%1:length(algs)
                for iter_ind = 1:length(max_iters)
                    cur_max_iter = max_iters{iter_ind};
                    cur_alg = algs{alg_ind};
                    cur_alg_name = algs_names{alg_ind};
                    fprintf('runing %s !!! m is %d , k is %d\n', cur_alg_name, cur_m, cur_k);
                    data_description = sprintf('%d, %d, %d, %s, %d', cur_m, cur_k, cur_seed, cur_alg_name, cur_max_iter);
                    max_iter_terminator = Algorithm.Terminator.MaxIterTerminator();
                    est_gap_teminator = Algorithm.Terminator.EstGapTerminator();

                    max_iter_terminator.MAXITERATION = cur_max_iter;
                    terminator = Algorithm.Terminator.CompositeTerminator({max_iter_terminator, est_gap_teminator});

                    alg = cur_alg(ref_problem_eu, terminator);
                    best_val = Inf;
                    best_gap = Inf;
                    best_ind = 2;
                    ind = 0
                    while alg.nextGridParam()
                        ind  = ind + 1;
                        param_str = alg.showGridParam(ind);

                        [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run();
                        res_str = sprintf('%d, %d, %s, %s, %s, %s', ind, num_iters, num2str(true_gap), num2str(est_gap), num2str(time_elapsed), num2str(cur_obj_val) );
                        fprintf(raw_f, '%s,%s, %s, Euclidean\n', data_description, param_str, res_str);
                        if cur_obj_val < best_val
                            best_val = cur_obj_val;
                            best_gap = true_gap;
                            best_ind = ind;
                        end
                    end

                    
                    alg.setGridParam(best_ind);
                    param_str = alg.showGridParam(best_ind);
                    [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run();

                    res_str = sprintf('%d, %d, %s, %s, %s, %s', best_ind, num_iters, num2str(true_gap), num2str(est_gap), num2str(time_elapsed), num2str(cur_obj_val) );
                    fprintf(best_f, '%s,%s, %s, Euclidean\n', data_description, param_str, res_str);

                    %%

                    % alg = cur_alg(ref_problem_en, terminator);
                    % best_val = Inf;
                    % best_gap = Inf;
                    % best_ind = 2;
                    % ind = 0
                    % while alg.nextGridParam()
                    %     ind  = ind + 1;
                    %     param_str = alg.showGridParam(ind);

                    %     [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run();
                    %     res_str = sprintf('%d, %d, %s, %s, %s, %s', ind, num_iters, num2str(true_gap), num2str(est_gap), num2str(time_elapsed), num2str(cur_obj_val) );
                    %     fprintf(raw_f, '%s,%s, %s, Entropy\n', data_description, param_str, res_str);
                    %     if cur_obj_val < best_val
                    %         best_val = cur_obj_val;
                    %         best_gap = true_gap;
                    %         best_ind = ind;
                    %     end
                    % end

                    
                    % alg.setGridParam(best_ind);
                    % param_str = alg.showGridParam(best_ind);
                    % [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run();

                    % res_str = sprintf('%d, %d, %s, %s, %s, %s', ind, num_iters, num2str(true_gap), num2str(est_gap), num2str(time_elapsed), num2str(cur_obj_val) );
                    % fprintf(best_f, '%s,%s, %s, Entropy\n', data_description, param_str, res_str);
                end
            end
            fclose(best_f);
            fclose(raw_f);
        end
    end
end









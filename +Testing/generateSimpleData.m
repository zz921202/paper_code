function problems =generateSimpleData()
    rng_seeds = 100;
    max_iters = [20];
    num_rows = [5];
    num_scenarios = [20];
    breg_dists = {'Euclidean', 'BoxEntropy'};
    problems = [];
    %%%%%%%%%%%%%%%
    for m_ind = 1 : length(num_rows)
        for k_ind = 1: length(num_scenarios)
            for seed_ind = 1: length(rng_seeds)
                for iter_ind = 1: length(max_iters)
                    for breg_ind = 1: length(breg_dists)
                        cur_problem = Testing.SingleRandomProblem();
                        cur_problem.rng_seed = rng_seeds(seed_ind); 
                        cur_problem.num_rows = num_rows(m_ind); 
                        cur_problem.num_scenarios = num_scenarios(k_ind);
                        cur_problem.breg_dist = breg_dists{breg_ind}; 
                        cur_problem.max_iter = max_iters(iter_ind);
                        problems = [problems, cur_problem];
                    end
                end
            end
        end
    end


end 
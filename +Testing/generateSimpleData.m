function problems =generateSimpleData()
    rng_seeds = 100:100:300;
    max_iters = [200];
    num_rows = [20];
    num_scenarios = [20000];
    breg_dists = {'Euclidean'};
    % beta_alpha_pairs = {[0, 10000]};% [0 10], [0 100], [0 1000]};
    beta_alpha_pairs = {[0 20000]}
    problems = [];
    %%%%%%%%%%%%%%%
    for prob_ind = 1: length(beta_alpha_pairs)
        for m_ind = 1 : length(num_rows)
            for k_ind = 1: length(num_scenarios)
                for seed_ind = 1: length(rng_seeds)
                    for iter_ind = 1: length(max_iters)
                        for breg_ind = 1: length(breg_dists)
                            cur_al_beta = beta_alpha_pairs{prob_ind};
  
                            cur_problem = Testing.SingleRandomProblem();
                            cur_problem.alpha = cur_al_beta(1);
                            cur_problem.beta = cur_al_beta(2);
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


end 
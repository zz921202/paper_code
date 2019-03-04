function problems =generateSimpleData()
    rng_seeds = 100;
    % max_tune_iters = [200];
    num_rows = [2];   
    num_scenarios = [2];
    am_handles = { @Problem.TransportUncertainty, @Problem.X2Uncertainty, @Problem.RatioAmbiguitySet};
    beta_alpha_pairs = {[1 1 0], [1/2 2 0.2], [1/10 10 5], [1/40 40 20], [0 400 100]};
    beta_alpha_pairs = {[0.9 1.1 0.1]};
    problems = [];
    %%%%%%%%%%%%%%%
    for prob_ind = 1: length(beta_alpha_pairs)
        for m_ind = 1 : length(num_rows)
            for k_ind = 1: length(num_scenarios)
                for seed_ind = 1: length(rng_seeds)
                    % for iter_ind = 1: length(max_iters)
                        for am_ind = 1: length(am_handles)

                            cur_num_rows = num_rows(m_ind); 
                            mat_gen = DataGenerator.ToyRecourseRandomDataGenerator(cur_num_rows, cur_num_rows*2);
                            cur_al_beta = beta_alpha_pairs{prob_ind};
                            problem_handle = @Problem.ToyLinearProblem;
                            am_handle = am_handles{am_ind};
                            cur_problem = Testing.SingleRandomProblem(problem_handle, mat_gen, am_handle);
                            cur_problem.alpha = cur_al_beta(1);
                            cur_problem.beta = cur_al_beta(2);
                            cur_problem.radius = cur_al_beta(3);
                            cur_problem.rng_seed = rng_seeds(seed_ind); 
                            
                            cur_problem.num_scenarios = num_scenarios(k_ind);
                            % cur_problem.breg_dist = breg_dists{breg_ind}; 
                            % cur_problem.max_iter = max_iters(iter_ind);
                            problems = [problems, cur_problem];
                        end
                    % end
                end
            end
        end
    end


end 
% this file tests all algorithms by iteration number 



k = 5;
rng_seed = 100;
m = 20;
n = 40;
dist = 'Euclidean';
alpha = 1;
beta =1;
radius = 5;
% ref_problem = Testing.getRefProblem(n,m, k, alpha, beta, rng_seed, dist)
% ref_problem = Testing.getRefTransportProblem(n, m, k, rng_seed, dist, radius);
% ref_problem = Testing.getRefX2(n, m, k, rng_seed, dist, radius);

ref_problem = Testing.getRefSSN(n,m, k, alpha, beta, rng_seed, dist)

terminator = Algorithm.Terminator.MaxIterTerminator();


dsl_type_alg = {@Algorithm.ExperimentDSLAlgorithm, @Algorithm.ABPAlgorithm};

terminator.MAXITERATION = 200
true_gap_history = [];
names = {};
for idx = 1: length(dsl_type_alg)
    cur_alg_handle = dsl_type_alg{idx};
    cur_alg = cur_alg_handle(ref_problem, terminator);
    [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = cur_alg.run();
    cur_gap_history = cur_alg.obj_history - ref_problem.opt_val;
    true_gap_history = [true_gap_history, cur_gap_history];
    names = [names, cur_alg.getName()];
end

fo_type_alg = { @Algorithm.PDDAlgorithm, @Algorithm.PDHGAlgorithm, @Algorithm.DynamicNestrov}

for idx = 1: length(fo_type_alg)
    terminator.MAXITERATION = 50;
    cur_alg_handle = fo_type_alg{idx};
    alg = cur_alg_handle(ref_problem, terminator);
    best_val = Inf;
    best_gap = Inf;
    best_ind = 0;
    ind = 0;
    % alg.setGridParam(ind+1)
    while alg.nextGridParam()
        ind  = ind + 1;
        disp(alg.showGridParam(ind));
        % alg.setGridParam(ind);
        [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run();
        fprintf('%s th parameter choice with obj_val%s and gap %s\n', num2str(ind), num2str(cur_obj_val), num2str(true_gap));
        if cur_obj_val < best_val
            best_val = cur_obj_val;
            best_gap = true_gap;
            best_ind = ind;
        end
    end
    terminator.MAXITERATION = 200;
    alg.showGridParam(best_ind);
    alg.setGridParam(best_ind);

    [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run()
    cur_gap_history = alg.obj_history - ref_problem.opt_val;
    true_gap_history = [true_gap_history, cur_gap_history];
    names = [names, alg.getName()];
end

figure 
semilogy(true_gap_history, 'LineWidth',8);
legend(names{:});
title(sprintf('m = %d, n = %d, k =%d, radius %s, objval=%s', m,n,k,num2str(radius), num2str(ref_problem.opt_val)));
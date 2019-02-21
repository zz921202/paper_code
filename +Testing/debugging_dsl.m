k = 10;
rng_seed = 100;
m = 20;
n = 40;
dist = 'BoxEntropy';
alpha = 0;
beta =k;
ref_problem = Testing.getRefProblem(n,m, k, alpha, beta, rng_seed, dist)

max_iter_terminator = Algorithm.Terminator.MaxIterTerminator();
est_gap_teminator = Algorithm.Terminator.EstGapTerminator();

max_iter_terminator.MAXITERATION = 200;
terminator = Algorithm.Terminator.CompositeTerminator({max_iter_terminator, est_gap_teminator});


max_iter_terminator1 = Algorithm.Terminator.MaxIterTerminator();
est_gap_teminator1 = Algorithm.Terminator.EstGapTerminator();

max_iter_terminator1.MAXITERATION = 200;
terminator1 = Algorithm.Terminator.CompositeTerminator({max_iter_terminator1, est_gap_teminator1});


alg = Algorithm.ExperimentDSLAlgorithm(ref_problem, terminator);
% ref_alg = Algorithm.
ref_alg = Algorithm.ABPAlgorithm(ref_problem, terminator1);

% alg.setGridParam(1);
best_val = Inf;
best_gap = Inf;
best_ind = 1;
ind = 0
% while alg.nextGridParam()
%     ind  = ind + 1;
%     [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run()
%     if cur_obj_val < best_val
%         best_val = cur_obj_val;
%         best_gap = true_gap;
%         best_ind = ind
%     end
% end

alg.showGridParam(best_ind);
alg.setGridParam(best_ind);
[cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run();
[cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = ref_alg.run();

figure

semilogy(alg.obj_history - ref_problem.opt_val);
hold on
semilogy(ref_alg.obj_history - ref_problem.opt_val);
legend('exp', 'ref')
% plot([0, length(alg.obj_history)], [ref_problem.opt_val, ref_problem.opt_val], 'r')
title('obj history: DSL')
figure
semilogy(alg.gap_history)
hold on
semilogy(ref_alg.gap_history)
legend('exp', 'ref')
title('gap history: DSL history')

% alg.bundle.showConstraints();
% figure
% plot(alg.xt_history');
% title('xt history');
% xt_history = alg.xt_history
% [A, b] = ref_problem.getXConstraint();
% A = full(A)
% b = b

% % figure;
% % hold on
% % plot(alg.fake_cost_history);
% % plot(alg.temp_cost_history, 'x');
% % legend('fake cost', 'true cost');
% figure

% plot(alg.dist_to_x_history);
% title('dist to true x');
% [gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = alg.getProxParams();
% fake_cost_for_opt_x = ref_problem.evalFakeCost(ref_problem.ref_x, mu_p_t, alg.smooth_p, mu_pi_t, alg.smooth_pis)
% fake_cost_for_sol =  ref_problem.evalFakeCost(cur_x, mu_p_t, alg.smooth_p, mu_pi_t, alg.smooth_pis)
% figure
% plot(alg.x_history(:, 100:200)', '-x');
% legend('1', '2', '3', '4', '5', '6');
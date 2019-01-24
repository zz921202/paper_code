k = 10;
n = 20;
m = n/2;

% pis = {};
% for i = 1: k
%     pis = [pis, ones(m,1)];
% end

indi_costs = zeros(k, 1);
indi_costs(1) = 1;
mat_gen = DataGenerator.SimpleCompleteRecourseRandomData(m,n);

ref_problem = Problem.LinearRatioUncertainty(mat_gen,  'Euclidean');
ref_problem.alpha = 0;%x
ref_problem.beta = k;
rng(100)
ref_problem.generateData(k);


max_iter_terminator = Algorithm.Terminator.MaxIterTerminator();
est_gap_teminator = Algorithm.Terminator.EstGapTerminator();

max_iter_terminator.MAXITERATION = 100;
terminator = Algorithm.Terminator.CompositeTerminator({max_iter_terminator, est_gap_teminator});



alg = Algorithm.DSLAlgorithm(ref_problem, terminator);

% alg = Algorithm.ABPAlgorithm(ref_problem, terminator);

% alg.setGridParam(1);
best_val = Inf;
best_gap = Inf;
best_ind = 2;
ind = 0
while alg.nextGridParam()
    ind  = ind + 1;
    [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run()
    if cur_obj_val < best_val
        best_val = cur_obj_val;
        best_gap = true_gap;
        best_ind = ind
    end
end

alg.showGridParam(best_ind);
alg.setGridParam(best_ind);
[cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run()

figure

plot(alg.obj_history);
hold on

plot([0, length(alg.obj_history)], [ref_problem.opt_val, ref_problem.opt_val], 'r')
title('obj history: DSL')
figure
plot(alg.gap_history)
title('gap history: DSL obj history')

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
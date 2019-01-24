k = 10;
n = 20;
m = 10;

pis = {};
for i = 1: k
    pis = [pis, ones(m,1)];
end

indi_costs = zeros(k, 1);
indi_costs(1) = 1;
mat_gen = DataGenerator.SimpleCompleteRecourseRandomData(m,n);

ref_problem = Problem.LinearRatioUncertainty(mat_gen,  'Euclidean');
ref_problem.alpha = 0.5;%x
ref_problem.beta = k/2;
rng(100)
ref_problem.generateData(k);
% tic
% [ref_p, ref_grad, ref_cost] = ref_problem.projectP(1, ones(k,1)/ k, indi_costs, pis)
% toc


% cus_problem = Problem.LinearRatioUncertainty(mat_gen, 'BoxEntropy');
% % cus_problem.alpha = 0;
% % cus_problem.beta = k;
% rng(0)
% cus_problem.generateData(k);
% tic
% [cus_p, cus_grad, cus_cost] = cus_problem.projectP(1, ones(k,1)/ k, indi_costs, pis)
% toc
% y =[   
%        0.5737
%          0];

terminator = Algorithm.Terminator.MaxIterTerminator();
terminator.MAXITERATION = 100;
alg = Algorithm.DynamicNestrov(ref_problem, terminator);
% alg.setGridParam(1);
best_val = Inf;
best_gap = Inf;
best_ind = 0;
ind = 0
while alg.nextGridParam()
    ind  = ind + 1;
    [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run();
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
title('obj history: fixed Nestrov')
% title('obj history')
% % figure;
% % hold on
% % plot(alg.fake_cost_history);
% % plot(alg.temp_cost_history, 'x');
% % legend('fake cost', 'true cost');
% figure

% plot(alg.dist_to_x_history);
% title('dist to true x');
[gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = alg.getProxParams();
fake_cost_for_opt_x = ref_problem.evalFakeCost(ref_problem.ref_x, mu_p_t, alg.smooth_p, mu_pi_t, alg.smooth_pis)
fake_cost_for_sol =  ref_problem.evalFakeCost(cur_x, mu_p_t, alg.smooth_p, mu_pi_t, alg.smooth_pis)
% figure
% plot(alg.x_history(:, 100:200)', '-x');
% legend('1', '2', '3', '4', '5', '6');
% testing euclidean projector 


% c = [1;1];
% A = [2, 1];
% b = [3];
% x_pre = [5; 5];

% eu_projector = Helper.EuclideanProjector();
% eu_projector.setConstraint(A, b);
% [val, soln]= eu_projector.solve(c);
% [soln, obj] = eu_projector.project(1, x_pre, c);

% A = [1];
% b = -3;
% x_pre = 1;
% c = 1;
% eu_projector = Helper.EntropyProjector();
% eu_projector.setConstraint(A, b);
% [soln, obj] = eu_projector.project(1, x_pre, c)
% model = eu_projector.model

k = 6;
n = 10;
m = 5;

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
y =[   
       0.5737
         0];

terminator = Algorithm.Terminator.MaxIterTerminator();
terminator.MAXITERATION = 100;
alg_sub = Algorithm.Subgradient(ref_problem, terminator);
alg_sub.setGridParam(1);
[x, obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg_sub.run()

figure

plot(alg_sub.obj_history);
hold on
title('obj history')

plot([0, length(alg_sub.obj_history)], [ref_problem.opt_val, ref_problem.opt_val], 'r')
% title('obj history: fixed Nestrov')
% figure;
% hold on
% plot(alg_sub.fake_cost_history);
% plot(alg_sub.temp_cost_history, 'x');
% legend('fake cost', 'true cost');
% figure

% plot(alg_sub.dist_to_x_history);
% title('dist to true x');
[gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = alg_sub.getProxParams();
fake_cost_for_opt_x = ref_problem.evalFakeCost(ref_problem.ref_x, mu_p_t, alg_sub.smooth_p, mu_pi_t, alg_sub.smooth_pis)
fake_cost_for_sol =  ref_problem.evalFakeCost(x, mu_p_t, alg_sub.smooth_p, mu_pi_t, alg_sub.smooth_pis)
% figure
% plot(alg_sub.x_history(:, 100:200)', '-x');
% legend('1', '2', '3', '4', '5', '6');
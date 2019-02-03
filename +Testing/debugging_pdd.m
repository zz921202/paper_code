k = 20;
n = 20;
m = n/2;

pis = {};
for i = 1: k
    pis = [pis, ones(m,1)];
end

indi_costs = zeros(k, 1);
indi_costs(1) = 1;
mat_gen = DataGenerator.SimpleCompleteRecourseRandomData(m,n);

ref_problem = Problem.LinearRatioUncertainty(mat_gen,  'BoxEntropy');
% ref_problem.alpha = 1;%x
% ref_problem.beta = 1;
rng(100)
ref_problem.generateData(k);


terminator = Algorithm.Terminator.MaxIterTerminator();
terminator.MAXITERATION = 400;
alg = Algorithm.PDDAlgorithm(ref_problem, terminator);
ref_problem.alpha = 0.99;
ref_problem.beta = 1.02;


best_val = Inf;
best_gap = Inf;
best_ind = 15;
ind = 0
% while alg.nextGridParam()
%     ind  = ind + 1;
    
%     [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run();
%     fprintf('%s th parameter choice with obj_val%s and gap %s\n', num2str(ind), num2str(cur_obj_val), num2str(true_gap));
%     if cur_obj_val < best_val
%         best_val = cur_obj_val;
%         best_gap = true_gap;
%         best_ind = ind;
%     end
% end

alg.showGridParam(best_ind);
alg.setGridParam(best_ind);
[cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run()
figure


semilogy(alg.obj_history - ref_problem.opt_val);
% hold on

% plot([0, length(alg.obj_history)], [ref_problem.opt_val, ref_problem.opt_val], 'r')
title('obj history: PDD')


figure
semilogy(alg.gap_history)
title('gap history: PDD')

% figure
% plot(alg.x_history')
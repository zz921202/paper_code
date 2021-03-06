k = 4;
rng_seed = 100;
m = 2;
n = 4;
dist = 'Euclidean';
alpha = 0;
beta =k;
ref_problem = Testing.getRefProblem(n,m, k, alpha, beta, rng_seed, dist)

terminator = Algorithm.Terminator.MaxIterTerminator();
terminator.MAXITERATION = 200;
alg = Algorithm.PDHGAlgorithm(ref_problem, terminator);



best_val = Inf;
best_gap = Inf;
% best_ind = 52;
ind = 0;
best_ind = 0;

while alg.nextGridParam()
    ind  = ind + 1;
    disp(alg.showGridParam(ind));
    % alg.setGridParam(ind);
    try
        [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run();

        fprintf('%s th parameter choice with obj_val%s and gap %s\n', num2str(ind), num2str(cur_obj_val), num2str(true_gap));
        if cur_obj_val < best_val
            best_val = cur_obj_val;
            best_gap = true_gap;
            best_ind = ind;
        end
    catch
        fprintf('%s, numerically unstable', alg.showGridParam(ind));
    end
end
terminator.MAXITERATION = 1000;
alg = Algorithm.PDHGAlgorithm(ref_problem, terminator);
alg.showGridParam(best_ind);
alg.setGridParam(best_ind);

[cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run()
figure


semilogy(alg.obj_history - ref_problem.opt_val);
% hold on

% plot([0, length(alg.obj_history)], [ref_problem.opt_val, ref_problem.opt_val], 'r')
title('tunting at 50 obj history: PDHG')

figure 
semilogy(alg.x_diff_history);
title('x diff')

figure
semilogy(alg.ws_history');
title('ws history')

% figure
% semilogy(alg.gap_history)
% title('gap history: PDD')

% figure
% plot(alg.x_history')
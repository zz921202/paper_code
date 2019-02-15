k = 20;
n = 40;
m = n/2;


mat_gen = DataGenerator.SimpleCompleteRecourseRandomData(m,n);

ref_problem = Problem.LinearRatioUncertainty(mat_gen,  'Euclidean');
% ref_problem.alpha = 1;%x
% ref_problem.beta = 1;
% ref_problem.alpha = 1;
% ref_problem.beta = 1;
rng(10)
ref_problem.generateData(k);


terminator = Algorithm.Terminator.MaxIterTerminator();
terminator.MAXITERATION = 20;
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
    % try
        [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run();

        fprintf('%s th parameter choice with obj_val%s and gap %s\n', num2str(ind), num2str(cur_obj_val), num2str(true_gap));
        if cur_obj_val < best_val
            best_val = cur_obj_val;
            best_gap = true_gap;
            best_ind = ind;
        end
    % catch
        % fprintf('%s, no good\n', alg.showGridParam(ind));
    % end
end
terminator.MAXITERATION = 200;
alg = Algorithm.PDHGAlgorithm(ref_problem, terminator);
alg.showGridParam(best_ind);
alg.setGridParam(best_ind);

[cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = alg.run()
figure


semilogy(alg.obj_history - ref_problem.opt_val);
% hold on

% plot([0, length(alg.obj_history)], [ref_problem.opt_val, ref_problem.opt_val], 'r')
title('obj history: PDHG')

figure 
plot(alg.x_diff_history);
title('x diff')

figure
plot(alg.ws_history');
title('ws history')

% figure
% semilogy(alg.gap_history)
% title('gap history: PDD')

% figure
% plot(alg.x_history')
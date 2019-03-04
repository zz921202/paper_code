alg_handles = {  @Algorithm.DynamicNestrov, @Algorithm.PDDAlgorithm, ...
        @Algorithm.ExperimentDSLAlgorithm, @Algorithm.ABPAlgorithm, @Algorithm.PDHGAlgorithm};
% alg_handles = { @Algorithm.PDDAlgorithm}% @Algorithm.PDDAlgorithm}
best_writer = Helper.Writer('../result/testing33.csv');
% raw_writer = Helper.Writer('../result/modify_pdd_best.csv');
best_writer.write('\n\n')
tuning_iters = 50;
% raw_writer.write('\n\n')
algs = [];

for alg_ind  = 1 : length(alg_handles)
    cur_alg = Testing.SingleTuningAlgorithm(alg_handles{alg_ind},  best_writer, tuning_iters);
    algs = [algs, cur_alg];
end


problems = Testing.generateSimpleData();

for ind = 1: length(problems)
    % raw_writer.write('\n')
    % best_writer.write('\n')
    cur_problem = problems(ind);
    fprintf('testing problem %s', cur_problem.getProblemStr);
    % ref_problem = cur_problem.generateData();
    for alg_ind  =1 : length(algs)
        cur_alg = algs(alg_ind);

        cur_alg.runAndRecord(cur_problem);
    end
end
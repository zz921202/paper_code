alg_handles = {@Algorithm.Subgradient, @Algorithm.FixedNestrov, @Algorithm.DynamicNestrov, @Algorithm.PDDAlgorithm, ...
        @Algorithm.DSLAlgorithm, @Algorithm.ABPAlgorithm};
alg_handles = { @Algorithm.DynamicNestrov}% @Algorithm.PDDAlgorithm}
best_writer = Helper.Writer('/Users/Zhe/Desktop/result/nes3_entropy_best.csv');
raw_writer = Helper.Writer('/Users/Zhe/Desktop/result/nes3_entropy_raw.csv');
best_writer.write('\n\n')
raw_writer.write('\n\n')
algs = [];
for alg_ind  = 1 : length(alg_handles)
    cur_alg = Testing.SingleTuningAlgorithm(alg_handles{alg_ind},  raw_writer, best_writer);
    algs = [algs, cur_alg];
end


problems = Testing.generateSimpleData();

for ind = 1: length(problems)
    raw_writer.write('\n')
    best_writer.write('\n')
    cur_problem = problems(ind);
    fprintf('testing problem %s', cur_problem.getProblemStr);
    ref_problem = cur_problem.generateData();
    for alg_ind  =1 : length(algs)
        cur_alg = algs(alg_ind);

        cur_alg.runAndRecord(cur_problem, ref_problem);
    end
end
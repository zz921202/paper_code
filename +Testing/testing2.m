am_handle = @Problem.RatioAmbiguitySet;


mat_gen = DataGenerator.RealDataGenerator('./+Data/SSN_50.mat');
cur_al_beta = [0.9 1.1 0.1]
problem_handle = @Problem.RealProblem;

cur_problem = Testing.SingleRandomProblem(problem_handle, mat_gen, am_handle);
cur_problem.alpha = cur_al_beta(1);
cur_problem.beta = cur_al_beta(2);
cur_problem.radius = cur_al_beta(3);
cur_problem.rng_seed = 100; 

cur_problem.num_scenarios = 50;


alg_handles = {  @Algorithm.DynamicNestrov, @Algorithm.PDDAlgorithm, ...
        @Algorithm.ExperimentDSLAlgorithm, @Algorithm.ABPAlgorithm};
% alg_handles = { @Algorithm.PDDAlgorithm}% @Algorithm.PDDAlgorithm}
best_writer = Helper.Writer('../result/testing_real_ratio34.csv');
% raw_writer = Helper.Writer('../result/modify_pdd_best.csv');
best_writer.write('\n\n')
tuning_iters = 50;%TODO
% raw_writer.write('\n\n')
algs = [];

for alg_ind  = 1 : length(alg_handles)
    cur_alg = Testing.SingleTuningAlgorithm(alg_handles{alg_ind},  best_writer, tuning_iters);
    algs = [algs, cur_alg];
end



fprintf('testing problem %s', cur_problem.getProblemStr());
% ref_problem = cur_problem.generateData();
for alg_ind  =1 : length(algs)
       cur_alg = algs(alg_ind);
       cur_alg.runAndRecord(cur_problem);
end

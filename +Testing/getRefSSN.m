function ref_problem = getRefTransportProblem(n,m, k, alpha, beta, rng_seed, dist)


    mat_gen = DataGenerator.RealDataGenerator('+Data/SSN_50.mat')
    ambiguity_set = Problem.RatioAmbiguitySet(dist);
    ambiguity_set.alpha = alpha;
    ambiguity_set.beta = beta;
    ref_problem = Problem.RealProblem(mat_gen,  ambiguity_set);
    rng(rng_seed);
    ref_problem.generateData(k);
end
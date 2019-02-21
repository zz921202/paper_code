function ref_problem = getRefProblem(n,m, k, alpha, beta, rng_seed, dist)


    mat_gen = DataGenerator.ToyRecourseRandomDataGenerator(m,n);
    ambiguity_set = Problem.RatioAmbiguitySet(dist);
    ambiguity_set.alpha = alpha;
    ambiguity_set.beta = beta;
    ref_problem = Problem.ToyLinearProblem(mat_gen,  ambiguity_set);
    rng(rng_seed);
    ref_problem.generateData(k);
end
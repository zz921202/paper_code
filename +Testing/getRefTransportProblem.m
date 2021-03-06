function ref_problem = getRefTransportProblem(n,m, k, rng_seed, dist, radius)


    mat_gen = DataGenerator.ToyRecourseRandomDataGenerator(m,n);
    ambiguity_set = Problem.TransportUncertainty(dist);
    ambiguity_set.radius = radius;
    
    ref_problem = Problem.ToyLinearProblem(mat_gen,  ambiguity_set);
    rng(rng_seed);
    ref_problem.generateData(k);
end
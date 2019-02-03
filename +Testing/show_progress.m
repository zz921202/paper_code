m = 10;
n = 20;
k = 10;
mat_gen = DataGenerator.SimpleCompleteRecourseRandomData(m,n);

ref_problem = Problem.LinearRatioUncertainty(mat_gen,  'Euclidean');


cur_seed = rng_seeds{seed_ind}; cur_m = ms{m_ind}; cur_k = ks{k_ind};
            
mat_gen = DataGenerator.SimpleCompleteRecourseRandomData(cur_m,cur_m*2);
ref_problem_eu = Problem.LinearRatioUncertainty(mat_gen,  'Euclidean');
ref_problem_en = Problem.LinearRatioUncertainty(mat_gen,  'Entropy');
rng(cur_seed);
ref_problem_en.generateData(cur_k);
rng(cur_seed);
ref_problem_eu.generateData(cur_k);
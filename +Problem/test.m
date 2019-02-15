% mat_gen = DataGenerator.SimpleCompleteRecourseRandomData(1, 2); 
% ref_problem = Problem.LinearRatioUncertainty(mat_gen,  'Euclidean');

% ref_problem.alpha =.99;
% ref_problem.beta = 1.02;

% ref_problem.generateData(10);
% [omega_p2, ratio, A] = ref_problem.getProbParamEstimate()

disp('hello world')
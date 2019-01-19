% testing euclidean projector 


% c = [1;1];
% A = [2, 1];
% b = [3];
% x_pre = [5; 5];

% eu_projector = Helper.EuclideanProjector();
% eu_projector.setConstraint(A, b);
% [val, soln]= eu_projector.solve(c);
% [soln, obj] = eu_projector.project(1, x_pre, c);

% A = [1];
% b = -3;
% x_pre = 1;
% c = 1;
% eu_projector = Helper.EntropyProjector();
% eu_projector.setConstraint(A, b);
% [soln, obj] = eu_projector.project(1, x_pre, c)
% model = eu_projector.model

k = 200;
n = 2;
m = 1;

pis = {};
for i = 1: k
    pis = [pis, ones(m,1)];
end

indi_costs = zeros(k, 1);
indi_costs(1) = 1;
mat_gen = DataGenerator.SimpleCompleteRecourseRandomData(m,n)

ref_problem = Problem.LinearRatioUncertainty(mat_gen,  'Entropy');
% ref_problem.alpha = 0;
% ref_problem.beta = k;
rng(0)
ref_problem.generateData(k);
tic
[ref_p, ref_grad, ref_cost] = ref_problem.projectP(1, ones(k,1)/ k, indi_costs, pis)
toc


cus_problem = Problem.LinearRatioUncertainty(mat_gen, 'BoxEntropy');
% cus_problem.alpha = 0;
% cus_problem.beta = k;
rng(0)
cus_problem.generateData(k);
tic
[cus_p, cus_grad, cus_cost] = cus_problem.projectP(1, ones(k,1)/ k, indi_costs, pis)
toc
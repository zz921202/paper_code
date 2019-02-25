% k = 000;
k = 5000
% distance_mat = [0 16.5476; 16.5476 0];
radius = 1;
% rng(200)
individual_cost = rand(k, 1) * 10;
% individual_cost = [8; 300];
ref_p = ones(k, 1) ./ k;
DELTA = 1e-16;
prox_center = rand(k, 1);
prox_center = prox_center ./ sum(prox_center);

prox_param = 10;
tic
cvx_begin quiet
    variable p(k)
    variable q(1)
    maximize(p' * individual_cost - q^2 * prox_param/2);
    subject to 
         norm((p - ref_p) ./ ref_p) <= sqrt(radius)

         norm(p - prox_center) <= q
        p' * ones(k, 1) == 1;
        p >= 0
cvx_end
toc

full(p);
% smoothing_cost = q^2 * prox_param/2
x2 = Problem.X2Uncertainty('Euclidean');

x2.radius = radius;
x2.generateData(k);
% x2.solveForP(individual_cost)
tic
% [total_cost, next_p, ~] = x2.solveForP(individual_cost);
[next_p, total_cost, p_radius] = x2.projectP(prox_param, prox_center, individual_cost);
toc
% Ptable
norm(next_p - p)
% norm(P - Ptable, 'fro')
total_cost - cvx_optval
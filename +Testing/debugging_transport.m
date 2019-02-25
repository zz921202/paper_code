% k = 5000;
k = 2
distance_mat = [0 16.5476; 16.5476 0];
epsilon = 0.1;
% rng(200)
% individual_cost = rand(k, 1) * 1;
individual_cost = [8; 300];
ref_p = ones(k, 1) ./ k;
DELTA = 1e-16;
prox_center = rand(k, k) + ones(k, k) .* 0.1;
diag(sum(prox_center, 2))
prox_center = prox_center * diag(1./sum(prox_center, 2))
prox_param = 1
tic
cvx_begin quiet
    variable P(k, k)
    variable q(1)
    maximize(sum(P) * individual_cost -  q^2 * k * prox_param/2);
    subject to 
         norm(P - prox_center, 'fro') <= q
        P * ones(k, 1) == ref_p;
        sum(sum(P .* distance_mat)) <= epsilon
        P >= 0
cvx_end
toc

% full(P)

transport = Problem.TransportUncertainty('Euclidean');
transport.distance_table = distance_mat;
transport.epsilon = epsilon;
transport.generateData(k);
% transport.solveForP(individual_cost)
tic
% [total_cost, p, ~] = transport.solveForP(individual_cost);
[next_p, total_cost, ~, Ptable] = transport.projectP(prox_param, prox_center, individual_cost);
toc
% Ptable
% norm(sum(P) - p')
norm(P - Ptable, 'fro')
total_cost - cvx_optval
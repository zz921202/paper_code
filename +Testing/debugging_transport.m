k = 5000;
distance_mat = ones(k, k) - eye(k);
epsilon = 1e-1;
individual_cost = (1:k)';
ref_p = ones(k, 1) ./ k;
tic
cvx_begin quiet
    variable P(k, k)
    maximize(sum(P) * individual_cost)
    subject to 
        P * ones(k, 1) == ref_p;
        sum(sum(P .* distance_mat)) <= epsilon
        P >= 0
cvx_end
toc



transport = Problem.TransportUncertainty('Entropy');
transport.distance_table = distance_mat;
transport.epsilon = epsilon;
transport.generateData(k);
% transport.solveForP(individual_cost)
tic
[total_cost, p] = transport.solveForP(individual_cost);
toc

norm(sum(P) - p')
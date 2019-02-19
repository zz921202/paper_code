% testing euclidean projector 




% eu_projector = Helper.EntropyProjector();
% eu_projector.setConstraint(A, b);
k = 10000

ref_p = rand(k, 1);
% ref_p = ones(k, 1);
c = (1:k)';
ref_p = ref_p ./ sum(ref_p);
alpha = 0;
beta = 100;
ref_projector = Helper.BoxPProjector.BoxPEntropyProjector();
ref_projector.setConstraint(ref_p * alpha, ref_p * beta);
% [val, soln]= eu_projector.solve(c);
tic
[ref_soln, ref_obj] = ref_projector.project(100, ref_p, c);
toc

mosek_projector = Helper.EntropyProjector();

A = [-ones(1, k); ones(1, k); -eye(k); eye(k); eye(k)];
b = [-1 ; 1; -beta * ref_p; alpha * ref_p; zeros(k, 1)];
mosek_projector.setConstraint(A, b);
tic
[mos_soln, mos_obj] = mosek_projector.project(100, ref_p, c);
toc


eu_projector = Helper.BoxPProjector.VectorBoxPEntropyProjector();
eu_projector.setUpperLowerConstraint(ref_p * alpha, ref_p * beta);


tic
[eu_soln, eu_obj] = eu_projector.project(100, ref_p, c);
toc
% sum(abs(eu_obj - ref_obj))
sum(abs(eu_soln - mos_soln))
sum(abs(ref_soln - mos_soln))
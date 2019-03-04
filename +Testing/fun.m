% self = ref_problem
% n1 = self.n1, n2 = self.n2
% cvx_begin quiet
%     variable  y(n2, self.k)
%     variable x(n1, 1)
%     variable t
%     variable f(self.k,1)
%     variable fplus(self.k, 1)
%     % minimize( self.c' * x)
%     minimize( (1 - self.alpha)*( t + beta_hat * self.reference_p' * fplus) + self.alpha * self.reference_p' * f + self.c'  * x)
%     subject to
% %                         self.A * x >= self.b;
%     for i =1 :self.k
%         f(i) >= (self.eks{i})' * y(:, i)
%         fplus(i) >= f(i) - t
%         fplus(i) >= 0
%         self.Wks{i} * y(:, i) == self.Tks{i} * x + self.dks{i}
%     end
%     x >= 0
%     y >= 0
% cvx_end  

a = 0;
b = rand(100);
ticBytes(gcp);
parfor i = 1:100
    a = a + sum(b(:, i));
end
tocBytes(gcp)
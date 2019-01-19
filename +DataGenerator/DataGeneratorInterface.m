classdef DataGeneratorInterface < handle

    properties

        EIGENMIN = 1;
    end
 % generator for linear 2 stage problem only
    methods(Abstract)
    [c, A, b, dks, eks, Wks, Tks, Mt, omega_x2, omega_pi2] = generateData(self,k) % k is the number of scenarios

    [n1, m1, n2, m2] = getProblemDimension(self);
    end

    methods 
        function [mat, largest_eig] =randMatrixEig(self, m, n, start_value, end_value)
            mat = rand(m, n) * (end_value - start_value) + start_value;
            [U,S,V] = svd(mat);
            l = min(m, n);
            add_eig = zeros(size(S));
            add_eig(1:l, 1:l) = self.EIGENMIN;
            new_S = S+ add_eig;
            mat = U * (new_S) * V';
            largest_eig  = new_S(1, 1);
        end 

        function [mat, largest_eig] = randMatrix(self, m, n, start_value, end_value)
            mat = rand(m, n) * (end_value - start_value) + start_value;
            largest_eig = norm(mat); 
        end

    
         function vec = genVector(self, m, start_value, end_value) 
            vec = rand(m, 1) * (end_value - start_value) + start_value;
        end 

        function vs = genVectors(self, n, k, start_value, end_value)
            vs = {};
            for ind = 1: k
                vs = [vs, self.genVector(n,  start_value, end_value)];
            end
        end
    end

end
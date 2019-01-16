classdef RandomMatrixGeneratorInterface

    methods(Abstract)
        [Mt, Tks] = generateTechMatrices(self, m, n, k)%Tks should be negative in our formulation 
        Wks = generateRecourseMatrices(self, m, n, k)
        A = generateAMatrix(self, m, n)
    end

    methods (Access = protected)
        function [mat] =randMatrix(obj, m, n, start_val, end_val)
            eig_lower_bound = 1;
            eig_upper_bound = 10;
            k = min(m, n);
            eig_vals = obj.helper_rand_vec(k, eig_lower_bound, eig_upper_bound);
            mat = zeros(m, n);
            for i = 1:k
                cur_vet = obj.helper_rand_vec(m, start_val, end_val);
                cur_hor = obj.helper_rand_vec(n, start_val, end_val);
                mat = mat + cur_vet * cur_hor' .* eig_vals(i);
            end
            
        end 


        function [mat, largest_eig] =randMatrixWithEig(obj, m, n, start_val, end_val)
            eig_lower_bound = 1;
            eig_upper_bound = 10;
            k = min(m, n);
            eig_vals = obj.helper_rand_vec(k, eig_lower_bound, eig_upper_bound);
            mat = zeros(m, n);
            largets_eig = 0;
            for i = 1:k
                cur_vet = obj.helper_rand_vec(m, start_val, end_val);
                cur_hor = obj.helper_rand_vec(n, start_val, end_val);
                mat = mat + cur_vet * cur_hor' .* eig_vals(i);
                largest_eig = max(largest_eig, norm(cur_vet) * norm(cur_hor) * eig_vals(i));
            end
        end 
    
    
    
         function vec = helper_rand_vec(obj, m,  start_val, end_val) 
            vec = rand(m, 1) * (end_val - start_val) + start_val;
        end 
    end

end
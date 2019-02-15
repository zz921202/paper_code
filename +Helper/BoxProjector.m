classdef BoxProjector < Helper.LinearProjectorInterface
    properties
        lower_bound;
        upper_bound;
%         n;
    end
    
    methods

        function [val, soln] = solve(self, grad)
            soln = zeros(self.n, 1);
            % init_pi = grad;
            pos_ind = (grad < 0);
            neg_ind = (grad >= 0);
            soln(pos_ind) = self.upper_bound(pos_ind);
            soln(neg_ind) = self.lower_bound(neg_ind);
            % soln
            val = grad' * soln;
        end

        function setUpperLowerConstraint(self, lower_bound, upper_bound)
            % Ax >= b
            self.lower_bound = lower_bound;
            self.upper_bound = upper_bound;
            self.n = length(lower_bound);
        end



        function [soln, val] = project(self,prox_param, prox_center, grad)
            % prox_param = prox_param
            % prox_center = prox_center
            % grad = grad
            if prox_param == 0
                [val, soln] = self.solve(grad);
                return
            end
            temp_next_x = prox_center - 1/prox_param * grad;
            n_zero = zeros(self.n, 1);
            next_x = temp_next_x + max([self.lower_bound - temp_next_x, n_zero], [], 2) - max([temp_next_x - self.upper_bound, n_zero], [], 2);
            val = grad' * next_x + prox_param * norm(next_x - prox_center)^2 / 2;
            soln = next_x;
        end

        function str =getDistanceName(self)
            str = 'CustomEuclidean';
        end

    end

end
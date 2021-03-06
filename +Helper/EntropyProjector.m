classdef EntropyProjector < Helper.LinearProjectorInterface
    
    

    methods
        
        function [soln, obj] = project(self, prox_param, prox_center, grad)
        % solves min <x, grad>  + U(prox_center, x) * 
                
        k = length(grad);
            cvx_begin quiet
                variable p(self.n)
                minimize( grad' * p + ones(1, k) * (kl_div(p, prox_center)) * prox_param)
                subject to 
                    self.model.A * p >= self.model.rhs
            cvx_end

            soln = p;
            obj = cvx_optval;
        end

        function str =getDistanceName(self)
            str = 'Entropy';
        end
    end
end
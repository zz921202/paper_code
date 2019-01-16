classdef EntropyProjector < LinearProjectorInterface
    


    methods
        
        function [soln, obj] = project(self, prox_param, prox_center, grad)
        % solves min <x, grad>  + U(prox_center, x) * 
                
        
            cvx_begin quiet
                variable p(self.n)
                minimize( grad' * p + ones(1, obj.k) * (kl_div(prev_p, p)) * prox_param)
                subject to 
                    self.model.A * p >= self.model.rhs
            cvx_end
            
            soln = p;
            obj = cvx_optval;
        end
end
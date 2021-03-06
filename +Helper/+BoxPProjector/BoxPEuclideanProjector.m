classdef BoxPEuclideanProjector < Helper.BoxPProjector.BoxPProjector
    
    properties
        
        prox_param;
    end
    methods


        function compute_base(self, prox_param, prox_center, grad)
            % update base for all components
            self.prox_param = prox_param;
%             self.p_components
            for ind = 1: self.k
                cur_p_component = self.p_components(ind);
                cur_p_component.base = prox_center(ind) - 1/prox_param * grad(ind);
            end
        end

        
        function p =  computeP(self, lambda) % compute the probility corresponding to the chosen lambda value
            p = zeros(self.k, 1); 
            for ind  = 1: self.k
                cur_p_component = self.p_components(ind);
                cur_p = cur_p_component.base + lambda / self.prox_param;
                p(ind) = cur_p_component.getPValue(cur_p);
            end
        end
        function lam = computeLambda(self, base, target, num_active_terms)
        % solves lambda such that sum(p(base, lambda)) == target
            
            lam = (target - base) * self.prox_param/num_active_terms;
        end

        function objective = computeCost(self, prox_param, prox_center, grad, p)
            objective = grad' * p + prox_param * 1/2 * norm(p - prox_center)^2;
        end

    end

end
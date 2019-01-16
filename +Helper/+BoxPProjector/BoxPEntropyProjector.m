classdef BoxPEntropyProjector < Helper.BoxPEntropyProjector.BoxPEntropyProjector
    
    properties
        DELTA = 1e-16;
        prox_param;
    end
    methods


        function compute_base(self, prox_param, prox_center, grad)
            % update base for all components
            self.prox_param = prox_param;
            for ind = 1: self.k
                cur_p_component = self.p_components{ind};
                cur_p_component.base = exp(log(prox_center(ind) + self.DELTA) -1 - 1/prox_parm*(grad(ind)));
            end
        end

        
        function p =  computeP(self, lambda) % compute the probility corresponding to the chosen lambda value
            p = zeros(self.k, 1); 
            for ind  = 1: self.k
                cur_p_component = self.p_components{ind};
                cur_p = cur_p_component.base * exp(1/self.prox_param * lambda) - self.DELTA;
                p(ind) = cur_p_component.getPValue(cur_p);
            end
        end
        function lam = computeLambda(self, base, target, num_active_terms)
        % solves lambda such that sum(p(base, lambda)) == target
            real_target = target + self.DELTA * num_active_terms;
            lam = log(real_target / base) * self.prox_param;
        end

        function objective = computeCost(self, prox_param, prox_center, grad, p)
            objective = grad' * p + prox_param * ((p + self.DELTA)' * log(p + self.DELTA) - (p + self.DELTA)' * log(prox_center + self.DELTA));
        end

    end

end
classdef VectorBoxPEntropyProjector < Helper.BoxPProjector.VectorBoxPProjector

    properties 
        DELTA = 1e-16;
    end

    methods 
        
        function cost = computeCost(self, prox_param, prox_center, grad, p)
            cost = grad' * p + prox_param * ((p + self.DELTA)' * log(p + self.DELTA) - (p + self.DELTA)' * log(prox_center + self.DELTA));
        end
        function lam = computeLambda(self, active_base_indicator, target_val)
            % active_base_indicator
            % self.base(active_base_indicator)
            lam = (log(target_val + sum(active_base_indicator) * self.DELTA) - Helper.logsumexp(self.base(active_base_indicator))) * self.prox_param;
        end

        function prob = computeActiveProb(self, active_base_indicator, lambda)
        % only update probility for the active base_idx
%             prob = zeros(self.k, 1);
            exponent = self.base + lambda / self.prox_param;
            prob = exp(exponent) - self.DELTA;
            prob(~active_base_indicator) = 0;
        end


        function combined_lambda = computeBase(self, prox_param, prox_center, grad)
            self.base = log(prox_center + self.DELTA) -1 - grad ./ prox_param;
%             base = self.base
%             self.lower_prob + self.DELTA
            lower_lambda = (log(self.lower_prob + self.DELTA) - self.base) * prox_param;
            upper_lambda = (log(self.upper_prob + self.DELTA) - self.base)  * prox_param;
            combined_lambda = [lower_lambda ; upper_lambda];
        end

        function str= getDistanceName(self)
            str = 'BoxEntropy2';
        end
    end
end
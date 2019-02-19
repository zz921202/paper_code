classdef VectorBoxPEuclideanProjector < Helper.BoxPProjector.VectorBoxPProjector



    methods 
        
        function cost = computeCost(self, prox_param, prox_center, grad, p)
            cost = grad' * p + prox_param * 1/2 * norm(p - prox_center)^2;
        end
        function lam = computeLambda(self, active_base_indicator, target_val)
            % active_base_indicator
            % self.base(active_base_indicator)
            active_base_val = sum(self.base(active_base_indicator));
            lam = (target_val - active_base_val)/sum(active_base_indicator) * self.prox_param;
        end

        function prob = computeActiveProb(self, active_base_indicator, lambda)
        % only update probility for the active base_idx
%             prob = zeros(self.k, 1);
            prob = self.base + lambda / self.prox_param;
            prob(~active_base_indicator) = 0;
        end


        function combined_lambda = computeBase(self, prox_param, prox_center, grad)
            self.base = prox_center - grad ./ prox_param;
%             base = self.base
%             self.lower_prob + self.DELTA
            lower_lambda = (self.lower_prob - self.base) * prox_param;
            upper_lambda = (self.upper_prob - self.base)  * prox_param;
            combined_lambda = [lower_lambda ; upper_lambda];
        end

        function str= getDistanceName(self)
            str = 'BoxEuclidean2';
        end
    end
end
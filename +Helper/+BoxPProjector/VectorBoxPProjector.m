classdef VectorBoxPProjector < Helper.LinearProjectorInterface
    properties
        lower_prob_aug;
        upper_prob_aug;
        lower_prob;
        upper_prob;
        k;
        combined_sorted_lambda_val;
        combined_sorted_lambda_upper_idx;
        combined_sorted_lambda_lower_idx;
        base;
        prox_param;
        prox_center;
        grad;
    end
    
    methods



        function setUpperLowerConstraint(self, lower_bound, upper_bound)
            % Ax >= b
            self.lower_prob_aug = [lower_bound; 0];
            self.upper_prob_aug = [upper_bound; 0];
            self.lower_prob = lower_bound;
            self.upper_prob = upper_bound;
            self.k = length(lower_bound);
        end



        function [soln, val] = project(self,prox_param, prox_center, grad)
            self.prox_param = prox_param;
            self.prox_center = prox_center;
            self. grad = grad;
            if prox_param == 0
                [val, soln] = self.solve(grad);
                return
            end
            combined_lambda = self.computeBase(prox_param, prox_center, grad);            
            lower_idx = [(1: self.k)'; (self.k + 1) * ones(self.k, 1)];
            upper_idx = [(self.k + 1) * ones(self.k, 1); (1: self.k)'];
            [self.combined_sorted_lambda_val, sorting_idx] = sort(combined_lambda);
            self.combined_sorted_lambda_lower_idx = lower_idx(sorting_idx);
            self.combined_sorted_lambda_upper_idx = upper_idx(sorting_idx);

            start_idx = 1; 
            end_idx = 2 * self.k;
            while true
                [prob, start_idx, end_idx] = self.search(start_idx, end_idx);
                if ~isempty(prob)
                    break
                end
            end
            soln = prob;
            val = self.computeCost(prox_param, prox_center, grad, prob);
        end


        function [prob, next_start_idx, next_end_idx] = search(self, start_idx, end_idx)
            % prob is empty if current search fails
            prob = [];
            next_start_idx = start_idx; next_end_idx = end_idx;
            search_idx = floor((start_idx + end_idx)/2);
            cur_lam = self.combined_sorted_lambda_val(search_idx);
            % fprintf('current search idx is %d \n', search_idx)
            if search_idx + 1 > self.k * 2
                next_lam = Inf;
            else
                next_lam = self.combined_sorted_lambda_val(search_idx  + 1);
            end
            fixed_lower_idx = self.combined_sorted_lambda_lower_idx(search_idx + 1 : end);
            fixed_upper_idx = self.combined_sorted_lambda_upper_idx(1 : search_idx);
            active_base_aug = ones(self.k+1, 1, 'logical');
            active_base_aug(fixed_lower_idx) = 0;
            active_base_aug(fixed_upper_idx) = 0;
            active_base_indicator = active_base_aug(1 : self.k);
            fixed_val = sum(self.lower_prob_aug(fixed_lower_idx)) + sum(self.upper_prob_aug(fixed_upper_idx));
            % active_base_indicator
            if fixed_val > 1
                next_end_idx = search_idx;
                return
            end

            if ~any(active_base_indicator)
                if fixed_val < 1
                    next_start_idx = search_idx;
                else
                    aug_prob = zeros(self.k + 1, 1);
                    aug_prob(fixed_lower_idx) = self.lower_prob_aug(fixed_lower_idx);
                    aug_prob(fixed_upper_idx) = self.upper_prob_aug(fixed_upper_idx);
                    prob = aug_prob(1 : self.k);
                end
                return
            end

            guess_lambda  = self.computeLambda(active_base_indicator, 1 - fixed_val);
            if guess_lambda < cur_lam - 1e-10
                % fprintf('updating upper %d\n', search_idx)
                next_end_idx = search_idx;
                % next_start_idx = start_idx;
            elseif guess_lambda > next_lam + 1e-10
                % fprintf('updating lower %d\n', search_idx+1)
                % next_end_idx = end_idx;
                next_start_idx = search_idx + 1;
            else 
                active_prob = self.computeActiveProb(active_base_indicator, guess_lambda);
                aug_prob = [active_prob; 0];
                aug_prob(fixed_lower_idx) = self.lower_prob_aug(fixed_lower_idx);
                aug_prob(fixed_upper_idx) = self.upper_prob_aug(fixed_upper_idx);
                prob = aug_prob(1 : self.k);
            end
        end

        function str =getDistanceName(self)
            str = 'CustomEuclidean';
        end


    end

    methods(Abstract)
        computeBase(self, prox_param, prox_center, grad)
        % update base,  return combined lambda := [lower lambda, upper lambda]
        cost = computeCost(self, prox_param, prox_center, grad, p)
        lam = computeLambda(self, active_base_indicator, target_val)
        prob = computeActiveProb(self, active_base_indicator, lambda)
        % only update probility for the active base_idx
        
    end

end
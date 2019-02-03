classdef BoxPProjector < handle
    properties
        p_components = {};
        k;
        lower_prob_sum;
    end
    
    methods
        function setConstraint(self, lower_probs, upper_probs)
            k = length(upper_probs);
            self.k = k;
            self.lower_prob_sum = sum(lower_probs);
            p_components = {};
            % lower_probs
            % upper_probs
            for i =1 : k
                cur_p_component = Helper.BoxPProjector.PComponent();
                cur_p_component.upper_p = upper_probs(i);
                cur_p_component.lower_p = lower_probs(i);
                p_components = [p_components, cur_p_component];
            end
            self.p_components = p_components;
        end
        
        function [p, objective] = project(self, prox_param, prox_center, grad)
            % solves min <grad, p> + prox_param * U(prox_cneter, p)
            
            % initialize p components
            self.compute_base(prox_param, prox_center, grad);
            upper_lambdas = [];
            lower_lambdas = [];
            for i = 1:self.k
                cur_p_component = self.p_components(i);
                cur_p_component.upper_lambda = self.computeLambda(cur_p_component.base, cur_p_component.upper_p, 1);
                cur_p_component.lower_lambda = self.computeLambda(cur_p_component.base, cur_p_component.lower_p, 1);
                upper_lambdas = [upper_lambdas, cur_p_component.upper_lambda];
                lower_lambdas = [lower_lambdas, cur_p_component.lower_lambda];
            end
            
            [~, lower_lambda_sort_ind] = sort(lower_lambdas);
            [~, upper_lambda_sort_ind] = sort(upper_lambdas);
            
            % initialize the sequential search alg
            fixed_value = self.lower_prob_sum;
            acc_active_base = 0;
            num_active_base = 0;
            old_lambda = Inf;
            cur_lower_ind = 1;
            cur_upper_ind = 1;
            % break point search algorithm
            while true
                 
                if cur_lower_ind <= self.k
                    cur_lower_p_component = self.p_components(lower_lambda_sort_ind(cur_lower_ind));
                end
%                 cur_upper_ind
                cur_upper_p_component = self.p_components(upper_lambda_sort_ind(cur_upper_ind));
                if cur_lower_ind <= self.k && cur_lower_p_component.lower_lambda < cur_upper_p_component.upper_lambda
                    if old_lambda <= cur_lower_p_component.lower_lambda
                        p = self.computeP(old_lambda);
                        break
                    else
                        fixed_value = fixed_value - cur_lower_p_component.lower_p;
                        acc_active_base = acc_active_base + cur_lower_p_component.base;
                        num_active_base = num_active_base + 1;
                        cur_lower_ind  = cur_lower_ind + 1;
                    end
                    
                    
                else
                    if old_lambda <= cur_upper_p_component.upper_lambda

                        p = self.computeP(old_lambda);
                        break
                    else
                        fixed_value = fixed_value + cur_upper_p_component.upper_p;
                        acc_active_base = acc_active_base - cur_upper_p_component.base;
                        num_active_base = num_active_base - 1;
                        cur_upper_ind  = cur_upper_ind + 1;
                    end
                end
                % cur_upper_ind
                % num_active_base
                % acc_active_base
                % cur_lower_ind
%                 acc_active_basex
                old_lambda = self.computeLambda(acc_active_base, 1 - fixed_value, num_active_base);
                % fixed_value
                % target = 1 - fixed_value
                % self.computeP(old_lambda)

                
            end
            objective = self.computeCost(prox_param, prox_center, grad, p);
        end
        
        
        
        
    end
    
    
    methods(Abstract)
        compute_base(self, prox_param, prox_center, grad)
        % update base for all components
        p =  computeP(self, lambda) % compute the probility corresponding to the chosen lambda value
        lam = computeLambda(self, base, target, num_active_terms)
        objective = computeCost(self, prox_param, prox_center, grad, p)
        % solves lambda such that sum(p(base, lambda)) == target
        
    end
end
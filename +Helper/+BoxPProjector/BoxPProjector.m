classdef BoxPProjector < handle
    properties
        p_components = {};
        k;
        lower_prob_sum;
    end
    
    methods
        function setConstraint(self, lower_probs, upper_probs)
            [~, k] = size(upper_probs);
            self.k = k;
            self.lower_prob_sum = sum(lower_probs);
            p_components = {};
            for i =1 : k
                cur_p_component = Helper.BoxPPropector.PComponent();
                cur_p_component.upper_p = upper_probs(i);
                cur_p_component.lower_p = lower_probs(i);
                p_components = [p_components, cur_p_component];
            end
            self.p_components = p_components;
        end
        
        function p = project(self, prox_param, prox_center, grad)
            % solves min <grad, p> + prox_param * U(prox_cneter, p)
            
            % initialize p components
            self.compute_base(prox_param, prox_center, grad);
            upper_lambdas = [];
            lower_lambdas = [];
            for i = 1:k
                cur_p_component = self.p_components{i};
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
                cur_lower_p_component = self.p_components{lower_lambda_sort_ind(cur_lower_ind)};
                cur_upper_p_component = self.p_components{upper_lambda_sort_ind(cur_upper_ind)};
                if cur_lower_p_component.lower_lambda < cur_upper_p_component.upper_lamda
                    if old_lambda <= cur_lower_p_component.lower_lambda
                        p = computeP(lambda);
                        break
                    else
                        fixed_value = fixed_value - cur_lower_p_component.lower_p;
                        acc_active_base = acc_active_base + cur_lower_p_component.base;
                        num_active_base = num_active_base + 1;
                        cur_lower_ind  = cur_lower_ind + 1;
                    end
                    
                    
                else
                    if old_lambda <= cur_upper_p_component.upper_lambda
                        p = computeP(lambda);
                        break
                    else
                        fixed_value = fixed_value + cur_lower_p_component.upper_p;
                        acc_active_base = acc_active_base - cur_upper_p_component.base;
                        num_active_base = num_active_base - 1;
                        cur_upper_ind  = cur_upper_ind + 1;
                    end
                end
                
                old_lambda = self.computeLambda(acc_active_base, 1 - fixed_value, num_active_base);
            end
        end
        
        
        
        
    end
    
    
    methods(Abstract)
        compute_base(self, prox_param, prox_center, grad)
        % update base for all components
        p =  computeP(self, lambda) % compute the probility corresponding to the chosen lambda value
        lam = computeLambda(self, base, target, num_active_terms)
        % solves lambda such that sum(p(base, lambda)) == target
        
    end
end
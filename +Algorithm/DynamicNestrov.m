classdef DynamicNestrov < Algorithm.NestrovTypeAlgorithm 
    properties
        % omega_p_ratios = [1, 1e-1, 1e-2];
        omega_p_ratios = 10 .^ (-2:0);
        % omega_pi_ratios = 3.^(-4:1);%[ 1e-2, 1e-3, 1e-1];
        omega_pi_ratios = 10 .^ (-4:0);
        % omega_x_ratios = [1];
        omega_x_ratios = 10.^(-1:1);%[1e-1, 1e-2, 1e-3];
        MAXITER;
    end
    methods
    
        function self = DynamicNestrov(problem_data, terminator)
            self@Algorithm.NestrovTypeAlgorithm(problem_data, terminator);
            
        end

        function [gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = getProxParams(self)

            t = self.num_iters;
            CONSTANT = self.cp * self.cur_M_pi * self.cur_omega_p + self.cur_omega_pi;
%             Mt = self.cur_omega_x
            c = self.cur_omega_x / (self.Mt * CONSTANT);
            mu_t = c / t;
            alpha_t = 1 / t;
            
            gamma_t_inv = 1 / c * self.omega_x_ratios(self.x_choice) ;
            
            mu_pi_t = mu_t * self.Mt^2 * CONSTANT / self.cur_omega_pi * self.omega_pi_ratios(self.pi_choice);
            mu_p_t = mu_pi_t * self.getOptimalSmoothingRatio() *self.omega_p_ratios(self.p_choice);
        end

        function str = getName(self)
            str = 'dynNes';
        end
    end
end
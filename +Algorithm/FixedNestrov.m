classdef FixedNestrov < Algorithm.NestrovTypeAlgorithm 
    properties
        omega_p_ratios = [1, 0.5, 0.1, 0.01];
        omega_pi_ratios = [1, 0.5, 0.1, 0.01];
        MAXITER;
    end
    methods


        function [gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = getProxParams(self)

            t = self.num_iters;
            CONSTANT = self.cp * self.cur_M_pi * self.cur_omega_p + self.cur_omega_pi;
            c = self.cur_omega_x / (self.Mt * CONSTANT);
            mu_t = c / t;
            alpha_t = 1 / t;
            gamma_t_inv = 1 / c;
            mu_p_t = mu_t * self.cp * self.cur_M_pi * self.curMt^2 * CONSTANT / self.cur_omega_p;
            mu_pi_t = mu_t * self.curMt^2 * CONSTANT / self.cur_omega_pi;
        end
    end
end
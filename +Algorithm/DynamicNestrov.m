classdef DynamicNestrov < Algorithm.NestrovTypeAlgorithm 
    properties
        omega_p_ratios = [1];
        omega_pi_ratios = [1];
        MAXITER;
    end
    methods
        function self = FixedNestrov(problem_data, terminator)
            self = self@Algorithm.NestrovTypeAlgorithm(problem_data, terminator);
            ter = Algorithm.Terminator.MaxIterTerminator();
            self.MAXITER = ter.MAXITERATION;
        end

        function [gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = getProxParams(self)
            t = self.num_iters;
            CONSTANT = self.cp * self.cur_M_pi * self.cur_omega_p + self.cur_omega_pi;
            mu = 2 * self.cur_omega_x / (self.MAXITER * self.Mt * CONSTANT);
            alpha_t = 2 / (t+1);
            gamma_t_inv = 2 / (t * mu);
            mu_p_t = mu * self.cp * self.cur_M_pi * self.curMt^2 * CONSTANT / self.cur_omega_p;
            mu_pi_t = mu * self.curMt^2 * CONSTANT / self.cur_omega_pi;
        end
    end
end
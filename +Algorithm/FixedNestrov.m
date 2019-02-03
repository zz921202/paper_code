classdef FixedNestrov < Algorithm.NestrovTypeAlgorithm 
    properties
        omega_p_ratios = [1e-1, 1e-2, 1e-3, 1e-4];
        omega_pi_ratios = [ 1e-1, 1e-2, 1e-3, 1e-4];
        omega_x_ratios = [ 1e-1, 1e-2, 1e-3, 1e-4]
        MAXITER;
    end
    methods
        function self = FixedNestrov(problem_data, terminator)
            self = self@Algorithm.NestrovTypeAlgorithm(problem_data, terminator);
            ter = Algorithm.Terminator.MaxIterTerminator();
            self.MAXITER = ter.MAXITERATION;
        end

        function [gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = getProxParams(self)
            T = self.MAXITER;

%             self.cur_omega_p = 0;%TODO
            CONSTANT = self.cp * self.cur_M_pi * self.cur_omega_p + self.cur_omega_pi;
            mu = 2 * self.cur_omega_x  / ( self.Mt * CONSTANT * T );
            alpha_t = 2 / (self.num_iters+1);
            
            gamma_t_inv = 2   / (self.num_iters * mu ) * self.omega_x_ratios(self.x_choice);
            % hacking * sqrt(self.Mt * self.cur_M_pi)
            mu1 = mu;
            % mu1 = mu / (  1e3);
            
            mu_p_t = mu1 * self.cp * self.cur_M_pi * self.Mt^2 * CONSTANT / self.cur_omega_p * self.omega_p_ratios(self.p_choice);
%             mu_p_t = 0;%TODO
            mu_pi_t = mu1 * self.Mt^2 * CONSTANT / self.cur_omega_pi * self.omega_pi_ratios(self.pi_choice);
        end

        function str = getName(self)
            str = 'fixNes';
        end
    end
end
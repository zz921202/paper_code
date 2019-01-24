classdef ABPAlgorithm < Algorithm.DSLAlgorithm 
    
    methods
            function init_phase(self)
%             v_lower = self.phase_lower;
            [v_upper, ~] = self.problem_data.evalX(self.phase_x0);
     
%             phase_l = self.BETA * v_lower + (1 - self.BETA) * v_upper
            self.phase_l = self.BETA * self.phase_lower + (1 - self.BETA) * v_upper;
            
            % mu = (self.THETA * (v_upper - self.phase_l)) / (2 * self.cur_omega_pi2_estimate *(1 + sqrt(2)) * self.cur_p_estimate * self.cp);
            % self.phase_mu_p = mu * self.cur_omega_pi2_estimate * self.cp * sqrt(2) / self.cur_p_estimate;
            % self.phase_mu_pi = mu;
            self.phase_mu_p = 0;
            self.phase_mu_pi = 0;

            self.phase_lower_threshold = self.phase_l - self.THETA * (self.phase_l - self.phase_lower);
%             self.phase_lower_threshold
            self.phase_upper_threshold = self.phase_l + self.THETA*(v_upper - self.phase_l);
        end

        function self  = ABPAlgorithm(problem_data, terminator)
            self@Algorithm.DSLAlgorithm(problem_data, terminator);
            self.omega_p_ratios = [1];
            self.num_param_choices =  1;
        end
    end
end
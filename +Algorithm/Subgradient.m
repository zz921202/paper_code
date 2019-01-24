classdef Subgradient < Algorithm.FixedNestrov
    properties
        % omega_p_ratios = [1, 0.5, 0.1, 0.01];
        % omega_pi_ratios = [1, 0.5, 0.1, 0.01];
        x_history = [];
        
    end

    methods




        function [x, obj_val, est_gap, true_gap, time_elapsed, num_iters] = run(self)
            % refer to page 5 of paper ALgorithm 1
            self.cleanUp();
            cur_x_temp = self.cur_x;
            smooth_p = self.smooth_p;
            smooth_pis = self.smooth_pis;
            [gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = self.getProxParams();
            mu_pi_t = 0;
            mu_p_t = 0;
            fake_cost = self.problem_data.evalFakeCost(self.problem_data.ref_x, mu_p_t, smooth_p, mu_pi_t, smooth_pis);
            temp_cost = self.problem_data.evalFakeCost(cur_x_temp, mu_p_t, smooth_p, mu_pi_t, smooth_pis);
            while ~self.terminator.terminate()
                
                self.mytimer.start()
                x_under = self.cur_x;
                [self.cur_pis, indi_costs] = self.problem_data.projectPis(mu_pi_t, smooth_pis, x_under);
%                 smooth_p
%                 indi_costs
                [self.cur_p, secon_grad, secon_cost] = self.problem_data.projectP(mu_p_t, smooth_p, indi_costs, self.cur_pis);
                x_grad = secon_grad + self.problem_data.c;
                cur_fake_gap = temp_cost - fake_cost ;
                self.cur_true_gap;
                stepsize = self.cur_true_gap / norm(x_grad)^2;
                x_root = self.cur_x  - stepsize * x_grad;

                [self.cur_x, ~] = self.problem_data.projectX(1, x_root, -self.problem_data.c);
                
                self.x_history = [self.x_history, self.cur_x];
                % self.cur_x = (1 - alpha_t) * self.cur_x + alpha_t * cur_x_temp;
                self.mytimer.pause();
                % fake_cost = self.problem_data.evalFakeCost(cur_x_temp, mu_p_t, smooth_p, mu_pi_t, smooth_pis);
                self.recordObjTime();
                [temp_cost, ~] = self.problem_data.evalX(x_under);

                self.temp_cost_history = [self.temp_cost_history, temp_cost];
                % self.fake_cost_history = [self.fake_cost_history, fake_cost];
            end
            x = self.cur_x; obj_val = self.cur_obj; est_gap = self.cur_est_gap;  true_gap=self.cur_true_gap; time_elapsed=self.time_elapsed; num_iters=self.num_iters;
        end

        function self = Subgradient(problem_data, terminator)
            self@Algorithm.FixedNestrov(problem_data, terminator);
            self.num_param_choices = length(self.omega_x_ratios) * length(self.omega_p_ratios) * length(self.omega_pi_ratios);
        end

        function getSubGradient(self, s)
            smooth_p = self.smooth_p;
            smooth_pis = self.smooth_pis;
            x_under = s;
            [gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = self.getProxParams();
            [cur_pis, indi_costs] = self.problem_data.projectPis(mu_pi_t, smooth_pis, x_under)
            pis = cur_pis{1}
            [cur_p, secon_grad, secon_cost] = self.problem_data.projectP(mu_p_t, smooth_p, indi_costs, cur_pis)

            x_grad = secon_grad + self.problem_data.c
        end

    end


end
classdef NestrovTypeAlgorithm < Algorithm.FOAlgorithm
    properties
        % omega_p_ratios = [1, 0.5, 0.1, 0.01];
        % omega_pi_ratios = [1, 0.5, 0.1, 0.01];
        omega_x_estimates = [10, 100, 1000, 10000];
        cur_omega_p;
        cur_omega_pi;
        cur_M_pi;
        cur_omega_x;
        cur_x_temp; %xt while cur_x refers to the accumulated a.k.a x_bar

    end

    methods



        function setGridParam(self, index)
            num_p_choices = length(self.omega_p_ratios);
            num_pi_choices = length(self.omega_pi_ratios);
            num_ppi = num_p_choices * num_pi_choices;
            x_choice = floor((index-0.1) / (num_ppi)) + 1;
            ppi_index = index - num_ppi * x_choice;
            p_choice = floor((ppi_index-0.1) / num_p_choices) + 1;
            pi_choice = ppi_index - p_choice * num_p_choices;
            self.cur_omega_x = self.omega_x_estimates(x_choice);
            self.cur_omega_p = self.omega_p_ratios(p_choice) * sqrt(self.problem_data.getConservativePBregDist());
            self.cur_omega_pi = self.omega_pi_ratios(pi_choice) * sqrt(self.problem_data.getConservativePiBregDist());
            self.cur_M_pi = sqrt(2) * self.cur_omega_pi;% CHEATING, valid only when smooth_pi = 0
        end

        function str = showGridParam(self, index)
            num_p_choices = length(self.omega_p_ratios);
            num_pi_choices = length(self.omega_pi_ratios);
            num_ppi = num_p_choices * num_pi_choices;
            x_choice = floor((index-0.1) / (num_ppi)) + 1;
            ppi_index = index - num_ppi * x_choice;
            p_choice = floor((ppi_index-0.1) / num_p_choices) + 1;
            pi_choice = ppi_index - p_choice * num_p_choices;
            str = sprintf('param choice x %s, p_ratio: %s, pi_ratio %s', num2str(self.omega_x_estimates(x_choice)), 
                num2str(self.omege_p_ratios(p_choice)), num2str(self.omega_pi_ratios(pi_choice)));
            disp(str)
        end

        function [x, obj_val, est_gap, true_gap, time_elapsed, num_iters] = run(self)
            % refer to page 5 of paper ALgorithm 1
            self.clean();
            cur_x_temp = self.cur_x;
            k = length(self.cur_p);
            smooth_p = 1/k * ones(k, 1);
            smooth_pis = self.cur_pis;
            while ~self.terminator.terminate()
                [gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = self.getProxParams()
                self.mytimer.start()
                x_under = (1 - alpha_t) * self.cur_x + alpha_t * cur_x_temp;
                self.cur_pis = self.problem_data.projectPis(mu_pi_t, smooth_pis, x_under);
                self.cur_p = self.problem_data.projectP(mu_p_t, smooth_p, x_under, self.cur_pis);
                cur_x_temp = self.problem_data.projectX(gamma_t_inv, cur_x_temp, self.cur_p, self.cur_pis);
                self.cur_x = (1 - alpha_t) * self.cur_x + alpha_t * cur_x_temp;
                self.mytimer.pause()

                %record down all info
                self.time_elapsed = self.mytimer.getTime();
                self.time_history = [self.timer_history, self.time_elasped];
                [self.cur_obj, self.cur_true_gap] = self.problem_data.evalX(self.cur_x);
                self.obj_history = [self.obj_history, self.cur_obj];
            end
        end

        function self = NestrovTypeAlgorithm(problem_data, terminator)
            self = self@Algorithm.FOAlgorithm(problem_data, terminator);
            self.num_param_choices = length(self.omega_x_estimates) * length(self.omega_p_ratios) * length(self.omega_pi_ratios);
        end

    end

    methods(Abstract)        
        [gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = getProxParams(self)
        
    end
end
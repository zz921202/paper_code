classdef NestrovTypeAlgorithm < Algorithm.FOAlgorithm
    properties
        % omega_p_ratios = [1, 0.5, 0.1, 0.01];
        % omega_pi_ratios = [1, 0.5, 0.1, 0.01];
        % omega_x_ratios = [1, 10, 0.1];
        cur_omega_p;
        cur_omega_pi;
        cur_M_pi;
        cur_omega_x;
        cur_x_temp; %xt while cur_x refers to the accumulated a.k.a x_bar
        fake_cost_history = [];
        temp_cost_history = [];
        x_choice;
        p_choice;
        pi_choice;
    end

    methods



        function setGridParam(self, index)
            num_p_choices = length(self.omega_p_ratios);
            num_pi_choices = length(self.omega_pi_ratios);
            num_ppi = num_p_choices * num_pi_choices;
            self.x_choice = floor((index-0.1) / (num_ppi)) + 1;
            ppi_index = index - num_ppi * (self.x_choice-1);
            self.p_choice = floor((ppi_index-0.1) / num_pi_choices) + 1;
            self.pi_choice = ppi_index - (self.p_choice-1) * (num_pi_choices);
            self.cur_omega_x =  sqrt(self.x_radius2_est); %self.omega_x_ratios(x_choice) *
            self.cur_omega_p =  sqrt(self.est_p_breg_dist); %self.omega_p_ratios(p_choice) *
            self.cur_omega_pi = sqrt(self.est_pi_breg_dist);%self.omega_pi_ratios(pi_choice) * 
            self.cur_M_pi = sqrt(2) * self.cur_omega_pi;% CHEATING, valid only when smooth_pi = 0
            self.showGridParam(index);
        end

        function str = showGridParam(self, index)
             num_p_choices = length(self.omega_p_ratios);
            num_pi_choices = length(self.omega_pi_ratios);
            num_ppi = num_p_choices * num_pi_choices;
            x_choice = floor((index-0.1) / (num_ppi)) + 1;
            ppi_index = index - num_ppi * (x_choice-1);
            p_choice = floor((ppi_index-0.1) / num_pi_choices) + 1;
            pi_choice = ppi_index - (p_choice-1) * (num_pi_choices);
            % fprintf('param choice x %s, p_ratio: %s, pi_ratio %s\n', num2str(self.omega_x_ratios(x_choice)),num2str(self.omega_p_ratios(p_choice)), num2str(self.omega_pi_ratios(pi_choice)));
            
            str = sprintf( '%s %s %s', num2str(self.omega_x_ratios(x_choice)),num2str(self.omega_p_ratios(p_choice)), num2str(self.omega_pi_ratios(pi_choice)));
        end

        function [x, obj_val, est_gap, true_gap, time_elapsed, num_iters] = run(self)
            % refer to page 5 of paper ALgorithm 1
            self.cleanUp();
            cur_x_temp = self.cur_x;
            smooth_p = self.smooth_p;
            smooth_pis = self.smooth_pis;
            best_val = Inf;
            while ~self.terminator.terminate()
                [gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = self.getProxParams();
                self.mytimer.start()

                x_under = (1 - alpha_t) * self.cur_x + alpha_t * cur_x_temp;
                % mu_pi = mu_pi_t
                % cur_pis = self.cur_pis{1}

                [self.cur_pis, indi_costs] = self.problem_data.projectPis(mu_pi_t, smooth_pis, x_under);
%                 smooth_p

%                 indi_costs
                [self.cur_p, secon_grad, secon_cost] = self.problem_data.projectP(mu_p_t, smooth_p, indi_costs, self.cur_pis);
                % fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!! iter %d', self.num_iters);
                [cur_x_temp, ~] = self.problem_data.projectX(gamma_t_inv, cur_x_temp, secon_grad); 
                cur_x = (1 - alpha_t) * self.cur_x + alpha_t * cur_x_temp;
                [cur_val, ~] =self.problem_data.evalX(cur_x);
                if cur_val < best_val
                    self.cur_x = cur_x;
                    best_val = cur_val;
                end
                self.mytimer.pause();
                % fake_cost = self.problem_data.evalFakeCost(cur_x_temp, mu_p_t, smooth_p, mu_pi_t, smooth_pis);
                self.recordObjTime();
                [temp_cost, ~] = self.problem_data.evalX(x_under);
                self.temp_cost_history = [self.temp_cost_history, temp_cost];
                % self.fake_cost_history = [self.fake_cost_history, fake_cost];
            end
            x = self.cur_x; obj_val = self.cur_obj; est_gap = self.cur_est_gap;  true_gap=self.cur_true_gap; time_elapsed=self.time_elapsed; num_iters=self.num_iters;
        end

        function self = NestrovTypeAlgorithm(problem_data, terminator)
            self@Algorithm.FOAlgorithm(problem_data, terminator);
            self.num_param_choices = length(self.omega_x_ratios) * length(self.omega_p_ratios) * length(self.omega_pi_ratios);
            % num_param_choices = self.num_param_choices
            % length(self.omega_p_ratios)
        end

    end

    methods(Abstract)        
        [gamma_t_inv, mu_pi_t, mu_p_t, alpha_t] = getProxParams(self)
        
    end
end
classdef PDDAlgorithm < Algorithm.FOAlgorithm
    properties
        omega_pi_ratios = [1, 1e1, 1e-1 ];
        omega_p_ratios = [1, 1e1,  1e-1];
        omega_x_ratios = [1,  1e-1, 1e-2];
        cur_omega_p;
        cur_omega_pi;
        cur_M_pi;
        cur_omega_x;
        cur_x_temp; %xt while cur_x refers to the accumulated a.k.a x_bar

        x_bar;
        cur_lower;
        cur_upper;
        p_bar;
        pis_bar;
        p_ratio;
        accum_p;

        x_choice;
        p_choice;
        pi_choice;

        x_history = [];
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
            self.cur_omega_p =  sqrt(self.conservative_p_breg_dist); %self.omega_p_ratios(p_choice) *
            self.cur_omega_pi = sqrt(self.conservative_pi_breg_dist);%self.omega_pi_ratios(pi_choice) * 
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
            str = sprintf('param choice x %s, p_ratio: %s, pi_ratio %s', num2str(self.omega_x_ratios(x_choice)),num2str(self.omega_p_ratios(p_choice)), num2str(self.omega_pi_ratios(pi_choice)));
            disp(str)
        end

        function self = PDDAlgorithm(problem_data, terminator)
            self = self@Algorithm.FOAlgorithm(problem_data, terminator);
            self.num_param_choices = length(self.omega_x_ratios) * length(self.omega_p_ratios) * length(self.omega_pi_ratios);
            
        end




        function [x, obj_val, est_gap, true_gap, time_elapsed, num_iters] = run(self)
            % refer to page 5 of paper ALgorithm 1
            self.cleanUp();
            
            self.cur_upper = Inf;
            self.cur_lower = -Inf;
            self.accum_p = zeros(self.problem_data.k, 1);
            self.x_bar = zeros(size(self.cur_x));
            self.p_bar = zeros(size(self.cur_p));
            self.pis_bar = {};
            for ind = 1: self.problem_data.k
                self.pis_bar = [self.pis_bar, zeros(size(self.cur_pis{1}))];
            end
            
            pre_x = self.cur_x;
            cur_x = pre_x;
            cur_p = self.cur_p;
            cur_pis = self.cur_pis;
            [lambda, sigma, eta, tau] = self.getProxParams();
            fprintf('param choices are pi: %s, x: %s\n', num2str(sigma), num2str(eta));
            momentum = cur_x - pre_x;
            while ~self.terminator.terminate()
                
                self.mytimer.start()
                x_under = cur_x + lambda *momentum;
                [next_pis, indi_costs] = self.problem_data.projectPis(sigma, cur_pis, x_under);
                [next_p, secon_grad, ~] = self.problem_data.projectP(tau, cur_p, indi_costs, next_pis);
                
                [next_x, total_cost] = self.problem_data.projectX(eta, cur_x, secon_grad);
                % self.cur_x = cur_x;
                
                self.setBarX(next_x, next_p, next_pis);
                momentum = next_x - cur_x;
                cur_x = next_x; cur_p = next_p; cur_pis = next_pis;
                self.mytimer.pause();
                self.recordObjTime();
            end
            x = self.cur_x; obj_val = self.cur_obj; est_gap = self.cur_est_gap;  true_gap=self.cur_true_gap; time_elapsed=self.time_elapsed; num_iters=self.num_iters;
        end
    end

    methods(Access = public)%
        function  [lambda, sigma, eta, tau] = getProxParams(self)
            lambda = 1;
            prob = self.problem_data; %TODO
            self.p_ratio = (prob.beta - prob.alpha) / prob.alpha; % ALSO Hardcoded
            A = prob.beta;
            Mpi = sqrt(2) * self.cur_omega_pi;

            sigma = self.cur_omega_x * self.Mt * (1 + sqrt(self.p_ratio)) / (Mpi) * self.omega_pi_ratios(self.pi_choice);
            eta = self.Mt * Mpi * (self.cp * self.cur_omega_p + (1 + sqrt(self.p_ratio)) *A) / self.cur_omega_x * self.omega_x_ratios(self.x_choice);
            tau = self.cp * self.Mt * Mpi * self.cur_omega_x / self.cur_omega_p * self.omega_p_ratios(self.p_choice);

            % lambda = 0;
            % sigma = 0;
            % eta = 1e2;
        end

        function setBarX(self, x, p, pis)
            % compute xbar, pbar, pibars and update the lower and upper bound and gap estimate
            self.accum_p = self.accum_p + p;
            
            % self.x_history = [self.x_history, x];
            % x_bar = mean(self.x_history, 2);
            % self.x_bar = self.x_bar + 1/self.num_iters * (x - self.x_bar);
            % x_bar = self.x_bar;
            [x_val, ~] = self.problem_data.evalX(x);
            if x_val < self.cur_obj
                self.cur_obj = x_val;
                self.cur_x = x;
            end
            
            self.cur_upper = self.cur_obj;
            

            self.p_bar = self.p_bar + 1/self.num_iters * (p - self.p_bar);
            for ind = 1: self.problem_data.k
                self.pis_bar{ind} = self.pis_bar{ind} + p(ind) / self.accum_p(ind) * (pis{ind} - self.pis_bar{ind});
            end
           
            [cur_lower, ~] = self.problem_data.solveForX(self.p_bar, self.pis_bar);

            if cur_lower > self.cur_lower
                self.cur_lower = cur_lower;
                self.cur_p = self.p_bar;
                self.cur_pis = self.pis_bar;
            end

            self.cur_est_gap = self.cur_upper - self.cur_lower;
            
        end
    end


    


     

end
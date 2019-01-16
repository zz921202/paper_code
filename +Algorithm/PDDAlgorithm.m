classdef PDDAlgorithm < Algorithm.FOAlgorithm
    properties
        omega_p_ratios = [1, 0.5, 0.1, 0.01];
        omega_pi_ratios = [1, 0.5, 0.1, 0.01];
        omega_x_estimates = [10, 100, 1000, 10000];
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

        function self = PDDAlgorithm(problem_data, terminator)
            self = self@Algorithm.FOAlgorithm(problem_data, terminator);
            self.num_param_choices = length(self.omega_x_estimates) * length(self.omega_p_ratios) * length(self.omega_pi_ratios);
            self.p_ratio = self.problem_data.get_p_ratio();
        end




        function [x, obj_val, est_gap, true_gap, time_elapsed, num_iters] = run(self)
            % refer to page 5 of paper ALgorithm 1
            self.clean();
            pre_x = self.cur_x;
            cur_x = pre_x;
            cur_p = self.cur_p;
            cur_pis = self.cur_pis;
            [lambda, sigma, eta, tau] = self.getProxParams();

            while ~self.terminator.terminate()
                self.mytimer.start()
                x_under = cur_x + lambda (cur_x - pre_x);
                [cur_pis, indi_costs] = self.problem_data.projectPis(sigma, cur_pis, x_under);
                [cur_p, ~, ~] = self.problem_data.projectP(tau, cur_p, indi_costs, cur_pis);
                x_pre = cur_x;
                [cur_x, ~] = self.problem_data.projectX(eta, cur_x, secon_grad);
                self.cur_x = (1 - alpha_t) * self.cur_x + alpha_t * cur_x_temp;
                self.mytimer.pause();
                self.setBarX(cur_x, cur_p, cur_pis);
                self.recordObjTime();
            end
            x, obj_val, est_gap, true_gap, time_elapsed, num_iters = self.x, self.cur_obj, self.cur_est_gap, self.cur_true_gap, self.time_elaspese, self.num_iters;
        end
    end

    methods(access = private)%
        function  [lambda, sigma, eta, tau] = getProxParams(self)
            lambda = 1;
            prob = Problem.LinearRatioUncertainty(None, None, "Entropy"); %TODO
            self.p_ratio = prob.beta / prob.alpha; % ALSO Hardcoded
            A = prob.beta;
            Mpi = sqrt(2) * self.cur_omega_pi;
            sigma = self.cur_omega_x * self.Mt * (1 + sqrt(self.p_ratio)) / (Mpi);
            eta = self.Mt * Mpi * (self.cp * self.cur_omega_p + (1 + self.p_ratio) *A) / self.cur_omega_x;
            tau = self.cp * self.Mt * Mpi * self.cur_omega_x / self.cur_omega_p;
        end

        function setBarX(self, x, p, pis)
            % compute xbar, pbar, pibars and update the lower and upper bound and gap estimate
            self.accum_p = self.accum_p + p;
            self.x_bar = self.x_bar + 1/self.numiters * (x - self.x_bar);
            self.p_bar = self.p_bar + 1/self.numiters * (p - self.p_bar);
            for ind = 1: self.problem_data.k
                self.pis_bar{ind} = self.pis_bar{ind} + p{ind} / self.accum{ind} * (pis{ind} _- self.pis_bar{ind});
            end
            [cur_upper, ~] = self.problem_data.evalX(self.x_bar);
            [cur_lower, ~] = self.problem_data.solveForX(self.p_bar, self.pis_bar);
            if cur_upper < self.cur_upper
                self.cur_upper = cur_upper;
                self.cur_x = self.x_bar;
            end
            if cur_lower > self.cur_lower
                self.cur_lower = cur_lower;
                self.cur_p = self.p_bar;
                self.cur_pis = self.pis_bar;
            end

            self.cur_est_gap = self.cur_upper - self.cur_lower;
        end
    end


    


     

end
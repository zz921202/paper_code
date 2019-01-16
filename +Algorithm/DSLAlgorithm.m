classdef DSLAlgorithm < FOAlgorithm 
    properties
        bundle;

        BETA = 0.5;
        THETA = 0.5;
        OEMGAPI = 10;

        p_esitmate_ratios = [1, 0.5, 0.1, 0.01];
        cur_p_estimate;

        phase_x0;
        phase_v0_lower;
        phase_v0_upper;
        cur_omega_pi2_estimate;

        phase_l;
        phase_mu_pi;
        phase_mu_p;
        lower_threshold;
        upper_threshold;
        v_lower = -Inf;
        v_upper = Inf;

        
    end

    methods
        
        function setGridParam(self, index)
            self.cur_p_estimate = sqrt(self.conservative_p_breg_dist) * self.p_estimate_ratios(index);
        end
        % disp parameter choice in words
        % just call showGridParam 
        function str = showGridParam(self, index)
            str = sprintf('p ratio is %s', num2str(self.p_estimate_ratios(index)));
            disp(str);
        end

        function [x, obj_val, est_gap, true_gap, time_elapsed, num_iters] = run(self)
            self.cleanUp()
            x0 = self.problem_data.getInitialX();
            [indi_costs, pi0] = self.problem_data.solveForX(x0);
            [p0, secon_grad, secon_cost] = self.problem_data.projectP(0, self.smooth_p, inidi_costs, pi0);
            linear_approx_grad = self.problem_data.c + second_grad;
            cst_offset = secon_costs - secon_grad' * x0;
            bundle = Algorithm.Bundle(self.problem_data);
            [x1, obj] = bundle.solve(linear_approx_grad);
            self.phase_v0_lower = obj + cst_offset;

            [obj_x0, ~] = obj.problem_data.eval(x0);
            [obj_x1, ~] = obj.problem_data.eval(x1);
            
            if obj_x0 < obj_x1
                self.phase_v0_upper = obj_x0;
                self.phase_x0 = x0;
            else
                self.phase_v0_upper = obj_x1;
                self.phase_x0 = x1;
            end

            while ~self.terminator.terminate()
                self.dslPhase()
            end
            x, obj_val, est_gap, true_gap, time_elapsed, num_iters = self.x, self.cur_obj, self.cur_est_gap, self.cur_true_gap, self.time_elaspese, self.num_iters;
        end

        function  dslPhase(self)
            self.init_phase();
            t = 1;
            self.cur_x = self.phase_x0;
            x_t = self.cur_x;            
            cur_lower;
            cur_upper;
            self.bundle = Algorithm.Bundle(self.problem_data);
            while ~self.terminator.terminate()
                self.timer.start()
                alpha_t = 2 / (t+1);
                x_tl = (1 - alpha_t) * self.cur_x + alpha_t * x_t;
                [self.cur_pis, indi_costs] = self.problem_data.projectPis(self.phase_mu_pi, self.smooth_pis, x_tl);
                [self.cur_p, secon_grad, secon_cost] = self.problem_data.projectP(self.mu_p, self.smooth_p, indi_costs, self.cur_pis);
                % [cur_x_temp, total_cost] = self.problem_data.projectX(gamma_t_inv, cur_x_temp, secon_grad);
                obj_vec = self.problem_data.c + secon_grad;%TODO change for more general f composition
                [temp_lower, ~] = Bundle.solve(obj_vec);
                cur_lower = max(temp_lower, cur_lower);
                if cur_lower > self.lower_threshold
                    self.terminate_phase(cur_lower, self.cur_omega_pi2_estimate);


                b = self.phase_l - secon_cost + secon_grad' * x_tl;
                self.bundle.addConstraint(-obj_vec, -b);
                [x_t, ~] = self.bundle.project(self.phase_x0); % Projection finding the next xt
                x_tmd = (1 - alpha_t) * self.cur_x + alpha_t * x_t;
                ft = self.problem_data.solve(x_tmd);
                if ft < cur_upper
                    self.cur_x = x_tmd;
                    cur_upper = ft;
                    self.cur_obj = cur_upper;
                end
                    
                if cur_upper < self.phase_upper_threshold 
                    self.terminate(cur_lower, self.cur_omega_pi2_estimate);
                end

                %MAYBE NOT NEEDED in practice
                [tmd_pis, indi_costs] = self.problem_data.projectPis(self.phase_mu_pi, self.smooth_pis, x_tmd);
                Q = max(Helper.computeMaxL2norm2(tmd_pis), Helper.computeMaxL2norm2(self.cur_pis));
                if Q > self.cur_omega_pi2_estimate
                    self.terminate(cur_lower, self.cur_omega_pi2_estimate *2 );
                end
                self.
                self.timer.pause()
                self.recordObjTime();
                t = t + 1;
                self.num_iters = num_iters + 1;
                self.cur_est_gap = cur_upper - cur_lower;
            end
        end

        function terminate_phase(self, cur_lower, Q)
            self.phase_l = cur_lower;
            self.cur_omega_pi2_estimate = Q;
            self.phase_x0 = self.cur_x;
        end

        function init_phase(self)
            v_lower = self.phase_l;
            [v_upper, ~] = self.problem_data.evalX(phase_x0);

            self.phase_l = self.BETA * v_lower + (1 - self.BETA) * v_upper;
            mu = (self.THETA * (v_upper - self.phase_l)) / (2 * self.cur_omega_pi2_estimate *(1 + sqrt(2)) * 
                self.cur_p_estimate * self.cp);
            self.phase_mu_p = mu * self.cur_omega_pi2_estimate * self.cp * sqrt(2) / self.cur_p_estimate;
            self.phase_mu_pi = mu;
            self.lower_threshold = phase_l - self.THETA * (phase_l - self.phase_v0_lower);
            self.upper_threshold = phase_l + self.THETA(self.phase_v0_upper - phase_l);
        end
    end
end
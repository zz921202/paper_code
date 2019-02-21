classdef ExperimentDSLAlgorithm < Algorithm.FOAlgorithm 
    properties
        bundle;

        BETA = 0.5;
        THETA = 0.5;
        OMEGAPIS = [100, 500, 1000] ;
        OMEGAPI ;

        omega_p_ratios = 3.^(0:-1:-2);
        cur_p_estimate;

        phase_x0;
%         phase_v0_lower;
%         phase_v0_upper;
        cur_omega_pi2_estimate =0;
%         cur_p_estimate;
        phase_l;
        phase_lower;
        phase_mu_pi;
        phase_mu_p;
        phase_lower_threshold;
        phase_upper_threshold;
        % v_lower = -Inf;
        % v_upper = Inf;
        xt_history = [];
        cst_ratio = 2^-5;
        
    end

    methods
        
        function setGridParam(self, index)
            num_p = length(self.omega_p_ratios);
            pi_index = floor((index-0.1)/num_p) + 1;
            self.OMEGAPI = self.OMEGAPIS(pi_index);
            p_index = index - (pi_index -1) * num_p;
            self.cur_p_estimate = sqrt(self.est_p_breg_dist) * self.omega_p_ratios(p_index);
        end
        % disp parameter choice in words
        % just call showGridParam 
        function str = showGridParam(self, index)
            % fprintf('p ratio is %s\n', num2str(self.omega_p_ratios(index)));
             num_p = length(self.omega_p_ratios);
            pi_index = floor((index-0.1)/num_p) + 1;
            self.OMEGAPI = self.OMEGAPIS(pi_index);
            p_index = index - (pi_index -1) * num_p;
            
            str = sprintf('_ %s %s', num2str(self.omega_p_ratios(p_index)), num2str(self.OMEGAPIS(pi_index)));
            
        end

        function [x, obj_val, est_gap, true_gap, time_elapsed, num_iters] = run(self)
            self.cleanUp()
            self.cur_omega_pi2_estimate = 0;
            self.cur_p_estimate = 0;
            self.cst_ratio = 2^-5;
            x0 = self.problem_data.getInitialX();
            [indi_costs, pi0] = self.problem_data.solveForPis(x0);
            [p0, secon_grad, secon_costs] = self.problem_data.projectP(0, self.smooth_p, indi_costs, pi0);
            linear_approx_grad = self.problem_data.c + secon_grad;
            cst_offset = secon_costs - secon_grad' * x0;
            self.bundle = Algorithm.Bundle(self.problem_data);
            [obj, x1] = self.bundle.solve(linear_approx_grad);
            self.phase_lower = obj + cst_offset;
            self.updateRadiusEstimate(p0, pi0);
            
            [obj_x0, ~] = self.problem_data.evalX(x0);
            [obj_x1, ~] = self.problem_data.evalX(x1);
            
            if obj_x0 < obj_x1
%                 self.phase_v0_upper = obj_x0;
                self.phase_x0 = x0;
                
            else
%                 self.phase_v0_upper = obj_x1;
                self.phase_x0 = x1;
            end

            while ~self.terminator.terminate()
                self.dslPhase()
            end
            x = self.cur_x; obj_val = self.cur_obj; est_gap = self.cur_est_gap;  true_gap=self.cur_true_gap; time_elapsed=self.time_elapsed; num_iters=self.num_iters;
        end

        function  dslPhase(self)
            self.init_phase();
            
            t = 1;
            self.cur_x = self.phase_x0;
            x_t = self.cur_x;            
            cur_lower = self.phase_lower;
            [cur_upper, ~] = self.problem_data.evalX(self.phase_x0);
            self.bundle = Algorithm.Bundle(self.problem_data);
            while ~self.terminator.terminate()
                self.mytimer.start()
                self.cur_est_gap = cur_upper - cur_lower;
                alpha_t = 2 / (t+1);
                x_tl = (1 - alpha_t) * self.cur_x + alpha_t * x_t;
%                 self.phase_mu_pi
                [pisl, indi_costs] = self.problem_data.projectPis(self.phase_mu_pi, self.smooth_pis, x_tl);
                [pl, secon_grad, secon_costs] = self.problem_data.projectP(self.phase_mu_p, self.smooth_p, indi_costs, pisl);
                % [cur_x_temp, total_cost] = self.problem_data.projectX(gamma_t_inv, cur_x_temp, secon_grad);
                obj_vec = self.problem_data.c + secon_grad;%TODO change for more general f composition
                cst_offset = secon_costs - secon_grad' * x_tl;
                [temp_lower, x_min_support] = self.bundle.solve(obj_vec);
                temp_lower = temp_lower + cst_offset;
                cur_lower = max(cur_lower, min(temp_lower, self.phase_l));
                if cur_lower > self.phase_lower_threshold
%                     self
%                     self.phase_lower_threshold
                    fprintf('iter %s : terminating in lower bound improvement %s -> %s\n', num2str(self.num_iters), num2str(self.phase_lower_threshold), num2str(cur_lower))
                    self.terminate_phase(cur_lower, self.cur_omega_pi2_estimate);
                    break;
                end
                b = self.phase_l - cst_offset;
                self.bundle.addConstraint(-obj_vec, -b);
                [x_t, ~] = self.bundle.project(self.phase_x0); % Projection finding the next xt
                self.xt_history = [self.xt_history, x_t];
                x_tmd = (1 - alpha_t) * self.cur_x + alpha_t * x_t;
                [ft , ~]= self.problem_data.evalX(x_tmd);

                if ft < cur_upper
                    self.cur_x = x_tmd;
                    cur_upper = ft;
                    self.cur_obj = cur_upper;
                end
                
                if cur_upper < self.phase_upper_threshold 
                    fprintf('iter %s : terminating in upper bound improvement %s -> %s\n', num2str(self.num_iters), num2str(self.phase_upper_threshold),num2str(cur_upper))
                    self.terminate_phase(cur_lower, self.cur_omega_pi2_estimate);
                    break;
                 
                end

                %MAYBE NOT NEEDED in practice
                [tmd_pis, indi_costs] = self.problem_data.projectPis(self.phase_mu_pi, self.smooth_pis, x_tmd);
                [p_tmd, secon_grad, secon_cost] = self.problem_data.projectP(self.phase_mu_p, self.smooth_p, indi_costs, tmd_pis);
                tmd_cost = secon_cost + x_tmd' * self.problem_data.c;
                if tmd_cost < (self.phase_upper_threshold -self.phase_l)*0.5 + self.phase_l
                    self.cst_ratio = self.cst_ratio * 2;
                    fprintf('iter: %d updating cst ratio to %s\n', self.num_iters, num2str(self.cst_ratio));
                    self.terminate_phase(cur_lower, self.cur_omega_pi2_estimate);
                    break

                end

                if self.updateRadiusEstimate(pl, pisl) || self.updateRadiusEstimate(p_tmd, tmd_pis)
                    self.terminate_phase(cur_lower, self.cur_omega_pi2_estimate);
                    break
                end
                % Q = max(Helper.computeMaxL2norm2(tmd_pis), Helper.computeMaxL2norm2(cur_pis));
                % fprintf('iter: %d current radius is %s \n', self.num_iters ,num2str(Q));
                % if Q > self.cur_omega_pi2_estimate
                %     fprintf('iter: %d gap enlarging to %s \n', self.num_iters, num2str(self.cur_omega_pi2_estimate *2) );
                %     self.terminate_phase(self.phase_lower, self.cur_omega_pi2_estimate *2 );
                %     break;
                % end
                
                self.mytimer.pause()
                self.recordObjTime();
                t = t + 1;
             
                
                
            end
        end

        function terminate_phase(self, cur_lower, Q)
            self.phase_lower = cur_lower;
%             update_lower = self.phase_lower
            % self.cur_omega_pi2_estimate = Q;
            self.phase_x0 = self.cur_x;
            self.mytimer.pause()
            self.recordObjTime();
        end

        function flag = updateRadiusEstimate(self, p, pis)
            % hard-coding 
            p_radius = norm(p - self.smooth_p, 2);
            if self.problem_data.isEntropy
                % p
                % p .* (log(p) - log(self.smooth_p))
                radius_sq = real(sum((p + 1e-16) .* (log(p + 1e-16) - log(self.smooth_p + 1e-16))));
                p_radius = sqrt(radius_sq);
            end
            pis_radius_sq = Helper.computeMaxL2norm2(pis);
            % pis_radius_sq_sum = 0;
            % for ind = 1: length(pis)
            %     pis_radius_sq_sum = pis_radius_sq_sum + norm(pis{ind})^2;
            % end
            % pi_radius_sq = pis_radius_sq_sum / (self.cp^2);
            flag = false;
            if p_radius > self.cur_p_estimate
                self.cur_p_estimate = p_radius * 1.2;
                fprintf('iter: %d p_estimate enlarging to %s \n', self.num_iters, num2str(self.cur_p_estimate));
                flag = true;
            end
            if pis_radius_sq > self.cur_omega_pi2_estimate
                self.cur_omega_pi2_estimate = pis_radius_sq * 1.2;
                fprintf('iter: %d pi_estimate enlarging to %s \n', self.num_iters, num2str(self.cur_omega_pi2_estimate) );
                flag = true;
            end

        end

        function init_phase(self)
%             v_lower = self.phase_lower;
            [v_upper, ~] = self.problem_data.evalX(self.phase_x0);
     
%             phase_l = self.BETA * v_lower + (1 - self.BETA) * v_upper
            self.phase_l = self.BETA * self.phase_lower + (1 - self.BETA) * v_upper;
            
            mu = (self.THETA * (v_upper - self.phase_l)) / (self.cst_ratio * 2 * self.cur_omega_pi2_estimate *(1 + sqrt(2) * self.cur_p_estimate * self.cp));
            if self.cur_p_estimate < 1e-14
                optimal_ratio = 0;
            else
                optimal_ratio = self.cur_omega_pi2_estimate * self.cp * sqrt(2) / self.cur_p_estimate;
            end

            self.phase_mu_p = mu * optimal_ratio;
            self.phase_mu_pi = mu;
            % self.phase_mu_p = 0;
            % self.phase_mu_pi = 0;

            self.phase_lower_threshold = self.phase_l - self.THETA * (self.phase_l - self.phase_lower);
%             self.phase_lower_threshold
            self.phase_upper_threshold = self.phase_l + self.THETA*(v_upper - self.phase_l);
        end

        function self  = ExperimentDSLAlgorithm(problem_data, terminator)
            self@Algorithm.FOAlgorithm(problem_data, terminator);
            self.num_param_choices =  1;%length(self.omega_p_ratios) * length(self.OMEGAPIS);
        end

        function str = getName(self)
            str = 'ExpDsl';
        end
    end
end
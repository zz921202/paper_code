classdef ABPAlgorithm < Algorithm.DSLAlgorithm 
    
    methods

          function  dslPhase(self)
            % fprintf('phase lower is %s\n', num2str(self.phase_lower));
            self.init_phase();
            
            t = 1;
            self.cur_x = self.phase_x0;
            x_t = self.cur_x;            
            cur_lower = self.phase_lower;
            [cur_upper, ~] = self.problem_data.evalX(self.phase_x0);
            % self.phase_x0
            % fprintf('^^^^^^^^cur_upper is %s', num2str(cur_upper));
            self.bundle = Algorithm.Bundle(self.problem_data);
            self.bundle.LIMIT = 100;
            while ~self.terminator.terminate()
                self.mytimer.start()
                alpha_t = 2 / (t+1);
                x_tl = (1 - alpha_t) * self.cur_x + alpha_t * x_t;
                % smooth_pis = self.smooth_pis
                % phase_mu_pi = self.phase_mu_pi
%                 self.phase_mu_pi
                [cur_pis, indi_costs] = self.problem_data.projectPis(self.phase_mu_pi, self.smooth_pis, x_tl);
                % phase_mu_p = self.phase_mu_p
                % smooth_p = self.smooth_p
                % cur_pis = cur_pis
                % indi_costs = indi_costs
                [~, secon_grad, secon_costs] = self.problem_data.projectP(self.phase_mu_p, self.smooth_p, indi_costs, cur_pis);
                % [cur_x_temp, total_cost] = self.problem_data.projectX(gamma_t_inv, cur_x_temp, secon_grad);
                obj_vec = self.problem_data.c + secon_grad;%TODO change for more general f composition
                cst_offset = secon_costs - secon_grad' * x_tl;
                [temp_lower, x_min_support] = self.bundle.solve(obj_vec);
                temp_lower = temp_lower + cst_offset;
                cur_lower = max(cur_lower, min(temp_lower, self.phase_l));
                if cur_lower > self.phase_lower_threshold
%                     self
%                     self.phase_lower_threshold
                    % fprintf('cur upper %s minus lower is %s  ', num2str(cur_upper), num2str(cur_upper - cur_lower))
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
                [tmd_pis, ~] = self.problem_data.projectPis(self.phase_mu_pi, self.smooth_pis, x_tmd);
                
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
             
                
                self.cur_est_gap = cur_upper - cur_lower;
            end
        end
        
        function init_phase(self)
%             v_lower = self.phase_lower;
            [v_upper, ~] = self.problem_data.evalX(self.phase_x0);
            % [v_upper1, ~] = self.problem_data.evalX(self.phase_x0)
            % [v_upper2, ~] = self.problem_data.evalX(self.phase_x0)
            % [v_upper3, ~] = self.problem_data.evalX(self.phase_x0)

            % self.phase_x0
            % fprintf('upper is %s', num2str(v_upper));
     
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
            % fprintf('!!!!!!!!! phase lower is %s<- %s->%s phase l \n', num2str(self.phase_lower), num2str(self.phase_lower_threshold), num2str(self.phase_l));
        end

        function self  = ABPAlgorithm(problem_data, terminator)
            self@Algorithm.DSLAlgorithm(problem_data, terminator);
            self.omega_p_ratios = [1];
            self.num_param_choices =  1;
        end

        function str = getName(self)
            str = 'abl';
        end
    end
end
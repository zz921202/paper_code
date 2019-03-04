classdef PDHGAlgorithm < Algorithm.FOAlgorithm 
properties
    etas = 10 .^ (-10:10);
    taus = [1];
    cur_eta;
    cur_tau;
    phi = 1.5;
    ws_history = [];
    x_diff_history = [];
end
    methods 
        function str = getName(self);
            str = 'PDHG';
        end

        function setGridParam(self, index)
            num_eta_choices = length(self.etas);
            cur_tau_choice = floor((index - 0.1)/num_eta_choices) + 1;
            cur_eta_choice = index - (cur_tau_choice - 1) * num_eta_choices;
            self.cur_eta = self.etas(cur_eta_choice);
            self.cur_tau = 1 ./ self.cur_eta;
        end

        function str = showGridParam(self, index)
            num_eta_choices = length(self.etas);
            cur_tau_choice = floor((index - 0.1)/num_eta_choices) + 1;
            cur_eta_choice = index - (cur_tau_choice - 1) * num_eta_choices;
            cur_eta = self.etas(cur_eta_choice);
            cur_tau = 1 ./ self.cur_eta;
            str = sprintf('%s %s', num2str(cur_tau), num2str(cur_eta));
        end

        function self = PDHGAlgorithm(problem_data, terminator)
            self@Algorithm.FOAlgorithm(problem_data, terminator);
            self.num_param_choices = length(self.etas) * length(self.taus);
        end

        function [x, obj_val, est_gap, true_gap, time_elapsed, num_iters] = run(self)
            self.cleanUp();
            % initialization 
            % disp('helooooooooooo')
            xs = {};
            p = self.cur_p;
            ys = {};
            k = self.problem_data.k;
            sigma = zeros(k, 1);
            fake_y = zeros(self.problem_data.n2, 1);
            ws = {};
%             gap0 = self.cur_est_gap

            for ind = 1:k                
                [cur_sigma, cur_x, cur_y]  = self.problem_data.projectZ(0, self.cur_x, fake_y, ind);
                xs = [xs, cur_x]; ys = [ys, cur_y]; sigma(ind) = cur_sigma;
                ws = [ws, zeros( self.problem_data.n1, 1)];
            end

            sigma_next = sigma; xs_next = xs; ys_next = ys; ws_next = ws; p_next = p;
            %finishing initialization

            while ~self.terminator.terminate()
                self.mytimer.start();
                %prox step with respect to x
                % fprintf('iter %d est gap is %s', self.num_iters, num2str(self.cur_est_gap))
                % disp('helooooooooooo')
                x_temps = {};
                y_temps = ys;
                sigma_temps = zeros(k, 1);
                for ind = 1 : k
                    pre_ind = ind -1;
                    if ind == 1
                        pre_ind = k;%                     xs{ind}
%                     ws{ind}
%                     ws{pre_ind}

                    end
                    x_temps{ind} = xs{ind} - 1/2 * self.cur_tau *(ws{ind} - ws{pre_ind});
                    sigma_temps(ind) = sigma(ind) - self.cur_tau * p(ind);
                end
                    % [sig_ind, x_ind, y_ind] = self.problem_data.projectZ(sigma_temp, x_temp, y_temp, ind);
                    % sigma_next(ind) = sig_ind; xs_next{ind} = x_ind; ys_next{ind} = y_ind;
                [sigma_next, xs_next, ys_next] = self.problem_data.projectZs(sigma_temps, x_temps, y_temps);

                % prox step with respect to p
                % momentum update
                p_temp = p + self.cur_eta * (2 * sigma_next - sigma);
                [p_next, ~, ~] = self.problem_data.projectP(1, p_temp, zeros(k, 1), self.smooth_pis);
                for ind = 1 : k
                    next_ind = ind + 1;
                    if ind == k 
                        next_ind = 1;
                    end
                    ws_next{ind} = ws{ind} + self.cur_eta * 1/2 *...
                     (( 2 * xs_next{ind} - xs{ind}) - (2 * xs_next{next_ind}- xs{next_ind}));
                end

                % phi update
                for ind = 1:k
                    xs{ind} = (1 - self.phi) * xs{ind} + self.phi * xs_next{ind};
                    ys{ind} = (1 - self.phi) * ys{ind} + self.phi * ys_next{ind};
                    ws{ind} = (1 - self.phi) * ws{ind} + self.phi * ws_next{ind};
                end
                p = (1 - self.phi) * p + self.phi * p_next;
                sigma = (1 - self.phi) * sigma + self.phi * sigma_next;
                self.updateXSoln(xs);
                self.mytimer.pause();
                self.debugging(xs, ws);
                self.recordObjTime();
            end
            %one iteration

            x = self.cur_x; obj_val = self.cur_obj; est_gap = self.cur_est_gap;  true_gap=self.cur_true_gap; time_elapsed=self.time_elapsed; num_iters=self.num_iters;
        end

        function flag = solveEntropy(self)
            flag = false;
        end
            % just call showGridParam 
            
    end


    methods(Access = public)

        function debugging(self, xs, ws)
            x_diff = 0;
            k = self.problem_data.k;
            ws_val = zeros(k, 1);
            for ind = 1: k-1
                x_diff = x_diff + norm(xs{ind} - xs{ind + 1});
                ws_val(ind) = norm(ws{ind});
            end
            self.x_diff_history = [self.x_diff_history; x_diff];
            self.ws_history = [self.ws_history, ws_val];

        end

        function updateXSoln(self, xs)
            k = self.problem_data.k;
            x_bar = zeros(self.problem_data.n1, 1);
            for ind = 1: k
                x_bar = x_bar + 1/k * xs{ind};
            end
            [cur_cost, ~] = self.problem_data.evalX(x_bar);
            if cur_cost < self.cur_obj
                self.cur_x = x_bar;
            end
        end
    end
end
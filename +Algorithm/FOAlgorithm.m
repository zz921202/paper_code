classdef FOAlgorithm < handle

    properties
       

        mytimer = Helper.Timer();
        terminator;
        problem_data;
        
        num_param_choices;
        param_counter = 0;

        cur_x;
        cur_p;
        cur_pis;

        cur_obj = Inf;        
        obj_history = [];
        gap_history = [];

        time_elapsed = 0;
        time_history = [];

        cur_est_gap = Inf; 
        cur_true_gap = Inf;

        num_iters = 1;

        Mt;
        cp;
        est_p_breg_dist; 
        est_pi_breg_dist;
        x_radius2_est;
        smooth_pis;
        smooth_p;
        A;
        prob_ratio;
        dist_to_x_history = [];

        tuning_iter = 50;
        best_ind = 1;
    end

    methods
        function cleanUp(self)
            % clean up running information for another run with set parameter choices
            self.cur_x = self.problem_data.getInitialX();
            self.cur_p = self.problem_data.getInitialP();
            self.cur_pis = self.problem_data.getInitialPis();
            self.smooth_p = self.cur_p;
            self.smooth_pis = self.cur_pis;
            self.mytimer.reset();

            self.num_iters = 1;
            self.cur_est_gap = Inf;
            self.cur_true_gap = Inf;
            self.cur_obj = Inf;
            
            self.obj_history = [];
            self.gap_history = [];
            self.time_history = [];
            self.time_elapsed = 0;            

        end

        function recordObjTime(self)
            %eval cur_x and record down gap, true_gap, obj, 
            self.time_elapsed = self.mytimer.getTime();
            self.time_history = [self.time_history; self.time_elapsed];
            [self.cur_obj, self.cur_true_gap] = self.problem_data.evalX(self.cur_x);
            self.obj_history = [self.obj_history; self.cur_obj];
%             self.dist_to_x_history = [self.dist_to_x_history, norm(self.cur_x - self.problem_data.ref_x)];
            self.gap_history = [self.gap_history, self.cur_est_gap];
            self.num_iters = self.num_iters + 1;
        end

        function self = FOAlgorithm( problem_data, terminator)
            terminator.setFOAlgorithm(self);
            self.terminator = terminator;
            self.problem_data = problem_data; % problem 
            [omega_x2, omega_pi2, Mt] = problem_data.getParamsEstimate();
            [omega_p2, ratio, A, cp] = problem_data.getProbParamEstimate();
            self.cp = cp;
            self.Mt = Mt;
            self.prob_ratio = ratio;
            self.A = A;
            self.x_radius2_est = omega_x2;
            self.est_p_breg_dist = omega_p2;
            self.est_pi_breg_dist = omega_pi2;
            % counting param choices will be initialized in the subclass
            self.smooth_pis = problem_data.getInitialPis();
            k = self.problem_data.k;
            smooth_p = 1/k * ones(k, 1);
        end

        function flag = nextGridParam(self)
            self.param_counter = self.param_counter + 1;
            flag = true;
            if self.param_counter > self.num_param_choices
                flag = false;
                self.param_counter = 1;% reset
                return 
            end
            self.setGridParam(self.param_counter);
            self.showGridParam(self.param_counter);
        end

        function p_to_pi_ratio = getOptimalSmoothingRatio(self)
            % be careful about the case when alpha = 0 and ratio = inf
            if self.est_p_breg_dist < 1e-14
                p_to_pi_ratio = 0;
            else
                p_to_pi_ratio = self.cp * sqrt(2) * self.est_pi_breg_dist / sqrt(self.est_p_breg_dist);
            end
        end

        function setTuningIter(self, tunning_iter)
            self.tuning_iter = tunning_iter;
        end

        function startTuning(self)
            
            true_terminator = self.terminator;
            fake_terminator = Algorithm.Terminator.MaxIterTerminator();
            fake_terminator.MAXITERATION = self.tuning_iter;
            self.terminator = fake_terminator;
            fake_terminator.setFOAlgorithm(self);
            ind = 0;
            best_val = Inf;
            while self.nextGridParam()
                ind  = ind + 1;
                [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = self.run();
                fprintf('tuning %s th parameter choice %s  with obj_val%s and gap %s\n', num2str(ind), self.showGridParam(ind), num2str(cur_obj_val), num2str(true_gap));
                if cur_obj_val < best_val
                    best_val = cur_obj_val;
                    best_gap = true_gap;
                    best_ind = ind;
                end
            end
            fprintf('finished tuning with best param %d, %s', best_ind, self.showGridParam(ind));
            self.best_ind = ind;
            self.setGridParam(best_ind);
            self.terminator = true_terminator;
        end

        function flag = solveEntropy(self)
            flag = true;
        end

    end

    methods(Abstract)
        str = getName(self);
        [x, obj_val, est_gap, true_gap, time_elapsed, num_iters] = run(self)
        setGridParam(self, index)
        % disp parameter choice in words
        % just call showGridParam 
        str = showGridParam(self, index)
    end
end
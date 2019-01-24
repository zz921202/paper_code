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
        conservative_p_breg_dist; 
        conservative_pi_breg_dist;
        x_radius2_est;
        smooth_pis;
        smooth_p;

        dist_to_x_history = [];
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
            self.dist_to_x_history = [self.dist_to_x_history, norm(self.cur_x - self.problem_data.ref_x)];
            self.gap_history = [self.gap_history, self.cur_est_gap];
            self.num_iters = self.num_iters + 1;
        end

        function self = FOAlgorithm( problem_data, terminator)
            terminator.setFOAlgorithm(self);
            self.terminator = terminator;
            self.problem_data = problem_data; % problem 
            [omega_x2, omega_p2, omega_pi2, cp, Mt] = problem_data.getParamsEstimate();
            self.cp = cp;
            self.Mt = Mt;
            self.x_radius2_est = omega_x2;
            self.conservative_p_breg_dist = omega_p2;
            self.conservative_pi_breg_dist = omega_pi2;
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

    end

    methods(Abstract)
        [x, obj_val, est_gap, true_gap, time_elapsed, num_iters] = run(self)
        setGridParam(self, index)
        % disp parameter choice in words
        % just call showGridParam 
        str = showGridParam(self, index)
    end
end
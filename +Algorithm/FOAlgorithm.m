classdef FOAlgoritm < handle

    properties
       

        mytimer = Helper.Timer();
        terminator;
        problem_data;
        
        num_param_choices;
        param_counter = 1;

        cur_x;
        cur_p;
        cur_pis;

        cur_obj = Inf;        
        obj_history = [];

        time_elapsed = 0;
        time_history = [];

        cur_est_gap = Inf; 
        cur_true_gap = Inf;

        num_iters = 1;

        Mt;
        cp;
        conservative_p_breg_dist; 
        conservative_pi_breg_dist;

    end

    methods
        function cleanUp(self)
            % clean up running information for another run with set parameter choices
            self.cur_x = self.problem_data.getInitialX();
            self.cur_p = self.problem_data.getInitialP();
            self.cur_pis = self.problem_data.getInitialPis();
            self.timer.reset();

            self.num_iters = 1;
            self.cur_est_gap = Inf;
            self.cur_true_gap = Inf;
            self.cur_obj = Inf;

            self.time_history = [];
            self.time_elapsed = 0;            

        end

        function self = FOAlgorithm(self, problem_data, terminator)
            terminator.setFOAlgorithm(terminator);
            self.terminator = terminator;
            self.problem_data = problem_data; % problem 
            self.cp = problem_data.getCP();
            self.Mt = problem_data.getMt();
            self.conservative_p_breg_dist = problem_data.getConservativePBregDist();
            self.conservative_pi_breg_dist = problem_data.getConservativePiBregDist();
            % counting param choices will be initialized in the subclass
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
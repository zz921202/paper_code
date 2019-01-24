classdef LinearProblem < Problem.ProblemDataInterface 
    
    properties
        data_generator ; 
        EuclideanProjectorHandle = @Helper.EuclideanProjector;
        x_projector = [];
        pi_projectors = {};
        p_projector = [];
        % opt_val = 0;
        Tks = {};
        Wks = {};
        eks = {};
        dks = {};
        A = [];
        b = [];
        % c
        % m = 0;
        % n = 0;
        % k = 0;
        Mt;
        cp;
        omega_pi2;
        omega_x2;
        BOXBOUNDX = 20;
        

    end

    methods 
        function self = LinearProblem(data_generator)
            
            self.data_generator = data_generator;
        end

        function [time_elapsed] = generateData(self, k)
            
            
            self.k= k; 
            
            [c, A, b, dks, eks, Wks, Tks, Mt, omega_x2, omega_pi2] = self.data_generator.generateData(k); % k is the number of scenarios
            
            self.Tks=Tks; self.Wks=Wks; self.dks=dks; self.eks=eks; self.A=A; self.b=b; self.c =c; 
            self.Mt = Mt; self.omega_pi2 = omega_pi2; self.Mt = Mt; self.omega_x2 = omega_x2;
            [n1, m1, n2, m2] = self.data_generator.getProblemDimension();
            self.n1 = n1; self.m1 = m1; self.n2 = n2; self.m2 = m2;
            % generate projector for each individual scenario
            pi_projectors = {};
            for scen_ind = 1 :k 
                cur_projector = self.EuclideanProjectorHandle();
                cur_projector.setConstraint(-(Wks{scen_ind})', -eks{scen_ind});
                pi_projectors = [pi_projectors, cur_projector];
            end
            self.pi_projectors = pi_projectors;
            self.x_projector = self.EuclideanProjectorHandle();
            n = length(c);
            Aplus = [A; eye(length(c)); -eye(n)];
            bplus = [b; zeros(length(c), 1); ones(n, 1) * -self.BOXBOUNDX];

            self.x_projector.setConstraint(Aplus, bplus);

            self.generatePData( k);

            myTimer = Helper.Timer();
            myTimer.start();
            self.opt_val = self.getReferenceObjective();
            myTimer.pause();
            time_elapsed = myTimer.getTime();
            disp(sprintf('Optimal Objective is %s computed in %s sceonds', num2str(self.opt_val), num2str(time_elapsed)));            
        end

        function [next_x, total_cost] = projectX(self, prox_param, prox_center, secon_grad)
            % solves the following problem 
            %  min    stepsize * <grad, x> + V(prox_center, x)
            %  x in X
            overall_grad = self.c + secon_grad;
            [next_x, total_cost] = self.x_projector.project(prox_param, prox_center, overall_grad);
            
        end

        function [next_pis, indi_costs] = projectPis(self, prox_param, prox_centers, x)
            %   returns all pis
            next_pis = {};
            indi_costs = [];
            for ind = 1: self.k
                cur_grad = self.Tks{ind} * x + self.dks{ind};
                cur_projector = self.pi_projectors(ind);
                [cur_next_pi, neg_cost] = cur_projector.project(prox_param, prox_centers{ind}, -cur_grad);
                next_pis = [next_pis, cur_next_pi];
                indi_costs = [indi_costs; -neg_cost];
            end
        end

        function init_pis = getInitialPis(self)
            init_pis = {};
            
            for ind = 1: self.k
                cur_projector = self.pi_projectors(ind);
                [~, cur_pi] = cur_projector.solve(zeros(self.m2, 1));
                init_pis = [init_pis, cur_pi];
            end
        end

        function init_x = getInitialX(self)
            [~, init_x] = self.x_projector.solve(zeros(self.n1, 1));
        end

        function [total_cost, x] = solveForX(self, p, pis)
            [~, n] = size(self.A);
            second_stage_grad = zeros(n, 1);
            additional_cost = 0;% cost associated with p * <pik, dk>
            for ind = 1: self.k
                second_stage_grad = second_stage_grad + p(ind) * (self.Tks{ind})' * pis{ind};
                additional_cost = additional_cost + (self.dks{ind})' * pis{ind} * p(ind);
            end
            overall_grad = self.c + second_stage_grad;
            [total_cost, x] = self.x_projector.solve(overall_grad);
            total_cost = total_cost + additional_cost;
        end

        function [indi_costs, pis] = solveForPis(self, x)
            indi_costs = zeros(self.k, 1);
            pis = {};
            for ind = 1: self.k
                cur_pi_projector = self.pi_projectors(ind);
                grad = self.Tks{ind} * x + self.dks{ind};
                [neg_cur_cost, cur_pi] = cur_pi_projector.solve(-grad);
                pis = [pis, cur_pi];
                indi_costs(ind) = -neg_cur_cost;
            end
        end

        function [objective_val, gap] = evalX(self, x)
            
            [indi_costs, ~] = self.solveForPis(x);%%%%%%%%%%%%%%%%%%%%%
            [second_stage_cost, p]= self.solveForP(indi_costs);%%%%%%%%%%%%%%%
            first_stage_cost = x' * self.c;
            objective_val = second_stage_cost + first_stage_cost;
            gap = objective_val - self.opt_val;
        end



        function [A, b] = getXConstraint(self)
            A = self.x_projector.model.A;
            b = self.x_projector.model.rhs;
        end

        function [omega_x2, omega_p2, omega_pi2, cp, Mt] = getParamsEstimate(self)
            omega_x2 = self.omega_x2;
            Mt = self.Mt;
            omega_pi2 = self.omega_pi2;
            cp = self.getCp();
            omega_p2 = self.getConservativePBregDist();
        end
        function cost = evalFakeCost(self, x, p_param, p_center, pi_param, pi_center)
            [pis, indi_costs] = self.projectPis(pi_param, pi_center, x);
            [~, ~, second_cost] = self.projectP(p_param, p_center, indi_costs, pis);
            cost = (self.c)' * x + second_cost;
        end
    end

    methods(Access= protected)

        function cp = getCp(self)
            if strcmp(self.p_projector.getDistanceName(), 'Entropy')
                cp = 1;
            elseif strcmp(self.p_projector.getDistanceName() ,"Euclidean")
                cp = sqrt(self.k);
            else
                error('Unknown p_projector, please compute cp beforehand');
            end
        end




        function p_dist_est = getConservativePBregDist(self)
            if strcmp(self.p_projector.getDistanceName(), 'Entropy')
                p_dist_est = log(self.k);
            elseif strcmp(self.p_projector.getDistanceName() ,"Euclidean")
                p_dist_est = 1;
            else
                error('Unknown p_projector, please compute cp beforehand');
            end
        end


    end


    
    methods(Abstract)
        % handles projection
        
            %   returns all pis
        % handles evaluation, quality of solution 
        
            % a copy will be stored as RefObjective
        

        % [indi_costs, pis] = solveForPis(self, x)
        % [total_cost, x] = solveForX(self, grad)
        % init_pis = getInitialPis(self)
        % init_x = getInitialX(self)

        init_p = getInitialP(self)
        [next_p, secon_cost, secon_grad] = projectP(self, prox_param, prox_center, grad)
        generatePData(self, k)
        [total_cost, p] = solveForP(self, individual_costs)
        optimal_val = getReferenceObjective(self)
        % generate both the P projector

    end
end
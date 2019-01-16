classdef LinearProblem < Problem.ProblemDataInterface 
    
    properties
        randomMatrixGenerator = [];
        randomVectorGenerator = [];
        EuclideanProjectorHandle = @Helper.EuclideanProjector;
        x_projector = [];
        pi_projectors = {};
        p_projector = [];
        opt_val = 0;
        Tks = {};
        Wks = {};
        eks = {};
        dks = {};
        A = [];
        b = [];
        
        m = 0;
        n = 0;
        k = 0;
        Mt;
        cp;
        

    end

    methods 
        function self = LinearProblem(randomMatrixGenerator, randomVectorGenerator)
            self.randomMatrixGenerator = randomMatrixGenerator;
            self.randomVectorGenerator = randomVectorGenerator;
        end

        function [para_str, time_elapsed] = generateData(self, k, m, n)
            param_str = sprintf("%d scenarios with %d rows and %d columns", k, m, n);
            self.k, self.m, self.n = k, m, n;
            disp("Generating Example with "+ param_str);
            [Tks, self.Mt] = self.randomMatrixGenerator.generateTechMatrices(m,n,k);
            Wks = self.randomMatrixGenerator.generateRecourseMatrices(m, n, k);
            A = self.randomMatrixGenerator.generateAMatrix(m, n);
            b = self.randomVectorGenerator.generateVector(m);
            c = self.randomVectorGenerator.generateVector(n);
            dks = self.randomVectorGenerator.generateVector(m, k);
            eks = self.randomVectorGenerator.generateVector(n, k);

            self.Tks, self.Wks, self.dks, self.eks, self.A, self.b, self.c = Tks, Wks, A, b, c, dks, eks;
            % generate projector for each individual scenario
            pi_projectors = {};
            for scen_ind = 1 :k 
                cur_projector = self.EuclideanProjectorHandle();
                cur_projector.setConstraint(-(Wks{scen_ind})', -eks{scne_ind});
                pi_projectors = [pi_projectors, cur_projector];
            end
            self.pi_projectors = pi_projectors;
            self.x_projector = self.EuclideanProjectorHandle();
            self.x_projector.setConstraint(A, b);

            self.generatePData(k);

            myTimer = Helper.timer();
            myTimer.start();
            self.opt_val = self.getReferenceObjective();
            myTimer.pasue();
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
                cur_projector = self.pi_projectors{ind};
                [cur_next_pi, neg_cost] = cur_projector.project(prox_param, prox_centers{ind}, -cur_grad);
                next_pis = [next_pis, cur_next_pi];
                indi_costs = [indi_costs; -neg_cost];
            end
        end

        function init_pis = getInitialPis(self)
            init_pis = {};
            for ind = 1: self.k
                cur_projector = self.pi_projectors{ind};
                [~, cur_pi] = cur_projector.solve{zeros(self.m, 1)};
                init_pis = [init_pis, cur_pi];
            end
        end

        function init_x = getInitialX(self)
            [~, init_x] = self.x_projector.solve(zeros(self.n, 1));
        end

        function [total_cost, x] = solveForX(self, p, pis)
            second_stage_grad = zeros(self.n, 1);
            additional_cost = 0;% cost associated with p * <pik, dk>
            for ind = 1: self.k
                second_stage_grad = second_stage_grad + p(ind) * (self.Tks{ind})' * pis{ind};
                addtional_cost = additional_cost + (self.dks{ind})' * pis{ind} * p(ind);
            end
            overall_grad = self.c + second_stage_grad;
            [x, total_cost] = self.x_projector.solve(overall_grad);
            
        end

        function [indi_costs, pis] = solveForPis(self, x)
            indi_costs = zeros(1, self.k);
            pis = {};
            for ind = 1: self.k
                cur_pi_projector = self.pi_projectors{ind};
                grad = self.Tks{ind} * x + self.dks{ind};
                [cur_cost, cur_pi] = cur_pi_projector.solve{-grad};
                pis = [pis, cur_pi];
                indi_costs(ind) = -cur_cost;
            end
        end

        function [objective_val, gap] = evalX(self, x)
            [indi_costs, pis] = self.solveForPis(x);
            [second_stage_cost, p]= self.solveForP(indi_costs);
            first_stage_cost = x' * self.c;
            objective_val = second_stage_cost + first_stage_cost;
            gap = objective_val - self.opt_val;
        end

        function Mt = getMt(self)
            self.Mt;
        end

        function [A, b] = getXConstraint(self)
            A = self.x_projector.model.A;
            b = self.x_projector.model.b;
        end

    end
    methods(Abstract)
        %% handles projection
        
            %   returns all pis
        %% handles evaluation, quality of solution 
        
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
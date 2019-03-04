classdef LinearProblem < Problem.ProblemDataInterface 
    
    properties
        data_generator ; 
        EuclideanProjectorHandle = @Helper.EuclideanProjector;
        x_projector = [];
        pi_projectors = {};
        p_projector = [];
        z_projectors = {};
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
        omega_p2;
        ambiguity_set;
        solver;
        % opt_val

    end

    methods 
        function self = LinearProblem(data_generator, ambiguity_set)
            
            self.data_generator = data_generator;
            self.ambiguity_set = ambiguity_set;
            ambiguity_set.setProblem(self);
        end

        function [param_str, time_elapsed] = generateData(self, k)                        
            self.k= k; 
            disp('starting to generate data')
            [c, A, b, dks, eks, Wks, Tks, Mt, omega_x2, omega_pi2] = self.data_generator.generateData(k); % k is the number of scenarios
            disp('finished generating data')
            self.Tks=Tks; self.Wks=Wks; self.dks=dks; self.eks=eks; self.A=A; self.b=b; self.c =c; 
            self.Mt = Mt; self.omega_pi2 = omega_pi2; self.Mt = Mt; self.omega_x2 = omega_x2;
            [n1, m1, n2, m2] = self.data_generator.getProblemDimension();
            self.n1 = n1; self.m1 = m1; self.n2 = n2; self.m2 = m2;

            % generate projector for each individual scenario

            % initialize z projectors
            self.setUpXPiProjectors(c, A, b, Wks, Tks, dks, eks, k);
            self.setUpZProjectors(c, A, b, Wks, Tks, dks, eks, k);            
            % finish z projector init

            self.ambiguity_set.generateData(k);
            myTimer = Helper.Timer();
            myTimer.start();
            self.opt_val = self.getReferenceObjective();
            myTimer.pause();
            time_elapsed = myTimer.getTime();
            param_str = sprintf('k=%d, n1=%d, n2=%d, m2=%d, %s', k, n1, n2, m2, self.ambiguity_set.getName());
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
            indi_costs = zeros(self.k, 1);
            
            pi_projectors = self.pi_projectors;
            Tks = self.Tks;
            dks = self.dks;
            grads = {};
            for ind = 1: self.k
                grads{ind} = self.Tks{ind} * x + self.dks{ind};
            end
%             cur_projector = pi_projectors(1);
            % tic
            % ticBytes(gcp)
            for ind = 1: self.k
                cur_projector = pi_projectors(ind);
                % tic
                [cur_next_pi, neg_cost] = cur_projector.project(prox_param, prox_centers{ind}, -grads{ind});
                % toc%%TODO
                next_pis{ind} = cur_next_pi;
                indi_costs(ind) =  -neg_cost;
            end
            % tocBytes(gcp)
            % toc
            
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
%             [~, n] = size(self.A);
            n = self.n1;
            
            second_stage_grad = zeros(n, 1);
            additional_cost = 0;% cost associated with p * <pik, dk>
            for ind = 1: self.k
%                 pis{ind}
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

            pi_projectors = self.pi_projectors;
            Tks = self.Tks;
            dks = self.dks;
            grads = {};
            for ind = 1: self.k
                grads{ind} = self.Tks{ind} * x + self.dks{ind};
            end
            for ind = 1: self.k
                cur_pi_projector = pi_projectors(ind);
                [neg_cur_cost, cur_pi] = cur_pi_projector.solve(-grads{ind});
                pis{ind} = cur_pi;
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

        function [omega_x2, omega_pi2, Mt] = getParamsEstimate(self)
            omega_x2 = self.omega_x2;
            Mt = self.Mt;
            omega_pi2 = self.omega_pi2;
        end

        function cost = evalFakeCost(self, x, p_param, p_center, pi_param, pi_center)
            [pis, indi_costs] = self.projectPis(pi_param, pi_center, x);
            [~, ~, second_cost] = self.projectP(p_param, p_center, indi_costs, pis);
            cost = (self.c)' * x + second_cost;
        end

        function [next_sigma, next_x, next_y] = projectZ(self, sigma, x, y, ind)
%             sigma
%             x
%             y
            z_center = [sigma; x; y];
            
            % self.z_projectors(1)
            
            cur_projector = self.z_projectors(ind);
            total_len = 1 + self.n1 + self.n2;
            
            [next_z, ~]= cur_projector.project(1, z_center, zeros(total_len, 1));
            next_sigma = next_z(1);
            next_x = next_z(2: (1+ self.n1));
            next_y = next_z((2+self.n1):total_len);
        end

        function [next_sigmas, next_xs, next_ys] = projectZs(self, sigmas, xs, ys)
            % note that sigmas are stored as a double vector
            z_centers = {};
            for idx = 1: self.k
                z_centers{idx} = [sigmas(idx); xs{idx}; ys{idx}];
            end
            next_sigmas = zeros(self.k, 1);
            next_xs ={}; next_ys = {};
            total_len = 1 + self.n1 + self.n2;
            z_projectors = self.z_projectors;

            for idx = 1: self.k
                cur_projector = z_projectors(idx);
                [next_z, ~]= cur_projector.project(1, z_centers{idx}, zeros(total_len, 1));
                next_sigmas(idx) = next_z(1);
                next_xs{idx} = next_z(2: (1+ self.n1));
                next_ys{idx} = next_z((2+self.n1):total_len);
            end
        end

        function table = getDistanceTable(self)
            table = zeros(self.k, self.k);
            for i = 1: self.k
                for j = 1: self.k
                    table(i, j) = norm(self.dks{i} - self.dks{j});
                end
            end
            
        end

        function optimal_val = getReferenceObjective(self)
                                    
            est_gap_terminator = Algorithm.Terminator.EstGapTerminator();
            est_gap_terminator.GAP = 1e-3;
            max_iter_terminator = Algorithm.Terminator.MaxIterTerminator();
            max_iter_terminator.MAXITERATION = 500;
            terminator = Algorithm.Terminator.CompositeTerminator({max_iter_terminator, est_gap_terminator});
            self.solver = Algorithm.ABPAlgorithm(self, terminator);
            self.solver.setGridParam(1);
            [cur_x, cur_obj_val, est_gap, true_gap, time_elapsed, num_iters] = self.solver.run();
            self.ref_x = cur_x;
            x = cur_x
            optimal_val = cur_obj_val;
            fprintf('ABP with obj_val%s in   %s sec\n', num2str(cur_obj_val), num2str(time_elapsed));
            omega_x2 = norm(self.getInitialX() - cur_x)^2 * 10;
            fprintf('omega_x2_ratio is %s', num2str(self.omega_x2 / omega_x2) )
            self.omega_x2 = omega_x2;                            
        end
               % Make a copy of a handle object.
        function new = copy(this)
            % Instantiate new object of the same class.
            new = Problem.LinearProblem(this.data_generator, this.ambiguity_set);
 
            % Copy all non-hidden properties.
            p = properties(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
        end

       function str = getClassName(self)
            str = 'GenLin'; 
       end



    end

    methods(Access = protected)

        function setUpXPiProjectors(self, c, A, b, Wks, Tks, dks, eks, k)
                
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
        end


        function setUpZProjectors(self, c, A, b, Wks, Tks, dks, eks, k)
            n1 = self.n1; n2 = self.n2; m1 = self.m1; m2 = self.m2;                
            for ind = 1:k
                cur_z_projector = self.EuclideanProjectorHandle();
                obj_cons = [1, -c', -(eks{ind})'];
                scen_cons_pos = [zeros(m2, 1), -Tks{ind}, Wks{ind}];
                scen_cons_neg = [zeros(m2, 1), Tks{ind}, -Wks{ind}];
                fir_cons = [zeros(m1, 1), A, zeros(m1, n2)];
                non_neg_x = [zeros(n1, 1), eye(n1), zeros(n1, n2)];
                non_neg_y = [zeros(n2, 1), zeros(n2, n1), eye(n2)];
                lhs_cons = [obj_cons; scen_cons_pos; scen_cons_neg; fir_cons; non_neg_x; non_neg_y];
                rhs_cons = [0; dks{ind}; -dks{ind}; b; zeros(n1+n2, 1)];
                cur_z_projector.setConstraint(lhs_cons, rhs_cons);
                self.z_projectors = [self.z_projectors; cur_z_projector];
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


    
    methods % trivial post processing of result from the P maxization stage

        function init_p = getInitialP(self)
            init_p = self.ambiguity_set.getInitialP();
        end
        function [next_p, secon_grad, secon_cost, dist_est] = projectP(self, prox_param, prox_center, indi_cost, pis)
            if prox_param == 0
                [secon_cost, next_p, dist_est] = self.solveForP(indi_cost);
            else
                [next_p, secon_cost, dist_est] = self.ambiguity_set.projectP(prox_param, prox_center, indi_cost);
            end
            

            n1 = length(self.c);
            second_stage_grad = zeros(n1, 1);
            for ind = 1: self.k
                second_stage_grad = second_stage_grad + next_p(ind) * (self.Tks{ind})' * pis{ind};
            end
            % secon_cost
            secon_grad = second_stage_grad;
        end

        
        function [total_cost, p, dist_est] = solveForP(self, individual_costs)
            [total_cost, p, dist_est] = self.ambiguity_set.solveForP(individual_costs);
        end
        
        function [omega_p2, ratio, A, cp] = getProbParamEstimate(self)
            [omega_p2, ratio, A, cp] = self.ambiguity_set.getProbParamEstimate();
        end
        function flag = isEntropy(self)
            flag = self.ambiguity_set.isEntropy();
        end
    end
end
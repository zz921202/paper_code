classdef LinearRatioUncertainty < Problem.LinearProblem
    properties
        alpha = 0.9;
        beta = 1.1;
        reference_p = [];
        box_projector = false;
        box_projector_flag = false;
        distance_handle = @(x, y) norm(x-y, 2)^2;
        solver;
    end

    methods

        function self = LinearRatioUncertainty(data_generator, projector_name)
            % projector name: Entropy,  Euclidean, BoxEntropy, BoxEuclidean
            self@Problem.LinearProblem(data_generator);
            
            if projector_name == "Entropy"
                self.p_projector = Helper.EntropyProjector();
            elseif projector_name == "Euclidean"
                self.p_projector = Helper.EuclideanProjector();
            elseif projector_name == "BoxEntropy"
                self.box_projector = Helper.BoxPProjector.BoxPEntropyProjector();
                self.box_projector_flag = true;
                self.p_projector = Helper.EntropyProjector();
                DELTA = 1e-16;
                self.distance_handle = @(x, y) sum((y + DELTA).* (log(y + DELTA) - log(x + DELTA)));
            elseif projector_name == "BoxEuclidean"
                self.box_projector = Helper.BoxPProjector.BoxPEuclideanProjector();
                self.p_projector = Helper.EuclideanProjector();
                self.box_projector_flag = true;
            else
                error('Please choose projector from " Euclidean| BoxEntropy| BoxEuclidean"');
            end

        end

        function generatePData(self, k)
            self.reference_p = ones(k, 1) ./ k;%TODO
            % rand_vec = rand(k, 1) * 10 + 1;
            % rand_vec = [1;4];
            % self.reference_p = rand_vec / sum(rand_vec);
            A = [-ones(1, k); ones(1, k); -eye(k); eye(k); eye(k)];
            b = [-1 ; 1; -self.beta * self.reference_p; self.alpha * self.reference_p; zeros(k, 1)];
            self.p_projector.setConstraint(A, b);
            if self.box_projector_flag
                self.box_projector.setConstraint(self.reference_p * self.alpha, self.reference_p * self.beta);
            end

        end

        function [total_cost, p] = solveForP(self, individual_costs)
            if self.alpha == 1
                p = self.reference_p;
                total_cost = p' * individual_costs;
                return
            end
            [total_cost, p] = self.p_projector.solve(-individual_costs);%%%%%%%%%%%%%%%%%%%%
            total_cost = - total_cost;
        end

        function optimal_val = getReferenceObjective(self)
            [n1, m1, n2, m2] = self.data_generator.getProblemDimension();

            if self.k * n2 > 10000
                est_gap_terminator = Algorithm.Terminator.EstGapTerminator();
                est_gap_terminator.GAP = 1e-1;
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
                return
            end


            if self.alpha == 1
                beta_hat = self.beta;
            else
                beta_hat = (self.beta - self.alpha) / (1 - self.alpha);
            end

            if isempty(self.b)
                cvx_begin quiet
                    variable  y(n2, self.k)
                    variable x(n1, 1)
                    variable t
                    variable f(self.k,1)
                    variable fplus(self.k, 1)
                    % minimize( self.c' * x)
                    minimize( (1 - self.alpha)*( t + beta_hat * self.reference_p' * fplus) + self.alpha * self.reference_p' * f + self.c'  * x)
                    subject to
%                         self.A * x >= self.b;
                    for i =1 :self.k
                        f(i) >= (self.eks{i})' * y(:, i)
                        fplus(i) >= f(i) - t
                        fplus(i) >= 0
                        self.Wks{i} * y(:, i) == self.Tks{i} * x + self.dks{i}
                    end
                    x >= 0
                    y >= 0
                cvx_end  
            else
                
                cvx_begin quiet
                    variable  y(n2, self.k)
                    variable x(n1, 1)
                    variable t
                    variable f(self.k,1)
                    variable fplus(self.k, 1)
                    % minimize( self.c' * x)
                    minimize( (1 - self.alpha)*( t + beta_hat * self.reference_p' * fplus) + self.alpha * self.reference_p' * f + self.c'  * x)
                    subject to
                        self.A * x >= self.b;
                    for i =1 :self.k
                        f(i) >= (self.eks{i})' * y(:, i)
                        fplus(i) >= f(i) - t
                        fplus(i) >= 0
                        self.Wks{i} * y(:, i) == self.Tks{i} * x + self.dks{i}
                    end
                    x >= 0
                    y >= 0 
                cvx_end  
            end
            x
            optimal_val  = cvx_optval;
            omega_x2 = norm(self.getInitialX() - x)^2 * 10;
            fprintf('omega_x2_ratio is %s', num2str(self.omega_x2 / omega_x2) );
            self.omega_x2 = omega_x2;
            self.ref_x = x;
        end
        

        function init_p = getInitialP(self)
            init_p = self.reference_p;
        end

        function [next_p, secon_grad, secon_cost] = projectP(self, prox_param, prox_center, indi_costs, pis)
            % prox_param = prox_param
            % prox_center = prox_center
            % indi_costs = indi_costs
%             indi_costs
            if self.alpha == 1
                next_p = self.reference_p;
                neg_secon_cost = -next_p' * indi_costs;
            else
                if self.box_projector_flag
                    [next_p, neg_secon_cost] = self.fast_PProject(prox_param, prox_center, -indi_costs);
                else
                    [next_p, neg_secon_cost] = self.p_projector.project(prox_param, prox_center, -indi_costs);
                end
            end
            secon_cost = - neg_secon_cost;
            n1 = length(self.c);
            second_stage_grad = zeros(n1, 1);
            for ind = 1: self.k
                second_stage_grad = second_stage_grad + next_p(ind) * (self.Tks{ind})' * pis{ind};
            end
            secon_grad = second_stage_grad;

        end

        
    

        function [omega_p2, ratio, A] = getProbParamEstimate(self)
            ratio = (self.beta - self.alpha) / self.alpha;
            A = min(self.beta, self.k);
            ascending_vec = (1:self.k)';
            [~, p1] = self.solveForP(ascending_vec);%TODO
            [~, p2] = self.solveForP(-ascending_vec);
            omega_p2 = max(self.distance_handle(p1, p2), self.distance_handle(p2, p1));
            self.omega_p2 = omega_p2;
        end
        
    end

    methods(Access = private)
        function [next_p, val] = fast_PProject(self, prox_param, prox_center, grad)
            try
                
                [next_p, val]  = self.box_projector.project(prox_param, prox_center, grad);
                if any(isnan(next_p))
                     % fprintf('hello\n')
                    [next_p, val] = self.p_projector.project(prox_param, prox_center, grad);

                end
            catch
                % fprintf('hello    000\n')
                [next_p, val] = self.p_projector.project(prox_param, prox_center, grad);
            end
        end
    end
end
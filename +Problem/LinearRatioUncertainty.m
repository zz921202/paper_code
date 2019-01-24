classdef LinearRatioUncertainty < Problem.LinearProblem
    properties
        alpha = 0.9;
        beta = 1.1;
        reference_p = [];
        box_projector = false;
        box_projector_flag = false;
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
            elseif projector_name == "BoxEuclidean"
                self.box_projector = Helper.BoxPProjector.BoxPEuclideanProjector();
                self.p_projector = Helper.EuclideanProjector();
                self.box_projector_flag = true;
            else
                error('Please choose projector from "Entropy| Euclidean| BoxEntropy| BoxEuclidean"');
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
            [total_cost, p] = self.p_projector.solve(-individual_costs);%%%%%%%%%%%%%%%%%%%%
            total_cost = - total_cost;
        end

        function optimal_val = getReferenceObjective(self)
            [n1, m1, n2, m2] = self.data_generator.getProblemDimension();

            if self.alpha == 1
                beta_hat = self.beta;
            else
                beta_hat = (self.beta - self.alpha) / (1 - self.alpha);
            end

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
            x
            optimal_val  = cvx_optval;

            self.ref_x = x;
        end
        

        function init_p = getInitialP(self)
            init_p = self.reference_p;
        end

        function [next_p, secon_grad, secon_cost] = projectP(self, prox_param, prox_center, indi_costs, pis)

            if self.box_projector_flag
                [next_p, neg_secon_cost] = self.fast_PProject(prox_param, prox_center, -indi_costs);
            else
                [next_p, neg_secon_cost] = self.p_projector.project(prox_param, prox_center, -indi_costs);
            end
            secon_cost = - neg_secon_cost;
            n1 = length(self.c);
            second_stage_grad = zeros(n1, 1);
            for ind = 1: self.k
                second_stage_grad = second_stage_grad + next_p(ind) * (self.Tks{ind})' * pis{ind};
            end
            secon_grad = second_stage_grad;

        end

        
    end

    methods(Access = private)
        function [next_p, val] = fast_PProject(self, prox_param, prox_center, grad)
            
            [next_p, val]  = self.box_projector.project(prox_param, prox_center, grad);
        end
    end
end
classdef LinearRatioUncertainty < Problem.LinearProblem
    properties
        alpha = 0.9;
        beta = 1.1;
        reference_p = [];
        box_projector = false;
    end

    methods

        function self = LinearRatioUncertainty(randomMatrixGenerator, randomVectorGenerator, projector_name)
            % projector name: Entropy,  Euclidean, BoxEntropy, BoxEuclidean
            self@Problem.LinearProblem(randomMatrixGenerator, randomVectorGenerator);
            
            if projector_name == "Entropy"
                self.p_projector = Helper.EntropyProjector();
            elseif projector_name == "Euclidean"
                self.p_projector = Helper.EuclideanProjector();
            elseif projector_name == "BoxEntropy"
                self.box_projector = Helper.BoxPProjector.BoxPEntropyProjector();
                self.p_projector = Helper.EntropyProjector();
            elseif projector_name == "BoxEuclidean"
                self.box_projector = Helper.BoxPProjector.BoxPEuclideanProjector();
                self.p_projector = Helper.EuclideanProjector();
            else
                error('Please choose projector from "Entropy| Euclidean| BoxEntropy| BoxEuclidean"');
            end

        end

        function generatePData(self, k)
            self.reference_p = ones(k, 1) ./ k;
            A = [-ones(1, k); ones(1, k); -eye(k); eye(k); eye(k)];
            b = [-1 ; 1; -self.beta * self.reference_p; self.alpha * self.reference_p; zeros(k, 1)];
            self.p_projector.setconstraint(A, b);
            if box_projector
                box_projector.setConstraint(self.reference_p * alpha, self.reference_p * beta);
            end

        end

        function [total_cost, p] = solveForP(self, individual_costs)
            [total_cost, p] = self.p_projector.solve(-individual_costs);
            total_cost = - total_cost;
        end

        function optimal_val = getReferenceObjective(self)

            cvx_begin 
                variable  y(self.n, self.k)
                variable x(self.n)
                variable t
                variable f(self.k)
                variable fplus(self.k)

                minimize( (1 - self.alpha)*( t + self.beta * self.reference_p' * fplus) + self.alpha * self.reference_p' * f + self.c'  * x)
                subject to
                
                for i =1 :self.k
                    f(i) >= (self.eks{i})' * y(:, i)
                    fplus(i) >= f(i) - t
                    fplus(i) >= 0
                    self.Wks{i} * y(:, i) == self.Tks{i} * x + self.dks{i}
                end
                x >= 0
                y >= 0
            cvx_end                                       
            optimal_val  = cvx_optval;
        end
        

        function init_p = getInitialP(self)
            init_p = self.reference_p;
        end

        function [next_p, secon_grad, secon_cost] = projectP(self, prox_param, prox_center, indi_costs, pis)

            if self.custom_projector
                [next_p, neg_secon_cost] = self.fast_PProject(prox_param, prox_center, -indi_costs);
            else
                [next_p, neg_secon_cost] = self.p_projector.project(prox_param, prox_center, -indi_costs);
            end
            secon_cost = - neg_secon_cost;

            second_stage_grad = zeros(self.n, 1);
            for ind = 1: self.k
                second_stage_grad = second_stage_grad + next_p * (self.Tks{ind})' * pis{ind};
            end

        end

        
    end

    methods(Access = private)
        function [next_p, val] = fast_PProject(self, prox_param, prox_center, grad)
            [next_p, val]  = self.box_projector.project(prox_param, prox_center, grad);
        end
    end
end
classdef LinearRatioUncertainty < Problem.LinearProblem
    properties
        alpha = 0.9;
        beta = 1.1;
        reference_p = [];
        box_projector = false;
    end

    methods

        function self = LinearRatioUncertainty(randomMatrixGenerator, randomVectorGenerator, projector_name)
            % projector name: entropy, custom, euclidean
            self@Problem.LinearProblem(randomMatrixGenerator, randomVectorGenerator);
            self.p_projector = Helper.EuclideanProjector();
            if projector_name == "Entropy"
                self.p_projector = Helper.EntropyProjector();
            elseif projector_name == "Euclidean"
                self.p_projector = Helper.EuclideanProjector();
            elseif projector_name == "BoxEntropy"
                self.box_projector = Helper.BoxPProjector.BoxPEntropyProjector();
                
            elseif projector_name == "BoxEuclidean"
                self.box_projector = Helper.BoxPProjector.BoxPEuclideanProjector();
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

        function next_p = projectP(self, prox_param, prox_center, grad)
            if self.custom_projector
                next_p = self.fast_PProject(stepsize, prox_center, -grad);
            else
                next_p = self.p_projector.project(prox_param, prox_center, -grad);
            end
        end

        
    end

    methods(Access = private)
        function next_p = fast_PProject(self, prox_param, prox_center, grad)
            next_p = self.box_projector.project(prox_param, prox_center, grad);
        end
    end
end
classdef X2Uncertainty < Problem.ProbabilityAmbiguitySet

    properties
        TERMINATIONERR = 1e-8;
        is_entropy;
        projection_searcher= Helper.BinarySearcher();
        lp_searcher = Helper.BinarySearcher(); 
        
        distance_handle = @(x, y) 1/2 *  norm(x-y, 2)^2;
        k;
        reference_p;
        CUTOFF = 5000;
        problem_data;
        ambiguity_name = '';
        radius;
        eu_projector;
    end

    methods
        function self = X2Uncertainty(projector_name)
            self.ambiguity_name = 'X2Euclidean';
            self.is_entropy = false;
            if ~strcmp(projector_name, 'Euclidean')
                warning('only Euclidean is implemented for X2 uncertainty, so %s will be ignored', projector_name);
            end
            self.eu_projector = Helper.BoxPProjector.VectorBoxPEuclideanProjector();

        end

        function [omega_p2, ratio, A, cp] = getProbParamEstimate(self)
            cp = 1;
            if self.is_entropy
                cp = sqrt(self.k);
            end
            [~, p1] = self.solveForP(ascending_vec);
            [~, p2] = self.solveForP(-ascending_vec);
            omega_p2 = self.distance_handle(p1 - p2);
            ratio = (p1(1) - p2(1))/ p2(1);
            ratio = min(ratio, seld.CUTOFF);
            A = min(sqrt(self.CUTOFF), p1(1) * self.k);
        end

        function self = setProblem(self, problemData)
            self.problem_data = problemData;
        end
        function self = generateData(self, k)
            self.k = k;
            self.reference_p = ones(k, 1) ./ k;
            self.eu_projector.setUpperLowerConstraint(zeros(k, 1), ones(k, 1));
            A = [-ones(1, k); ones(1, k);eye(k)];
            b = [-1 ; 1; zeros(k, 1)];
            self.eu_projector.setConstraint(A, b);
        end
        % set up the p projector
        function str = getName(self)
            str = self.ambiguity_name;
        end

        function flag = isEntropy(self)
            flag = self.is_entropy;
        end

        function init_p = getInitialP(self)
            init_p = self.reference_p;
        end

        function [total_cost, p, p_radius2] = solveForP(self, individual_costs)
            closed_form_solver_handle = @(lambda) self.computeLP(lambda, -individual_costs);
            [p, ~] = self.searchLambda(self.lp_searcher, closed_form_solver_handle);
            total_cost = p' * individual_costs;
            p_radius2 = norm(p - self.reference_p)^2;

        end
        
        function [next_p, scen_cost, p_radius2] = projectP(self, prox_param, prox_center, indi_cost)
            closed_form_solver_handle = @(lambda) self.computeEuclidean(-indi_cost, lambda, prox_center, prox_param);
            [next_p, ~] = self.searchLambda(self.projection_searcher, closed_form_solver_handle);
             p_radius2 = norm(next_p - self.reference_p)^2;
             smoothing_cost = norm(next_p - prox_center)^2 * prox_param /2;
            scen_cost = next_p' * indi_cost - smoothing_cost;
        end


    end

    methods 
        function [next_p, lambda] = searchLambda(self, binary_searcher, closed_form_solver_handle)
            % will use binary searcher to search for lambda such that <P, M> <= epsilon 
            % consider epsilon - <P, M> , which is a monotone non-decreasing function of lambda 
            p0 = closed_form_solver_handle(0);
            epsilon = self.radius / self.k^2;
            if norm(p0 - self.reference_p)^2 <= epsilon 
                next_p = p0;
                lambda = 0;
            else
                while true 
                    % binary_searcher;
                    cur_lam = binary_searcher.nextSearchPoint();
                    
                    cur_p = closed_form_solver_handle(cur_lam);
                    x2 = norm(cur_p - self.reference_p)^2;
                    cur_err = epsilon - x2;

                    % fprintf('cur_lam %s cur_err is %s\n', num2str(cur_lam), num2str(cur_err));
                    % fprintf('cur_bd %s  <-> %s\n', num2str(self.lp_searcher.cur_lb), num2str(self.lp_searcher.cur_ub));
                    % self.distance_table
            

                    if cur_err < - self.TERMINATIONERR                        
                        binary_searcher.updateLowerBound(cur_lam);
                        % fprintf('updating lowerBound to %s :=%s \n', num2str(cur_lam), num2str(binary_searcher.cur_lb))
                    elseif cur_err > self.TERMINATIONERR

                        binary_searcher.updateUpperBound(cur_lam);
                        % fprintf('updating upperBound to %s :=%s \n', num2str(cur_lam), num2str(binary_searcher.cur_ub))
                    else
                        next_p = cur_p;
                        binary_searcher.terminate(cur_lam);
                        lambda = cur_lam;
                        break;
                    end
                end
            end
        end


        function next_p = computeEuclidean(self, grad, lambda, prox_center, prox_param)
            combined_prox_param = 2 * lambda  + prox_param;
            combined_prox_center = (2 * lambda * self.reference_p + prox_param * prox_center) ./ combined_prox_param;
            [next_p, ~] = self.eu_projector.project(combined_prox_param, combined_prox_center, grad);
        end

        function next_p = computeLP(self, lambda, grad)
            [next_p, ~] = self.eu_projector.project(2 * lambda, self.reference_p, grad);
        end
    end
end
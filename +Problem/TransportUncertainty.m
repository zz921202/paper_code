classdef TransportUncertainty < Problem.ProbabilityAmbiguitySet
    properties 
        distance_table;
        epsilon; 
        TERMINATIONERR = 1e-4;
        is_entropy;
        projection_searcher= Helper.BinarySearcher();
        lp_searcher = Helper.BinarySearcher(); 
        projection_solvers;
        distance_handle = @(x, y) 1/2 * norm(x-y, 'fro')^2;
        vec_distance_handle = @(x, y) 1/2 *  norm(x-y, 2)^2;
        k;
        reference_p;
        CUTOFF = 5000;
        problem_data;
        ambiguity_name = '';
        prev_p_table = [];
    end

    methods 
        function self = TransportUncertainty(projector_name)
            self.ambiguity_name = ['Transport ', projector_name];
            if projector_name == "Entropy"
                self.is_entropy = true;
                DELTA = 1e-16;
                self.distance_handle = @(x, y) sum(sum((y + DELTA).* (log(y + DELTA) - log(x + DELTA))));
                self.vec_distance_handle = @(x, y) sum((y + DELTA).* (log(y + DELTA) - log(x + DELTA)));
            elseif projector_name == "Euclidean"
                self.p_projector = Helper.EuclideanProjector();
            else
                error('Please choose projector from " Euclidean| BoxEntropy| BoxEuclidean"');
            end

        end
        flag = isEntropy(self)

        
        
        % set up the p projector
        str = getName(self)
        init_p = getInitialP(self)
        
        function [next_p, scen_cost, Ptable] = projectP(self, prox_param, prox_center, indi_cost)

            if (sum(self.prev_p_table))' == prox_center
                prox_center = self.prev_p_table;
            else
                prox_center = diag(self.reference_p);
            end
            grad = -indi_cost;
            if self.is_entropy
                closed_form_solver_handle = @(lambda) self.computeEntropy(grad, lambda, prox_center, prox_param);
                [Ptable, lambda] = self.searchLambda(self.lp_searcher, closed_form_solver_handle);
                next_p = (sum(Ptable))';
                scen_cost = indi_cost' * next_p - prox_center * self.distance_handle(prox_center, Ptable);
            end
        end

        % will reload prox_center if it is different from reference P
        

        function [omega_p2, ratio, A, cp] = getProbParamEstimate(self)
            cp = sqrt(self.k);
            ascending_vec = (1:self.k)';
            [~, P1] = self.solveForP(ascending_vec);
            [~, P2] = self.solveForP(-ascending_vec);
            p1 = sum(P1, 2); p2 = sum(P2, 2);
            omega_p2 = max(self.vec_distance_handle(p1, p2), self.vec_distance_handle(p2, p1));
            max_ratio = max(max(p1 ./ p2), max(p2 ./ p1));
            ratio = min(max_ratio, self.CUTOFF);
            max_A = max(max(p1), max(p2));
            A = min(max_A, sqrt(self.CUTOFF));
        end

        function self = setProblem(self, problem_data)
            self.problem_data = problem_data;

        end

        function self = generateData(self, k)
            self.reference_p = ones(k, 1) ./ k;
            self.k = k;
            % self.distance_table = self.problem_data.getDistanceTable();
            self.prev_p_table = eye(self.k) ./k;
        end

        function [total_cost, p] = solveForP(self, individual_costs)
            % closed_form_solver_handle = @(lam) self.computeLP(lam, -individual_costs);
            % [Ptable, ~] = self.searchLambda(self.lp_searcher, closed_form_solver_handle);
            % p = (sum(Ptable))';
            [next_p, scen_cost, p_table] = self.projectP(1e-10, diag(self.reference_p), individual_costs);
            p = next_p;
            total_cost = next_p' * individual_costs;
        end
    end

    methods 
        function [Ptable, lambda] = searchLambda(self, binary_searcher, closed_form_solver_handle)
            % will use binary searcher to search for lambda such that <P, M> <= epsilon 
            % consider epsilon - <P, M> , which is a monotone non-decreasing function of lambda 
            P0 = closed_form_solver_handle(0);
            epsilon = self.epsilon;
            if sum(dot(P0, self.distance_table)) <= self.epsilon 
                Ptable = P0;
                lambda = 0;
            else
                while true 
                    binary_searcher;
                    cur_lam = binary_searcher.nextSearchPoint();
                    fprintf('searching cur_lam %s\n', num2str(cur_lam));
                    cur_P = closed_form_solver_handle(cur_lam);
                    transport_cost = sum(dot(cur_P, self.distance_table));
                    cur_err = epsilon - sum(dot(cur_P, self.distance_table));
                    if cur_err < - self.TERMINATIONERR
                        
                        binary_searcher.updateLowerBound(cur_lam);
                        % fprintf('updating lowerBound to %s :=', num2str(cur_lam), num2str(binary_searcher.cur_lb))
                    elseif cur_err > self.TERMINATIONERR

                        binary_searcher.updateUpperBound(cur_lam);
                        % fprintf('updating upperBound to %s :=%s', num2str(cur_lam), num2str(binary_searcher.cur_ub))
                    else
                        Ptable = cur_P;
                        binary_searcher.terminate(cur_lam);
                        lambda = cur_lam;
                        break;
                    end
                end
            end
        end

        function Ptable = computeEntropy(self, grad, lambda, prox_center, prox_param)
            DELTA = 1e-16;
            if prox_param == 0
                prox_param = 1e-10;
            end
            aug_table = - 1./ prox_param * (lambda * self.distance_table +  ones(self.k, 1) * grad') + log(prox_center + DELTA);
            subtract_max = aug_table - max(aug_table, [], 2) * ones(1, self.k);
            exp_table = exp(subtract_max);
            normalized_mat = exp_table ./ (sum(exp_table, 2) * ones(1, self.k));
            mat = diag(self.reference_p + DELTA * self.k) * normalized_mat;
            Ptable = max(mat - DELTA, zeros(self.k, self.k));
        end

        function Ptable = computeEuclidean(self, grad, lambda, prox_center, prox_param)
        end

        function Ptable = computeLP(self, lambda, grad)
            % solves min_{Pi,. = refer_p_i} <lam * distance + c, p>
            aug_table = lambda .* self.distance_table + ones(self.k, 1) * grad';
            [~, min_col_idx] = min(aug_table, [], 2);
            Ptable = zeros(self.k, self.k);
            min_idx = min_col_idx + (0: (self.k-1))' .* self.k;
            Ptable(min_idx) = self.reference_p;
            Ptable = Ptable';
        end
    end
end
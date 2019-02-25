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
        eu_projector;
        amplify_prox_param_ratio = 1;
    end

    methods 
        function self = TransportUncertainty(projector_name)
            self.ambiguity_name = ['Transport ', projector_name];
            
            if strcmp(projector_name, "Entropy")
                self.is_entropy = true;
                DELTA = 1e-16;
                self.distance_handle = @(x, y) sum(sum((y + DELTA).* (log(y + DELTA) - log(x + DELTA))));
                self.vec_distance_handle = @(x, y) sum((y + DELTA).* (log(y + DELTA) - log(x + DELTA)));
            elseif strcmp(projector_name, "Euclidean")
                self.eu_projector = Helper.BoxPProjector.VectorBoxPEuclideanProjector();
                self.is_entropy = false;
            else
                error('Please choose projector from " Euclidean| BoxEntropy| BoxEuclidean"');
            end

        end

        function flag = isEntropy(self)
            flag = self.is_entropy;
        end

        
        
        % set up the p projector
        function str = getName(self)
            str = self.ambiguity_name;
        end

        function init_p = getInitialP(self)
            init_p = self.reference_p;
        end
        
        function [next_p, scen_cost, dist_est, Ptable] = projectP(self, prox_param, prox_center, indi_cost)
            [n,m] = size(prox_center);
            if m ~= n

                if (sum(self.prev_p_table))' == prox_center
                    prox_center = self.prev_p_table;
                else
                    prox_center = diag(self.reference_p);
                end
            end
            % prox_center
            grad = -indi_cost;
            if self.is_entropy

                closed_form_solver_handle = @(lambda) self.computeEntropy(grad, lambda, prox_center, prox_param);
            else
                closed_form_solver_handle = @(lambda) self.computeEuclidean(grad, lambda, prox_center, prox_param);
            end
            [Ptable, lambda] = self.searchLambda(self.lp_searcher, closed_form_solver_handle);
            next_p = (sum(Ptable))';
            scen_cost = indi_cost' * next_p - prox_param * self.distance_handle(prox_center, Ptable) * self.amplify_prox_param_ratio;
            dist_est = self.distance_handle(prox_center, Ptable) * self.amplify_prox_param_ratio;
        end

        % will reload prox_center if it is different from reference P
        

        function [omega_p2, ratio, A, cp] = getProbParamEstimate(self)
            cp = sqrt(self.k);
            ascending_vec = (1:self.k)';
            [~, P1] = self.solveForP(ascending_vec);
            [~, P2] = self.solveForP(-ascending_vec);
            p1 = sum(P1, 2); p2 = sum(P2, 2);
            omega_p2 = max(self.distance_handle(P1, P2), self.vec_distance_handle(P2, P1)) * self.amplify_prox_param_ratio;
            max_ratio = max(max((p1 - p2) ./ p2), max((p2 -p1) ./ p1));
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
%             self.distance_table = self.problem_data.getDistanceTable();
            self.prev_p_table = ones(self.k, self.k) ./k^2;
            

            if ~self.is_entropy
                self.eu_projector.setUpperLowerConstraint(zeros(k, 1), ones(k, 1));
                self.amplify_prox_param_ratio = k;
                A = [-ones(1, k); ones(1, k);eye(k)];
                b = [-1 ; 1; zeros(k, 1)];
                self.eu_projector.setConstraint(A, b);
            end
        end

        function [total_cost, p, dist_est] = solveForP(self, individual_costs)
            % closed_form_solver_handle = @(lam) self.computeLP(lam, -individual_costs);
            % [Ptable, ~] = self.searchLambda(self.lp_searcher, closed_form_solver_handle);
            % p = (sum(Ptable))';
            is_entropy_flag = self.is_entropy;
            self.is_entropy = true;
            % individual_costs
            fake_individual_cost = individual_costs ./ max(abs(individual_costs));
            [next_p, scen_cost, dist_est, p_table] = self.projectP(1e-10, ones(self.k, self.k) ./ (self.k)^2, fake_individual_cost);
            self.is_entropy = is_entropy_flag;
            p = next_p;
            total_cost = next_p' * individual_costs;
            % fprintf('total_cost is %s, scen_cost is %s, gap is %s\n', num2str(total_cost), num2str(scen_cost), num2str(total_cost - scen_cost));
        end
    end

    methods(Access = private)
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
                    % binary_searcher;
                    cur_lam = binary_searcher.nextSearchPoint();
                    
                    cur_P = closed_form_solver_handle(cur_lam);
                    transport_cost = sum(dot(cur_P, self.distance_table));
                    cur_err = epsilon - transport_cost;

                    % fprintf('cur_lam %s cur_err is %s\n', num2str(cur_lam), num2str(cur_err));
                    % fprintf('cur_bd %s  <-> %s\n', num2str(self.lp_searcher.cur_lb), num2str(self.lp_searcher.cur_ub));
                    % self.distance_table

                    % if self.lp_searcher.cur_ub- self.lp_searcher.cur_lb < 1e-8
                    %     Ptable = cur_P;
                    %     binary_searcher.terminate(cur_lam);
                    %     lambda = cur_lam;
                    %     break
                    % end

                    if cur_err < - self.TERMINATIONERR
                        
                        binary_searcher.updateLowerBound(cur_lam);
                        % fprintf('updating lowerBound to %s :=%s \n', num2str(cur_lam), num2str(binary_searcher.cur_lb))
                    elseif cur_err > self.TERMINATIONERR

                        binary_searcher.updateUpperBound(cur_lam);
                        % fprintf('updating upperBound to %s :=%s \n', num2str(cur_lam), num2str(binary_searcher.cur_ub))
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
            Ptable = zeros(self.k, self.k);
%             lambda 
%             self.distance_table
            grad_table = self.distance_table .* lambda + ones(self.k, 1) * grad';
            prox_center_aug =  prox_center .* self.k;
            for i = 1: self.k 
                [cur_p, ~] = self.eu_projector.project(prox_param, (prox_center_aug(i, :))', (grad_table(i, :))');
                
                Ptable(i, :) = cur_p ./ self.k;
            end
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
classdef RatioAmbiguitySet < Problem.ProbabilityAmbiguitySet
    properties
        % alpha = 0.9;
        % beta = 1.1;
        reference_p = [];
        omega_p2;
        box_projector = false;
        box_projector_flag = false;
        distance_handle = @(x, y) norm(x-y, 2)^2;
%         problem_data;        
        ambiguity_name = '';
        p_projector;        
        k;
        is_entropy=false;

    end

    methods
        function self = RatioAmbiguitySet(projector_name)
            self.ambiguity_name = ['Ratio ', projector_name];
            % if projector_name == "Entropy"
            %     self.p_projector = Helper.EntropyProjector();
            %     self.is_entropy = true;
            % elseif projector_name == "Euclidean"
            %     self.p_projector = Helper.EuclideanProjector();
            if strcmp(projector_name, "Entropy")
                self.box_projector = Helper.BoxPProjector.VectorBoxPEntropyProjector();
                self.box_projector_flag = true;
                self.is_entropy = true;
                self.p_projector = Helper.EntropyProjector();
                DELTA = 1e-16;
                self.distance_handle = @(x, y) sum((y + DELTA).* (log(y + DELTA) - log(x + DELTA)));
            elseif strcmp(projector_name ,"Euclidean")
                self.box_projector = Helper.BoxPProjector.VectorBoxPEuclideanProjector();
                self.p_projector = Helper.EuclideanProjector();
                self.box_projector_flag = true;
            else
                error('Please choose projector from " Euclidean| Entropy"');
            end
        end

        function self = setProblem(self, problem_data)
            self.problem_data = problem_data;
        end


        
        function self = generateData(self, k)
            self.reference_p = ones(k, 1) ./ k;
            A = [-ones(1, k); ones(1, k); -eye(k); eye(k); eye(k)];
            b = [-1 ; 1; -self.beta * self.reference_p; self.alpha * self.reference_p; zeros(k, 1)];
            self.p_projector.setConstraint(A, b);
            if self.box_projector_flag
                self.box_projector.setUpperLowerConstraint(self.reference_p * self.alpha, self.reference_p * self.beta);
                self.box_projector.setConstraint(A, b);
            end
            self.k = k;
        end
        % set up the p projector
        function str = getName(self)
            str = self.ambiguity_name;

       end

       function str = getInfo(self)
            str = sprintf('Ratio al %s be %s', num2str(self.alpha), num2str(self.radius));
       end

        function [omega_p2, ratio, A, cp] = getProbParamEstimate(self)
            if strcmp(self.p_projector.getDistanceName(), 'Entropy')
                cp = 1;
            elseif strcmp(self.p_projector.getDistanceName() ,"Euclidean")
                cp = sqrt(self.k);
            else
                error('Unknown p_projector, please compute cp beforehand');
            end

            ratio = (self.beta - self.alpha) / max(self.alpha, 0.01);
            ratio = min(ratio, 5000);
            A = min(self.beta, self.k);
            A = min(A, sqrt(5000));
            ascending_vec = (1:self.k)';
            [~, p1] = self.solveForP(ascending_vec);
            [~, p2] = self.solveForP(-ascending_vec);
            omega_p2 = max(self.distance_handle(p1, p2), self.distance_handle(p2, p1));

            self.omega_p2 = omega_p2;
        end

        function init_p = getInitialP(self)
            init_p = self.reference_p;
        end

        function [total_cost, p, dist_est] = solveForP(self, individual_costs)
            if self.alpha == 1
                p = self.reference_p;
                total_cost = p' * individual_costs;
                dist_est = 0;
                return
            end
            % individual_costs
            % full(self.p_projector.model.A)%TODO
%             individual_costs
            [total_cost, p] = self.p_projector.solve(-individual_costs);%%%%%%%%%%%%%%%%%%%%
            total_cost = - total_cost;
            dist_est = self.distance_handle(self.reference_p, p);
        end
        
        function [next_p, secon_cost, dist_est] = projectP(self, prox_param, prox_center, indi_costs)
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
            dist_est = self.distance_handle(self.reference_p, next_p);
        end
        function flag = isEntropy(self)
            flag = self.is_entropy;
        end


    end

     methods(Access = private)
        function [next_p, val] = fast_PProject(self, prox_param, prox_center, grad)
            % try
                % prox_center
                % grad
                [next_p, val]  = self.box_projector.project(prox_param, prox_center, grad);
                % next_p
                if any(isnan(next_p))
                     fprintf('hello\n')
                    [next_p, val] = self.p_projector.project(prox_param, prox_center, grad);

                end
            % catch
            %     fprintf('hello    000\n')
            %     [next_p, val] = self.p_projector.project(prox_param, prox_center, grad);
            % end
        end
    end


end
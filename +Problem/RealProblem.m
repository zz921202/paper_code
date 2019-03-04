classdef RealProblem < Problem.LinearProblem

    methods 
        function self = RealProblem(data_generator, ambiguity_set)
            self@Problem.LinearProblem(data_generator, ambiguity_set);
            self.BOXBOUNDX = 2000;
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
        parfor ind = 1: self.k
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
    
        function str = getClassName(self)
            str = 'SSN';
        end




    end

    methods(Access = protected)



        function setUpXPiProjectors(self, c, A, b, Wks, Tks, dks, eks, k)
                
            pi_projectors = {};
            ones_vec = ones(self.m2, 1);
            for scen_ind = 1 :k 
                cur_projector = self.EuclideanProjectorHandle();
                cur_projector.setConstraint(-(Wks{scen_ind})', -eks{scen_ind});
                cur_projector.setUpperLowerBound(-ones_vec * 1e6, ones_vec * 1e6);
                pi_projectors = [pi_projectors, cur_projector];

            end
            self.pi_projectors = pi_projectors;

            self.x_projector = self.EuclideanProjectorHandle();
            n = length(c);


            self.x_projector.setConstraint(A, b);
        end


        function setUpZProjectors(self, c, A, b, Wks, Tks, dks, eks, k)
            n1 = self.n1; n2 = self.n2; m1 = self.m1; m2 = self.m2; 
            z_projectors = {};
            ub_vec = [Inf; ones(n1, 1) * Inf; ones(n2, 1) * 1e6];               
            for ind = 1:k
                cur_z_projector = self.EuclideanProjectorHandle();
                obj_cons = [1, -c', -(eks{ind})'];
                scen_cons_pos = [zeros(m2, 1), -Tks{ind}, Wks{ind}];
                scen_cons_neg = [zeros(m2, 1), Tks{ind}, -Wks{ind}];
                % fir_cons = [zeros(m1, 1), A, zeros(m1, n2)];
                % non_neg_x = [zeros(n1, 1), eye(n1), zeros(n1, n2)];
                non_neg_y = [zeros(n2, 1), zeros(n2, n1), eye(n2)];
                % copied from Linear Problem, ignoring non_neg_x
                lhs_cons = [obj_cons; scen_cons_pos; scen_cons_neg;  non_neg_y];
                rhs_cons = [0; dks{ind}; -dks{ind};  zeros(n2, 1)];
                cur_z_projector.setConstraint(lhs_cons, rhs_cons);
                cur_z_projector.setUpperLowerBound(-ub_vec, ub_vec);
                self.z_projectors = [self.z_projectors; cur_z_projector];

            end
        end




    end

end
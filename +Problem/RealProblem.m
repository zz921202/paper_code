classdef RealProblem < Problem.LinearProblem

    methods 
        function self = RealProblem(data_generator, ambiguity_set)
            self@Problem.LinearProblem(data_generator, ambiguity_set);
            self.BOXBOUNDX = 2000;
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


            self.x_projector.setConstraint(A, b);
        end


        function setUpZProjectors(self, c, A, b, Wks, Tks, dks, eks, k)
            n1 = self.n1; n2 = self.n2; m1 = self.m1; m2 = self.m2; 
            z_projectors = {};               
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
                self.z_projectors = [self.z_projectors; cur_z_projector];
            end
        end

    end

end
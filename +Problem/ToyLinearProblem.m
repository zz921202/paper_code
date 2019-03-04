classdef ToyLinearProblem < Problem.LinearProblem

    methods 
        function self = ToyLinearProblem(data_generator, ambiguity_set)
            self@Problem.LinearProblem(data_generator, ambiguity_set);
        end
        
        function str = getClassName(self)
            str = 'ToyLin';
        end


    end

    methods(Access = protected)

        function setUpXPiProjectors(self, c, A, b, Wks, Tks, dks, eks, k)
                
            pi_projectors = {};
            for scen_ind = 1 :k 
                cur_projector = Helper.BoxProjector();
                cur_ek = eks{scen_ind};
                cur_projector.setUpperLowerConstraint(-cur_ek(self.m2+1 : 2 * self.m2), cur_ek(1: self.m2));
                pi_projectors = [pi_projectors, cur_projector];
            end
            self.pi_projectors = pi_projectors;
            self.x_projector = Helper.BoxProjector();
            self.x_projector.setUpperLowerConstraint(zeros(self.n1, 1), ones(self.n1, 1) * self.BOXBOUNDX);
            
            n = length(c);
            Aplus = [A; eye(length(c)); -eye(n)];
            bplus = [b; zeros(length(c), 1); ones(n, 1) * -self.BOXBOUNDX];

            self.x_projector.setConstraint(Aplus, bplus);
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
                non_neg_x = [zeros(n1, 1), eye(n1), zeros(n1, n2)];
                non_neg_y = [zeros(n2, 1), zeros(n2, n1), eye(n2)];
                lhs_cons = [obj_cons; scen_cons_pos; scen_cons_neg; non_neg_x; non_neg_y];
                rhs_cons = [0; dks{ind}; -dks{ind};  zeros(n1+n2, 1)];
                cur_z_projector.setConstraint(lhs_cons, rhs_cons);
                self.z_projectors = [self.z_projectors; cur_z_projector];
            end
        end



    end

end
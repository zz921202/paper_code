classdef ToyLinearProblem < Problem.LinearProblem

    methods 
        function self = ToyLinearProblem(data_generator, ambiguity_set)
            self@Problem.LinearProblem(data_generator, ambiguity_set);
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
    end

end
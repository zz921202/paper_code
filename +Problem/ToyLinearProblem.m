classdef ToyLinearProblem < Problem.LinearRatioUncertainty

    methods 
        function self = ToyLinearProblem(data_generator, projector_name)
            self@Problem.LinearRatioUncertainty(data_generator, projector_name);
        end

        function [time_elapsed] = generateData(self, k)
            
            %WARNING to be used with ToyRecourseRandomDataGenerator ONLY!
            self.k= k; 
            
            [c, A, b, dks, eks, Wks, Tks, Mt, omega_x2, omega_pi2] = self.data_generator.generateData(k); % k is the number of scenarios
            
            self.Tks=Tks; self.Wks=Wks; self.dks=dks; self.eks=eks; self.A=A; self.b=b; self.c =c; 
            self.Mt = Mt; self.omega_pi2 = omega_pi2; self.Mt = Mt; self.omega_x2 = omega_x2;
            [n1, m1, n2, m2] = self.data_generator.getProblemDimension();
            self.n1 = n1; self.m1 = m1; self.n2 = n2; self.m2 = m2;
            % generate projector for each individual scenario
            pi_projectors = {};
            for scen_ind = 1 :k 
                cur_projector = Helper.BoxProjector();
                cur_ek = eks{scen_ind};
                cur_projector.setUpperLowerConstraint(-cur_ek(m2+1 : 2 * m2), cur_ek(1: m2));
                pi_projectors = [pi_projectors, cur_projector];
            end
            self.pi_projectors = pi_projectors;
            self.x_projector = Helper.BoxProjector();
            self.x_projector.setUpperLowerConstraint(zeros(n1, 1), ones(n1, 1) * self.BOXBOUNDX);
            
            n = length(c);
            Aplus = [A; eye(length(c)); -eye(n)];
            bplus = [b; zeros(length(c), 1); ones(n, 1) * -self.BOXBOUNDX];

            self.x_projector.setConstraint(Aplus, bplus);
            


            
            % initialize z projectors
            % for ind = 1:k
            %     cur_z_projector = self.EuclideanProjectorHandle();
            %     obj_cons = [1, -c', -(eks{ind})'];
            %     scen_cons_pos = [zeros(m2, 1), -Tks{ind}, Wks{ind}];
            %     scen_cons_neg = [zeros(m2, 1), Tks{ind}, -Wks{ind}];
            %     fir_cons = [zeros(m1, 1), A, zeros(m1, n2)];
            %     non_neg_x = [zeros(n1, 1), eye(n1), zeros(n1, n2)];
            %     non_neg_y = [zeros(n2, 1), zeros(n2, n1), eye(n2)];
            %     lhs_cons = [obj_cons; scen_cons_pos; scen_cons_neg; fir_cons; non_neg_x; non_neg_y];
            %     rhs_cons = [0; dks{ind}; -dks{ind}; b; zeros(n1+n2, 1)];
            %     cur_z_projector.setConstraint(lhs_cons, rhs_cons);
            %     self.z_projectors = [self.z_projectors; cur_z_projector];
            % end
            % finish z projector init
            disp('generating data')
            self.generatePData(k);
            disp('finishing generating data')
            myTimer = Helper.Timer();
            myTimer.start();
            self.opt_val = self.getReferenceObjective();
            myTimer.pause();
            time_elapsed = myTimer.getTime();
            disp(sprintf('Optimal Objective is %s computed in %s sceonds', num2str(self.opt_val), num2str(time_elapsed)));            
        end
    end

end
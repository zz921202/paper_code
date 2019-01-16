classdef ProblemDataInterface 
    properties
        RefObjective;
        m = 0;
        n = 0;
        k = 0;
        c = [];
    end

    methods(Abstract)
        %% handles projection grad 
        [next_x, total_cost] = projectX(self, prox_param, prox_center, secon_grad)
            % solves the following problem 
            %  min    stepsize * <grad, x> + prox_param * V(prox_center, x)
            %  x in X
        [next_p, secon_grad, secon_cost] = projectP(self, prox_param, prox_center, indi_costs, pis)
        [next_pis, indi_costs] = projectPis(self, prox_param, prox_centers, x)
            %   returns all pis CONVERSION to min problem carried out inside the function

        %% handles evaluation, quality of solution 
        optimal_val = getReferenceObjective(self)
            % a copy will be stored as RefObjective
        [objective_val, gap] = evalX(self, x)
        [total_cost, x] = solveForX(self, p, pis)
        [total_cost, p] = solveForP(self, individual_costs)
        [indi_costs, pis] = solveForPis(self, x)

        %%initialization 
        [para_str, time] = generateData(self, k, m, n)
        %  optimal soln will be computed, prefereably extensive LP solver as a baseline comparison by calling getReferenceObjective
        %  will print "generating data", "computing optimal Soln" and "time taken to solve the problem"

        init_x = getInitialX(self)
        init_p = getInitialP(self)
        init_pis = getInitialPis(self)
        [A, b] = getXConstraint(self)

        %TODO 
        getCP(self) % get it from p projector
        getConservativePBregDist(self)
        getConservativePiBregDist(self)
        getMt(self)
        % getPRatio(self)

    end
end
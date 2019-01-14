classdef ProblemDataInterface 
    properties(Access = protected)
        RefObjective
    end

    methods(Abstract)
        %% handles projection grad 
        next_x = projectX(self, prox_param, prox_center, p, pis)
            % solves the following problem 
            %  min    stepsize * <grad, x> + prox_param * V(prox_center, x)
            %  x in X
        next_p = projectP(self, prox_param, prox_center, x, pis)
        next_pis = projectPis(self, prox_param, prox_centers, x)
            %   returns all pis

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

        %TODO 
        getCP(self) % get it from p projector
        getConservativePBregDist(self)
        getConservativePiBregDist(self)
        getMt(self)
    end
end
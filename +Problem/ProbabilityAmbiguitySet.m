classdef ProbabilityAmbiguitySet < handle 
    methods(Abstract)
        [omega_p2, ratio, A, cp] = getProbParamEstimate(self)
        self = setProblem(self, problemData)
        self = generateData(self, k)
        % set up the p projector
        str = getName(self)
        init_p = getInitialP(self)
        [total_cost, p, p_radius2] = solveForP(self, individual_costs)
        [next_p, scen_cost, p_radius2] = projectP(self, prox_param, prox_center, indi_cost)
        flag = isEntropy(self)
    end
end
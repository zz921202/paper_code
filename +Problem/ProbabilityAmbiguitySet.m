classdef ProbabilityAmbiguitySet < handle 
    properties
        radius = 0.1;
        alpha = 0.9;
        beta = 1.1;
        problem_data;
    end
    methods(Abstract)
        [omega_p2, ratio, A, cp] = getProbParamEstimate(self)
%         self = setProblem(self, problemData)
        self = generateData(self, k)
        % set up the p projector
        str = getName(self)
        init_p = getInitialP(self)
        [total_cost, p, p_radius2] = solveForP(self, individual_costs)
        [next_p, scen_cost, p_radius2] = projectP(self, prox_param, prox_center, indi_cost)
        flag = isEntropy(self)
        str = getInfo(self)
    end

    methods
        function self = setProblem(self, problem_data)
            self.problem_data = problem_data;

        end

    end
end
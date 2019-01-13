classdef LinearProjectorInterface < handle
    
    properties
        model=[];
        m;
        n;
    end
    methods
        % solves min <grad, x>
        function [val, soln] = solve(self, grad)
            self.model.obj = grad;
            params.outputflag = 0;
            result = gurobi(self.model, params);
            val = result.objval;
            soln = result.x;
        end

        function setConstraint(self, A, b)
            % Ax >= b
            [self.m, self.n] = size(A);
            self.model.A = sparse(A);
            self.model.rhs = b;
            self.model.sense = '>';
        end
    end


    methods(Abstract)
        
        soln = project(prox_param, prox_center, grad)
        % solves min <x, grad>  + U(prox_center, x) * prox_param
        end
end
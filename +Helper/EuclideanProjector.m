classdef EuclideanProjector < Helper.LinearProjectorInterface
    properties
        modelq;
    end
    
    methods
        function [soln, obj] = project(self, prox_param, prox_center, grad)
%                 prox_center
%                 prox_param
%                 grad
% %                 full(self.model.A)
%                 self.model.rhs
                self.modelq = self.model;
                self.modelq.obj = grad - prox_param * prox_center;
                self.modelq.Q =  sparse( 1/2 * prox_param * eye(self.n) );
                params.outputflag = 0; 
                
                res = gurobi(self.modelq, params);
                soln = res.x;
                obj = res.objval + 1/2 * prox_param * (prox_center)' * prox_center;

        end

        function str =getDistanceName(self)
            str = 'Euclidean';
        end
    end

end
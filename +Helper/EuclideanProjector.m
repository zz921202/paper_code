classdef EuclideanProjector < Helper.LinearProjectorInterface
    
    methods
        function [soln, obj] = project(self, prox_param, prox_center, grad)
                
                self.model.obj = grad - prox_param * prox_center;
                self.model.Q =  sparse( 1/2 * prox_param * eye(self.n) );
                params.outputflag = 0; 
                res = gurobi(self.model, params);
                soln = res.x;
                obj = res.objval + 1/2 * prox_param * (prox_center)' * prox_center;

        end

        function str =getDistanceName(self)
            str = 'Euclidean';
        end
    end

end
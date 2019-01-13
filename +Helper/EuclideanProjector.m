classdef EuclideanProjector < LinearProjectorInterface
    
    methods
        function soln = project(self, prox_param, prox_center, grad)

                self.model.obj = grad - prox_param * prox_center;
                self.model.Q =  sparse( 1/2 * prox_param * eye(self.n) * stepsize);
                params.outputflag = 0; 
                res = gurobi(self.model, params);
                x = res.x;

        end
    end

end
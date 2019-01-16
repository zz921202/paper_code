classdef Bundle < handle 

    properties
        LIMIT = 30;
        projector = Helper.EuclideanProjector();
        A;
        b;
        A_support = [];
        b_support = [];
    end

    methods
        function self = Bundle(problem_data)
            [self.A, self.b] = problem_data.getXConstraint()
        end

        function addConstraint(self, a, b)
            m = length(self.b_support);                                          
            if length(self.b_support) < self.LIMIT
                    self.A_support = [self.A_support; a'];
                    self.b_support = [self.b_support; b];
                else
                    self.A_support(1:m-1, :) = self.A_support(2:m, :);
                    self.b_support(1:m-1) = self.b_support(2:m);
                    self.A_support(m, :) = a';
                    self.b_support(m) = b;
                end
            end
        end

        function [next_x, val] = project(self, x_0)            
            self.loadConstraint();
            n = length(x_0);
            [next_x, val] = self.projector.project(1, x_0, zeros(n, 1));
        end

        function [obj, soln] = solve(self, grad)
            self.loadConstraint();
            [obj, soln] = self.projector.solve(grad);
        end
    end

    methods(access = private)

        function [A,b] =loadConstraint(self)
            if length(self.b_support) > 0
                A = [self.A; self.A_support];
                b = [self.b; self.b_support];
            else
                A = self.A;
                b = self.b;
            end
            self.projector.setConstraint(A,b);
        end
    end
end
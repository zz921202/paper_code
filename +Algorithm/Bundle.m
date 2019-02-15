classdef Bundle < handle 

    properties
        LIMIT = 40;
        projector = Helper.EuclideanProjector();
        A;
        b;
        A_support = [];
        b_support = [];
        BOXBOUND = 20;
    end

    methods
        function self = Bundle(problem_data)
            [A, b] = problem_data.getXConstraint();
            
            [~, n] = size(A);
            % self.A = [A; eye(n); -eye(n)]; 
            % self.b = [b ; -self.BOXBOUND * ones(n, 1); -self.BOXBOUND * ones(n, 1)];
            self.A = A;
            self.b = b;
        end
        
        function addConstraint(self, a, b)
            m = length(self.b_support);
            % a = a
            % b = b
            if length(self.b_support) < self.LIMIT
                % disp('adding in support: Bundle');
                self.A_support = [self.A_support; a'];
                self.b_support = [self.b_support; b];
            else
                self.A_support(1:m-1, :) = self.A_support(2:m, :);
                self.b_support(1:m-1) = self.b_support(2:m);
                self.A_support(m, :) = a';
                self.b_support(m) = b;
            end
            % self.showConstraints();
            % self.A_support 
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

        function showConstraints(self)

            for ind = 1: length(self.b_support)
                fprintf('%s : %s >= %s \n' , num2str(ind), num2str(self.A_support(ind, :)), num2str(self.b_support(ind)));
            end
        end
    end

    methods(Access = private)

        function [A,b] =loadConstraint(self)
            if ~isempty(self.b_support)
                A = [self.A; self.A_support];
                b = [self.b; self.b_support];
            else
                A = self.A;
                b = self.b;
            end
            % A_ = full(A)
            % b = b
            % fprintf('additional constrs %s', num2str(length(b)) );
            self.projector.setConstraint(A,b);
        end
    end
end
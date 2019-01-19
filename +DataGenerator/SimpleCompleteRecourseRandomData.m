classdef SimpleCompleteRecourseRandomData < DataGenerator.DataGeneratorInterface
properties
    m;
    n;
    X_BOX = 10;
    UPPERBOUND = 20;
    LOWERBOUND = 1;
end

methods
    function self = SimpleCompleteRecourseRandomData(m,n)
        self.m = m;
        self.n = n;
        if n < 2 * m
            error('SimpleCompleteRecourseRandomData: n must be larger than 2m');
        end
    end

    function  [c, A, b, dks, eks, Wks, Tks, Mt, omega_x2, omega_pi2] = generateData(self,k) % k is the number of scenarios
        param_str = sprintf("%d scenarios with %d rows and %d columns", k, self.m, self.n);
        disp("Generating Example with "+ param_str);
        c = self.genVector(self.n, self.LOWERBOUND, self.UPPERBOUND);
        b = self.genVector(self.n, self.LOWERBOUND, self.UPPERBOUND);
        omega_x2 = self.n * self.X_BOX;
        A = self.generateAMatrix( self.m, self.n);
        eks = self.genVectors(self.n, k, self.LOWERBOUND, self.UPPERBOUND);
        dks = self.genVectors(self.m, k, self.LOWERBOUND, self.UPPERBOUND);
        [Tks, Mt] = self.generateTechMatrices( self.m, self.n, k);
        Wks = self.generateRecourseMatrices(self.m, self.n, k);

        ext_vecs = {};
        for ind = 1: k
            cur_ek = eks{ind};
            cur_ext_vec = cur_ek(1: self.m) - cur_ek(self.m + 1: 2 * self.m);
            ext_vecs = [ext_vecs, cur_ext_vec];
        end
        
        omega_pi2 = Helper.computeMaxL2norm2(ext_vecs);
    end

    function [n1, m1, n2, m2] = getProblemDimension(self)
        n1 = self.n; n2 = self.n; m1 = self.m; m2 = self.m;
    end

end

methods(Access = protected)
        function [Tks, Mt] = generateTechMatrices(self, m, n, k)
            Tks = {};
            Mt = 0;

            for ind = 1: k
                [cur_mat, cur_mt] = self.randMatrixEig(m, n, 0, self.UPPERBOUND);
                Tks = [Tks, -cur_mat];
                Mt = max(cur_mt, Mt);
            end
        end



        function Wks = generateRecourseMatrices(self, m, n, k)
            if n < 2 * m
                error('PositiveCompleteRandomMatrix : cannot generate complete recourse when n < 2 * m')
            end 
            Wks = {};
            Wk = [eye(m), - eye(m), zeros(m, n - 2* m)];
            for ind = 1: k
                Wks = [Wks, Wk];
            end
        end

        function A = generateAMatrix(self, m, n)
            [A, ~] = self.randMatrixEig(m, n, self.LOWERBOUND, self.UPPERBOUND);
        end



end
end
classdef PositiveCompleteRandomMatrix < Helper.RandomMatrixGeneratorInterface
    properties
        UPPERBOUND = 100
    end
    methods
        function [Tks, Mt] = generateTechMatrices(self, m, n, k)
            Tks = {};
            Mt = 0;

            for ind = 1: k
                [cur_mat, cur_mt] = self.randMatrixWithEig(m, n, 0, self.UPPERBOUND);
                TKs = [Tks, -cur_mat];
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
                Wks = [Wks, Wk]
            end
        end

        function A = generateAMatrix(self, m, n)
            A = self.randMatrix(m, n, 0, self.UPPERBOUND);
        end
    end
end
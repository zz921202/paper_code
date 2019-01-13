classdef PositiveRandomVector < RandomVectorGeneratorInterface
    properties
        UPPERBOUND = 100;
    end
    methods
        function v = generateVector(self, n)
            v = self.randMatrix(n, 1, 0, self.UPPERBOUND);
        end
        function vs = generateVectors(self, n, k)
            vs = {}
            for ind = 1: k
                vs = [vs, self.randMatrix(n, 1, 0, self.UPPERBOUND)];
        end
    end
end
classdef CompositeTerminator < Algorithm.Terminator.TerminatorInterface
    properties
        all_terminators;
        N;
    end

    methods
        function self = CompositeTerminator(all_terminators)
            self.all_terminators = all_terminators;
            [~, self.N] = size(all_terminators);
        end

        function setFOAlgorithm(self, fo_algorithm)
            for ind = 1 : self.N
                cur_terminator = self.all_terminators{ind};
                cur_terminator.setFOAlgorithm(fo_algorithm);
            end
        end

        function flag = terminate(self)
            flag = false;
            for ind = 1 : self.N
                cur_terminator = self.all_terminators{ind};
                if cur_terminator.terminate()
                    flag = true;
                    break;
                end

            end

        end
    end
end
classdef TerminatorInterface < handle
    properties
        fo_algorithm;
    end

    methods
        function setFOAlgorithm(self, fo_algortihm) 
            self.fo_algorithm = fo_algorithm;
        end
    end

    methods(Abstract)
        flag = terminate(self) %decides if one should terminate it 
    end
end
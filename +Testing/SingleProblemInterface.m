classdef SingleProblemInterface < handle


    methods(Abstract)
        problem_data = generateData(self) % generate a Problem Data to put inside 
        % terminator = getTerminator(self) 
        problem_str = getProblemStr(self)
    end


end
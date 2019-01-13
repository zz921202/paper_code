classdef PComponent < handle 
    properties 
     upper_p;
     lower_p;
     upper_lambda;
     lower_lambda;
     base;
    end

    methods
        function p = getPValue(self, unconstrained_p)
            if unconstrained_p < self.lower_p
                p = self.lower_p;
            elseif unconstrained_p > self.upper_p
                p = self.upper_p;
            else
                p = unconstrained_p;
            end

        end
    end
end
classdef BinarySearcher < handle 
% searches the root for some monotone 1 dimension problem 
    properties 
        start_point = 1; 
        valid_lb = 0;
        valid_ub = Inf;
        cur_lb = [];
        cur_ub = [];
    end

    methods 
        function search_point = nextSearchPoint(self)

            if isempty(self.cur_lb) && isempty(self.cur_ub)
                search_point = self.start_point;
            elseif isempty(self.cur_lb)
                search_point = (self.valid_lb + self.cur_ub) / 2;
            elseif isempty(self.cur_ub)
                search_point = self.cur_lb * 2;
            else
                search_point = (self.cur_lb + self.cur_ub) / 2;
            end
        end

        function updateLowerBound(self, lb)
            self.cur_lb = lb;
        end

        function updateUpperBound(self, ub)
            self.cur_ub = ub;
        end

        function terminate(self, point)
            self.start_point = point;
            self.cur_lb = [];
            self.cur_ub = [];
        end
    end
end
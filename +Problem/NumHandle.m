classdef NumHandle < handle
   properties
      
   end

   methods
        function self = NumHandle()
            disp('hello world')
        end

        function hi(self)
            self.Number
        end
   end
end
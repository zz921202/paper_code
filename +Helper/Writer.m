classdef Writer < handle
    properties
        target_file_path = '';
    end

    methods
        function self = Writer(file_path)
            self.target_file_path = file_path;
        end
        function write(self, str)
            best_f = fopen(self.target_file_path, 'a');
            fprintf(best_f, [str, '\n']);
            fclose(best_f);
        end
    end
end


classdef SingleRandomProblem < Testing.SingleProblemInterface
    properties
        rng_seed ;
        max_iter ;
        % num_rows ;
        % assume numCols  =  2* numRows
        num_scenarios ;
        % breg_dist;
        radius;
        alpha = 0.9;
        beta = 1.1;
        ambiguity_set_handle;
        problem_handle;
        data_gen;
        problem_str;
        data_generator;
        eu_problem = [];
        en_problem;
    end

    methods
        function [eu_problem, en_problem] = generateData(self) % generate a Problem Data to put inside 
            % mat_gen = DataGenerator.SimpleCompleteRecourseRandomData(self.num_rows, self.num_rows * 2); 
            
            if isempty(self.eu_problem)
                eu_set = self.ambiguity_set_handle('Euclidean');
                en_set = self.ambiguity_set_handle('Entropy');

                eu_problem = self.problem_handle(self.data_generator,  eu_set);
                eu_set.setProblem(eu_problem);
                % ref_problem.alpha = self.alpha;
                % ref_problem.beta = self.beta;
                eu_set.alpha = self.alpha; eu_set.beta= self.beta; eu_set.radius = self.radius;
                rng(self.rng_seed);
                eu_problem.generateData(self.num_scenarios);
                en_problem = eu_problem.copy();
                en_set.alpha = self.alpha; en_set.beta= self.beta; en_set.radius = self.radius;
                en_set.setProblem(en_problem);
                en_set.generateData(self.num_scenarios);
                en_problem.ambiguity_set = en_set;
                self.eu_problem = eu_problem;
                self.en_problem = en_problem;
                self.problem_str = [eu_problem.getInfo(),',', eu_set.getInfo()];
                eu_problem%TODO
                eu_problem.getInfo()
            else
                eu_problem = self.eu_problem; en_problem = self.en_problem;
            end
        end

        function self =  SingleRandomProblem(problem_handle, data_generator, ambiguity_set_handle)
            %data generator should be initialized in the control file
            self.ambiguity_set_handle = ambiguity_set_handle;
            self.data_generator = data_generator;
            self.problem_handle = problem_handle;
        end
        


        % function terminator = getTerminator(self) 
        %     max_iter_terminator = Algorithm.Terminator.MaxIterTerminator();
        %     est_gap_teminator = Algorithm.Terminator.EstGapTerminator();
        %     max_iter_terminator.MAXITERATION = self.max_iter;
        %     terminator = Algorithm.Terminator.CompositeTerminator({max_iter_terminator, est_gap_teminator});
        % end

        function problem_str = getProblemStr(self)
            problem_str = sprintf('rng %d, %s', self.rng_seed, self.problem_str);
        end

    end
end
classdef SingleRandomProblem < Testing.SingleProblemInterface
    properties
        rng_seed ;
        max_iter ;
        num_rows ;
        % assume numCols  =  2* numRows
        num_scenarios ;
        breg_dist;
        alpha = 0.9;
        beta = 1.1;

    end

    methods
        function ref_problem = generateData(self) % generate a Problem Data to put inside 
            % mat_gen = DataGenerator.SimpleCompleteRecourseRandomData(self.num_rows, self.num_rows * 2); 
            mat_gen = DataGenerator.ToyRecourseRandomDataGenerator(self.num_rows, self.num_rows * 2); 
            ambiguity_set = Problem.RatioAmbiguitySet(self.breg_dist);
            ambiguity_set.alpha = self.alpha;
            ambiguity_set.beta = self.beta;
            ref_problem = Problem.ToyLinearProblem(mat_gen,  ambiguity_set);
            
            % ref_problem.alpha = self.alpha;
            % ref_problem.beta = self.beta;
            rng(self.rng_seed)
            ref_problem.generateData(self.num_scenarios);


            % ref_problem
        end

        function terminator = getTerminator(self) 
            max_iter_terminator = Algorithm.Terminator.MaxIterTerminator();
            est_gap_teminator = Algorithm.Terminator.EstGapTerminator();
            max_iter_terminator.MAXITERATION = self.max_iter;
            terminator = Algorithm.Terminator.CompositeTerminator({max_iter_terminator, est_gap_teminator});
        end

        function problem_str = getProblemStr(self)
            problem_str = sprintf('rng %d, rows %d, scenarios %d, %s, %d, alpha %s, beta %s',self.rng_seed, self.num_rows, self.num_scenarios, self.breg_dist, self.max_iter, num2str(self.alpha), num2str(self.beta));
        end

    end
end
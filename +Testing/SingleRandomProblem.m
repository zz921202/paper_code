classdef SingleRandomProblem < Testing.SingleProblemInterface
    properties
        rng_seed ;
        max_iter ;
        num_rows ;
        % assume numCols  =  2* numRows
        num_scenarios ;
        breg_dist;

    end

    methods
        function ref_problem = generateData(self) % generate a Problem Data to put inside 
            mat_gen = DataGenerator.SimpleCompleteRecourseRandomData(self.num_rows, self.num_rows * 2); 
            ref_problem = Problem.LinearRatioUncertainty(mat_gen,  self.breg_dist);
            rng(self.rng_seed);
            ref_problem.alpha =.99;
            ref_problem.beta = 1.02;
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
            problem_str = sprintf('rng %d, rows %d, scenarios %d, %s, %d',self.rng_seed, self.num_rows, self.num_scenarios, self.breg_dist, self.max_iter);
        end

    end
end
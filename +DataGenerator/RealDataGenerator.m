classdef RealDataGenerator < DataGenerator.DataGeneratorInterface
% all examples come form George Lan's APL-SP code base, conversion file is in Dropbox

    properties
        m1;
        n1;
        n2;
        m2;
        A; b; c;
        T; W; e; dks;
        k;
        x2; pi2;
        Mt;
    end

    methods
        function self = RealDataGenerator(file_path)
            % file path to a mat file with mydata
            fprintf('loading data from:%s', file_path);
            load(file_path);
            self.A = full(myData.A); self.b = full(myData.b); self.c= full(myData.c); 
            self.T = full(myData.T); self.W = full(myData.W); self.e = full(myData.e); self.dks = myData.dks;
            self.k = length(self.dks); [self.m1, self.n1] = size(self.A); [self.m2, self.n2] = size(self.W);

            self.x2 = self.getRadius2Estimate(self.A, self.b); 
            % full(self.e)
            self.pi2 = self.getPi2Estimate(self.W, self.e);
            self.Mt = self.getEigEstimate(self.T);
            fprintf('finished loading data');
        end

        function [n1, m1, n2, m2] = getProblemDimension(self)
            n1 = self.n1; m1 = self.m1;
            n2 = self.n2; m2 = self.m2;
        end


        function [c, A, b, dks, eks, Wks, Tks, Mt, omega_x2, omega_pi2] = generateData(self,k)
             % k is the number of scenarios

             if k > self.k
                error('requesting more scenario %d than existing in the data set %d', k, self.k);
            end
            chosen_idx = randsample(self.k, k);
            dks = self.dks(chosen_idx);
            eks = {}; Wks = {}; Tks = {};
            for i = 1: k
                eks = [eks, self.e];
                Wks = [Wks, self.W];
                Tks = [Tks, self.T];
                dks{i} = full(dks{i});
            end

            c = self.c; A = self.A; b = self.b;
            Mt = self.Mt; omega_x2 = self.x2; omega_pi2 = self.pi2; 

        end

        function pi2 = getPi2Estimate(self, W, e)
            cur_c = ones(self.m2, 1) * 100;
            cur_c(self.m2 - self.n1 + 1: self.m2) = 0.1;
            solver = Helper.EuclideanProjector();
            solver.setConstraint(-W', -e);
            % figure
            % subplot(1,2,1)
            % spy(W')
            % subplot(1,2,2)
            % spy(e)
            [~, ext] = solver.solve(-cur_c);
            pi2 = 4 * norm(ext)^2;
        end


    end
end
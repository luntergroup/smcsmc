from __future__ import print_function

import unittest
import populationmodels
from test_generic import TestGeneric


#
# Test various basic constant population size scenarios
#


class TestConstPopSize(TestGeneric):

    def setUp(self):
        TestGeneric.setUp(self, "testdata/constpopsize")
        self.seqlen = 1e7
        self.pop = populationmodels.Population( sequence_length = self.seqlen,
                                                scrmpath=self.scrmpath,
                                                change_points = [0, 0.01, 0.25, 0.5, 1, 1.5],
                                                num_populations = 1,
                                                population_sizes = [[1], [1], [1], [1], [1], [1]])
        
        # set default inference parameters
        self.popt = None
        self.smcsmc_initial_pop_sizes = self.pop.population_sizes
        self.lag = 2
        self.em = 0
        self.np = 1000
        self.debug = True
        self.bias_heights = [400]
        self.bias_strengths = [3,1]
        self.tmax = 4
        self.seed = (1,)

        # set targets
        self.targets = []
        self.targets.append({'type':"Recomb", 'min':9.6e-9, 'max':1.04e-8, 'truth':1e-8, 'ess':20})
        for idx in range(6):
            self.targets.append({'type':"Coal", 'pop':0, 'epoch':idx, 'min':9600, 'max':10400, 'truth':10000, 'ess':[1,15,40,60,60,35][idx]})
        self.targets[1]['min'] = 6000
        self.targets[1]['max'] = 15000
        self.max_out_of_range = 0



class TestConstPopSize_FourEpochs(TestConstPopSize):

    def setUp(self, name = "testdata/constpopsize_4epochs"):
        TestGeneric.setUp(self, name)
        self.seqlen = 1e7
        self.missing_leaves = []
        self.popt = None                               # don't use epoch pattern for inference
        self.smcsmc_change_points = [0, .02, .1, .5]   # but rather use four epochs

        self.pop = populationmodels.Pop2( sequence_length = self.seqlen,
                                          change_points = [0, .02, .1, .5],
                                          population_sizes = [[1], [1], [1], [1]],
                                          scrmpath=self.scrmpath )
        
        # set default inference parameters
        self.em = 0
        self.np = 1000
        self.bias_heights = [800]
        self.bias_strengths = [1.5,1]
        self.seed = (1,)
        self.debug = True

        # set targets
        self.targets = [{'type':"Recomb", 'min':0.98e-8, 'max':1.02e-8, 'truth':1e-8, 'ess':15},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':4000,  'max':18000, 'truth':10000, 'ess':1},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':8500,  'max':11500, 'truth':10000, 'ess':5},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':9700,  'max':10300, 'truth':10000, 'ess':20},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':9800,  'max':10200, 'truth':10000, 'ess':30}]

        self.max_out_of_range = 0


class TestConstPopSize_MissingData(TestConstPopSize):

    def setUp(self):
        TestConstPopSize.setUp(self)

        # identical to TestConstPopSize, except for name, and for running inference without any data
        self.prefix = "testdata/constpopsize_missingdata"
        self.missing_leaves = [0,1]

        self.np = 100
        self.em = 0
        self.seed=(4,)
        
        # set targets
        self.targets = [{'type':"Recomb", 'truth':1e-8, 'min':0.98e-8, 'max':1.02e-8, 'ess': 45},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'truth':10000, 'min': 6000, 'max': 12000, 'ess': 12},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'truth':10000, 'min': 9800, 'max': 10200, 'ess': 45},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'truth':10000, 'min': 9800, 'max': 10200, 'ess':50},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'truth':10000, 'min': 9800, 'max': 10200, 'ess':50},
                        {'type':"Coal", 'pop':0, 'epoch':4, 'truth':10000, 'min': 9800, 'max': 10200, 'ess':50}]
        self.max_out_of_range = 0   # set to -1 to always fail (and keep intermediate files)

        


class TestConstPopSize_FourEpochs_MissingData(TestConstPopSize_FourEpochs):

    def setUp(self):
        TestConstPopSize_FourEpochs.setUp(self)

        # identical to TestConstPopSize, except for name, and for running inference without any data
        self.prefix = "testdata/constpopsize_4epochs_missingdata"
        self.missing_leaves = [0,1]

        self.np = 100
        self.em = 0
        self.bias_heights = [800]
        self.bias_strengths = [3,1]
        self.seed=(4,)
        
        # set targets
        self.targets = [{'type':"Recomb", 'truth':1e-8, 'min':0.99e-8, 'max':1.01e-8, 'ess': 30},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'truth':10000, 'min': 9000, 'max': 11000, 'ess': 6},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'truth':10000, 'min': 9400, 'max': 10600, 'ess': 20},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'truth':10000, 'min': 9800, 'max': 10200, 'ess':35},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'truth':10000, 'min': 9900, 'max': 10100, 'ess':40}]
        self.max_out_of_range = 0   # set to -1 to always fail (and keep intermediate files)


class TestConstPopSize_FourEpochs_EightSamples(TestConstPopSize_FourEpochs):

    def setUp(self, name = "testdata/constpopsize_4epochs_8samples"):
        TestConstPopSize_FourEpochs.setUp(self, name)

        # relatively short sequence, so expect biased inferences due to lack of data.  Instead, testing log likelihood
        # The estimated log likelihood is about -23917 with np=3000 and bias=2.0
        # The estimated log likelihood is about -23871 with np=6000 and bias=2.0
        # The true log likelihood may be around -23814 (np=12000 bias=2.0)
        self.pop.sequence_length = 1e6
        self.pop.num_samples =  8
        
        # set default inference parameters
        self.em = 0
        self.np = 3000
        self.bias_heights = [800]
        self.bias_strengths = [2,1]
        self.seed = (2,)
        self.debug = True

        # set targets
        self.targets = [{'type':"Recomb", 'min':0.85e-8, 'max':1.15e-8, 'truth':1e-8, 'ess':1.3},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':4000,  'max':18000, 'truth':10000, 'ess':1},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':8500,  'max':11500, 'truth':10000, 'ess':1},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':9700,  'max':10300, 'truth':10000, 'ess':1.1},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':9600,  'max':10400, 'truth':10000, 'ess':1.7},
                        #{'type':"LogL", 'min':-23920, 'max':0, 'truth': -23814}]
                        {'type':"LogL", 'min':-26560, 'max':0, 'truth': -23814}]

        self.max_out_of_range = 0




        
if __name__ == "__main__":
    unittest.main()

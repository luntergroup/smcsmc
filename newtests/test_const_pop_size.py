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
        self.em = 4
        self.np = 200
        self.debug = True
        self.bias_heights = [400]
        self.bias_strengths = [3,1]
        self.tmax = 4
        self.seed = (1,)

        # set targets
        self.targets = []
        self.targets.append({'type':"Recomb", 'min':9e-9, 'max':1.1e-8, 'truth':1e-8})
        for idx in range(5):
            self.targets.append({'type':"Coal", 'pop':0, 'epoch':idx, 'min':9000, 'max':11000, 'truth':10000})
        self.targets[4]['min'] = 8000
        self.targets[4]['max'] = 12000
        self.max_out_of_range = 0


"""
class TestConstPopSize_ThreeEpochs(TestConstPopSize):

    def setUp(self, name = "testdata/constpopsize_3epochs"):
        TestGeneric.setUp(self, name)
        self.seqlen = 1e7
        self.missing_leaves = []
        self.popt = None                            # don't use epoch pattern for inference
        self.smcsmc_change_points = [0,1,1.99]         # but rather use three large epochs

        self.pop = populationmodels.Pop2( sequence_length = self.seqlen,
                                          population_sizes = [1, 1, 1],
                                          change_points = [0,1,1.99],
                                          scrmpath=self.scrmpath )
        
        # set default inference parameters
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':2e-9, 'max':7e-8, 'truth':5e-9},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':0,     'max':100000, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':0,     'max':100000, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':0,     'max':100000, 'truth':10000},]

        self.max_out_of_range = -1
"""     


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
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':5.62e-9, 'max':5.96e-8, 'truth':1e-8},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':4287,  'max':84767, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':15060, 'max':21934, 'truth':10000},  # it's worrying that the current version does not contain truth
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':10506, 'max':11696, 'truth':10000},  # it's worrying that the current version does not contain truth
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':8268,  'max':9136,  'truth':10000}]  # it's worrying that the current version does not contain truth

        self.max_out_of_range = -1


class TestConstPopSize_MissingData(TestConstPopSize):

    def setUp(self):
        TestConstPopSize.setUp(self)

        # identical to TestConstPopSize, except for name, and for running inference without any data
        self.prefix = "testdata/constpopsize_missingdata"
        self.missing_leaves = [0,1]

        self.np = 100
        self.em = 0
        self.seed=(3,)
        
        # set targets
        self.targets = [{'type':"Recomb", 'truth':1e-8, 'min':0.96e-8, 'max':1.04e-8},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'truth':10000, 'min': 5000, 'max': 12000},  # not happy with this...
                        {'type':"Coal", 'pop':0, 'epoch':1, 'truth':10000, 'min': 9800, 'max': 10200},                        
                        {'type':"Coal", 'pop':0, 'epoch':2, 'truth':10000, 'min': 9800, 'max': 10200},                        
                        {'type':"Coal", 'pop':0, 'epoch':3, 'truth':10000, 'min': 9800, 'max': 10200},                        
                        {'type':"Coal", 'pop':0, 'epoch':4, 'truth':10000, 'min': 9800, 'max': 10200}]
        self.max_out_of_range = 0   # set to -1 to always fail (and keep intermediate files)

        


class TestConstPopSize_FourEpochs_MissingData(TestConstPopSize_FourEpochs):

    def setUp(self):
        TestConstPopSize_FourEpochs.setUp(self)

        # identical to TestConstPopSize, except for name, and for running inference without any data
        self.prefix = "testdata/constpopsize_4epochs_missingdata"
        self.missing_leaves = [0,1]

        self.np = 100
        self.em = 0
        

if __name__ == "__main__":
    unittest.main()

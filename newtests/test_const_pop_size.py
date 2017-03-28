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
        self.pop = populationmodels.Pop2( sequence_length = self.seqlen,
                                          population_sizes = [1, 1, 1, 1, 1],
                                          scrmpath=self.scrmpath )
        
        # set default inference parameters
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = []
        self.targets.append({'type':"Recomb", 'min':1e-9, 'max':1e-8, 'truth':5e-9})
        for idx in range(17):
            self.targets.append({'type':"Coal", 'pop':0, 'epoch':idx, 'min':0, 'max':100000, 'truth':10000})
        self.max_out_of_range = -1



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
        


class TestConstPopSize_FourEpochs(TestConstPopSize):

    def setUp(self, name = "testdata/constpopsize_4epochs"):
        TestGeneric.setUp(self, name)
        self.seqlen = 1e7
        self.missing_leaves = []
        self.popt = None                               # don't use epoch pattern for inference
        self.smcsmc_change_points = [0, .02, .1, .5]   # but rather use four epochs

        self.pop = populationmodels.Pop2( sequence_length = self.seqlen,
                                          change_points = [0, .02, .1, .5],
                                          population_sizes = [1, 1, 1, 1],
                                          scrmpath=self.scrmpath )
        
        # set default inference parameters
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':5.62e-9, 'max':5.96e-8, 'truth':5e-9},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':4287,  'max':84767, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':15060, 'max':21934, 'truth':10000},  # it's worrying that the current version does not contain truth
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':10506, 'max':11696, 'truth':10000},  # it's worrying that the current version does not contain truth
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':8268,  'max':9136,  'truth':10000}]  # it's worrying that the current version does not contain truth

        self.max_out_of_range = -1


# Setting the initial smcsmc parameters requires a format which is inconsistent with the original pop_sizes format for one population (i.e. [1, 1, 1, 1])
# Consider cleaning this up?
class TestConstPopSize_FourEpochs_FalseStart(TestConstPopSize_FourEpochs):

    def setUp(self, name="testdata/constpopsize_4epochs_falsestart"):
        TestConstPopSize_FourEpochs.setUp(self, name)

        self.smcsmc_initial_pop_sizes = [ [1.2], [1.2], [1.2], [1.2] ]

class TestConstPopSize_MissingData(TestConstPopSize):

    def setUp(self):
        TestConstPopSize.setUp(self)

        # identical to TestConstPopSize, except for name, and for running inference without any data
        self.prefix = "testdata/constpopsize_missingdata"
        self.missing = [0,1]



class TestConstPopSize_FourEpochs_MissingData(TestConstPopSize_FourEpochs):

    def setUp(self):
        TestConstPopSize_FourEpochs.setUp(self)

        # identical to TestConstPopSize, except for name, and for running inference without any data
        self.prefix = "testdata/constpopsize_4epochs_missingdata"
        self.missing = [0,1]



if __name__ == "__main__":
    unittest.main()

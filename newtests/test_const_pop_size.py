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
        self.popt = None                     # don't use epoch pattern for inference
        self.smcsmc_change_points = [1, 2]   # but rather use three large epochs

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
        for idx in range(3):
            self.targets.append({'type':"Coal", 'pop':0, 'epoch':idx, 'min':9000, 'max':11000, 'truth':10000})
        self.max_out_of_range = -1

class TestConstPopSize_MissingData(TestConstPopSize):

    def setUp(self):
        TestConstPopSize.setUp(self)

        # identical to TestConstPopSize, except for name, and for running inference without any data
        self.prefix = "testdata/constpopsize_missingdata"
        self.missing = [0,1]

class TestConstPopSize_ThreeEpochs_MissingData(TestConstPopSize_ThreeEpochs):

    def setUp(self):
        TestConstPopSize_ThreeEpochs.setUp(self)

        # identical to TestConstPopSize, except for name, and for running inference without any data
        self.prefix = "testdata/constpopsize_3epochs_missingdata"
        self.missing = [0,1]



if __name__ == "__main__":
    unittest.main()

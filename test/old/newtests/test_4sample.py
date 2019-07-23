from __future__ import print_function

import unittest
import os
import subprocess

import populationmodels
from test_generic import *


#
# Test if we can properly handle missing data
#

class TestMissingDataOne(TestGeneric):

    def setUp(self):
        TestGeneric.setUp(self, "testdata/4sample_missingone")
        self.seqlen = 1e7
        self.missing_leaves = [0]
        self.pop = populationmodels.Pop4( sequence_length = self.seqlen,
                                          scrmpath=self.scrmpath )
        
        # set default inference parameters
        self.em = 4
        self.np = 1000
        self.seed = (1,)

        # set targets
        self.targets = []
        self.targets.append({'type':"Recomb", 'min':1e-9, 'max':1e-8, 'truth':5e-9})
        for idx in range(17):
            self.targets.append({'type':"Coal", 'pop':0, 'epoch':idx, 'min':9000, 'max':11000, 'truth':10000})
        self.max_out_of_range = 0
        

class TestMissingDataAll(TestGeneric):

    def setUp(self):
        TestGeneric.setUp(self, "testdata/4sample_missingall")
        self.seqlen = 1e7
        self.missing_leaves = [0,1,2,3]
        self.pop = populationmodels.Pop4( sequence_length = self.seqlen,
                                          scrmpath=self.scrmpath )
        
        # set default inference parameters
        self.em = 4
        self.np = 1000
        self.seed = (1,)

        # set targets
        self.targets = []
        self.targets.append({'type':"Recomb", 'min':1e-9, 'max':1e-8, 'truth':5e-9})
        for idx in range(3):
            self.targets.append({'type':"Coal", 'pop':0, 'epoch':idx, 'min':9000, 'max':11000, 'truth':10000})
        self.max_out_of_range = 3
        


if __name__ == "__main__":
    unittest.main()

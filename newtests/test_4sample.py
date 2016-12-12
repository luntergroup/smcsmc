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

        # set inference targets
        self.target_min = [0,     1000, 1000, 2000, 4000, 4500, 7000, 6000, 7000, 6000, 5000, 5000, 6000, 7000, 9000, 9500, 11000]
        self.target_max = [1e+10, 2000, 2000, 4000, 5000, 5500, 9000, 7000, 9000, 8000, 7000, 7000, 7500, 8500, 11000, 11500, 13000]
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

        # set inference targets
        self.target_min = [9800,  9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800,  9800, 9800]
        self.target_max = [10200,  10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200,  10200, 10200]
        self.max_out_of_range = 3
        


if __name__ == "__main__":
    unittest.main()

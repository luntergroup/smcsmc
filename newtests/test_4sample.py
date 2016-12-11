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
        
    # actual test
    def test_inference(self):
        self.infer( case = 1 )

        # Ne target ranges for the varous epochs
        target_min = [0,     1000, 1000, 2000, 4000, 4500, 7000, 6000, 7000, 6000, 5000, 5000, 6000, 7000, 9000, 9500, 11000]
        target_max = [1e+10, 2000, 2000, 4000, 5000, 5500, 9000, 7000, 9000, 8000, 7000, 7000, 7500, 8500, 11000, 11500, 13000]
        results = self.readResults()

        out_of_range = 0
        for epoch, emresults in enumerate(results['Coal']):
            # get the last EM estimate
            result = emresults['estimates'][0][-1]

            # see if it is within range
            print ("Checking epoch",epoch,result,target_min[epoch],target_max[epoch])
            #self.assertTrue( result >= target_min[epoch] )
            #self.assertTrue( result <= target_max[epoch] )
            if result < target_min[epoch] or result > target_max[epoch]:
                print ("Out of range!")
                out_of_range += 1

        self.assertTrue( out_of_range < 3 )
        self.success = True



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

        
    # actual test
    def test_inference(self):
        self.infer( case = 1 )

        # Ne target ranges for the varous epochs
        target_min = [9800,  9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800, 9800,  9800, 9800]
        target_max = [10200,  10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200, 10200,  10200, 10200]
        results = self.readResults()

        out_of_range = 0
        for epoch, emresults in enumerate(results['Coal']):
            # get the last EM estimate
            result = emresults['estimates'][0][-1]

            # see if it is within range
            print ("Checking epoch",epoch,result,target_min[epoch],target_max[epoch])
            #self.assertTrue( result >= target_min[epoch] )
            #self.assertTrue( result <= target_max[epoch] )
            if result < target_min[epoch] or result > target_max[epoch]:
                print ("Out of range!")
                out_of_range += 1

        self.assertTrue( out_of_range < 3 )
        self.success = True


if __name__ == "__main__":
    unittest.main()

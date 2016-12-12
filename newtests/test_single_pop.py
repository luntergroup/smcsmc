from __future__ import print_function

import unittest
import os
import subprocess

from test_generic import *
import populationmodels

class TestGenericSinglePop(TestGeneric):

    def setUp(self):
        TestGeneric.setUp(self, "testdata/2sampleSingleConst")
        self.seqlen = 1e7
        self.pop = populationmodels.PopSingleConst( sequence_length = self.seqlen,
                                                    scrmpath=self.scrmpath )

        # set default inference parameters
        self.popt = None        # use the model's change points for inference
        self.em = 5
        self.np = 200
        self.seed = (1,)

        # set default expectations for population sizes
        self.target_min = [0] * 5        # 5 epochs for default model
        self.target_max = [1e10] * 5
        self.max_out_of_range = 0


@unittest.skip("Runs too slow")
class TestTwoSampleSingleExpand(TestGenericSinglePop):
    def setUp(self):
        TestGeneric.setUp(self, "testdata/2sampleSingleExpand")
        self.seqlen = 1e7
        self.pop = populationmodels.PopSingleExpand( sequence_length = self.seqlen,
                                                     scrmpath=self.scrmpath )
        # set default inference parameters
        self.popt = None        # use the model's change points for inference
        self.em = 5
        self.np = 200
        self.seed = (1,)

        # set default expectations for population sizes
        self.target_min = [0] * 5        # 5 epochs for default model
        self.target_max = [1e10] * 5
        self.max_out_of_range = 0



if __name__ == "__main__":
    unittest.main()

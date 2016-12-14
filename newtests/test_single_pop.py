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

        # set targets
        self.targets = []
        self.targets.append({'type':"Recomb", 'min':1e-9, 'max':1e-8, 'truth':5e-9})
        for idx in range(5):
            self.targets.append({'type':"Coal", 'pop':0, 'epoch':idx, 'min':9000, 'max':11000, 'truth':10000})
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

        # set targets
        self.targets = []
        self.targets.append({'type':"Recomb", 'min':1e-9, 'max':1e-8, 'truth':5e-9})
        for idx in range(5):
            self.targets.append({'type':"Coal", 'pop':0, 'epoch':idx, 'min':9000, 'max':11000, 'truth':10000})
        self.max_out_of_range = 0



if __name__ == "__main__":
    unittest.main()

from __future__ import print_function

import unittest
import os
import subprocess

from test_generic import *
import populationmodels

# NB: all targets were set using commit 34d52e

class TestGenericSinglePop(TestGeneric):

    def setUp(self):
        TestGeneric.setUp(self, "testdata/2sampleSingleConst")
        self.seqlen = 1e7
        self.pop = populationmodels.PopSingleConst( sequence_length = self.seqlen,
                                                    scrmpath=self.scrmpath )

        # set default inference parameters
        self.popt = None        # use the model's change points for inference
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = []
        self.targets.append({'type':"Recomb", 'min':1e-9, 'max':1e-8, 'truth':5e-9})
        for idx in range(5):
            self.targets.append({'type':"Coal", 'pop':0, 'epoch':idx, 'min':9000, 'max':11000, 'truth':10000})
        self.max_out_of_range = 0


class TestTwoSampleSingleExpand(TestGenericSinglePop):
    def setUp(self):
        TestGeneric.setUp(self, "testdata/2sampleSingleExpand")
        self.seqlen = 1e7
        self.num_samples = 2
        self.pop = populationmodels.PopSingleExpand( sequence_length = self.seqlen,
                                                     scrmpath=self.scrmpath,
                                                     num_samples=self.num_samples )
        # set default inference parameters
        self.popt = None        # use the model's change points for inference
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':5.13e-9, 'max':5.51e-9, 'truth':5e-9},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':0, 'max':3.5e9, 'truth':20000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':9252, 'max':10035, 'truth':10000}]

        self.max_out_of_range = 0


class TestFourSampleSingleExpand(TestGenericSinglePop):
    def setUp(self):
        TestGeneric.setUp(self, "testdata/4sampleSingleExpand")
        self.seqlen = 1e7
        self.num_samples = 4
        self.pop = populationmodels.PopSingleExpand( sequence_length = self.seqlen,
                                                     scrmpath=self.scrmpath,
                                                     num_samples=self.num_samples )
        # set default inference parameters
        self.popt = None        # use the model's change points for inference
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':5.66e-9, 'max':5.96e-9, 'truth':5e-9},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':0, 'max':173274, 'truth':20000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':10291, 'max':10930, 'truth':10000}]
        # these targets are not ideal, without pilots has a much wider range for epoch 0, is there a bug?

        self.max_out_of_range = 0


class TestEightSampleSingleExpand(TestGenericSinglePop):
    def setUp(self):
        TestGeneric.setUp(self, "testdata/8sampleSingleExpand")
        self.seqlen = 1e7
        self.num_samples = 8
        self.pop = populationmodels.PopSingleExpand( sequence_length = self.seqlen,
                                                     scrmpath=self.scrmpath,
                                                     num_samples=self.num_samples )
        # set default inference parameters
        self.popt = None        # use the model's change points for inference
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':6.57e-9, 'max':6.87e-9, 'truth':5e-9},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':20288, 'max':40042, 'truth':20000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':10608, 'max':11290, 'truth':10000}]

        self.max_out_of_range = 0


class TestSixteenSampleSingleExpand(TestGenericSinglePop):
    def setUp(self):
        TestGeneric.setUp(self, "testdata/16sampleSingleExpand")
        self.seqlen = 1e7
        self.num_samples = 16
        self.pop = populationmodels.PopSingleExpand( sequence_length = self.seqlen,
                                                     scrmpath=self.scrmpath,
                                                     num_samples=self.num_samples )
        # set default inference parameters
        self.popt = None        # use the model's change points for inference
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':8.42e-9, 'max':8.86e-9, 'truth':5e-9},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':16487, 'max':21128, 'truth':20000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':11123, 'max':11811, 'truth':10000}]

        self.max_out_of_range = 0


class TestFourSampleSingleShrink(TestGenericSinglePop):
    def setUp(self):
        TestGeneric.setUp(self, "testdata/4sampleSingleShrink")
        self.seqlen = 1e7
        self.num_samples = 4
        self.pop = populationmodels.PopSingleShrink( sequence_length = self.seqlen,
                                                     scrmpath=self.scrmpath,
                                                     num_samples=self.num_samples )
        # set default inference parameters
        self.popt = None        # use the model's change points for inference
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':5.40e-9, 'max':5.77e-9, 'truth':5e-9},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':3500, 'max':34252, 'truth':5000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':9612, 'max':10314, 'truth':10000}]

        self.max_out_of_range = 0


class TestFourSampleSingleBottleneck(TestGenericSinglePop):
    def setUp(self):
        TestGeneric.setUp(self, "testdata/4sampleSingleBottleneck")
        self.seqlen = 1e7
        self.num_samples = 4
        self.pop = populationmodels.PopSingleBottleneck( sequence_length = self.seqlen,
                                                     scrmpath=self.scrmpath,
                                                     num_samples=self.num_samples )
        # set default inference parameters
        self.popt = None        # use the model's change points for inference
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':5.57e-9, 'max':5.89e-9, 'truth':5e-9},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':4606, 'max':60033, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':9604, 'max':26360, 'truth':5000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':9780, 'max':10516, 'truth':10000}]

        self.max_out_of_range = 0

if __name__ == "__main__":
    unittest.main()

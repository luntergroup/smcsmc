from __future__ import print_function

import unittest
import os
import subprocess

from test_generic import *
import populationmodels

class TestGenericSinglePop(TestGeneric):

    # Called by unittest.main when TestTwoSample is created.
    # Responsible for per-class creation of simulated data.
    # It is implemented as a class method to avoid simulating the same data multiple times
    @classmethod
    def setUpClass(cls):
        cls.prefix = "2sampleSingleConst"
        cls.segfile = cls.prefix + ".seg"
        cls.seqlen = 1e7
        cls.scrmpath = "../scrm"
        cls.pop = populationmodels.PopSingleConst( filename = cls.segfile,
                                                   sequence_length = cls.seqlen,
                                                   scrmpath=cls.scrmpath )
        cls.success = True   # set ourselves up for success

        print ("simulating for",cls.prefix,"...")
        cls.pop.simulate()


    def setUp(self):
        print ("called TestGenericSinglePop.setUp()")
        TestGeneric.setUp(self)
        self.popt = None        # use the model's change points for inference
        self.npop = 1


    # actual test
    def test_inference(self):
        self.em = 5
        self.np = 200
        self.seed = (3647837471,)
        self.infer( case = 1 )

        # Ne target ranges for the varous epochs
        target_min = self.pop.target_min
        target_max = self.pop.target_max
        results = self.readResults()

        out_of_range = 0
        for epoch, emresults in enumerate(results['Coal']):
            print(`epoch` + `emresults` )
            # get the last EM estimate
            result = emresults['estimates'][-1]

            # see if it is within range
            print ("Checking epoch",epoch,result,target_min[epoch],target_max[epoch])
            #self.assertTrue( result >= target_min[epoch] )
            #self.assertTrue( result <= target_max[epoch] )
            if result < target_min[epoch] or result > target_max[epoch]:
                print ("Out of range!")
                out_of_range += 1

        self.assertTrue( out_of_range < 3 )
        self.success = True


@unittest.skip("Runs too slow")
class TestTwoSampleSingleExpand(TestGenericSinglePop):
    @classmethod
    def setUpClass(cls):
        cls.prefix = "2sampleSingleExpand"
        cls.segfile = cls.prefix + ".seg"
        cls.seqlen = 1e7
        cls.scrmpath = "../scrm"
        cls.pop = populationmodels.PopSingleExpand( filename = cls.segfile,
                                                    sequence_length = cls.seqlen,
                                                    scrmpath=cls.scrmpath)
        cls.success = True   # set ourselves up for success

        print ("simulating for",cls.prefix,"...")
        cls.pop.simulate()


@unittest.skip("Runs too slow")
class TestTwoSampleSingleShrink(TestGenericSinglePop):
    @classmethod
    def setUpClass(cls):
        cls.prefix = "2sampleSingleShrink"
        cls.segfile = cls.prefix + ".seg"
        cls.seqlen = 1e7
        cls.scrmpath = "../scrm"
        cls.pop = populationmodels.PopSingleShrink( filename = cls.segfile,
                                                    sequence_length = cls.seqlen,
                                                    scrmpath=cls.scrmpath )
        cls.success = True   # set ourselves up for success

        print ("simulating for",cls.prefix,"...")
        cls.pop.simulate()


@unittest.skip("Runs too slow")
class TestThreeSampleSingleConst(TestGenericSinglePop):
    @classmethod
    def setUpClass(cls):
        cls.prefix = "3sampleSingleConst"
        cls.segfile = cls.prefix + ".seg"
        cls.seqlen = 1e7
        cls.scrmpath = "../scrm"
        cls.pop = populationmodels.PopSingleConst( filename = cls.segfile,
                                                   sequence_length = cls.seqlen,
                                                   scrmpath=cls.scrmpath,
                                                   num_samples=3 )
        cls.success = True   # set ourselves up for success

        print ("simulating for",cls.prefix,"...")
        cls.pop.simulate()


@unittest.skip("Runs too slow")
class TestThreeSampleSingleExpand(TestGenericSinglePop):
    @classmethod
    def setUpClass(cls):
        cls.prefix = "3sampleSingleExpand"
        cls.segfile = cls.prefix + ".seg"
        cls.seqlen = 1e7
        cls.scrmpath = "../scrm"
        cls.pop = populationmodels.PopSingleExpand( filename = cls.segfile,
                                                    sequence_length = cls.seqlen,
                                                    scrmpath=cls.scrmpath,
                                                    num_samples=3 )
        cls.success = True   # set ourselves up for success

        print ("simulating for",cls.prefix,"...")
        cls.pop.simulate()


@unittest.skip("Runs too slow")
class TestThreeSampleSingleShrink(TestGenericSinglePop):
    @classmethod
    def setUpClass(cls):
        cls.prefix = "3sampleSingleShrink"
        cls.segfile = cls.prefix + ".seg"
        cls.seqlen = 1e7
        cls.scrmpath = "../scrm"
        cls.pop = populationmodels.PopSingleShrink( filename = cls.segfile,
                                                    sequence_length = cls.seqlen,
                                                    scrmpath=cls.scrmpath,
                                                    num_samples=3 )
        cls.success = True   # set ourselves up for success

        print ("simulating for",cls.prefix,"...")
        cls.pop.simulate()

@unittest.skip("Runs too slow")
class TestFourSampleSingleConst(TestGenericSinglePop):
    @classmethod
    def setUpClass(cls):
        cls.prefix = "4sampleSingleConst"
        cls.segfile = cls.prefix + ".seg"
        cls.seqlen = 1e7
        cls.scrmpath = "../scrm"
        cls.pop = populationmodels.PopSingleConst( filename = cls.segfile,
                                                   sequence_length = cls.seqlen,
                                                   scrmpath=cls.scrmpath,
                                                   num_samples=4 )
        cls.success = True   # set ourselves up for success

        print ("simulating for",cls.prefix,"...")
        cls.pop.simulate()


@unittest.skip("Runs too slow")
class TestFourSampleSingleExpand(TestGenericSinglePop):
    @classmethod
    def setUpClass(cls):
        cls.prefix = "4sampleSingleExpand"
        cls.segfile = cls.prefix + ".seg"
        cls.seqlen = 1e7
        cls.scrmpath = "../scrm"
        cls.pop = populationmodels.PopSingleExpand( filename = cls.segfile,
                                                    sequence_length = cls.seqlen,
                                                    scrmpath=cls.scrmpath,
                                                    num_samples=4 )
        cls.success = True   # set ourselves up for success

        print ("simulating for",cls.prefix,"...")
        cls.pop.simulate()


@unittest.skip("Runs too slow")
class TestFourSampleSingleShrink(TestGenericSinglePop):
    @classmethod
    def setUpClass(cls):
        cls.prefix = "4sampleSingleShrink"
        cls.segfile = cls.prefix + ".seg"
        cls.seqlen = 1e7
        cls.scrmpath = "../scrm"
        cls.pop = populationmodels.PopSingleShrink( filename = cls.segfile,
                                                    sequence_length = cls.seqlen,
                                                    scrmpath=cls.scrmpath,
                                                    num_samples=4 )
        cls.success = True   # set ourselves up for success

        print ("simulating for",cls.prefix,"...")
        cls.pop.simulate()


if __name__ == "__main__":
    unittest.main()

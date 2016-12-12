from __future__ import print_function

import unittest
import populationmodels
from test_generic import TestGeneric


#
# Test various basic constant population size scenarios
#

class TestMigration(TestGeneric):

    def setUp(self):
        TestGeneric.setUp(self, "testdata/migration")
        self.seqlen = 1e7
        self.pop = populationmodels.TwoPopUniDirMigr( sequence_length = self.seqlen,
                                                      scrmpath=self.scrmpath )
        
        # set default inference parameters
        self.em = 4
        self.np = 1000
        self.seed = (1,)

        # set targets
        self.target_min = [45000, 45000, 9500]
        self.target_max = [55000, 55000, 10500]
        self.max_out_of_range = 1

        self.migr_targets = [{'from':0, 'to': 1, 'min':[0,0,0], 'max': [3e-5,3e-5,3e-5]}]
        self.migr_max_out_of_range = 1

    # actual test
    def test_inference(self):

        # Call test_inference on parent class to generate data, perform inference, and check population sizes
        TestGeneric.test_inference(self)
        self.success = False

        # Read results from .out file (again)
        results = self.readResults()

        # Check migration results and count out-of-range epochs
        out_of_range = 0
        for epoch, emresults in enumerate(results['Migr']):

            # TODO:
            # this needs work: currently migration results are not properly stored in the results[] array.
            # should probably store this differently; the results array is a bit of a hack.
            # perhaps have 3 results arrays, for recombination, coalescence, and migration rates?

        self.assertTrue( out_of_range < self.max_out_of_range )
        self.success = True


if __name__ == "__main__":
    unittest.main()

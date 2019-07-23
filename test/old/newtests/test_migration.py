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
        self.num_samples=8
        self.pop = populationmodels.TwoPopUniDirMigr( sequence_length = self.seqlen,
                                                      scrmpath=self.scrmpath,
                                                      num_samples=self.num_samples )
        
        # set default inference parameters
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':5.71e-9, 'max':5.92e-9, 'truth':5e-9},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':13866, 'max':19228, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':12688, 'max':14835, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':8943,  'max':9937,  'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':0, 'min':12569, 'max':17354, 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':1, 'min':9242,  'max':10604, 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':2, 'min':7821,  'max':8671, 'truth':10000},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0, 'min':1.42e-6, 'max':8.42e-6, 'truth':2.5e-5},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':1, 'min':1.64e-5, 'max':2.44e-5, 'truth':2.5e-5},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':2, 'min':2.62-5,  'max':3.06e-5, 'truth':2.5e-5},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':0, 'min':0,       'max':3.9e-6,  'truth':0},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':1, 'min':6.40e-6, 'max':1.35e-5, 'truth':0},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':2, 'min':2.43e-5, 'max':2.87e-5, 'truth':0}]

        self.max_out_of_range = 1


class TestTwoPopUniDirMigr(TestMigration):

    def setUp(self):
        TestGeneric.setUp(self, "testdata/uni_dir_migration")
        self.seqlen = 1e7
        self.num_samples=8
        self.pop = populationmodels.TwoPopUniDirMigr( sequence_length = self.seqlen,
                                                      scrmpath=self.scrmpath,
                                                      num_samples=self.num_samples )
        
        # set default inference parameters
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':5.71e-9, 'max':5.92e-9, 'truth':5e-9},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':13866, 'max':19228, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':12688, 'max':14835, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':8943,  'max':9937,  'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':0, 'min':12569, 'max':17354, 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':1, 'min':9242,  'max':10604, 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':2, 'min':7821,  'max':8671, 'truth':10000},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0, 'min':1.42e-6, 'max':8.42e-6, 'truth':2.5e-5},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':1, 'min':1.64e-5, 'max':2.44e-5, 'truth':2.5e-5},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':2, 'min':2.62-5,  'max':3.06e-5, 'truth':2.5e-5},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':0, 'min':0,       'max':3.9e-6,  'truth':0},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':1, 'min':6.40e-6, 'max':1.35e-5, 'truth':0},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':2, 'min':2.43e-5, 'max':2.87e-5, 'truth':0}]

        self.max_out_of_range = 1


class TestTwoPopBiDirMigr(TestMigration):

    def setUp(self):
        TestGeneric.setUp(self, "testdata/bi_dir_migration")
        self.seqlen = 1e7
        self.num_samples=8
        self.pop = populationmodels.TwoPopBiDirMigr( sequence_length = self.seqlen,
                                                      scrmpath=self.scrmpath,
                                                      num_samples=self.num_samples )
        
        # set default inference parameters
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':5.85e-9, 'max':6.08e-9, 'truth':5e-9},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':16000, 'max':23132, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':13929, 'max':16529, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':9774,  'max':10819, 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':0, 'min':16564, 'max':22886, 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':1, 'min':14485, 'max':16783, 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':2, 'min':9930,  'max':10765, 'truth':10000},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0, 'min':1.10e-5, 'max':2.68e-5, 'truth':2.5e-5},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':1, 'min':2.63e-5, 'max':3.58e-5, 'truth':2.5e-5},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':2, 'min':2.54e-5, 'max':2.91e-5, 'truth':2.5e-5},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':0, 'min':9.77e-6, 'max':2.93e-5, 'truth':2.5e-5},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':1, 'min':2.80e-5, 'max':3.69e-5, 'truth':2.5e-5},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':2, 'min':2.58e-5, 'max':2.91e-5, 'truth':2.5e-5}]

        self.max_out_of_range = 1


class TestTwoPopSplitNoMigr(TestMigration):

    def setUp(self):
        TestGeneric.setUp(self, "testdata/split_no_migration")
        self.seqlen = 1e7
        self.num_samples=8
        self.pop = populationmodels.TwoPopSplitNoMigr( sequence_length = self.seqlen,
                                                      scrmpath=self.scrmpath,
                                                      num_samples=self.num_samples )
        
        # set default inference parameters
        self.em = 4
        self.np = 100
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':6.08e-9, 'max':6.32e-9, 'truth':5e-9},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':13616, 'max':19257, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':9981,  'max':11617, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':8819,  'max':9665,  'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':0, 'min':12186, 'max':16499, 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':1, 'min':9294,  'max':10710, 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':2, 'min':9999,  'max':10001, 'truth':"NA"},  #defaults: pop no longer exists
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0, 'min':0,       'max':2.58e-6, 'truth':0},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':1, 'min':8.88e-6, 'max':1.49e-5, 'truth':0},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':2, 'min':0,       'max':0,       'truth':"NA"}, #defaults: to_pop no longer exists
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':0, 'min':0,       'max':2.28e-6, 'truth':0},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':1, 'min':8.71e-6, 'max':1.47e-5, 'truth':0},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':2, 'min':2.5e-5,  'max':2.51e-5, 'truth':"NA"}] #defaults: to_pop no longer exists

        self.max_out_of_range = 1

if __name__ == "__main__":
    unittest.main()

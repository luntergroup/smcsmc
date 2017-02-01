from __future__ import print_function

import unittest
import populationmodels
from test_generic import TestGeneric


#
# Test simple models under no-data conditions, but under various pilot biasing conditions
#
# - all tests run in <60 s
# - tests have been checked for robustness (seed 1-10)
#

class TestBias_00_NoBias(TestGeneric):

    def setUp(self, fn="testdata/bias_nobias"):
        TestGeneric.setUp(self, fn)
        self.seqlen = 1e6
        self.num_samples = 8
        self.pop = populationmodels.Population( sequence_length = self.seqlen,
                                                scrmpath = self.scrmpath,
                                                num_samples = self.num_samples,
                                                num_populations = 1,
                                                change_points = [0, 0.5, 1],
                                                population_sizes = [[1.2], [0.8], [1.0]] )

        # set default inference parameters
        self.missing_leaves = range(self.num_samples)
        self.smcsmc_initial_pop_sizes = self.pop.population_sizes
        self.lag = 2
        self.popt = None
        self.em = 0
        self.np = 100
        self.seed = (8,)
        self.bias_heights = None
        self.bias_strengths = None

        # set targets
        self.targets = [{'type':"Recomb",
                         'min':9.9e-9, 'max':1.01e-8, 'truth':1e-8, 'ess':self.np-1},
                        {'type':"Coal", 'pop':0, 'epoch':0,
                         'min':11800, 'max':12200, 'truth':12000, 'ess':self.np-1},
                        {'type':"Coal", 'pop':0, 'epoch':1,
                         'min':7800,  'max':8200,  'truth': 8000, 'ess':self.np-1},
                        {'type':"Coal", 'pop':0, 'epoch':2,
                         'min':9800,  'max':10200, 'truth':10000, 'ess':self.np-1}]
        self.max_out_of_range = 0


        
class TestBias_01_Bias2(TestBias_00_NoBias):

    def setUp(self):
        TestBias_00_NoBias.setUp(self, "testdata/bias_2bias")

        self.np = 200
        self.bias_heights = [400]
        self.bias_strengths = [2, 1]

        for idx in range(len(self.targets)):
            self.targets[idx]['ess'] = self.np / 3



class TestBias_02_Bias5(TestBias_00_NoBias):

    def setUp(self):
        TestBias_00_NoBias.setUp(self, "testdata/bias_5bias")

        self.np = 500
        self.bias_heights = [400]
        self.bias_strengths = [5, 1]

        # recombination; and epochs 1-3
        for idx,ess in enumerate( [30, 30, 50, 100] ):
            self.targets[idx]['ess'] = ess


class TestBias_0_Migr_NoBias(TestBias_00_NoBias):

    def setUp(self, fn="testdata/bias_migr"):
        TestBias_00_NoBias.setUp(self, fn)

        self.np = 100
        self.pop.num_populations = 2
        self.pop.sample_populations = [1,1,1,1,2,2,2,2]
        self.pop.population_sizes = [[1.2, 0.8], [0.8, 1.0], [1.0, 1.2]]
        self.pop.migration_rates = [ [ [0, 0.5],
                                       [0.25, 0] ],
                                     [ [0, 0.25],
                                       [0.5, 0] ],
                                     [ [0, 0.75],
                                       [0, 0] ] ]
        # make sure smc^2 is passed the right population sizes
        self.smcsmc_initial_pop_sizes = self.pop.population_sizes

        scaling = 1.0 / (4.0 * self.pop.N0)

        self.targets.append( {'type':"Coal", 'pop':1, 'epoch':0,
                              'min':7800,  'max':8200, 'truth':8000,  'ess':self.np-1} )
        self.targets.append( {'type':"Coal", 'pop':1, 'epoch':1,
                              'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )
        self.targets.append( {'type':"Coal", 'pop':1, 'epoch':2,
                              'min':11800, 'max':12200,'truth':12000, 'ess':self.np-1} )
        
        self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                              'min':0.45*scaling, 'max': 0.55*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':0,
                              'min':0.21*scaling, 'max': 0.29*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':1,
                              'min':0.23*scaling, 'max': 0.27*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':1,
                              'min':0.45*scaling, 'max': 0.55*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':2,
                              'min':0.72*scaling, 'max': 0.78*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':2,
                              'min':0.00*scaling, 'max': 0.00*scaling} )

        for t in self.targets:
            if t['type'] == "Migr":
                t['truth'] = scaling * self.pop.migration_rates[t['epoch']][t['from_pop']][t['to_pop']]



class TestBias_1_Migr_Bias2(TestBias_0_Migr_NoBias):

    def setUp(self):
        TestBias_0_Migr_NoBias.setUp(self, "testdata/bias2_migr")

        self.np = 200
        self.bias_heights = [400]
        self.bias_strengths = [2, 1]
        
        for idx in range(len(self.targets)):
            self.targets[idx]['ess'] = self.np / 3
        

            
        
if __name__ == "__main__":
    unittest.main()

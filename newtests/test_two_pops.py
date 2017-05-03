from __future__ import print_function
from context import populationmodels
import unittest
from test_generic import TestGeneric


#
# Test various two-population uni-directional migration scenarios
#

class TestTwoPops(TestGeneric):

    def setUp(self, fn="testdata/twopops"):
        TestGeneric.setUp(self, fn)
        self.seqlen = 1e6
        self.num_samples = 8
        self.pop = populationmodels.Population( sequence_length = self.seqlen,
                                                scrmpath = self.scrmpath,
                                                num_samples = self.num_samples,
                                                num_populations = 2,
                                                change_points = [0, 0.5, 1],
                                                population_sizes = [ [1,1], [1,1], [1,1] ],
                                                sample_populations = [1,1,1,1,2,2,2,2] )

        # set default inference parameters
        self.smcsmc_initial_pop_sizes = self.pop.population_sizes

        self.lag = 1
        self.popt = None
        self.em = 4
        self.np = 100
        self.seed = (8,)
        self.bias_heights = None
        self.bias_strengths = None

        # set targets
        self.targets = [{'type':"Recomb",
                         'min':9.9e-9, 'max':1.01e-8, 'truth':1e-8, 'ess':self.np-1},
                        {'type':"Coal", 'pop':0, 'epoch':0,
                         'min':9800, 'max':10200, 'truth':10000, 'ess':self.np-1},
                        {'type':"Coal", 'pop':0, 'epoch':1,
                         'min':9800,  'max':10200,  'truth': 10000, 'ess':self.np-1},
                        {'type':"Coal", 'pop':0, 'epoch':2,
                         'min':9800,  'max':10200, 'truth':10000, 'ess':self.np-1}]
        self.max_out_of_range = 0

class TestTwoPopsSplit(TestTwoPops):
    
    def setUp(self, fn="testdata/twopopssplit"):
        TestTwoPops.setUp(self, fn)
        
        self.pop.split_command = "-ej .5 2 1"

# This is an important regression test,
# Note: targets were calculated based on 100 runs of 2b3a_wo_apply_immediately_hack 
class TestTwoPopsSplitUniDirMigr(TestTwoPopsSplit):

    def setUp(self, fn="testdata/twopopssplit_unidirmigr"):
        TestTwoPopsSplit.setUp(self, fn)

        self.np = 1000
        self.seqlen = 30e6
        self.pop.sequence_length = self.seqlen
        self.pop.change_points = [0,0.1,0.5]
        self.pop.migration_rates = [ [ [0, 0.2],[0, 0] ],
                                     [ [0, 0.2],[0, 0] ],
                                     [ [0, 0]  ,[0, 0] ] ]

        self.smcsmc_initial_migr_rates = [ [[0,.2],[.2,0]],
                                           [[0,.2],[.2,0]],
                                           [[0,0] ,[0,0] ] ]

        # make sure smc^2 is passed the right population sizes
        self.smcsmc_initial_pop_sizes = self.pop.population_sizes

        scaling = 1.0 / (4.0 * self.pop.N0)

        self.targets = [{'type':"Coal", 'pop':0, 'epoch':0,
                         'min':8209, 'max':9457, 'truth':10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':1, 'epoch':0,
                         'min':8612,  'max':10077,'truth':10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':0, 'epoch':1,
                         'min':10333,  'max':10878,  'truth': 10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':1, 'epoch':1,
                         'min':9944,  'max':10504,  'truth': 10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':0, 'epoch':2,
                         'min':9374,  'max':9584, 'truth':10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':1, 'epoch':2,
                         'min':10000,  'max':10000, 'truth':10000, 'ess':self.np-1},
                         
                        {'type':"Recomb",
                         'min':9.91e-9, 'max':1.016e-8, 'truth':1e-8, 'ess':self.np-1},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':0, 'max':1.40e-6},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':0, 'max':5.32e-7},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':2.64e-6, 'max':3.98e-6},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':2.29e-6, 'max':3.47e-6},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':0, 'max':0},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':5e-6, 'max':5e-6} ]

        for t in self.targets:
            if t['type'] == "Migr":
                t['truth'] = scaling * self.pop.migration_rates[t['epoch']][t['from_pop']][t['to_pop']]


# This is an important regression test, 
# Note: targets were calculated based on 100 runs of 2b3a_wo_apply_immediately_hack 
class TestTwoPopsSplitUniDirMigrInRecentEpoch(TestTwoPopsSplit):

    def setUp(self, fn="testdata/twopopssplit_unidirmigrinrecentepoch"):
        TestTwoPopsSplit.setUp(self, fn)

        self.np = 1000
        self.seqlen = 30e6
        self.pop.change_points = [0,0.1,0.5]
        self.pop.migration_rates = [ [ [0, 0.2],[0, 0] ],
                                     [ [0, 0]  ,[0, 0] ],
                                     [ [0, 0]  ,[0, 0] ] ]

        self.smcsmc_initial_migr_rates = [ [[0,.2],[.2,0]],
                                           [[0,.2],[.2,0]],
                                           [[0,0] ,[0,0] ] ]

        # make sure smc^2 is passed the right population sizes
        self.smcsmc_initial_pop_sizes = self.pop.population_sizes

        scaling = 1.0 / (4.0 * self.pop.N0)

        self.targets = [{'type':"Coal", 'pop':0, 'epoch':0,
                         'min':8539, 'max':9749, 'truth':10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':1, 'epoch':0,
                         'min':8453,  'max':9765,'truth':10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':0, 'epoch':1,
                         'min':9710,  'max':10165,  'truth': 10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':1, 'epoch':1,
                         'min':9427,  'max':9827,  'truth': 10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':0, 'epoch':2,
                         'min':9548,  'max':9755, 'truth':10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':1, 'epoch':2,
                         'min':10000,  'max':10000, 'truth':10000, 'ess':self.np-1},
                         
                        {'type':"Recomb",
                         'min':9.87e-9, 'max':1.015e-8, 'truth':1e-8, 'ess':self.np-1},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':0, 'max':1.10e-6},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':0, 'max':3.33e-7},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':2.02e-6, 'max':2.98e-6},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':1.72e-6, 'max':2.77e-6},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':0, 'max':0},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':5e-6, 'max':5e-6} ]

        for t in self.targets:
            if t['type'] == "Migr":
                t['truth'] = scaling * self.pop.migration_rates[t['epoch']][t['from_pop']][t['to_pop']]

# This is an important regression test, 
# Note: targets were calculated based on 100 runs of 2b3a_wo_apply_immediately_hack                
class TestTwoPopsSplitUniDirMigrInMidEpoch(TestTwoPopsSplit):

    def setUp(self, fn="testdata/twopopssplit_unidirmigrinmidepoch"):
        TestTwoPopsSplit.setUp(self, fn)

        self.np = 1000
        self.seqlen = 30e6
        self.pop.change_points = [0,0.1,0.5]
        self.pop.migration_rates = [ [ [0, 0]  ,[0, 0] ],
                                     [ [0, 0.2],[0, 0] ],
                                     [ [0, 0]  ,[0, 0] ] ]

        self.smcsmc_initial_migr_rates = [ [[0,.2],[.2,0]],
                                           [[0,.2],[.2,0]],
                                           [[0,0] ,[0,0] ] ]

        # make sure smc^2 is passed the right population sizes
        self.smcsmc_initial_pop_sizes = self.pop.population_sizes

        scaling = 1.0 / (4.0 * self.pop.N0)

        self.targets = [{'type':"Coal", 'pop':0, 'epoch':0,
                         'min':8653, 'max':10198, 'truth':10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':1, 'epoch':0,
                         'min':7897,  'max':9085,'truth':10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':0, 'epoch':1,
                         'min':10057,  'max':10550,  'truth': 10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':1, 'epoch':1,
                         'min':9866,  'max':10294,  'truth': 10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':0, 'epoch':2,
                         'min':9195,  'max':9410, 'truth':10000, 'ess':self.np-1},
                         
                        {'type':"Coal", 'pop':1, 'epoch':2,
                         'min':10000,  'max':10000, 'truth':10000, 'ess':self.np-1},
                         
                        {'type':"Recomb",
                         'min':9.88e-9, 'max':1.014e-8, 'truth':1e-8, 'ess':self.np-1},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':0, 'max':8.49e-7},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':0, 'max':4.47e-7},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':2.64e-6, 'max':3.88e-6},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':2.29e-6, 'max':3.33e-6},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':0, 'max':0},
                         
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                         'min':5e-6, 'max':5e-6} ]

        for t in self.targets:
            if t['type'] == "Migr":
                t['truth'] = scaling * self.pop.migration_rates[t['epoch']][t['from_pop']][t['to_pop']]



#class TestTwoPopsBiDirMigr(TestTwoPops):

    #def setUp(self):
        #TestTwoPops.setUp(self, "testdata/twopops_bidirmigr")

        #sself.np = 100
        #self.pop.migration_rates = [ [ [0, 0.2], [0.2, 0] ],
                                     #[ [0, 0.2], [0.2, 0] ],
                                     #[ [0, 0.2], [0.2, 0] ] ]

        ## make sure smc^2 is passed the right population sizes
        #self.smcsmc_initial_pop_sizes = self.pop.population_sizes

        #scaling = 1.0 / (4.0 * self.pop.N0)

        #self.targets.append( {'type':"Coal", 'pop':1, 'epoch':0,
                              #'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )
        #self.targets.append( {'type':"Coal", 'pop':1, 'epoch':1,
                              #'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )
        #self.targets.append( {'type':"Coal", 'pop':1, 'epoch':2,
                              #'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )

        #self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                              #'min':0.15*scaling, 'max': 0.25*scaling} )
        #self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':0,
                              #'min':0.15*scaling, 'max': 0.25*scaling} )
        #self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':1,
                              #'min':0.15*scaling, 'max': 0.25*scaling} )
        #self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':1,
                              #'min':0.15*scaling, 'max': 0.25*scaling} )
        #self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':2,
                              #'min':0.15*scaling, 'max': 0.25*scaling} )
        #self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':2,
                              #'min':0.15*scaling, 'max': 0.25*scaling} )

        #for t in self.targets:
            #if t['type'] == "Migr":
                #t['truth'] = scaling * self.pop.migration_rates[t['epoch']][t['from_pop']][t['to_pop']]




if __name__ == "__main__":
    unittest.main()

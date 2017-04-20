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
        self.em = 0
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
# TO DO: calculate targets based on earlier stable verions 
#        old targets used 34d52e63b6e7, lets update to a59a96d439c52e as it is less likely to hang in calculate_median_survival_distance
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

        self.targets.append( {'type':"Coal", 'pop':1, 'epoch':0,
                              'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )
        self.targets.append( {'type':"Coal", 'pop':1, 'epoch':1,
                              'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )
        self.targets.append( {'type':"Coal", 'pop':1, 'epoch':2,
                              'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )

        self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                              'min':0.15*scaling, 'max': 0.25*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':0,
                              'min':0*scaling, 'max': 0.0001*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':1,
                              'min':0.15*scaling, 'max': 0.25*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':1,
                              'min':0*scaling, 'max': 0.0001*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':2,
                              'min':0.15*scaling, 'max': 0.25*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':2,
                              'min':0*scaling, 'max': 0.0001*scaling} )

        for t in self.targets:
            if t['type'] == "Migr":
                t['truth'] = scaling * self.pop.migration_rates[t['epoch']][t['from_pop']][t['to_pop']]


# This is an important regression test, 
# TO DO: calculate targets based on earlier stable verions 
#        old targets used 34d52e63b6e7, lets update to a59a96d439c52e as it is less likely to hang in calculate_median_survival_distance
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

        self.targets.append( {'type':"Coal", 'pop':1, 'epoch':0,
                              'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )
        self.targets.append( {'type':"Coal", 'pop':1, 'epoch':1,
                              'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )
        self.targets.append( {'type':"Coal", 'pop':1, 'epoch':2,
                              'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )

        self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                              'min':0.15*scaling, 'max': 0.25*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':0,
                              'min':0*scaling, 'max': 0.0001*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':1,
                              'min':0.15*scaling, 'max': 0.25*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':1,
                              'min':0*scaling, 'max': 0.0001*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':2,
                              'min':0.15*scaling, 'max': 0.25*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':2,
                              'min':0*scaling, 'max': 0.0001*scaling} )

        for t in self.targets:
            if t['type'] == "Migr":
                t['truth'] = scaling * self.pop.migration_rates[t['epoch']][t['from_pop']][t['to_pop']]

# This is an important regression test, 
# TO DO: calculate targets based on earlier stable verions 
#        old targets used 34d52e63b6e7, lets update to a59a96d439c52e as it is less likely to hang in calculate_median_survival_distance                
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

        self.targets.append( {'type':"Coal", 'pop':1, 'epoch':0,
                              'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )
        self.targets.append( {'type':"Coal", 'pop':1, 'epoch':1,
                              'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )
        self.targets.append( {'type':"Coal", 'pop':1, 'epoch':2,
                              'min':9800,  'max':10200,'truth':10000, 'ess':self.np-1} )

        self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0,
                              'min':0.15*scaling, 'max': 0.25*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':0,
                              'min':0*scaling, 'max': 0.0001*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':1,
                              'min':0.15*scaling, 'max': 0.25*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':1,
                              'min':0*scaling, 'max': 0.0001*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':2,
                              'min':0.15*scaling, 'max': 0.25*scaling} )
        self.targets.append( {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':2,
                              'min':0*scaling, 'max': 0.0001*scaling} )

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

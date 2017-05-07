from __future__ import print_function
from context import populationmodels
import unittest
import itertools
from test_generic import TestGeneric


#
# Test various basic constant population size scenarios
#

# NOTE: it seems that a lag of 4 gives more accurate parameter estimates, despite a lower ESS

class TestConstPopSize(TestGeneric):

    def setUp(self, name="testdata/constpopsize"):
        TestGeneric.setUp(self, name)
        self.seqlen = 1e7
        self.pop = populationmodels.Population( sequence_length = self.seqlen,
                                                scrmpath=self.scrmpath,
                                                change_points = [0, 0.01, 0.25, 0.5, 1, 1.5],
                                                num_populations = 1,
                                                population_sizes = [[1], [1], [1], [1], [1], [1]])

        # use python front-end
        self.smcsmcpath = "../python/smcsmc.py"
        
        # set default inference parameters
        self.popt = None
        self.smcsmc_initial_pop_sizes = self.pop.population_sizes
        self.lag = 2
        self.em = 1
        self.np = 1000
        self.debug = True
        self.bias_heights = [400]
        self.bias_strengths = [3,1]
        self.tmax = 4
        self.seed = (1,)

        self.alpha = 0.5  #TEST

        # set targets
        self.targets = []
        self.targets.append({'type':"Recomb", 'min':9.6e-9, 'max':1.04e-8, 'truth':1e-8, 'ess':20})
        for idx in range(6):
            self.targets.append({'type':"Coal", 'pop':0, 'epoch':idx, 'min':9600, 'max':10400,
                                 'truth':10000, 'ess':[1,15,40,60,60,35][idx]})
        self.targets[1]['min'] = 6000
        self.targets[1]['max'] = 15000
        self.max_out_of_range = -1



class TestConstPopSize_ThreeEpochs(TestGeneric):

    def setUp(self, name="testdata/constpopsize_3epochs"):
        # base class for a number of experiments.
        TestGeneric.setUp(self, name)
        self.seqlen = 10000   # dummy length, for test
        self.pop = populationmodels.Population( sequence_length = self.seqlen,
                                                scrmpath=self.scrmpath,
                                                change_points = [0, 1, 2],   # three epochs
                                                num_populations = 1,
                                                population_sizes = [[1], [1], [1]] )
        self.em = 5
        self.np = 100
        self.bias_heights = None
        self.bias_strenghts = None
        self.alpha = 0
        self.debug = True

        # no tests here
        self.targets = []
        self.max_out_of_range = 0



class TestConstPopSize_Migration(TestGeneric):

    def setUp(self):
        TestGeneric.setUp(self, "testdata/constpopsize_migration")
        self.seqlen = 1e6
        self.pop = populationmodels.Population( sequence_length = self.seqlen,
                                                scrmpath=self.scrmpath,
                                                change_points = [0, 0.01, 0.25, 0.5, 1, 1.5],
                                                num_populations = 2,
                                                num_samples = 8,
                                                sample_populations = [1,1,1,1,2,2,2,2],
                                                population_sizes = [[1,1]] * 6,
                                                migration_rates = [ [[0,0],[1,0]], # -em 0 2 1 1
                                                                    [[0,0],[1,0]], # em 0.01 2 1 1
                                                                    [[0,0],[1,0]],
                                                                    [[0,0],[1,0]],
                                                                    [[0,0],[1,0]],
                                                                    [[0,0],[1,0]] ] )
        # use python front-end
        self.smcsmcpath = "../python/smcsmc.py"
        
        # set default inference parameters
        self.popt = None
        self.smcsmc_initial_pop_sizes = self.pop.population_sizes
        self.lag = 2
        self.em = 1
        self.np = 1000
        self.debug = True
        self.bias_heights = [800]
        self.bias_strengths = [2,1]
        self.tmax = 4
        self.seed = (1,)

        # set targets
        self.targets = [{'type':"Recomb", 'min':9.6e-9, 'max':1.04e-8, 'truth':1e-8, 'ess':1.5}]
        for idx, pop in itertools.product(range(6), range(2)):
            self.targets.append({'type':"Coal", 'pop':pop, 'epoch':idx,
                                 'min':9300  if idx>1 else 6000,
                                 'max': [15000,11500,10600,10600,10600,11500][idx],
                                 'truth':10000, 'ess':[1,1,1.3,1.5,2,3][idx]})
        self.targets += [ {'type':'Migr', 'from_pop':fp, 'to_pop':tp, 'epoch':ep,
                           'min':mi, 'max':ma, 'truth':tr}
                          for (fp,tp,ep,mi,ma,tr) in
                          [ (0,1,0, 0.0, 1e-5, 0),
                            (0,1,1, 0.0, 1e-5, 0),
                            (0,1,2, 0.0, 1e-5, 0),
                            (0,1,3, 0.0, 1e-5, 0),
                            (0,1,4, 0.0, 1e-5, 0),
                            (0,1,5, 0.0, 1e-5, 0),
                            (1,0,0, 0.0, 3e-5, 2.5e-5),
                            (1,0,1, 0.0, 3e-5, 2.5e-5),
                            (1,0,2, 0.0, 3e-5, 2.5e-5),
                            (1,0,3, 0.0, 3e-5, 2.5e-5),
                            (1,0,4, 0.0, 3e-5, 2.5e-5),
                            (1,0,5, 0.0, 3e-5, 2.5e-5) ] ]
                               
        self.max_out_of_range = -1
        


class TestConstPopSize_FourEpochs(TestConstPopSize):

    def setUp(self, name = "testdata/constpopsize_4epochs"):
        TestConstPopSize.setUp(self, name)
        self.seqlen = 1e7
        self.missing_leaves = []
        self.popt = None                               # don't use epoch pattern for inference
        self.smcsmc_change_points = [0, .02, .1, .5]   # but rather use four epochs

        self.pop = populationmodels.Pop2( sequence_length = self.seqlen,
                                          change_points = [0, .02, .1, .5],
                                          population_sizes = [[1], [1], [1], [1]],
                                          scrmpath=self.scrmpath )
        
        # use python front-end
        self.smcsmcpath = "../python/smcsmc.py"

        # set default inference parameters
        self.em = 1
        self.np = 1000
        self.bias_heights = [800]
        self.bias_strengths = [1.5,1]
        self.seed = (1,)
        self.debug = True

        self.alpha = 0.5  #TEST

        # set targets
        self.targets = [{'type':"Recomb", 'min':0.97e-8, 'max':1.03e-8, 'truth':1e-8, 'ess':15},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':4000,  'max':18000, 'truth':10000, 'ess':1},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':8500,  'max':11500, 'truth':10000, 'ess':5},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':9700,  'max':10300, 'truth':10000, 'ess':20},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':9800,  'max':10200, 'truth':10000, 'ess':30}]

        self.max_out_of_range = -1


        
class TestConstPopSize_MissingData(TestConstPopSize):

    def setUp(self):
        TestConstPopSize.setUp(self)

        # identical to TestConstPopSize, except for name, and for running inference without any data
        self.prefix = "testdata/constpopsize_missingdata"
        self.missing_leaves = [0,1]

        self.np = 100
        self.em = 1
        self.seed=(1,)
        
        # set targets
        self.targets = [{'type':"Recomb", 'truth':1e-8, 'min':0.98e-8, 'max':1.02e-8, 'ess': 45},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'truth':10000, 'min': 6000, 'max': 12000, 'ess': 12},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'truth':10000, 'min': 9800, 'max': 10200, 'ess': 45},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'truth':10000, 'min': 9800, 'max': 10200, 'ess':50},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'truth':10000, 'min': 9800, 'max': 10200, 'ess':50},
                        {'type':"Coal", 'pop':0, 'epoch':4, 'truth':10000, 'min': 9800, 'max': 10200, 'ess':50}]
        self.max_out_of_range = -1   # set to -1 to always fail (and keep intermediate files)

        


class TestConstPopSize_FourEpochs_MissingData(TestConstPopSize_FourEpochs):

    def setUp(self):
        TestConstPopSize_FourEpochs.setUp(self)

        # identical to TestConstPopSize, except for name, and for running inference without any data
        self.prefix = "testdata/constpopsize_4epochs_missingdata"
        self.missing_leaves = [0,1]

        self.np = 100
        self.em = 1
        self.bias_heights = [800]
        self.bias_strengths = [3,1]
        self.seed=(1,)
        
        # set targets
        self.targets = [{'type':"Recomb", 'truth':1e-8, 'min':0.99e-8, 'max':1.01e-8, 'ess': 30},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'truth':10000, 'min': 9000, 'max': 11000, 'ess': 6},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'truth':10000, 'min': 9400, 'max': 10600, 'ess': 20},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'truth':10000, 'min': 9800, 'max': 10200, 'ess':35},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'truth':10000, 'min': 9900, 'max': 10100, 'ess':40}]
        self.max_out_of_range = -1   # set to -1 to always fail (and keep intermediate files)


class TestConstPopSize_FourEpochs_EightSamples(TestConstPopSize_FourEpochs):

    def setUp(self, name = "testdata/constpopsize_4epochs_8samples"):
        TestConstPopSize_FourEpochs.setUp(self, name)

        # relatively short sequence, so expect biased inferences due to lack of data.  Instead, testing log likelihood
        # The estimated log likelihood is about -23917 with np=3000 and bias=2.0
        # The estimated log likelihood is about -23871 with np=6000 and bias=2.0
        # The true log likelihood may be around -23814 (np=12000 bias=2.0)
        self.pop.sequence_length = 1e6
        self.pop.num_samples =  8
        
        # set default inference parameters
        self.em = 1
        self.np = 3000
        self.bias_heights = [800]
        self.bias_strengths = [2,1]
        self.seed = (1,)
        self.debug = True

        # set targets
        self.targets = [{'type':"Recomb", 'min':0.85e-8, 'max':1.15e-8, 'truth':1e-8, 'ess':1.3},
                        {'type':"Coal", 'pop':0, 'epoch':0, 'min':4000,  'max':18000, 'truth':10000, 'ess':1},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':8500,  'max':11500, 'truth':10000, 'ess':1},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':9700,  'max':10300, 'truth':10000, 'ess':1.1},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':9600,  'max':10400, 'truth':10000, 'ess':1.7},
                        #{'type':"LogL", 'min':-23920, 'max':0, 'truth': -23814}]
                        {'type':"LogL", 'min':-26560, 'max':0, 'truth': -23814}]

        self.max_out_of_range = -1


class TestConstPopSize_Runtime(TestGeneric):

    def setUp(self, name="testdata/runtime"):
        TestGeneric.setUp(self, name)
        self.seqlen = 1e6
        self.pop = populationmodels.Population( sequence_length = self.seqlen,
                                                scrmpath=self.scrmpath,
                                                change_points = [0],
                                                num_populations = 1,
                                                num_samples = 8,
                                                population_sizes = [[1]])

        # use smcsmc front-end
        self.smcsmcpath = "../smcsmc"
        
        # use pattern for inference (default)
        #self.popt = None

        self.em = 0
        self.np = 1000
        self.bias_heights = [800]
        self.bias_strengths = [1,1]
        self.seed = (1,)
        self.debug = True
        self.targets = []
        self.max_out_of_range = -1
        
        
if __name__ == "__main__":
    unittest.main()

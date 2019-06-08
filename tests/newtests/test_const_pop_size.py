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
        self.em = 0
        self.np = 1000
        self.debug = True
        self.bias_heights = [400]
        self.bias_strengths = [3,1]
        self.tmax = 4
        self.seed = (1,)
        
        #self.guided_recomb_alpha = 0.5 #TEST

        # set targets
        self.targets = [{'type':"Coal", 'pop':0, 'epoch':0, 'min':0    , 'max':124573, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':10202, 'max':10471 , 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':9927 , 'max':10072 , 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':9963 , 'max':10066 , 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':4, 'min':9934 , 'max':10086 , 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':5, 'min':9862 , 'max':10098 , 'truth':10000},
                        {'type':"Recomb", 'min':9.77e-9, 'max':9.89e-9, 'truth':1e-8}]
        self.max_out_of_range = 0



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
        self.em = 0
        self.np = 1000
        self.bias_heights = [800]
        self.bias_strengths = [1.5,1]
        self.seed = (1,)
        self.debug = True

        # set targets
        self.targets = [{'type':"Coal", 'pop':0, 'epoch':0, 'min':4769 , 'max':15911, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':9164 , 'max':10081, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':9828 , 'max':9973 , 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':9807 , 'max':9896 , 'truth':10000},
                        {'type':"Recomb", 'min':9.78e-9, 'max':9.91e-9, 'truth':1e-8}]
        self.max_out_of_range = 0




        
class TestConstPopSize_MissingData(TestConstPopSize):

    def setUp(self):
        TestConstPopSize.setUp(self)

        # identical to TestConstPopSize, except for name, and for running inference without any data
        self.prefix = "testdata/constpopsize_missingdata"
        self.missing_leaves = [0,1]

        self.np = 100
        self.em = 0
        self.seed=(4,)
        
        # set targets
        self.targets = [{'type':"Coal", 'pop':0, 'epoch':0, 'min':5912 , 'max':14591, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':9879 , 'max':10116, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':9913 , 'max':10091, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':9930 , 'max':10070, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':4, 'min':9909 , 'max':10097, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':5, 'min':9885 , 'max':10105, 'truth':10000},
                        {'type':"Recomb", 'min':9.94e-9, 'max':1.01e-8, 'truth':1e-8}]
        self.max_out_of_range = 0   # set to -1 to always fail (and keep intermediate files)

        


class TestConstPopSize_FourEpochs_MissingData(TestConstPopSize_FourEpochs):

    def setUp(self):
        TestConstPopSize_FourEpochs.setUp(self)

        # identical to TestConstPopSize, except for name, and for running inference without any data
        self.prefix = "testdata/constpopsize_4epochs_missingdata"
        self.missing_leaves = [0,1]

        self.np = 100
        self.em = 0
        self.bias_heights = [800]
        self.bias_strengths = [3,1]
        self.seed=(4,)
        
        # set targets
        self.targets = [{'type':"Coal", 'pop':0, 'epoch':0, 'min':8303 , 'max':11882, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':9699 , 'max':10262, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':9910 , 'max':10089, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':9952 , 'max':10054, 'truth':10000},
                        {'type':"Recomb", 'min':9.91e-9, 'max':1.01e-8, 'truth':1e-8}]
        self.max_out_of_range = 0   # set to -1 to always fail (and keep intermediate files)


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
        self.em = 0
        self.np = 3000
        self.bias_heights = [800]
        self.bias_strengths = [2,1]
        self.seed = (2,)
        self.debug = True

        # set targets
        self.targets = [{'type':"Coal", 'pop':0, 'epoch':0, 'min':6132 , 'max':20605, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':9789 , 'max':12770, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':9434 , 'max':10716, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':9436 , 'max':10364, 'truth':10000},
                        {'type':"Recomb", 'min':9.99e-9, 'max':1.10e-8, 'truth':1e-8}]
        self.max_out_of_range = 0


class TestConstPopSize_FourEpochs_FalseStart(TestConstPopSize_FourEpochs):

    def setUp(self, name="testdata/constpopsize_4epochs_falsestart"):
        TestConstPopSize_FourEpochs.setUp(self, name)

        self.smcsmc_initial_pop_sizes = [[1.2], [1.2], [1.2], [1.2]]
        
        # set targets
        self.targets = [{'type':"Coal", 'pop':0, 'epoch':0, 'min':5382 , 'max':18567, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':10345, 'max':11398, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':11154, 'max':11298, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':11259, 'max':11354, 'truth':10000},
                        {'type':"Recomb", 'min':9.71e-9, 'max':9.84e-9, 'truth':1e-8}]
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
        self.em = 0
        self.np = 1000
        self.debug = True
        self.bias_heights = [800]
        self.bias_strengths = [2,1]
        self.tmax = 4
        self.seed = (1,)

        # set targets
        scaling = 1.0 / (4.0 * self.pop.N0)
        self.targets = [{'type':"Coal", 'pop':0, 'epoch':0, 'min':0    , 'max':1.97e8, 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':0, 'min':0    , 'max':2.14e8, 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':1, 'min':9211 , 'max':11671 , 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':1, 'min':9833 , 'max':12253 , 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':2, 'min':8920 , 'max':11205 , 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':2, 'min':8783 , 'max':12353 , 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':3, 'min':9130 , 'max':10790 , 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':3, 'min':8695 , 'max':11327 , 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':4, 'min':9084 , 'max':10583 , 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':4, 'min':8128 , 'max':12132 , 'truth':10000},
                        {'type':"Coal", 'pop':0, 'epoch':5, 'min':9460 , 'max':10360 , 'truth':10000},
                        {'type':"Coal", 'pop':1, 'epoch':5, 'min':8913 , 'max':11451 , 'truth':10000},
                        {'type':"Recomb", 'min':9.89e-9, 'max':1.06e-8, 'truth':1e-8},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':0, 'min':0      , 'max':0},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':0, 'min':0      , 'max':1.47e-4},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':1, 'min':0      , 'max':0},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':1, 'min':1.66e-5, 'max':3.21e-5},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':2, 'min':0      , 'max':0},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':2, 'min':2.00e-5, 'max':3.36e-5},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':3, 'min':0      , 'max':0},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':3, 'min':2.24e-5, 'max':3.28e-5},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':4, 'min':0      , 'max':0},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':4, 'min':2.32e-5, 'max':3.42e-5},
                        {'type':"Migr", 'from_pop':0, 'to_pop':1, 'epoch':5, 'min':0      , 'max':0},
                        {'type':"Migr", 'from_pop':1, 'to_pop':0, 'epoch':5, 'min':2.39e-5, 'max':3.00e-5} ]
        for t in self.targets:
            if t['type'] == "Migr":
                t['truth'] = scaling * self.pop.migration_rates[t['epoch']][t['from_pop']][t['to_pop']]
        self.max_out_of_range = 0
        


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

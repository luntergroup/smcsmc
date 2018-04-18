import itertools
import math
import sys

import experiment_base

###################################################################################
#
# experiment:
#  investigate effects of guiding, biasing, and number of particles
#  on ability to infer in population of 8 samples
#
###################################################################################




# parameters for this experiment.  Estimated to take 2.5 days for 100 chunks on the cluster
missingdata = False
inferences = range(12)
particles = [30000]
emiters = 30
seqlen_infpar = 2000e6
simseed = 1
chunks = 100  # don't forget to use -c -C "-P lunter.prjb -q long.qb"

mu = 1.25e-8  # from msmc paper (pg. 920)
rho = 3.5e-9  # from msmc paper (supp pg. 11):  zigzag simulation uses -t 7156 -r 2000 => mu/rho = 3.578
N0 = 14312    # from msmc paper (supp pg. 11):  7156 = 4 N0 mu L, L = 10e6

# name of this experiment
experiment_name = "zigzag_vb"

# class implementing the experiment outline
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize_FourEpochs


# main experiment - different clumped resampling strategies, keyed by smcseed  (see particle.cpp - hardcoded)
experiment_pars = [{'length':length, 'smcseed':smcseed, 'np':np, 'em':em, 'guide':guide, 'bias':bias, 'vb':vb }
                   for ( length, smcseed, np, em, guide, bias, vb) in itertools.product(
                           [seqlen_infpar],
                           inferences,
                           particles,
                           [emiters],
                           [0],
                           [2],
                           [True] ) ]


if missingdata:
    # for running without data;
    experiment_pars = [{'length':seqlen_infpar,
                        'smcseed':12345,
                        'np':10,
                        'em':0,
                        'guide':0,
                        'bias':1,
                        'vb':True}]
    inference_reps = 1



def schiffels_ne( t, n0=5 ):
    T = [0.000582262,0.00232905,0.00931619,0.0372648,0.149059,0.596236,1e99]
    e = [0,          1318.18,   -329.546,  82.3865,  -20.5966,5.14916, 0]
    ne = n0
    i = 0
    t0 = 0
    while t0 < t:
        t1 = min( T[i], t )
        ne = ne * math.exp( -e[i] * (t1 - t0) )
        t0 = t1
        i = i+1
    return ne


# run an experiment; keyword parameters agree with those in experiment_pars
def run_experiment( length, smcseed, np, em, guide, bias, vb ):
    simseed = smcseed
    if missingdata:
        missing_leaves = [0,1,2,3,4,5,6,7]
        label = "L{}_S{}_N{}_missing".format(int(length),smcseed,np)
    else:
        missing_leaves = []
        label = "L{}_S{}_N{}_G{}_B{}_VB{}".format(int(length),smcseed,np,guide,bias,vb)
    strpar = "guide{}_bias{}_vb{}".format(guide,bias,vb)
    if experiment_base.have_result( name=experiment_name,
                                    sequence_length=length,
                                    dataseed=simseed,
                                    infseed=smcseed,
                                    np=np,
                                    str_parameter=strpar):
        print "Skipping " + label
        return

    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( experiment_base.datapath + experiment_name )
    e.pop.num_samples = 8
    e.str_parameter = strpar

    # set simulation parameters
    rel_ne = 1         # Schiffel's paper uses 5
    #rel_ne = 5         # TEST: use Schiffel's original model
    e.pop.N0 = N0
    e.pop.mutation_rate = mu
    e.pop.recombination_rate = rho
    e.pop.change_points = [0, 0.596236]
    e.pop.population_sizes = [[rel_ne], [rel_ne/10.0]]
    e.missing_leaves = missing_leaves
    # use "migration commands" to set exponential growth parameters for zigzag model (see msmc paper, supp info, pg. 11)
    e.pop.migration_commands = ["-eG 0.000582262 1318.18 -eG 0.00232905 -329.546 -eG 0.00931619 82.3865 -eG 0.0372648 -20.5966 -eG 0.149059 5.14916",""]
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = experiment_base.scrmpath
    e.filename_disambiguator = label

    # set inference parameters
    e.popt = None
    nT = 256.0     # number of segments before joining; msmc paper, online methods "parameter inference"
    scale = 8
    e.smcsmc_change_points = [-math.log( 1-i/nT ) / scale for i in 
                               range(5) +         # up to ~0.0023; 1+4 epochs
                               range(6,18,2) +    # up to ~0.009;  6 epochs
                               range(18,58,5) +   # up to ~0.037;  8 epochs
                               range(58,162,13) + # up to ~0.149;  8 epochs
                               range(164,253,11)] # up to ~0.596 and beyond; 8+1 epochs
    e.smcsmc_change_points.append(0.693)

    # need to set pop sizes, migr rates and commands explicitly when change pts differ between simulation and inference:
    if missingdata:
        e.smcsmc_initial_pop_sizes = [ [ schiffels_ne(cp, rel_ne) ] for cp in e.smcsmc_change_points ]
    else:
        e.smcsmc_initial_pop_sizes = [ [rel_ne] for cp in e.smcsmc_change_points ]

    e.smcsmc_initial_migr_rates = [ [[0]] for cp in e.smcsmc_change_points ]
    e.smcsmc_migration_commands = [ None for cp in e.smcsmc_change_points ]

    ## testing the simulation
    if missingdata:
        e.pop.change_points = e.smcsmc_change_points          ## TEST
        e.pop.population_sizes = e.smcsmc_initial_pop_sizes   ## TEST
        e.pop.migration_commands = None                       ## TEST
    
    e.maxNE = 1e7   # was 5e4, but changed to show that EM converges to infty
    e.seed = (smcseed,)
    e.seqlen = length
    e.smcsmcpath = experiment_base.smcsmcpath
    e.infer_recombination = False
    e.np = np
    e.em = em
    e.vb = vb
    e.chunks = chunks
    e.submit_chunks = True
    e.aux_part_filt = 2       # always use apf=2
    e.delay = 0               # constpopsize_aux_part_filt experiments suggests delay never helps, and causes bias
    if bias > 0:
        e.bias_heights = [ e.smcsmc_change_points[i] * 4 * N0 for i in [1] ]  # 32
        e.bias_strengths = [bias,1]
    else:
        e.bias_heights = [ e.smcsmc_change_points[i] * 4 * N0 for i in [1,3,5] ] # 32, 96, 193
        e.bias_strengths = [1.0-(bias+1), 1-(bias+1)/2.0, 1-(bias+1)/10.0, 1]
    e.guided_recomb_alpha = guide
    e.m_step = True

    # perform inference and store results
    e.infer( case = smcseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label


def plot_experiment():
    import pandas as pd
    from rpy2.robjects import pandas2ri, r
    pandas2ri.activate()
    records = experiment_base.get_data_from_database( experiment_name )
    
    df = pd.DataFrame( data = records )
    r.source("zigzag_inference_guiding.R")
    g = 30.0
    rel_ne = 1  # must be synced with code above!
    truth = pd.DataFrame( data = [ (t*4*N0, N0*schiffels_ne(t, n0=rel_ne))
                                   for t in [ math.pow(10,j/100.0)/(g*4*N0)
                                              for j in range(250,650)]],
                          columns = ["t","Ne"] )
    r('plot.smcsmc')( df, truth, g )
        
if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.plot_experiment = plot_experiment
    experiment_base.main( experiment_name, experiment_pars )

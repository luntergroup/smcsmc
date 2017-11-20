import experiment_base
import itertools
import math

###################################################################################
#
# experiment:
#  investigate interaction between bias, delay, and aux particle filter, focusing
#  on the first epoch (<400 generations, 0.2% of coalescences)
# 
#  initial run (5/10/17) shows that bias=2 has benefits over vanilla; delay 0.5 has
#  limited benefits, and increasing bias to 5 with delay results in positive bias in
#  coalescent opportunity, and positive bias in recombination rate.  Using aux part
#  filt mode 2 (singletons and doubletons) improves likelihood, ess with minimal bias
#  in recombination rate inference.
#
#  next experiment: use fixed data to reduce variance in log likelihood; try delay
#  0.1 to see if that mitigates bias; try bias 10 to look at more extreme case.
#  More significant digits in LogL reporting.
#
###################################################################################


# parameters for this experiment
inference_reps = 10
seqlen = 50e6
particles = [1000, 3000, 10000, 30000]
emiters = 15
lagfactor = 2
nsam = 8
chunks = 12
submit_jobs = True          # if submit_jobs=False, don't forget to use -C "-P lunter.prjb -q long.qb -pe shmem <chunks>"
initial_ne_factor = 0.75
N0 = 10000

# name of this experiment
experiment_name = "constpopsize_auxpartfilt"


# class defining default population parameters
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize


# define the experiments
# HERE BE DRAGONS: 'delay' = -1 for unphased
experiment_pars = [{'length':seqlen,             # global
                    'simseed':0,                 # repetitions with the same data
                    'infseed':seed,
                    'numparticles':numparticles, 
                    'lag':lagfactor,             # global
                    'nsam':nsam,                 # global
                    'bias':bias,
                    'delay':delay,
                    'apf':apf}
                   for (seed, numparticles, (bias, delay, apf)) in itertools.chain(
                           itertools.product(
                               range(inference_reps),
                               particles,
                               [(1,0,0),
                                (1,0,2),
                                (1,-1,0),
                                (1,-1,2)]
                           )
                   )]

# run an experiment.  keyword parameters must match those in experiment_pars
def run_experiment( length, simseed, infseed, numparticles, lag, nsam, bias, delay, apf ):
    label = "P{}_N{}_S{}_B{}_D{}_APF{}".format(numparticles,nsam,simseed,bias,delay,apf)
    if experiment_base.have_result( name = experiment_name,
                                    infseed = infseed,
                                    np = numparticles,
                                    num_samples = nsam,
                                    bias_strengths = ' '.join(map(str,[bias,1])),
                                    int_parameter = int(delay),
                                    aux_part_filt = apf ):
        print "Skipping " + label
        return
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( experiment_base.datapath + experiment_name )

    # set simulation parameters
    e.pop.num_samples = nsam
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.change_points = [0, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.5, 1, 1.5]
    e.pop.population_sizes = [[1], [1], [1], [1], [1], [1], [1], [1], [1], [1]]
    e.pop.scrmpath = experiment_base.scrmpath
    e.filename_disambiguator = label

    # set inference parameters
    # need to set pop sizes, migr rates and commands explicitly when change pts
    # differ between simulation and inference:
    e.seqlen = length
    e.seed = (infseed,)
    e.smcsmc_change_points = [0, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.5, 1, 1.5]
    e.smcsmc_initial_pop_sizes = [ [initial_ne_factor * ne]
                                   for ne in [1,1,1,1,1,1,1,1,1,1] ]
    e.smcsmc_initial_migr_rates = [ [[0]] for cp in e.smcsmc_change_points ]
    e.smcsmc_migration_commands = [ None for cp in e.smcsmc_change_points ]
    e.smcsmcpath = experiment_base.smcsmcpath
    e.np = numparticles
    e.lag = lag
    e.em = emiters
    e.bias_heights = [400]
    e.bias_strengths = [bias, 1]
    e.delay = 0                    # I came off the idea of using delay, in favour of APF.  However for unphased data it may still have merit
    e.delay_type = "coal"
    e.int_parameter = int(delay)   # database does not store delay, so store in spare slot
    if delay == -1:
        e.phased = False
    e.aux_part_filt = apf
    e.chunks = chunks
    e.submit_chunks = submit_jobs  # submit chunks using qsub, rather than submitting toplevel script using qrsh & threading
    e.popt = None

    # perform inference and store results
    e.infer( case = infseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label
    return


def plot_experiment():
    import pandas as pd
    from rpy2.robjects import pandas2ri, r
    pandas2ri.activate()
    records = experiment_base.get_data_from_database( experiment_name )
    
    df = pd.DataFrame( data = records )
    r.source("constpopsize_auxpartfilt_experiment.R")
    g = 30.0
    truth = pd.DataFrame( data = [ (t*4*N0, N0*1)
                                   for t in [ math.pow(10,j/10.0)/(g*4*N0)
                                              for j in range(25,65)]],
                          columns = ["t","Ne"] )
    r('plot.smcsmc')( df, truth, g )
        
if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.plot_experiment = plot_experiment
    experiment_base.main( experiment_name, experiment_pars )


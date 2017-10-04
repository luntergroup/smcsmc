import experiment_base
import itertools

###################################################################################
#
# experiment:
#  estimate bias and variance, dependent on lag
#
###################################################################################


# parameters for this experiment
inference_reps = 5
seqlen = 50e6
particles = [3000, 10000, 30000]
emiters = 1
lagfactor = 2
nsam = 8
chunks = 4    # don't forget to use -C "-P lunter.prjb -q long.qb -pe shmem <chunks>"

# name of this experiment
experiment_name = "constpopsize_auxpartfilt"


# class defining default population parameters
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize


# define the experiments
experiment_pars = [{'length':seqlen,             # global
                    'simseed':simseed,
                    'infseed':simseed,           # repetitions with unique data
                    'numparticles':numparticles, 
                    'lag':lagfactor,             # global
                    'nsam':nsam,                 # global
                    'bias':bias,
                    'delay':delay,
                    'apf':apf}
                   for (simseed, numparticles, (bias, delay, apf)) in itertools.product(
                           range(inference_reps),
                           particles,
                           [(0,0,0),    # vanilla
                            (2,0,0),    # bias 2, no delay
                            (5,0,0),    # bias 5, no delay
                            (2,0.5,0),  # bias 2, delay
                            (5,0.5,0),  # bias 5, delay
                            (2,0.5,1),  # bias 2, delay, apf 1
                            (5,0.5,1),  # bias 5, delay, apf 1
                            (2,0,1),    # bias 2, apf 1
                            (5,0,1),    # bias 5, apf 1
                            (2,0,2),    # bias 2, apf 2
                            (5,0,2),    # bias 5, apf 2
                            (2,0,3),    # bias 2, apf 3
                            (5,0,3),    # bias 5, apf 3
                            (2,0,4),    # bias 2, apf 4
                            (5,0,4)     # bias 5, apf 4
                            ])]

# run an experiment.  keyword parameters must match those in experiment_pars
def run_experiment( length, simseed, infseed, numparticles, lag, nsam, bias, delay, apf ):
    label = "P{}_S{}_B{}_D{}_APF{}".format(numparticles,simseed,bias,delay,apf)
    if experiment_base.have_result( name = experiment_name,
                                    infseed = infseed,
                                    np = numparticles,
                                    bias_strengths = ' '.join(map(str,[bias,1])),
                                    int_parameter = int(2 * delay),
                                    aux_part_filt = apf ):
        print "Skipping " + label
        return
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( experiment_base.datapath + experiment_name )

    # set simulation parameters
    e.pop.num_samples = nsam
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = experiment_base.scrmpath
    e.filename_disambiguator = label
    e.pop.change_points = [0, 0.01, .25, 0.5, 1, 1.5]
    e.pop.population_sizes = [[1], [1], [1], [1], [1], [1]]

    # set inference parameters
    e.seqlen = length
    e.seed = (infseed,)
    e.np = numparticles
    e.lag = lag
    e.em = emiters-1
    e.bias_heights = [400]
    e.bias_strengths = [bias, 1]
    e.delay = delay
    e.int_parameter = int(2*delay)   # database does not store delay, so store in spare slot
    e.aux_part_filt = apf
    e.smcsmcpath = experiment_base.smcsmcpath
    e.chunks = chunks
    e.popt = None
    e.smcsmc_change_points = [0, 0.01, 0.25, 0.5, 1, 1.5]
    # need to set pop sizes, migr rates and commands explicitly when change pts
    # differ between simulation and inference:
    e.smcsmc_initial_pop_sizes = [ [1], [1], [1], [1], [1], [1] ]
    e.smcsmc_initial_migr_rates = [ [[0]] for cp in e.smcsmc_change_points ]
    e.smcsmc_migration_commands = [ None for cp in e.smcsmc_change_points ]

    # perform inference and store results
    e.infer( case = simseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label
    return


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )

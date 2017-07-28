import experiment_base
import itertools

###################################################################################
#
# experiment:
#  estimate bias and variance, dependent on lag
#
###################################################################################


# parameters for this experiment
inference_reps = 10
seqlen = 50e6
particles = 10000
emiters = 30
lagfactors = [0, 0.25, 0.5, 1.0, 2.0, 4.0]
nsam = [2,4,8]
chunks = 4    # don't forget to use -C "-P lunter.prjb -q long.qb -pe shmem <chunks>"

# name of this experiment
experiment_name = "bottleneck_bylag_falsestart"


# class defining default population parameters
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize_FourEpochs


# define the experiments
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'infseed':infseed,
                    'numparticles':numparticles, 'lag':lag, 'nsam':nsam}
                   for (seqlen, simseed, infseed, numparticles, lag, nsam) in (
                        # repetitions with unique data
                        [(seqlen, 100+rep, 100+rep, particles, lag, ns)
                         for rep, lag, ns in itertools.product(range(inference_reps), lagfactors, nsam)] )]


# run an experiment.  keyword parameters must match those in experiment_pars
def run_experiment( length, simseed, infseed, numparticles, lag, nsam ):
    label = "L{}_S{}_I{}_P{}_G{}_NSAM{}".format(int(length),simseed,infseed,numparticles, lag, nsam)
    if experiment_base.have_result( name = experiment_name,
                                    sequence_length = length,
                                    dataseed = simseed,
                                    infseed = infseed,
                                    np = numparticles,
                                    lag = lag,
                                    num_samples = nsam):
        print "Skipping " + label
        return
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( experiment_base.datapath + experiment_name )
    
    # set simulation parameters
    e.pop.change_points = [0, .01, 0.06, 0.2,  1, 2]
    e.pop.population_sizes = [[1], [0.1], [1], [0.5], [1], [2]]
    e.pop.num_samples = nsam
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = experiment_base.scrmpath
    e.filename_disambiguator = label
    
    # set inference parameters
    e.seqlen = length
    e.popt = None
    e.smcsmc_change_points = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
                              0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20,
                              0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.2, 1.4, 1.6, 1.8, 1.9, 2.0]
    # need to set pop sizes, migr rates and commands explicitly when change pts
    # differ between simulation and inference:
    e.smcsmc_initial_pop_sizes = [ [1], [1], [1], [1], [1], [1], [1], 
                                   [1], [1], [1], [1], [1], [1], [1],
                                   [1], [1], [1], [1], [1], [1], [1], [1],
                                   [1], [1], [1], [1], [1], [1] ]
    e.smcsmc_initial_migr_rates = [ [[0]] for cp in e.smcsmc_change_points ]
    e.smcsmc_migration_commands = [ None for cp in e.smcsmc_change_points ]
    e.seed = (infseed,)
    e.np = numparticles
    e.lag = lag
    e.em = emiters-1
    e.bias_heights = None
    e.smcsmcpath = experiment_base.smcsmcpath
    e.chunks = chunks
    
    # perform inference and store results
    e.infer( case = simseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label
    return


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )

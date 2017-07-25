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
particles = [1000,10000]
emiters = 10
lagfactors = [0, 0.25, 0.5, 1.0, 2.0, 4.0]
nsam = [2,4,8]
chunks = 4    # don't forget to use -C "-P lunter.prjb -q long.qb -pe shmem <chunks>"

# name of this experiment
experiment_name = "constpopsize_bylag_truestart"


# class defining default population parameters
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize_FourEpochs


# define the experiments
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'infseed':infseed,
                    'numparticles':numparticles, 'lag':lag, 'nsam':nsam}
                   for (seqlen, simseed, infseed, numparticles, lag, nsam) in (
                        # repetitions with unique data
                        [(seqlen, 100+rep, 100+rep, np, lag, ns)
                         for rep, np, lag, ns in itertools.product(range(inference_reps), particles, lagfactors, nsam)] )]


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
    e.pop.num_samples = nsam
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = experiment_base.scrmpath
    e.filename_disambiguator = label
    e.pop.change_points = [0, .02, .1, .5]
    e.pop.population_sizes = [[1], [1], [1], [1]]

    # set inference parameters
    e.seqlen = length
    e.seed = (infseed,)
    e.np = numparticles
    e.lag = lag
    e.em = emiters-1
    e.bias_heights = None
    e.smcsmcpath = experiment_base.smcsmcpath
    e.chunks = chunks
    e.popt = None
    e.smcsmc_change_points = [0, 0.02, .1, .5]
    # need to set pop sizes, migr rates and commands explicitly when change pts
    # differ between simulation and inference:
    e.smcsmc_initial_pop_sizes = [ [1], [1], [1], [1] ]
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


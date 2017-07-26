import experiment_base
import itertools

###################################################################################
#
# experiment:
#  uni directional migration, try inferring without a focused sampler
#
###################################################################################


# parameters for this experiment
inference_reps = 10
seqlen = 100e6
particles = [1000,10000]
emiters = 30
lagfactor = 1.0
nsam = [4,8]
chunks = 4    # don't forget to use -C "-P lunter.prjb -q long.qb -pe shmem <chunks>"
focus_heights = None
focus_strengths = None

# name of this experiment
experiment_name = "unidirmigr_nofocus"


# class defining default population parameters
import test_two_pops
experiment_class = test_two_pops.TestTwoPopsSplitUniDirMigr


# define the experiments
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'infseed':infseed,
                    'numparticles':numparticles, 'lag':lag, 'nsam':nsam}
                   for (seqlen, simseed, infseed, numparticles, lag, nsam) in (
                        # repetitions with unique data
                        [(seqlen, 100+rep, 100+rep, np, lagfactor, ns)
                         for rep, np, ns in itertools.product(range(inference_reps), particles, nsam)] )]


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
    if nsam == 4:
        e.pop.sample_populations = [1,1,2,2]
    e.pop.change_points = [0, .1, .5]
    #e.pop.population_sizes = [[1,1], [1,1], [1,1]]
    #e.pop.migration_rates = [ [ [0, 0.2],[0, 0] ],
    #                          [ [0, 0.2],[0, 0] ],
    #                          [ [0, 0]  ,[0, 0] ] ]

    # set inference parameters
    e.seqlen = length
    e.seed = (infseed,)
    e.np = numparticles
    e.lag = lag
    e.em = emiters-1
    e.bias_heights = focus_heights
    e.bias_strengths = focus_strengths
    e.smcsmcpath = experiment_base.smcsmcpath
    e.chunks = chunks
    e.popt = None
    #e.smcsmc_change_points = [0, .1, .5]
    ## I shouldn't have to specify this here as I want to use the values in newtests/
    ## but if I don't specify this here, it claims [0.0] is the only change point....
    #e.smcsmc_initial_pop_sizes = [ [1,1], [1,1], [1,1] ]
    #e.smcsmc_initial_migr_rates =  [[0,.2],[.2,0]],
    #                               [[0,.2],[.2,0]],
    #                               [[0,0] ,[0,0] ] ]
    #e.smcsmc_migration_commands = [ ????? ]

    # perform inference and store results
    e.infer( case = simseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label
    return


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )


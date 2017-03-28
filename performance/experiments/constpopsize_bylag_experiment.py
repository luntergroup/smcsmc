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
particles = 3000
emiters = 30
lagfactors = [0, 0.25, 0.5, 1.0, 2.0, 4.0]
nsam = [2,4,8]

# name of this experiment
experiment_name = "constpopsize_4epochs_falsestart_lagdependence"


# class defining default population parameters
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize_FourEpochs_FalseStart


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
    e.pop.num_samples = nsam
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = experiment_base.scrmpath
    e.filename_disambiguator = label
    # set inference parameters
    e.seqlen = length
    e.seed = (infseed,)
    e.np = numparticles
    e.lag = lag
    e.em = emiters-1
    e.bias_heights = None
    e.smcsmcpath = experiment_base.smcsmcpath
    # perform inference and store results
    e.infer( case = simseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    e.tearDown()
    print "Done " + label
    return


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )


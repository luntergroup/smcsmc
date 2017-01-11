import experiment_base
import itertools

###################################################################################
#
# experiment:
#  estimate bias and variance, dependent on lag
#
###################################################################################


# parameters for this experiment
inference_reps = 25
seqlen = 50e6
particles = 500
emiters = 30
lagfactors = [0.25, 0.5, 1.0, 2.0, 4.0]


# name of this experiment
experiment_name = "constpopsize_3epochs_lagdependence"


# class defining default population parameters
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize_ThreeEpochs


# define the experiments
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'infseed':infseed,
                    'numparticles':numparticles, 'lag':lag}
                   for (seqlen, simseed, infseed, numparticles, lag) in (
                        # repetitions with unique data
                        [(seqlen, 100+rep, 100+rep, particles, lag)
                         for rep, lag in itertools.product(range(inference_reps), lagfactors)] +
                        # repetitions with fixed data
                        [(seqlen, 100,     100+rep, particles, lag)
                         for rep, lag in itertools.product(range(1,inference_reps), lagfactors)])]


# run an experiment.  keyword parameters must match those in experiment_pars
def run_experiment( length, simseed, infseed, numparticles, lag ):
    label = "L{}_S{}_I{}_P{}_G{}".format(int(length),simseed,infseed,numparticles, lag)
    if experiment_base.have_result( name = experiment_name,
                                    sequence_length = length,
                                    dataseed = simseed,
                                    infseed = infseed,
                                    np = numparticles,
                                    lag = lag ):
        print "Skipping " + label
        return
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( experiment_base.datapath + experiment_name )
    # set simulation parameters
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


import experiment_base
import itertools

###################################################################################
#
# experiment:
#  investigate influence of number of particles on bias in inference once more,
#  but now using more EM iterations
#
###################################################################################


# parameters for this experiment
inference_reps = 15
seqlen = 50e6
particles = [500, 1000, 2000, 5000, 10000, 20000]
emiters = 30
lagfactor = 1.0


# name of this experiment
experiment_name = "constpopsize_3epochs_particles2"

# class
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize_ThreeEpochs


# define the experiments
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'infseed':infseed,
                    'numparticles':numparticles, 'lag':lagfactor}
                   for (simseed, infseed, numparticles) in (
                        # repetitions with unique data
                        [(100+rep, 100+rep, numparticles)
                         for rep, numparticles in itertools.product(range(inference_reps),
                                                                    particles)])]


# run a single experiment
def run_experiment( length, simseed, infseed, numparticles, lag ):
    label = "S{}_I{}_P{}".format(simseed,infseed,numparticles)
    if have_result( name = experiment_name,
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
    e.bias_heights = [400]
    e.bias_strengths = [10,1]
    e.smcsmcpath = experiment_base.smcsmcpath
    # perform inference and store results
    e.infer( case = simseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    e.tearDown()
    print "Done " + label


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )

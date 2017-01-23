import experiment_base
import itertools

###################################################################################
#
# experiment:
#  Test inference on non-contemporaneous models
#
###################################################################################


# parameters for this experiment
inference_reps = 5
seqlen = 50e6
particles = [500,1000,3000]
biases = [1.001, 1.5, 2, 3, 5, 10]
emiters = 5
lagfactor = 4.0
leaves = 8


# name of this experiment
experiment_name = "constpopsize_sampletimes"


# class
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize_ThreeEpochs


# define the experiments
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'infseed':infseed,
                    'numparticles':numparticles, 'bias': bias, 'missing': missing}
                   for (simseed, infseed, numparticles, bias, missing) in (
                           # repetitions without data; vary bias and number of particles
                           [(rep, rep, numparticles, bias, True)
                            for rep, numparticles, bias in itertools.product(range(inference_reps),
                                                                             particles,
                                                                             biases)])]


# run a single experiment.  keyword parameters must match those in experiment_pars
def run_experiment( length, simseed, infseed, numparticles, bias, missing ):
    
    label = "L{}_S{}_I{}_P{}_B{}_M{}".format(length,simseed,infseed,numparticles,bias,missing)
    missingness = list(range(leaves)) if missing else []
    if experiment_base.have_result( name = experiment_name,
                                    sequence_length = length,
                                    dataseed = simseed,
                                    infseed = infseed,
                                    np = numparticles,
                                    int_parameter = int(bias * 1000),
                                    missing_leaves = str(missingness) ):
        print "Skipping " + label
        return
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( experiment_base.datapath + experiment_name )

    # set simulation parameters
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.change_points = [0, 0.01, 0.05, .1]
    e.pop.population_sizes = [1, 1, 1, 1]
    e.pop.num_samples = leaves
    e.pop.sample_times = ([0] * (leaves/2)) + ([0.01] * (leaves/2))
    e.pop.populations = [1] * leaves
    e.pop.scrmpath = experiment_base.scrmpath
    e.filename_disambiguator = label
    e.missing_leaves = missingness
    e.seqlen = length

    # set inference parameters
    e.seed = (infseed,)
    e.np = numparticles
    e.em = emiters-1
    e.smcsmc_change_times = [0, 0.01, 0.05, .1]  # for now initial pop sizes are hardcoded as 1
    e.smcsmcpath = experiment_base.smcsmcpath
    e.bias_heights = [400]
    e.bias_strenghts = [bias, 1]
    # perform inference and store results
    e.infer( case = simseed )
    e.int_parameter = int(bias * 1000)
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    e.tearDown()
    print "Done " + label


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )

import experiment_base
import itertools

###################################################################################
#
# experiment:
#  estimate bias and variance, dependent on length of input sequence
#  also look at contribution of stochasticity of coalescent itself
#
###################################################################################


# parameters for this experiment
inference_reps = 10
seqlens = [1e6, 2e6, 5e6, 10e6, 20e6, 50e6, 100e6]


# name of this experiment
experiment_name = "constpopsize_3epochs_lengthdependence"

# class
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize_ThreeEpochs

# define the experiments
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'smcseed':smcseed,'missing':missing}
                   for (seqlen, simseed, smcseed, missing) in (
                           # repetitions with unique data
                           [(seqln, 1+rep, 1+rep, False)
                            for seqln, rep in
                            itertools.product(seqlens, range(inference_reps))] +
                           # repetitions with fixed data (avoiding duplications)
                           [(seqln, 1,     1+rep, False)
                            for seqln, rep in
                            itertools.product(seqlens, range(1, inference_reps))] +
                           # repetitions with missing data.  Use seed to separate filenames
                           [(seqln, 1000+rep,  1000+rep, True)
                            for seqln, rep in
                            itertools.product(seqlens, range(inference_reps))])]


# run an experiment; keyword parameters agree with those in experiment_pars
def run_experiment( length, simseed, smcseed, missing ):
    missing_leaves = str( [0,1] if missing else [] )
    label = "L{}_S{}_I{}_M{}".format(int(length),simseed,smcseed,missing)
    if experiment_base.have_result( name=experiment_name,
                                    sequence_length=length,
                                    dataseed=simseed,
                                    infseed=smcseed,
                                    missing_leaves=missing_leaves ):
        print "Skipping " + label
        return
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( experiment_base.datapath + experiment_name )
    # set simulation parameters
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = experiment_base.scrmpath
    if missing: e.missing_leaves = [0,1]
    e.filename_disambiguator = label
    # set inference parameters
    e.seqlen = length
    e.seed = (smcseed,)
    e.bias_heights = None
    e.smcsmcpath = experiment_base.smcsmcpath
    # perform inference and store results
    e.infer( case = smcseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )

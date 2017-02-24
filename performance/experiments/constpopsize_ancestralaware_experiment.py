import experiment_base
import itertools

###################################################################################
#
# experiment:
#  estimate bias and variance, dependent on known or unknown ancestral information
#
###################################################################################


# parameters for this experiment
inference_reps = 25
seqlen = 50e6
emiters = 30
particles = [100,500,1000]
anc_aware = [True, False]


# name of this experiment
experiment_name = "constpopsize_4epochs_ancestralawaredependence"


# class defining default population parameters
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize_FourEpochs


# define the experiments
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'infseed':infseed, 'numparticles':numparticles, 'ancestral_aware':anc_aware}
                   for (seqlen, simseed, infseed, numparticles, anc_aware) in (
                           # repetitions with unique data
                           [(seqlen, 100+rep, 100+rep, np, aa)
                            for rep, np, aa in itertools.product(range(inference_reps), particles, anc_aware)] +
                           # repetitions with fixed data
                           [(seqlen, 100,     100+rep, np, aa)
                            for rep, np, aa in itertools.product(range(1,inference_reps), particles, anc_aware)])]


# run an experiment.  keyword parameters must match those in experiment_pars
def run_experiment( length, simseed, infseed, numparticles, ancestral_aware ):
    label = "L{}_S{}_I{}_P{}_AA{}".format(int(length),simseed,infseed,numparticles,ancestral_aware)
    if experiment_base.have_result( name=experiment_name,
                                    sequence_length=length,
                                    dataseed=simseed,
                                    infseed=infseed,
                                    np=numparticles,
                                    ancestral_aware=ancestral_aware ): # should I create a column for ancestral aware? (ALTER TABLE table_name; ADD column_name datatype; UPDATE table_name SET newcolumn = <something>)
        print "Skipping " + label
        return
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( experiment_base.datapath + experiment_name )
    # set simulation parameters
    e.pop.num_samples = 8
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = experiment_base.scrmpath
    e.filename_disambiguator = label
    # set inference parameters
    e.seqlen = length
    e.seed = (infseed,)
    e.np = numparticles
    e.bias_heights = None
    e.smcsmcpath = experiment_base.smcsmcpath
    e.ancestral_aware = ancestral_aware
    # perform inference and store results
    e.infer( case = simseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label
    return


    
if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )

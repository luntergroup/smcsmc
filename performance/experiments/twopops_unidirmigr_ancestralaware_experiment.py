import experiment_base
import itertools

###################################################################################
#
# experiment:
#  estimate bias and variance, dependent on known or unknown ancestral information
#
###################################################################################


# parameters for this experiment
inference_reps = 10
seqlen = 50e6
emiters = 4
particles = [100,500,1000,5000]
anc_aware = [True, False]
nsam = [4,8]


# name of this experiment
experiment_name = "twopops_unidirmigr_falsestart_ancestralawaredependence"


# class defining default population parameters
import test_two_pops
experiment_class = test_two_pops.TestTwoPopsUniDirMigr


# define the experiments
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'infseed':infseed, 'numparticles':numparticles,
                    'ancestral_aware':anc_aware, 'nsam':nsam}
                   for (seqlen, simseed, infseed, numparticles, anc_aware, nsam) in (
                           # repetitions with unique data
                           [(seqlen, 100+rep, 100+rep, np, aa, ns)
                            for rep, np, aa, ns in itertools.product(range(inference_reps), particles, anc_aware, nsam)] )]


# run an experiment.  keyword parameters must match those in experiment_pars
def run_experiment( length, simseed, infseed, numparticles, ancestral_aware, nsam ):
    label = "L{}_S{}_I{}_P{}_AA{}_NSAM{}".format(int(length),simseed,infseed,numparticles,ancestral_aware,nsam)
    if experiment_base.have_result( name=experiment_name,
                                    sequence_length=length,
                                    dataseed=simseed,
                                    infseed=infseed,
                                    np=numparticles,
                                    ancestral_aware=ancestral_aware,
                                    num_samples=nsam ):
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
    if nsam==4:
        # need to change from default population assignment which assumes 8 samples
        e.pop.sample_populations = [1,1,2,2]
    # set inference parameters
    e.seqlen = length
    e.em = emiters
    e.seed = (infseed,)
    e.np = numparticles
    e.bias_heights = None
    e.lag = 1
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

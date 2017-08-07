import experiment_base
import itertools

###################################################################################
#
# experiment:
#  split no migration, try inferring with a focused sampler, with or without delay
#
###################################################################################


# parameters for this experiment
inference_reps = 10
seqlen = 100e6
particles = [10000]
emiters = 30
lagfactor = 1.0
nsam = [8]
chunks = 4    # don't forget to use -C "-P lunter.prjb -q long.qb -pe shmem <chunks>"
focus_heights = [[400]]
focus_strengths = [[2,1],[3,1]]
delays = [0,0.5]

# name of this experiment
experiment_name = "splitnomigr_withfocus"


# class defining default population parameters
import test_two_pops
experiment_class = test_two_pops.TestTwoPopsSplitNoMigr


# define the experiments
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'infseed':infseed,
                    'numparticles':numparticles, 'lag':lag, 'nsam':nsam,
                    'biasheights':biasheights, 'biasstrengths':biasstrengths, 'delay':delay}
                   for (seqlen, simseed, infseed, numparticles, lag, nsam, biasheights, biasstrengths, delay) in (
                        # repetitions with unique data
                        [(seqlen, 100+rep, 100+rep, np, lagfactor, ns, bh, bs, d)
                         for rep, np, ns, bh, bs, d in itertools.product(range(inference_reps), particles, nsam, focus_heights, focus_strengths, delays)] )]


# run an experiment.  keyword parameters must match those in experiment_pars
def run_experiment( length, simseed, infseed, numparticles, lag, nsam, biasheights, biasstrengths, delay ):
    label = "L{}_S{}_I{}_P{}_NSAM{}_BH{}_BS{}_D{}".format(int(length),simseed,infseed,numparticles, nsam, biasheights, biasstrengths, delay)
    label = label.replace(" ","")
    if experiment_base.have_result( name = experiment_name,
                                    sequence_length = length,
                                    dataseed = simseed,
                                    infseed = infseed,
                                    np = numparticles,
                                    num_samples = nsam,
                                    bias_heights = str(biasheights),
                                    bias_strengths = str(biasstrengths),
                                    str_parameter = str(delay)):
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

    # set inference parameters
    e.seqlen = length
    e.seed = (infseed,)
    e.np = numparticles
    e.lag = lag
    e.em = emiters-1
    e.bias_heights = biasheights
    e.bias_strengths = biasstrengths
    e.smcsmcpath = experiment_base.smcsmcpath
    e.chunks = chunks
    e.popt = None
    e.delay = delay
    e.delay_type = "coal"

    # perform inference and store results
    e.infer( case = simseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label
    return


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )


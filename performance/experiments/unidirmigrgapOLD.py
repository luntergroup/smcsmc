import experiment_base
import itertools

###################################################################################
#
# experiment:
#  uni directional migration with gap OLD (parameters of model the same as the first experiments)
#  try inferring with many focusing choices
#
#
###################################################################################


# parameters for this experiment
inference_reps = 10
seqlen = 100e6
particles = [3000,10000]
emiters = 30
lagfactor = 1.0
nsam = [8]
chunks = 4    # don't forget to use -C "-P lunter.prjb -q long.qb -pe shmem <chunks>"
focus_heights = [[400],[800]]
focus_strengths = [[2,1],[3,1],[1,1],[2,.9999],[3,.9999]]
delays = [0.5]

# name of this experiment
experiment_name = "unidirmigrgapOLD"


# class defining default population parameters
import test_two_pops
experiment_class = test_two_pops.TestTwoPopsSplitUniDirMigrGapOLD


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
    e.pop.change_points = [0, .058, .116, .4]

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
    e.smcsmc_change_points = [0, 0.01, 0.03, 0.058,
                              0.08, 0.116,
                              0.17, 0.24, 0.31, 0.4,
                              0.6, 0.9, 1.2, 1.5, 2.0]
    # need to set pop sizes, migr rates and commands explicitly when change pts
    # differ between simulation and inference:
    e.smcsmc_initial_pop_sizes = [ [1,1] for cp in e.smcsmc_change_points ]
    # the true rate in this model is 16. This is huge, but I think I got it from an example in MSMC...
    e.smcsmc_initial_migr_rates = [ [[0,.1],[.1,0]],
                                    [[0,.1],[.1,0]],
                                    [[0,.1],[.1,0]],
                               	    [[0,.1],[.1,0]],
                                    [[0,.1],[.1,0]],
                               	    [[0,.1],[.1,0]],
                                    [[0,.1],[.1,0]],
                               	    [[0,.1],[.1,0]],
                                    [[0,.1],[.1,0]],
                                    [[0,0] ,[0,0] ],
                                    [[0,0] ,[0,0] ],
                                    [[0,0] ,[0,0] ],
                                    [[0,0] ,[0,0] ],
                                    [[0,0] ,[0,0] ],
                                    [[0,0] ,[0,0] ] ]
    e.smcsmc_migration_commands = [ None,
                                    None,
       	       	       	       	    None,
       	       	       	       	    None,
       	       	       	       	    None,
       	       	       	       	    None,
       	       	       	       	    None,
       	       	       	       	    None,
       	       	       	       	    None,
                                    "-ej .4 2 1",
       	       	       	       	    None,
       	       	       	       	    None,
       	       	       	       	    None,
       	       	       	       	    None,
                                    None ]


    # perform inference and store results
    e.infer( case = simseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label
    return


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )


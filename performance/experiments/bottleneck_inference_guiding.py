import experiment_base
import itertools

###################################################################################
#
# experiment:
#  investigate effects of guiding, biasing, and number of particles
#  on ability to infer in population of 8 samples
#
###################################################################################


# parameters for this experiment
inference_reps = 10
particles = [500, 1000, 2000, 5000, 10000, 20000, 40000]
particles2 = [500, 1000, 2000, 5000, 10000]
emiters = 20
seqlen_infpar = 50e6
simseed = 1
chunks = 4    # don't forget to use -C "-P lunter.prjb -q long.qb -pe shmem <chunks>"

# name of this experiment
experiment_name = "bottleneck_inference_guiding"

# class
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize_FourEpochs

# define the experiments
experiment_pars = [{'length':length, 'smcseed':smcseed, 'np':np, 'em':em, 'guide':guide, 'bias':bias, 'mstep':mstep }
                   for ( length, smcseed, simseed, np, em, guide, bias, mstep) in (
                           # unguided, unbiased, various numbers of particles, inferring parameters; variable data
                           [(seqlen_infpar, 1+rep, 1+rep, np, emiters, 0.0, 1.0, True)    for np, rep in itertools.product(particles, range(inference_reps))] +
                           # unguided, biased, various numbers of particles, inferring parameters; variable data
                           [(seqlen_infpar, 1+rep, 1+rep, np, emiters, 0.0, 2.0, True)    for np, rep in itertools.product(particles, range(inference_reps))] +
                           # guided, unbiased, smaller range of particles, inferring parameters; variable data
                           [(seqlen_infpar, 1+rep, 1+rep, np, emiters, 0.5, 1.0, True)    for np, rep in itertools.product(particles2, range(inference_reps))] +
                           # guided, biased, smaller range of particles, inferring parameters; variable data
                           [(seqlen_infpar, 1+rep, 1+rep, np, emiters, 0.5, 2.0, True)    for np, rep in itertools.product(particles2, range(inference_reps))])]




# run an experiment; keyword parameters agree with those in experiment_pars
def run_experiment( length, smcseed, np, em, guide, bias, mstep ):
    missing_leaves = []
    label = "L{}_S{}_N{}_G{}_B{}_M{}".format(int(length),smcseed,np,guide,bias,mstep)
    strpar = "guide{}_bias{}_mstep{}".format(guide,bias,mstep)
    if experiment_base.have_result( name=experiment_name,
                                    sequence_length=length,
                                    dataseed=simseed,
                                    infseed=smcseed,
                                    np=np,
                                    str_parameter=strpar):
        print "Skipping " + label
        return
    
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( experiment_base.datapath + experiment_name )
    e.pop.num_samples = 8
    e.str_parameter = strpar

    # set simulation parameters
    e.pop.change_points = [0, .01, 0.06, 0.2,  1, 2]
    e.pop.population_sizes = [[1], [0.1], [1], [0.5], [1], [2]]
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = experiment_base.scrmpath
    e.filename_disambiguator = label

    # set inference parameters
    e.popt = None
    e.smcsmc_change_points = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
                              0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20,
                              0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                              1.2, 1.4, 1.6, 1.8, 1.9, 2.0]
    # need to set pop sizes, migr rates and commands explicitly when change pts
    # differ between simulation and inference:
    e.smcsmc_initial_pop_sizes = [ [1] for cp in e.smcsmc_change_points ]
    e.smcsmc_initial_migr_rates = [ [[0]] for cp in e.smcsmc_change_points ]
    e.smcsmc_migration_commands = [ None for cp in e.smcsmc_change_points ]
    e.seed = (smcseed,)
    e.seqlen = length
    e.smcsmcpath = experiment_base.smcsmcpath
    e.np = np
    e.em = em
    e.maxNE = 1e5
    e.bias_heights = [800]
    e.bias_strengths = [bias,1]
    e.guided_recomb_alpha = guide
    e.m_step = mstep
    e.chunks = chunks

    # testing
    #e.np = 1
    #e.em = 0
    #e.missing_leaves = range(8)

    # perform inference and store results
    e.infer( case = smcseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )

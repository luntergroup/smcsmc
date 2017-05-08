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
emiters = 5
seqlen = 1e6
simseed = 1


# name of this experiment
experiment_name = "constpopsize_4epochs_guiding"

# class
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize_FourEpochs

# define the experiments
experiment_pars = [{'length':seqlen, 'smcseed':smcseed, 'np':np, 'em':em, 'guide':guide, 'bias':bias, 'mstep':mstep }
                   for ( length, smcseed, simseed, np, em, guide, bias, mstep) in (
                           # unguided, unbiased, various numbers of particles
                           [(seqlen, 1+rep, 1, np, 1, 0.0, 1.0, False)         for np, rep in itertools.product(particles, range(inference_reps))] +
                           # unguided, biased, various numbers of particles
                           [(seqlen, 1+rep, 1, np, 1, 0.0, 2.0, False)         for np, rep in itertools.product(particles, range(inference_reps))] +
                           # guided, unbiased, smaller range of particles
                           [(seqlen, 1+rep, 1, np, emiters, 0.5, 1.0, False)   for np, rep in itertools.product(particles2, range(inference_reps))] +
                           # guided, biased, smaller range of particles
                           [(seqlen, 1+rep, 1, np, emiters, 0.5, 2.0, False)   for np, rep in itertools.product(particles2, range(inference_reps))] +
                           # unguided, unbiased, various numbers of particles, inferring parameters
                           [(seqlen, 1+rep, 1, np, emiters, 0.0, 1.0, True)    for np, rep in itertools.product(particles, range(inference_reps))] +
                           # unguided, biased, various numbers of particles, inferring parameters
                           [(seqlen, 1+rep, 1, np, emiters, 0.0, 2.0, True)    for np, rep in itertools.product(particles, range(inference_reps))] +
                           # guided, unbiased, smaller range of particles, inferring parameters
                           [(seqlen, 1+rep, 1, np, emiters, 0.5, 1.0, True)    for np, rep in itertools.product(particles2, range(inference_reps))] +
                           # guided, biased, smaller range of particles, inferring parameters
                           [(seqlen, 1+rep, 1, np, emiters, 0.5, 2.0, True)    for np, rep in itertools.product(particles2, range(inference_reps))])]




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
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = experiment_base.scrmpath
    e.filename_disambiguator = label

    # set inference parameters
    e.seed = (smcseed,)
    e.seqlen = length
    e.smcsmcpath = experiment_base.smcsmcpath
    e.np = np
    e.em = em
    e.bias_heights = [800]
    e.bias_strengths = [bias,1]
    e.alpha = guide
    e.m_step = mstep

    # perform inference and store results
    e.infer( case = smcseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )

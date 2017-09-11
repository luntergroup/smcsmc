import experiment_base
import itertools
import math

###################################################################################
#
# experiment:
#  zigzag model, {naive, FSDR} x {guide=.25, guide=.5}
#
###################################################################################


# parameters for this experiment
inference_reps = 10
seqlen = 100e6
particles = 10000
emiters = 30
lagfactors = [1.0]
nsam = [8]
chunks = 4    # don't forget to use -C "-P lunter.prjb -q long.qb -pe shmem <chunks>"
focus_heights=[[800]]
focus_strengths=[[1,1],[3,0.9999]]
delays=[0.5]
guides=[0.25,0.5]

# name of this experiment
experiment_name = "zigzag_focusing_with_guiding"

mu = 1.25e-8  # from msmc paper (pg. 920)
rho = 3.5e-9  # from msmc paper (supp pg. 11):  zigzag simulation uses -t 7156 -r 2000 => mu/rho = 3.578
N0 = 14312    # from msmc paper (supp pg. 11):  7156 = 4 N0 mu L, L = 10e6

# class defining default population parameters
import test_const_pop_size
experiment_class = test_const_pop_size.TestConstPopSize_FourEpochs

# define the experiments
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'infseed':infseed,
                    'numparticles':numparticles, 'lag':lag, 'nsam':nsam,
                    'biasheights':biasheights, 'biasstrengths':biasstrengths, 'guide':guide}
                   for (seqlen, simseed, infseed, numparticles, lag, nsam, biasheights, biasstrengths, guide) in (
                        # repetitions with unique data
                        [(seqlen, 100+rep, 100+rep, particles, lag, ns, bh, bs, guide)
                         for rep, lag, ns, bh, bs, guide in itertools.product(range(inference_reps), lagfactors, nsam, focus_heights, focus_strengths, guides)] )]


# run an experiment.  keyword parameters must match those in experiment_pars
def run_experiment( length, simseed, infseed, numparticles, lag, nsam, biasheights, biasstrengths, guide ):
    label = "L{}_S{}_I{}_P{}_G{}_NSAM{}_BH{}_BS{}_GUIDE{}".format(int(length),simseed,infseed,numparticles, lag, nsam, biasheights, biasstrengths, guide)
    label = label.replace(" ","")
    if experiment_base.have_result( name = experiment_name,
                                    sequence_length = length,
                                    dataseed = simseed,
                                    infseed = infseed,
                                    np = numparticles,
                                    lag = lag,
                                    num_samples = nsam,
                                    bias_heights = str(biasheights),
                                    bias_strengths = str(biasstrengths),
                                    guided_recomb_alpha = guide ):
        print "Skipping " + label
        return
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( experiment_base.datapath + experiment_name )
    
    # set simulation parameters
    e.pop.N0 = N0
    e.pop.mutation_rate = mu
    e.pop.recombination_rate = rho
    e.pop.change_points = [0, 0.596236]
    e.pop.population_sizes = [[5], [0.5]]
    # use "migration commands" to set exponential growth parameters for zigzag model (see msmc paper, supp info, pg. 11)
    e.pop.migration_commands = ["-eG 0.000582262 1318.18 -eG 0.00232905 -329.546 -eG 0.00931619 82.3865 -eG 0.0372648 -20.5966 -eG 0.149059 5.14916",""]
    e.pop.num_samples = nsam
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = experiment_base.scrmpath
    e.filename_disambiguator = label
    
    # set inference parameters
    e.popt = None
    nT = 256.0     # number of segments before joining; msmc paper, online methods "parameter inference"
    scale = 7
    e.smcsmc_change_points = [-math.log( 1-i/nT ) / scale for i in 
                               range(5) +         # up to ~0.0023; 1+4 epochs
                               range(10,18,2) +   # up to ~0.009;  4 epochs
                               range(18,58,5) +   # up to ~0.037;  8 epochs
                               range(58,162,13) + # up to ~0.149;  8 epochs
                               range(162,256,10)] # up to ~0.596 and beyond; 9+1 epochs
    # need to set pop sizes, migr rates and commands explicitly when change pts
    # differ between simulation and inference:
    e.smcsmc_initial_pop_sizes = [ [1] for cp in e.smcsmc_change_points ]
    e.smcsmc_initial_migr_rates = [ [[0]] for cp in e.smcsmc_change_points ]
    e.smcsmc_migration_commands = [ None for cp in e.smcsmc_change_points ]
    e.maxNE = 1e5

    e.seqlen = length
    e.seed = (infseed,)
    e.np = numparticles
    e.lag = lag
    e.em = emiters-1
    e.bias_heights = biasheights
    e.bias_strengths = biasstrengths
    e.smcsmcpath = experiment_base.smcsmcpath
    e.chunks = chunks
    e.delay_type = "recomb"
    e.guided_recomb_alpha = guide
    
    # perform inference and store results
    e.infer( case = simseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label
    return


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )

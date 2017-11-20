import experiment_base
import itertools
import math

###################################################################################
#
# experiment:
#  uni directional migration, test with apf and FSDR
#
###################################################################################


# parameters for this experiment
inference_reps = range(5)
focus_heights_strengths_delays = [ (None, None, 0) ] # , ([2000],[3,.9999],0.5) ]
aux_part_filt = [2] # [0,2]
seqlen = 1000e6
particles = [30000]
emiters = 30
lagfactor = 1.0
nsam = 8
chunks = 50   # don't forget to use -C "-P lunter.prjb -q long.qb -pe shmem <chunks>"

mu = 1.25e-8  # from msmc paper (pg. 920)
rho = 3.5e-9  # from msmc paper (supp pg. 11):  zigzag simulation uses -t 7156 -r 2000 => mu/rho = 3.578
N0 = 14312    # from msmc paper (supp pg. 11):  7156 = 4 N0 mu L, L = 10e6
g = 30.0      # from msmc paper
migr_time = 50e3
split_time = 100e3

# name of this experiment
experiment_name = "unidirmigr_apf"


# class defining default population parameters
import test_two_pops
experiment_class = test_two_pops.TestTwoPops


# define the experiments
experiment_pars = [{'length':seqlen,
                    'simseed':simseed,
                    'infseed':infseed,
                    'numparticles':numparticles,
                    'lag':lag,
                    'nsam':nsam,
                    'biasheights':biasheights,
                    'biasstrengths':biasstrengths,
                    'delay':delay,
                    'apf':apf}
                   for (seqlen, simseed, infseed, numparticles, lag, nsam, biasheights, biasstrengths, delay, apf) in (
                           # repetitions with unique data
                           [(seqlen, 100+rep, 100+rep, np, lagfactor, nsam, bh, bs, d, apf)
                            for rep, np, (bh, bs, d), apf in
                            itertools.product(inference_reps,
                                              particles,
                                              focus_heights_strengths_delays,
                                              aux_part_filt
                                              )
                            ]
                           )
                   ]


def quantile( k, n, M=4 ):
    if k == 0: return 0.0
    return -math.log(1-k/(n+0.0)) / (M * (M-1)/2.0)


# run an experiment.  keyword parameters must match those in experiment_pars
def run_experiment( length, simseed, infseed, numparticles, lag, nsam, biasheights, biasstrengths, delay, apf ):
    label = "L{}_S{}_I{}_P{}_NSAM{}_BH{}_BS{}_D{}_APF{}".format(int(length),simseed,infseed,numparticles, nsam, biasheights, biasstrengths, delay, apf)
    label = label.replace(" ","").replace("[","").replace("]","")
    if experiment_base.have_result( name = experiment_name,
                                    sequence_length = length,
                                    dataseed = simseed,
                                    infseed = infseed,
                                    np = numparticles,
                                    num_samples = nsam,
                                    bias_heights = str(biasheights),
                                    bias_strengths = str(biasstrengths),
                                    aux_part_filt = apf,
                                    str_parameter = str(delay)):
        print "Skipping " + label
        return
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( experiment_base.datapath + experiment_name )

    # set simulation parameters
    e.pop.N0 = N0
    e.pop.mutation_rate = mu
    e.pop.recombination_rate = rho
    e.pop.num_samples = nsam
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = experiment_base.scrmpath
    e.filename_disambiguator = label
    e.pop.sample_populations = [1]*(nsam/2) + [2]*(nsam/2)
    e.pop.change_points = [ 0, migr_time / (4*N0*g), split_time / (4*N0*g) ]
    e.pop.migration_rates = [ [[0,0],[0,0]],
                              [[0,1],[0,0]],
                              [[0,0],[0,0]] ]
    e.pop.migration_commands = [ None, None, "-ej {} 2 1".format( split_time / (4*N0*g) ) ]

    quantiles = 30
    e.smcsmc_change_points = [ quantile(k,quantiles) for k in range(quantiles) ]
    # need to set pop sizes, migr rates and commands explicitly when change pts
    # differ between simulation and inference:
    e.smcsmc_initial_pop_sizes = [ [1,1] for cp in e.smcsmc_change_points ]
    e.smcsmc_initial_migr_rates = [[[0,.5],[.5,0]]] * 10 + [[[0,0],[0,0]]] * (quantiles-10)
    e.smcsmc_migration_commands = [ None ] * 10 + [ "-ej {} 2 1".format( quantile(10, quantiles) ) ] + [None] * (quantiles-11)
    
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
    e.submit_chunks = True
    e.popt = None
    e.delay = delay
    e.delay_type = "coal"
    e.aux_part_filt = apf

    # perform inference and store results
    e.infer( case = simseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label
    return


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.main( experiment_name, experiment_pars )


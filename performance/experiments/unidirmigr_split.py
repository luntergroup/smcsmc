import experiment_base
import itertools
import math

#########################################################################################
#
# experiment:
#  uni directional migration, test placement of population split
#  modified model to be more asymmetric, and have continuous population size across split
#
#########################################################################################


# parameters for this experiment
inference_reps = range(100,115)
split_epochs = [12] # [11,12,13]
seqlen = 1000e6
particles = [30000]
models = [0,1,2,3]
emiters = 30
lagfactor = 1.0
nsam = 8
chunks = 20   # don't forget to use -C "-P lunter.prjb -q long.qb -pe shmem <chunks>"

mu = 1.25e-8  # from msmc paper (pg. 920)
rho = 3.5e-9  # from msmc paper (supp pg. 11):  zigzag simulation uses -t 7156 -r 2000 => mu/rho = 3.578
N0 = 10000    # msmc paper has 14312 (supp pg. 11):  7156 = 4 N0 mu L, L = 10e6
g = 30.0      # from msmc paper
migr_time = 50e3
split_time = 100e3

# name of this experiment
experiment_name = "unidirmigr_split"


# class defining default population parameters
import test_two_pops
experiment_class = test_two_pops.TestTwoPops

def uniq_chain(*args, **kwargs):
    seen = set()
    for x in itertools.chain(*args, **kwargs):
        if x in seen: continue
        seen.add(x)
        yield x


# define the experiments - repetitions with fixed data, varying split point
experiment_pars = [{'length':seqlen, 'simseed':102, 'infseed':infseed, 'numparticles':numparticles, 'lag':lagfactor,
                    'nsam':nsam, 'biasheights':None, 'biasstrengths':None, 'delay':0, 'apf':2, 'splitepoch':splitepoch, 'model':model}
                   for (infseed, numparticles, splitepoch, model) in uniq_chain(
                           itertools.product(
                               inference_reps,
                               particles,
                               split_epochs,
                               models ))]

def quantile( k, n, M=4 ):
    if k == 0: return 0.0
    return -math.log(1-k/(n+0.0)) / (M * (M-1)/2.0)


# run an experiment.  keyword parameters must match those in experiment_pars
def run_experiment( length, simseed, infseed, numparticles, lag, nsam, biasheights, biasstrengths, delay, apf, splitepoch, model ):
    label = "L{}_I{}_P{}_NSAM{}_SPLT{}_M{}".format(int(length),infseed,numparticles,nsam,splitepoch,model)
    label = label.replace(" ","").replace("[","").replace("]","")
    str_parameter = "{}".format(model)
    int_parameter = splitepoch
    if experiment_base.have_result( name = experiment_name,
                                    sequence_length = length,
                                    dataseed = simseed,
                                    infseed = infseed,
                                    np = numparticles,
                                    num_samples = nsam,
                                    bias_heights = str(biasheights),
                                    bias_strengths = str(biasstrengths),
                                    aux_part_filt = apf,
                                    int_parameter = int_parameter,
                                    str_parameter = str_parameter):
        print "Skipping " + label
        return
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( experiment_base.datapath + experiment_name )
    e.filename_disambiguator = label

    # set simulation parameters
    e.str_parameter = str_parameter
    e.int_parameter = int_parameter
    e.pop.N0 = N0
    e.pop.mutation_rate = mu
    e.pop.recombination_rate = rho
    e.pop.num_samples = nsam
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = experiment_base.scrmpath
    e.pop.sample_populations = [1]*(nsam/2) + [2]*(nsam/2)
    e.pop.population_sizes = [ [1.0, 0.5], [0.95, 0.05] ,[1,1] ]         # Using 0.95/0.05 here works with 1000pp and 150mbp.  0.9/0.1 should work with 10k/1000mbp
    e.pop.change_points = [ 0, migr_time / (4*N0*g), split_time / (4*N0*g) ]
    e.pop.migration_commands = [ None, None, "-ej {} 2 1".format( split_time / (4*N0*g) ) ]

    # need to set pop sizes, migr rates and commands explicitly when change pts
    # differ between simulation and inference.

    quantiles = 30
    e.smcsmc_change_points = [ quantile(k,quantiles) for k in range(quantiles) ]
    e.smcsmc_migration_commands = [ None ] * splitepoch + \
                                  [ "-ej {} 2 1".format( quantile(splitepoch, quantiles) ) ] + \
                                  [None] * (quantiles-splitepoch-1)
    e.smcsmc_initial_pop_sizes = ([ [1.0,1.0] ] * (splitepoch//2)) + \
                                 ([ [1.0,1.0] ] * (splitepoch - splitepoch//2)) + \
                                 ([ [1.0,1.0] ] * (quantiles-splitepoch))
    e.smcsmc_initial_migr_rates = [[[0,1],[1,0]]] * (splitepoch // 2) + \
                                  [[[0,1],[1,0]]] * (splitepoch - splitepoch//2) + \
                                  [[[0,0],[0,0]]] * (quantiles-splitepoch)

    if model == 0:
        e.pop.migration_rates = [ [[0,0],[0,0]],
                                  [[0,5],[0,0]],    # [[x, M12],[M21, x]] -- migration rate 5 from pop 2 to 1
                                  [[0,0],[0,0]] ]
    elif model == 1:
        e.pop.migration_rates = [ [[0,0],[0,0]],
                                  [[0,0],[5,0]],    # [[x, M12],[M21, x]] -- migration rate 5 from pop 1 to 2
                                  [[0,0],[0,0]] ]
    elif model == 2:
        e.pop.migration_rates = [ [[0,0],[0,0]],
                                  [[0,0],[0,0]],    # [[x, M12],[M21, x]] -- migration rate 5 from pop 1 to 2
                                  [[0,0],[0,0]] ]
    elif model == 3:
        e.pop.migration_rates = [ [[0,0],[0,0]],
                                  [[0,5],[5,0]],    # [[x, M12],[M21, x]] -- migration rate 5 from pop 1 to 2
                                  [[0,0],[0,0]] ]
    else:
        raise ValueError("bad model")

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
    e.guided_recomb_alpha = -1  # remove .recomb.gz files

    # perform inference and store results
    e.infer( case = simseed )
    e.resultsToMySQL( db = experiment_base.db )
    e.success = True   # remove files
    print "Done " + label
    return


def plot_experiment():
    import pandas as pd
    from rpy2.robjects import pandas2ri, r
    pandas2ri.activate()
    records = experiment_base.get_data_from_database( experiment_name )

    df = pd.DataFrame( data = records )
    r.source("unidirmigr_split.R")
    g = 30.0

    r('plot.smcsmc')( df, g )


if __name__ == "__main__":
    experiment_base.run_experiment = run_experiment
    experiment_base.plot_experiment = plot_experiment
    experiment_base.main( experiment_name, experiment_pars )


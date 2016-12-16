import sys
import os
import inspect
import argparse
import itertools
from multiprocessing import Pool
from sqlalchemy import Table, MetaData, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

Base = declarative_base()

# various files and directory paths
db = "experimentsdb"
scrmpath = "../../scrm"
smcsmcpath = "../../smcsmc"
datapath = "data/"

###################################################################################
#
# experiment:
#  estimate bias and variance, dependent on length of input sequence
#  also look at contribution of stochasticity of coalescent itself
#
###################################################################################


sys.path.extend( ["../../newtests","../newtests"] )
import test_const_pop_size


# name of this experiment
experiment_name = "constpopsize_3epochs_particledependence"

# class
experiment_class = test_const_pop_size.TestConstPopSize_ThreeEpochs

# parameters for this experiment
inference_reps = 10
seqlen = 50e6
particles = [50, 100, 200, 500, 1000]
experiment_pars = [{'length':seqlen, 'seed':seed, 'numparticles':numparticles}
                   for (seqlen, seed, numparticles) in (
                           # repetitions with unique data
                           [(seqlen, 100+rep, np)
                            for rep, np in itertools.product(range(inference_reps), particles)])]

# run a single experiment
def run_experiment_map( pars ):
    return run_experiment( **pars )


def run_experiment( length, seed, numparticles ):
    if have_result( length, seed, numparticles ):
        return "np={} seed={} -- skipped".format(numparticles, seed)
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( datapath + experiment_name )
    # set simulation parameters
    e.pop.sequence_length = length
    e.pop.seed = (seed,)
    e.pop.scrmpath = scrmpath
    e.filename_disambiguator = "_L{}_S{}_P{}".format(int(length),seed,numparticles)
    # set inference parameters
    e.seqlen = length
    e.seed = (seed,)
    e.np = numparticles
    e.smcsmcpath = smcsmcpath
    # perform inference and store results
    e.infer( case = seed )
    e.resultsToMySQL( db = db )
    #e.success = True   # remove files
    return "np={} seed={}".format(numparticles, seed)


def have_result( length, seed, numparticles ):
    """ see if the database already contains the required result """
    engine = create_engine("sqlite:///" + db)
    if not engine.dialect.has_table(engine, "experiment"):
        return False
    metadata = MetaData(bind=engine)
    class Experiment(Base):
        __table__ = Table('experiment', metadata, autoload=True)
    Session = sessionmaker(bind=engine)
    session = Session()
    result = session.query(Experiment).filter_by(name = experiment_name,
                                                 sequence_length = length,
                                                 dataseed = seed,
                                                 infseed = seed,
                                                 np = numparticles).first()
    session.commit()
    session.close()
    return result != None


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Run experiment ' + experiment_name)
    parser.add_argument('-j', '--jobs', type=int, default=1, help='set number of threads')
    args = parser.parse_args()

    if args.jobs == 1:
        for i in map( run_experiment_map, experiment_pars ):
            pass
    else:
        pool = Pool(processes = args.jobs)
        p = pool.map_async(run_experiment_map, experiment_pars )
        try:
            results = p.get(100000)
            pool.terminate()
        except KeyboardInterrupt:
            print "Received keyboard interrupt -- stopping"

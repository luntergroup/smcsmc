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
experiment_name = "constpopsize_3epochs_lagdependence"

# class
experiment_class = test_const_pop_size.TestConstPopSize_ThreeEpochs

# parameters for this experiment
inference_reps = 10
seqlen = 50e6
particles = 500
lagfactors = [0.5, 1.0, 2.0, 4.0, 8.0]
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'infseed':infseed, 'numparticles':numparticles, 'lag':lag}
                   for (seqlen, simseed, infseed, numparticles, lag) in (
                           # repetitions with unique data
                           [(seqlen, 100+rep, 100+rep, particles, lag)
                            for rep, lag in itertools.product(range(inference_reps), lagfactors)] +
                           # repetitions with fixed data
                           [(seqlen, 100,     100+rep, particles, lag)
                            for rep, lag in itertools.product(range(1,inference_reps), lagfactors)])]

# run a single experiment
def run_experiment_map( pars ):
    return run_experiment( **pars )


def run_experiment( length, simseed, infseed, numparticles, lag ):
    if have_result( length, simseed, infseed, numparticles, lag ):
        return
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( datapath + experiment_name )
    # set simulation parameters
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = scrmpath
    e.filename_disambiguator = "_L{}_S{}_I{}_P{}_G{}".format(int(length),simseed,infseed,numparticles, lag)
    # set inference parameters
    e.seqlen = length
    e.seed = (infseed,)
    e.np = numparticles
    e.lag = lag
    e.smcsmcpath = smcsmcpath
    # perform inference and store results
    e.infer( case = simseed )
    e.resultsToMySQL( db = db )
    #e.success = True   # remove files
    return


def have_result( length, simseed, infseed, numparticles, lag ):
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
                                                 dataseed = simseed,
                                                 infseed = infseed,
                                                 np = numparticles,
                                                 lag = lag).first()
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

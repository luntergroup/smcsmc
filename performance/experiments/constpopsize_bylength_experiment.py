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


sys.path.extend( ["../../newtests"] )
import test_const_pop_size


# name of this experiment
experiment_name = "constpopsize_3epochs_lengthdependence"

# class
experiment_class = test_const_pop_size.TestConstPopSize_ThreeEpochs

# parameters for this experiment
inference_reps = 5
seqlens = [1e6, 2e6, 5e6, 10e6, 20e6]
experiment_pars = [{'length':seqlen, 'simseed':simseed, 'smcseed':smcseed,'missing':missing}
                   for (seqlen, simseed, smcseed, missing) in (
                           # repetitions with unique data
                           [(seqln, 1+rep, 1+rep, False)
                            for seqln, rep in
                            itertools.product(seqlens, range(inference_reps))] +
                           # repetitions with fixed data (avoiding duplications)
                           [(seqln, 1,     1+rep, False)
                            for seqln, rep in
                            itertools.product(seqlens, range(1, inference_reps))] +
                           # repetitions with missing data.  Use seed to separate filenames
                           [(seqln, 1000+rep,  1000+rep, True)
                            for seqln, rep in
                            itertools.product(seqlens, range(inference_reps))])]


# run a single experiment
def run_experiment_map( pars ):
    return run_experiment( **pars )


def run_experiment( length, simseed, smcseed, missing ):
    if have_result( length, simseed, smcseed, missing ):
        return "L={} SeqSeed={} SmcSeed={} -- skipped".format(length, simseed, smcseed)
    e = experiment_class( 'setUp' )  # fake test fn to keep TestCase.__init__ happy
    e.setUp( datapath + experiment_name )
    # set simulation parameters
    e.pop.sequence_length = length
    e.pop.seed = (simseed,)
    e.pop.scrmpath = scrmpath
    if missing:
        e.missing_leaves = [0,1]
    e.filename_disambiguator = "_L{}_S{}_I{}_M{}".format(int(length),simseed,smcseed,missing)
    # set inference parameters
    e.seqlen = length
    e.seed = (smcseed,)
    e.smcsmcpath = smcsmcpath
    # perform inference and store results
    e.infer( case = smcseed )
    e.resultsToMySQL( db = db )
    #e.success = True   # remove files
    return "L={} SeqSeed={} SmcSeed={} missing={}".format(length, simseed, smcseed, missing)


def have_result( length, simseed, smcseed, missing ):
    """ see if the database already contains the required result """
    engine = create_engine("sqlite:///" + db)
    if not engine.dialect.has_table(engine, "experiment"):
        return False
    metadata = MetaData(bind=engine)
    class Experiment(Base):
        __table__ = Table('experiment', metadata, autoload=True)
    Session = sessionmaker(bind=engine)
    session = Session()
    missing_leaves = str( [0,1] if missing else [] )
    result = session.query(Experiment).filter_by(name = experiment_name,
                                                 sequence_length = length,
                                                 dataseed = simseed,
                                                 infseed = smcseed,
                                                 missing_leaves = missing_leaves).first()
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

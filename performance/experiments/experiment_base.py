import sys
import os
import inspect
import argparse
import itertools
from multiprocessing import Pool
from sqlalchemy import Table, MetaData, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

###################################################################################
#
# stuff shared by all experiments
#
###################################################################################

Base = declarative_base()

# various files and directory paths
db = "experimentsdb"
scrmpath = "../../scrm"
smcsmcpath = "../../python/smcsmc.py"
datapath = "data/"

sys.path.extend( ["../../newtests","../newtests"] )


# run a single experiment
def run_experiment_map( pars ):
    return run_experiment( **pars )


# check if this experiment has already been run
def have_result( **kwargs ):
    """ see if the database already contains the required result """
    engine = create_engine("sqlite:///" + db)
    if not engine.dialect.has_table(engine, "experiment"): return False
    metadata = MetaData(bind=engine)
    class Experiment(Base):
        __table__ = Table('experiment', metadata, autoload=True)
    Session = sessionmaker(bind=engine)
    session = Session()
    result = session.query(Experiment).filter_by( **kwargs ).first()
    session.commit()
    session.close()
    return result != None


# remove all data for thsi experiment
def remove( experiment_name ):
    """ remove any previous results for this experiment """
    engine = create_engine("sqlite:///" + db)
    if not engine.dialect.has_table(engine, "experiment"): return False
    metadata = MetaData(bind=engine)
    class Experiment(Base):
        __table__ = Table('experiment', metadata, autoload=True)
    class Result(Base):
        __table__ = Table('result', metadata, autoload=True)
    Session = sessionmaker(bind=engine)
    session = Session()
    session.execute("DELETE FROM result WHERE exp_id IN ( SELECT id FROM experiment WHERE name = '{}')".format(
        experiment_name
    ))
    session.execute("DELETE FROM experiment WHERE name = '{}'".format(
        experiment_name
    ))
    session.commit()
    session.close()


#
# main code.  Requires one function to be defined:
#
# - run_experiment: function that runs an experiment; takes values in parameter dictionary as keyword args
#
def main( experiment_name, experiment_pars ):

    parser = argparse.ArgumentParser(description='Run experiment ' + experiment_name)
    parser.add_argument('-j', '--jobs', type=int, default=1, help='set number of threads')
    parser.add_argument('-f', '--force', action='store_true', help='overwrite existing results in db')
    args = parser.parse_args()

    if args.force:
        remove( experiment_name )
    if args.jobs == 1:
        for i in map( run_experiment_map, experiment_pars ):
            pass
    else:
        pool = Pool(processes = args.jobs)
        p = pool.map_async(run_experiment_map, experiment_pars )
        try:
            results = p.get(200000)
            pool.terminate()
        except KeyboardInterrupt:
            print "Received keyboard interrupt -- stopping"

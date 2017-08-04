from __future__ import print_function
import sys
import os
import inspect
import argparse
import itertools
import threading
import traceback
import Queue
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

sys.path.extend( ["../../newtests","../newtests","../../python"] )

import execute


# run a single experiment
def run_experiment_map( pars ):
    global bucket
    try:
        return run_experiment( **pars )
    except Exception:
        bucket.put(sys.exc_info())
        raise


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

    global db
    global bucket

    parser = argparse.ArgumentParser(description='Run experiment ' + experiment_name)
    parser.add_argument('-j', '--jobs', type=int, default=1, help='set number of threads')
    parser.add_argument('-f', '--force', action='store_true', help='overwrite existing results in db; reuse existing iteration files')
    parser.add_argument('-c', '--cluster', action='store_true', help='use qsub to submit job(s) to cluster')
    parser.add_argument('-C', '--cconfig', help='qsub config parameter(s); override ./qsub.conf')
    parser.add_argument('-D', '--db', help='set database to write result into (default: {})'.format(db))
    args = parser.parse_args()

    if args.cluster:
        if args.cconfig is None:
            try:
                args.cconfig = open("./qsub.conf").read().replace('\n',' ')
            except:
                args.cconfig = ""
        execute.qsub_config = args.cconfig   # assign to global variable within module 'execute'
    if args.db:
        db = args.db
    if args.force:
        remove( experiment_name )
    bucket = Queue.Queue()                   # to hold exceptions that occur in threads
    if args.jobs == 1:
        for i in map( run_experiment_map, experiment_pars ):
            pass
    else:
        threads = []
        for par in experiment_pars:
            t = threading.Thread( target=run_experiment_map, args=(par,) )
            t.start()
            threads.append( t )

        for t in threads:
            t.join()

        try:
            exc = bucket.get(block=False)
        except Queue.Empty:
            pass
        else:
            exc_type, exc_obj, exc_trace = exc
            print("Caught exception %s %s" % (exc_type, exc_obj), file=sys.stderr)
            traceback.print_tb( exc_trace )
            raise RuntimeError("Error in thread")

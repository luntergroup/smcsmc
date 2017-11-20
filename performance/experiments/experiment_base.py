from __future__ import print_function
from collections import namedtuple
import sys
import os
import inspect
import argparse
import itertools
import threading
import traceback
import time
import random
import Queue

from sqlalchemy import Table, MetaData, create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, scoped_session

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

sys.path.extend( ["../../newtests","../../python"] )
import execute
import test_generic


# stub
def plot_experiment():
    print ("No plots defined - nothing to do")


# run a single experiment
def run_experiment_queue():
    global error_bucket
    global task_queue
    global idle_queue
    time.sleep( 3 * random.random() )    # avoid issues with sqlalchemy...
    while True:
        try:
            pars = task_queue.get(False)
            run_experiment( **pars )
        except Exception as e:
            if isinstance(e, Queue.Empty):
                # signal that thread is done, and kill it
                idle_queue.put( None )
                return
            # store error on queue, and raise error (which will kill thread, and nothing else)
            error_bucket.put(sys.exc_info())
            raise
        task_queue.task_done()
        


# check if this experiment has already been run
def have_result( **kwargs ):
    """ see if the database already contains the required result """
    if db != None:
        engine = create_engine("sqlite:///" + db)
    else:
        engine = create_engine("sqlite://")     # dummy, e.g. for cleanup action
    if not engine.dialect.has_table(engine, "experiment"): return False
    metadata = MetaData(bind=engine)
    Session = sessionmaker(bind=engine)
    session = Session()
    class Experiment(Base):
        __table__ = Table('experiment', metadata, autoload=True)
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
    Session = sessionmaker(bind=engine)
    session = Session()
    class Experiment(Base):
        __table__ = Table('experiment', metadata, autoload=True)
    class Result(Base):
        __table__ = Table('result', metadata, autoload=True)
    session.execute("DELETE FROM result WHERE exp_id IN ( SELECT id FROM experiment WHERE name = '{}')".format(
        experiment_name
    ))
    session.execute("DELETE FROM experiment WHERE name = '{}'".format(
        experiment_name
    ))
    session.commit()
    session.close()


# action: clean up
def clean_up( experiment_class ):
    prefix = experiment_class.caseprefix
    print("Removing non-critical output files for case {}".format(prefix))
    os.system("rm {}.stderr".format(prefix))
    os.system("rm {}.stdout".format(prefix))
    os.system("rm {}/emiter*/*.log".format(prefix))
    os.system("rm {}/emiter*/*.stderr".format(prefix))
    os.system("rm {}/emiter*/*.stdout".format(prefix))
    os.system("rm {}/emiter*/*.recomb.gz".format(prefix))
    os.system("for f in {}/emiter* ; do pth=`pwd`; cd $f; tar --remove-files -czf chunkNN.tar.gz chunk[0123456789]*.out ; cd $pth ; done".format(prefix))


# helper: extract data from database as named tuples
def get_data_from_database( experiment_name ):
    if Base == None: raise ValueError("SQLAlchemy not installed - bailing out")
    if db == None: raise ValueError("Have no database - ?")
    dbtype = "sqlite:///"
    engine = create_engine(dbtype + db)
    if not engine.dialect.has_table(engine, "experiment"): raise ValueError("No Experiment table in database")
    Session = scoped_session(sessionmaker(bind=engine))
    session = Session()
    dbvars = ['result.'+v for v in ['id','exp_id','epoch','start','end','iter','type','ne','rate','frm']] + \
             ['experiment.'+v for v in ['np','num_samples','sequence_length','dataseed','infseed','bias_heights','bias_strengths','str_parameter','int_parameter','aux_part_filt']]
    result = session.execute('SELECT {} FROM experiment INNER JOIN result ON experiment.id = result.exp_id WHERE experiment.name = \'{}\''
                             .format( ', '.join(dbvars), experiment_name ))
    Record = namedtuple('Record', result.keys())
    records = [Record(*r) for r in result.fetchall()]
    return records
    
    
#
# main code.  Requires one function to be defined:
#
# - run_experiment: function that runs an experiment; takes values in parameter dictionary as keyword args
#
def main( experiment_name, experiment_pars ):

    global db
    global datapath
    global error_bucket
    global task_queue
    global idle_queue

    parser = argparse.ArgumentParser(description='Run experiment ' + experiment_name)
    parser.add_argument('-j', '--jobs', type=int, default=1, help='set number experiments to run in parallel')
    parser.add_argument('-f', '--force', action='store_true', help='overwrite existing results in db; reuse existing iteration files')
    parser.add_argument('-c', '--cluster', action='store_true', help='use qsub to submit job(s) to cluster')
    parser.add_argument('-C', '--cconfig', help='qsub config parameter(s); override ./qsub.conf')
    parser.add_argument('-D', '--db', help='set database to write result into (default: {})'.format(db))
    parser.add_argument('-p', '--datapath', help='set data directionry (default: {})'.format(datapath))
    parser.add_argument('-X', '--clean', action='store_true', help='remove non-critical output files; do not run experiment')
    parser.add_argument('-P', '--plot', action='store_true', help='make experiment-specific plots')
    args = parser.parse_args()

    if args.clean:
        db = None                         # force a 'not found', so that experiment is called back
        test_generic.inference_callback = clean_up
        for pars in experiment_pars:
            run_experiment( **pars )
        sys.exit(0)
    if args.cluster:
        if args.cconfig is None:
            try:
                args.cconfig = open("./qsub.conf").read().replace('\n',' ')
            except:
                args.cconfig = ""
        execute.qsub_config = args.cconfig   # assign to global variable within module 'execute'
    if args.db:
        db = args.db
    if args.datapath:
        datapath = args.datapath + "/"
    if args.plot:
        plot_experiment()
        sys.exit(0)
    if args.force:
        remove( experiment_name )
    error_bucket = Queue.Queue()                   # to hold exceptions that occur in threads
    task_queue = Queue.Queue()
    idle_queue = Queue.Queue()
    for par in experiment_pars:
        task_queue.put( par )

    if args.jobs == 1:
        run_experiment_queue()
    else:
        for i in range(args.jobs):
            t = threading.Thread( target=run_experiment_queue )
            t.daemon = True
            t.start()

        # wait for all tasks to have been processed, and check for errors
        while idle_queue.qsize() < args.jobs:
            if not error_bucket.empty():
                break
            time.sleep( 10 )

        try:
            exc = error_bucket.get(block=False)
        except Queue.Empty:
            pass
        else:
            exc_type, exc_obj, exc_trace = exc
            print("Caught exception %s %s" % (exc_type, exc_obj), file=sys.stderr)
            traceback.print_tb( exc_trace )
            raise RuntimeError("Error in thread")

        # terminate program.  Since no error was produced, and args.jobs workers are idle,
        # they must all have finished.  But just to be sure...
        task_queue.join()

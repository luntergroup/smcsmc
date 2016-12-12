from __future__ import print_function

import unittest
import os
import subprocess
import itertools

try:
    from sqlalchemy import Column, Integer, Float, String, DateTime, ForeignKey, Table, MetaData, create_engine
    from sqlalchemy.ext.declarative import declarative_base
    from sqlalchemy.orm import sessionmaker
    from sqlalchemy.sql import func
    Base = declarative_base()
    Session = sessionmaker()
except:
    Base = None

import populationmodels


# All tests are derived from TestGeneric.  A class implementing a test should
#
# 1. implement a setUp(self) method that
#    - calls TestGeneric.setUp( "prefix" )
#    - optionally modifies or sets further parameters
#    - assigns a Population to self.pop
#
# 2. implement at least one test_xxx method that calls infer()
#
# infer() can be called multiple times, but should each have their own unique case
# number. These all use the same generated data, and should be run in series.
# Temporary files will be cleaned up when the class instance goes out of scope.


class TestGeneric(unittest.TestCase):

    # called every time an instance of TestGeneric is made -- set inference defaults
    def setUp(self, prefix = None):
        # these two must be set by deriving class:
        self.prefix = prefix
        self.pop = None
        # the following are all default values, and may be changed:
        self.scrmpath = "../scrm"
        self.np = 1000
        self.em = 0
        self.seed = (3647837471,)
        self.popt = "-p 1*3+15*4+1"
        self.tmax = 4
        self.smcsmc_change_points = None
        self.missing_leaves = []
        # state:
        self.simulated = False
        self.success = False
        self.cases = []

    # called every time an instance of TestGeneric is destroyed -- remove output file
    def tearDown(self):
        if self.success and self.prefix != None:
            toplevel = itertools.product( [self.prefix],
                                          [".seg",".seg.recomb"] )
            percase = itertools.product( [self.prefix + str(case) for case in self.cases ],
                                         ['.out','.log','.stdout','.stderr','.recomb','.resample'] )
            for prefix, suffix in itertools.chain( toplevel, percase ):
                try:
                    os.unlink( prefix + suffix )
                except OSError:
                    print ("Warning: file ",prefix + suffix," expected but not found")

    # method to build smcsmc command corresponding to the simulated data
    def build_command(self):
        nsamopt = "-nsam {}".format(self.pop.num_samples)
        topt = "-t {}".format(self.pop.mutations)
        ropt = "-r {} {}".format(self.pop.recombinations,self.pop.sequence_length)
        particlesopt = "-Np {np}".format(np=self.np)
        emopt = "-EM {em}".format(em=self.em)
        seedopt = "-seed {seed}".format(seed=' '.join(map(str,self.seed)))
        segopt = "-seg {}".format( self.prefix + ".seg" )
        addopt = self.pop.additional_commands
        migropt = " ".join( [ cmd for cmd in self.pop.migration_commands if cmd != None ] )

        if self.popt:
            # use -p / -tmax pattern to specify epochs for inference
            epochopt = "{popt} -tmax {tmax}".format(popt = self.popt,
                                                    tmax = self.tmax)
        else:
            if self.smcsmc_change_points != None:
                # use explicit epochs for inference
                epochs = self.smcsmc_change_points
            else:
                # use model epochs for inference
                epochs = self.pop.change_points
            epochopt = " ".join(["-eN {time} 1".format(time=time)
                                 for time in epochs])
        self.inference_command = "../smcsmc {nsam} {t} {add} {r} {np} {em} {epochs} {seed} {seg} {migr}".format(
            nsam = nsamopt,
            t = topt,
            r = ropt,
            np = particlesopt,
            em = emopt,
            epochs = epochopt,
            seed = seedopt,
            seg = segopt,
            add = addopt,
            migr = migropt)
        return self.inference_command

    # helper -- generate simulated data, if this has not been done yet
    def simulate(self):
        if self.simulated: return
        if not self.pop: raise ValueError("Must define population before simulating")
        if not self.prefix: raise ValueError("Must define prefix before simulating")
        print ("simulating for",self.__class__.__name__,"...")
        self.pop.filename = self.prefix + ".seg"
        pathname = os.path.dirname(self.pop.filename)
        os.makedirs(pathname, exist_ok=True)
        self.pop.simulate( self.missing_leaves )
        self.simulated = True
    
    # helper -- run smcsmc
    def infer(self, case = 0):
        self.simulate()
        if case in self.cases: raise ValueError("Must run case " + str(case) + " only once")
        self.cases.append(case)
        self.caseprefix = self.prefix + str(case)
        print (" running smcsmc for",self.__class__.__name__,", case",case,"...")
        self.outfile = self.caseprefix + ".out"
        cmd = "{cmd} -o {caseprefix} > {caseprefix}.stdout 2> {caseprefix}.stderr".format(
            cmd = self.build_command(),
            caseprefix = self.caseprefix )
        print (" command:",cmd)
        returnvalue = subprocess.check_call(cmd, shell = True)
        self.assertTrue( returnvalue == 0 )

    # actual test.  However, should not be run on the TestGeneric class itself, as there is 
    # nothing to test
    def test_inference(self):

        if self.__class__.__name__ == "TestGeneric": 
            return

        # Generate data, and perform inference
        self.infer( case = 1 )

        # Store results in db
        self.resultsToMySQL()

        # Read results from .out file
        results = self.readResults()

        # Check results and count out-of-range epochs
        out_of_range = 0
        for epoch, emresults in enumerate(results['Coal']):

            # get the last EM estimate (for 'population' 0)
            results = emresults['estimates'][0]
            result = results[-1]
            print ("Checking epoch",epoch,results,self.target_min[epoch],self.target_max[epoch])

            # see if it is within range
            if result < self.target_min[epoch] or result > self.target_max[epoch]:
                print ("Out of range!")
                out_of_range += 1

        self.assertTrue( out_of_range <= self.max_out_of_range )
        self.success = True

    # helper -- parse results of smcsmc run
    def readResults(self):
        results = {'Coal':[],    # type -> epoch -> {'start' -> # ,'end' -> # ,'estimates' -> from_pop -> iteration -> # }
                   'Recomb':[],
                   'Migr':[]}
        with open(self.outfile) as f:
            for line in f:
                elts = line.strip().split()
                if elts[0] == "Iter": continue
                it, epoch, start, end, typ, from_pop, to_pop, rate, ne = elts[0],elts[1],elts[2],elts[3],elts[4],elts[5],elts[6],elts[9],elts[10]
                if int(epoch) == -1:
                # recombination rate is assigned to epoch -1
                    epoch = 0
                    frop_pop = 0
                if len(results[typ]) <= int(epoch):
                    # there is one estimate of population size or migration per population
                    # (TODO: in fact, there are more for migration, except for two-population models)
                    # (TODO: not storing recombination estimates for now)
                    results[typ].append( {'start': float(start),
                                          'end': float(end),
                                          'estimates': [[] for i in range(self.pop.npop)]} )  
                if typ == "Coal":
                    result = float(ne)
                elif typ == "Migr" or typ == "Recomb":
                    result = float(rate)
                results[typ][int(epoch)]['estimates'][int(from_pop)].append( result )
        return results

    # helper -- store results in mysql database
    def resultsToMySQL(self, dbtype = "sqlite:///", db = "resultsdb"):

        # bail out of SQLAlchemy not installed
        if Base == None: return
        
        # make a connection
        engine = create_engine(dbtype + db)

        # load or create tables
        if engine.dialect.has_table(engine, "experiment"):

            # reflect the database tables we use, using metadata
            metadata = MetaData(bind=engine)
            class Experiment(Base):
                __table__ = Table('experiment', metadata, autoload=True)
            class Result(Base):
                __table__ = Table('result', metadata, autoload=True)

        else:

            # create tables.  Just declaring them does the trick
            class Experiment(Base):
                __tablename__ = "experiment"
                keep_existing = True
                id = Column(Integer, primary_key=True)
                date = Column(DateTime, default=func.now())
                simulate_command = Column(String)
                inference_command = Column(String)
                np = Column(Integer)
                num_samples = Column(Integer)
                sequence_length = Column(Float)
                recombination_rate = Column(Float)
                mutation_rate = Column(Float)
                dataseed = Column(Integer)
                infseed = Column(Integer)
            class Result(Base):
                __tablename__ = "result"
                keep_existing = True
                id = Column(Integer, primary_key=True)
                exp_id = Column(Integer, ForeignKey("experiment.id"))
                iter = Column(Integer)
                epoch = Column(Integer)
                start = Column(Float)
                end = Column(Float)
                type = Column(String)
                frm = Column(Integer)
                to = Column(Integer)
                opp = Column(Float)
                count = Column(Float)
                rate = Column(Float)
                ne = Column(Float)
                ess = Column(Float)

            # commit tables
            Base.metadata.create_all(engine, checkfirst=True)

        # add data for this experiment
        Session.configure(bind=engine)
        session = Session()
        this_exp = Experiment( np = self.np, num_samples=self.pop.num_samples, sequence_length=self.pop.sequence_length,
                               dataseed=self.pop.seed[0], infseed=self.seed[0], simulate_command=self.pop.simulate_command,
                               recombination_rate = self.pop.recombination_rate, mutation_rate = self.pop.mutation_rate,
                               inference_command = self.inference_command )
        session.add( this_exp )
        session.commit()         # this also sets this_exp.id

        # add results
        session = Session()
        with open(self.outfile) as f:
            for line in f:
                elts = line.strip().split()
                if elts[0] == "Iter": continue                     # skip header line
                assert len(elts) == 12
                session.add( Result( exp_id = this_exp.id,
                                     iter = int(elts[0]), epoch = int(elts[1]), start=float(elts[2]), end=float(elts[3]),
                                     type = elts[4], frm = int(elts[5]), to = int(elts[6]), opp=float(elts[7]), count=float(elts[8]),
                                     rate = float(elts[9]), ne = float(elts[10]), ess = float(elts[11]) ) )
        session.commit()



from __future__ import print_function

import unittest
import os
import subprocess

try:
    from sqlalchemy import Column, Integer, Float, String, DateTime, ForeignKey, create_engine
    from sqlalchemy.ext.declarative import declarative_base
    from sqlalchemy.orm import sessionmaker
    from sqlalchemy.sql import func
    Base = declarative_base()
except:
    Base = None

import populationmodels


class TestGeneric(unittest.TestCase):

    # remove simulated data
    @classmethod
    def tearDownClass(cls):
        if cls.segfile != None and cls.success:
            os.unlink( cls.segfile )
            cls.filename = None

    # called every time an instance of TestGeneric is made -- set inference defaults
    def setUp(self):
        self.np = 1000
        self.em = 0
        self.popt = "-p 1*3+15*4+1"
        self.tmax = 4
        self.seed = (3647837471,)
        self.smcsmc_change_points = None
        self.success = False
        self.npop = 1

    # called every time an instance of TestGeneric is destroyed -- remove output file
    def tearDown(self):
        if self.success and 'caseprefix' in self.__dict__:
            if self.caseprefix != None:
                for suffix in ['Resample','.out','.log','.stdout','.stderr']:
                    try:
                        os.unlink( self.caseprefix + suffix )
                    except OSError:
                        print ("Warning: file ",self.caseprefix + suffix," expected but not found")
                        pass
                self.caseprefix = None
        if not self.success:
            # tell derived class that test failed, so don't delete intermediate files
            self.__class__.success = False


    # method to build smcsmc command corresponding to the simulated data
    def build_command(self):
        nsamopt = "-nsam {}".format(self.pop.num_samples)
        topt = "-t {}".format(self.pop.mutations)
        ropt = "-r {} {}".format(self.pop.recombinations,self.pop.sequence_length)
        npopt = "-Np {np}".format(np=self.np)
        emopt = "-EM {em}".format(em=self.em)
        seedopt = "-seed {seed}".format(seed=' '.join(map(str,self.seed)))
        segopt = "-seg {}".format( self.segfile )
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
        return "../smcsmc {nsam} {t} {add} {r} {np} {em} {epochs} {seed} {seg} {migr}".format(
            nsam = nsamopt,
            t = topt,
            r = ropt,
            np = npopt,
            em = emopt,
            epochs = epochopt,
            seed = seedopt,
            seg = segopt,
            add = addopt,
            migr = migropt)

    # helper -- run smcsmc
    def infer(self, case = 0):
        print (" running smcsmc for case",case,"...")
        self.caseprefix = self.prefix + str(case)
        self.outfile = self.caseprefix + ".out"
        cmd = "{cmd} -o {caseprefix} > {caseprefix}.stdout 2> {caseprefix}.stderr".format(
            cmd = self.build_command(),
            caseprefix = self.caseprefix )
        print (" command:",cmd)
        returnvalue = subprocess.check_call(cmd, shell = True)
        self.assertTrue( returnvalue == 0 )

    # helper -- parse results of smcsmc run
    def readResults(self):
        results = {'Coal':[],    # type -> epoch -> {'start' -> # ,'end' -> # ,'estimates' -> iteration -> # }
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
                    # the length of the estimates should be number of populations
                    results[typ].append( {'start': float(start),
                                          'end': float(end),
                                          'estimates': [[] for i in range(self.npop)]} )  
                if typ == "Coal":
                    results[typ][int(epoch)]['estimates'][int(from_pop)].append( float(ne) )
                elif typ == "Migr":
                    results[typ][int(epoch)]['estimates'][int(from_pop)].append( float(rate) )
        return results

    # helper -- store results in mysql database
    def resultsToMySQL(self, db = "sqlite:///resultsdb"):
        # bail out of SQLAlchemy not installed
        if Base == None: return
        
        # make a connection
        engine = create_engine(db)

        # create tables, if they don't already exist
        class Experiment(Base):
            __tablename__ = "experiment"
            id = Column(Integer, primary_key=True)
            np = Column(Integer)
            num_samples = Column(Integer)
            sequence_length = Column(Float)
            dataseed = Column(Integer)
            infseed = Column(Integer)
            date = Column(DateTime, default=func.now())
        class Result(Base):
            __tablename__ = "result"
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
        Base.metadata.create_all(engine, checkfirst=True)
        Session = sessionmaker(bind=engine)

        # add data for this experiment
        session = Session()
        this_exp = Experiment( np = self.np, num_samples=self.pop.num_samples, sequence_length=self.pop.sequence_length,
                               dataseed=self.pop.seed[0], infseed=self.seed[0] )
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



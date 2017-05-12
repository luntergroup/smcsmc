from __future__ import print_function

import unittest
import os
import time
import subprocess
import itertools
import shutil

try:
    from sqlalchemy import Column, Integer, Float, Boolean, String, DateTime, ForeignKey, Table, MetaData, create_engine
    from sqlalchemy.ext.declarative import declarative_base
    from sqlalchemy.orm import sessionmaker
    from sqlalchemy.sql import func
    Base = declarative_base()
    Session = sessionmaker()
except:
    Base = None

import populationmodels
import execute


# All tests are derived from TestGeneric.  A class implementing a test should
# implement a setUp(self) method that
#    - calls TestGeneric.setUp( "prefix" )
#    - sets inference targets for population sizes (target_min, target_max, max_out_of_range)
#    - optionally modifies or sets further parameters
#    - assigns a Population to self.pop
# It can optionally overload test_inference() or add more test_xxx methods, to perform
# tests using the population model.
#
# infer() can be called multiple times, but should each have their own unique case
# number. These all use the same generated data, and should be run in series.
# Temporary files will be cleaned up when the class instance goes out of scope.


class TestGeneric(unittest.TestCase):

    # called every time an instance of TestGeneric is made -- set inference defaults
    def setUp(self, prefix = None, name = None):
        # these must be set by deriving class:
        self.prefix = prefix
        self.name = name
        self.pop = None
        # the following are all default values, and may be changed:
        self.scrmpath = "../scrm"
        self.smcsmcpath = "../smcsmc"
        self.np = 100
        self.em = 0
        self.seed = (1,)
        self.popt = "-p 1*3+15*4+1"
        self.tmax = 4
        self.infer_recombination = True
        self.m_step = True
        self.ancestral_aware = False
        self.phased = True
        self.chunks = 1
        self.lag = 1.0
        self.smcsmc_change_points = None
        self.smcsmc_initial_pop_sizes = None
        self.smcsmc_initial_migr_rates = None
        self.missing_leaves = []            # list of 0-based missing leaves
        self.guided_recomb_alpha = 0      # proportion of posterior recombination mixed in; default 0 (none)
        self.guided_recomb_beta = 4       # extent of smoothing of the recombination change points
        self.bias_heights = [400]
        self.bias_strengths = [2,1]
        self.filename_disambiguator = ""
        self.int_parameter = 0
        self.str_parameter = ""
        # state:
        self.simulated = False
        self.success = False
        self.debug = False
        self.store_in_db = False
        self.cases = []
        self.smcsmc_runtime = -1
        self.smcsmc_version = ""
        self.scrm_version = ""


    # called every time an instance of TestGeneric is destroyed -- remove output file
    def tearDown(self):
        if self.success and self.prefix != None:
            toplevel = itertools.product( [self.pop.filename],
                                          ["",".recomb",".seg.scrm"] )
            percase = itertools.product( ["{}_C{}".format(self.pop.filename[:-4], case)
                                          for case in self.cases ],
                                         ['.out','.log','.stdout','.stderr','.recomb.gz'] )
            for prefix, suffix in itertools.chain( toplevel, percase ):
                try:
                    os.unlink( prefix + suffix )
                except OSError:
                    #print ("Warning: file ",prefix + suffix," expected but not found")
                    # ignore non-existent files; smcsmc.py does not produce .log or .recomb.gz files
                    pass
            shutil.rmtree( os.path.dirname(self.pop.filename), ignore_errors = True )

    # method to build smcsmc command corresponding to the simulated data
    def build_command(self):
        num_epochs = len(self.pop.population_sizes)
        num_populations = len(self.pop.population_sizes[0])

        if self.smcsmc_initial_pop_sizes is not None:
            inf_popsizes = self.smcsmc_initial_pop_sizes
        else:
            inf_popsizes = [ [1] * num_populations
                             for j in range(num_epochs) ]     # starting values for inference

        if self.smcsmc_initial_migr_rates is not None:
            inf_migrationrates = self.smcsmc_initial_migr_rates
        else:
            inf_migrationrates = self.pop.migration_rates         # starting values for inference

        core_cmd = self.pop.core_command_line( inference_popsizes = inf_popsizes,
                                               inference_migrationrates = inf_migrationrates )

        nsamopt = "-nsam {}".format(self.pop.num_samples)

        lagopt = "-calibrate_lag {}".format(self.lag)
        particlesopt = "-Np {np}".format(np=self.np)
        emopt = "-EM {em}".format(em=self.em)
        seedopt = "-seed {seed}".format(seed=' '.join(map(str,self.seed)))
        segopt = "-seg {}".format( self.pop.filename )
        guideopt = "-alpha {}".format( self.guided_recomb_alpha ) if self.guided_recomb_alpha>0 else ""

        ### TODO: I suspect this doesn't work anymore, as self.pop.core_command_line already
        ###       specifies the epochs to infer.  That's probably what we want to do anyway,
        ###       but in that case the test code shouldn't suggest that it supports this!
        if self.popt:
            # use -p / -tmax pattern to specify epochs for inference
            epochopt = "{popt} -tmax {tmax}".format(popt = self.popt,
                                                    tmax = self.tmax)
            # extract the number of epochs from the pattern string:
            # add up the Ns in the 'N*M' and 'N' patterns, which are separated with '+'es
            num_epochs = sum( [ int(elt.split('*')[0])
                                for elt in self.popt.split(' ')[1].split('+') ] )
        else:
            if self.smcsmc_change_points != None:
                # use explicit epochs for inference
                epochs = self.smcsmc_change_points
            else:
                # use model epochs for inference
                epochs = self.pop.change_points
            epochopt = "-tmax {tmax}".format(tmax = self.tmax)
            num_epochs = len(epochs)

        pilotsopt = ""
        recinfopt = ""
        ancawareopt = ""
        mstepopt = ""
            
        if self.bias_heights != None:
            pilotsopt = "-bias_heights {} -bias_strengths {}".format(
                ' '.join(map(str,self.bias_heights)),
                ' '.join(map(str,self.bias_strengths))
            )

        if not self.infer_recombination:
            recinfopt = "-no_infer_recomb" if ".py" in self.inference_command else "-xr 1-{}".format(num_epochs)

        if self.ancestral_aware:
            ancawareopt = "-ancestral_aware"

        if not self.m_step:
            mstepopt = "-no_m_step"
            
        self.inference_command = "{smcsmc} {core} {nsam} {recinf} {np} {em} {guide} " \
                                 "{lag} {epochs} {seed} {seg} {pilots} {ancestral_aware} {mstep}".format(
                                     smcsmc = self.smcsmcpath,
                                     core = core_cmd,
                                     nsam = nsamopt,
                                     recinf = recinfopt,
                                     np = particlesopt,
                                     em = emopt,
                                     guide = guideopt,
                                     lag = lagopt,
                                     epochs = epochopt,
                                     seed = seedopt,
                                     seg = segopt,
                                     pilots = pilotsopt,
                                     ancestral_aware = ancawareopt,
                                     mstep = mstepopt)
        if self.debug:
            print (self.inference_command)
        return self.inference_command

    # helper -- generate simulated data, if this has not been done yet
    def simulate(self):
        if self.simulated: return
        if not self.pop: raise ValueError("Must define population before simulating")
        if not self.prefix: raise ValueError("Must define prefix before simulating")
        self.pop.filename = self.prefix + self.filename_disambiguator + ".seg"
        print ("\n simulating for",self.__class__.__name__," to ",self.pop.filename)
        pathname = os.path.dirname(self.pop.filename)
        if pathname != '' and not os.path.exists(pathname):
            try:
                os.makedirs(pathname)
            except:
                # to cope with race conditions when multiple tests are run in parallel
                pass
        self.pop.simulate( self.missing_leaves, self.phased, self.debug )
        self.simulated = True

    # helper -- run smcsmc
    def infer(self, case = 0):
        self.simulate()
        if case in self.cases: raise ValueError("Must run case " + str(case) + " only once")
        self.cases.append(case)
        self.caseprefix = "{}_C{}".format(self.pop.filename[:-4], case)
        if os.path.exists(self.caseprefix):
            raise ValueError("File {} already exists".format(self.caseprefix))
        print (" running smcsmc for",self.__class__.__name__,", case",case,"...")
        self.outfile = self.caseprefix + ".out"
        cmd = "{cmd} -o {caseprefix} > {caseprefix}.stdout 2> {caseprefix}.stderr".format(
            cmd = self.build_command(),
            caseprefix = self.caseprefix )
        versioncmd = [self.smcsmcpath,"-v"]
        try:
            versiondata = subprocess.Popen(versioncmd,
                                           stdout=subprocess.PIPE).communicate()[0].split('\n')
        except:
            raise ValueError("Failed to execute {}".format(" ".join(versioncmd)))
        if len(versiondata) >= 3:
            self.smcsmc_version = versiondata[1].strip().split()[-1]
            self.scrm_version   = versiondata[2].strip().split()[-1]
        start = time.time()
        # run, or submit to the cluster; note that the execution time is meaningless in the latter case...
        returnvalue = execute.check_call(cmd,
                                         outputdir=os.path.dirname(self.caseprefix),
                                         name=os.path.basename(self.caseprefix))
        end = time.time()
        self.assertTrue( returnvalue == 0 )
        self.smcsmc_runtime = end-start

    # actual test.  However, should not be run on the TestGeneric class itself, as there is
    # nothing to test
    def test_inference(self):

        if self.__class__.__name__ == "TestGeneric":
            return

        # Generate data, and perform inference
        self.infer( case = 1 )

        # Store results in db
        if self.store_in_db:
            self.resultsToMySQL()

        # Read results from .out file
        results = self.readResults()

        # Check results and count out-of-range epochs
        out_of_range = 0
        for result in results:

            if result['type'] == "Coal":

                estimate = result['ne'][-1]
                ess = result['ess'][-1]
                epoch = result['epoch']
                pop = result['pop']
                # extract the target range from the population model
                try:
                    this_parameter = (item
                                      for item in self.targets
                                      if (item["type"] == "Coal"
                                          and item["pop"] == pop
                                          and item["epoch"] == epoch)).next()
                except:
                    continue
                print ("  checking Ne of pop",pop,"epoch",epoch,end=" ")

            elif result['type'] == "Recomb":

                estimate = result['rate'][-1]
                ess = result['ess'][-1]
                # extract the target range
                try:
                    this_parameter = (item for item in self.targets if item["type"] == "Recomb").next()
                except:
                    continue
                print ("  checking recomb rate        ",end=" ")


            elif result['type'] == "Migr":

                estimate = result['rate'][-1]
                ess = result['ess'][-1]
                epoch = result['epoch']
                from_pop = result['from_pop']
                to_pop = result['to_pop']
                # extract the target range
                try:
                    this_parameter = (item
                                      for item in self.targets
                                      if (item["type"] == "Migr"
                                          and item["from_pop"] == from_pop
                                          and item["to_pop"] == to_pop
                                          and item["epoch"] == epoch)).next()
                except:
                    continue
                print ("  checking migr",from_pop,"->",to_pop,"epoch",result['epoch'],end=" ")

            elif result['type'] == "LogL":
                estimate = result['rate'][-1]
                ess = 1.0
                try:
                    this_parameter = (item for item in self.targets if item["type"] == "LogL").next()
                except:
                    continue
                print ("  checking log likelihood     ",end=" ")                
                
            target_min = this_parameter['min']
            target_max = this_parameter['max']
            truth      = this_parameter['truth']
            min_ess    = this_parameter.get('ess',0.0)

            msg = ""
            if estimate < target_min or estimate > target_max:
                msg += "  ** Out of range! **"
                out_of_range += 1
            if ess < min_ess:
                msg += "  ** ESS too low! **"
                out_of_range += 1
            print (" True {:9.4g} Est {:9.4g} Range {:9.4g} - {:9.4g}; ESS {:7.3g} Min {:7.3g}{}".format(
                truth, estimate, target_min, target_max, ess, min_ess, msg))

        self.assertTrue( out_of_range <= self.max_out_of_range )
        self.success = True

    # helper -- parse results of smcsmc run
    def readResults(self):
        results = []

        with open(self.outfile) as f:
            for line in f:
                elts = line.strip().split()
                if elts[0] == "Iter": continue
                it, epoch, start, end, typ, from_pop, to_pop, _, _, rate, ne, ess = elts[:12]

                if int(it) == 0:
                # this is the first time we've seen this parameter, create a dictionary
                    if typ == "Coal":
                        results.append( {'type' :'Coal', 'pop': int(from_pop),
                                         'epoch': int(epoch), 'start': float(start),
                                         'end': float(end), 'ne': [float(ne)], 'ess': [float(ess)]} )
                    elif typ == "Recomb":
                        results.append( {'type' :'Recomb', 'rate': [float(rate)], 'ess': [float(ess)]} )
                    elif typ == "Migr":
                        results.append( {'type' :'Migr',
                                         'from_pop': int(from_pop), 'to_pop': int(to_pop),
                                         'epoch': int(epoch), 'start': float(start),
                                         'end': float(end), 'rate': [float(rate)], 'ess': [float(ess)]} )
                    elif typ == "LogL":
                        results.append( {'type': 'LogL', 'rate': [float(rate)] } )
                else:
                # append this parameter estimate to the list within the correct dictionary
                # first extract the correct dictionary, then append it
                    if typ == "Coal":
                        result = (item for item in results
                                  if (item['type'] == "Coal" and
                                      item['pop'] == int(from_pop) and
                                      item['epoch'] == int(epoch))).next()
                        result['ne'].append( float(ne) )
                        result['ess'].append( float(ess) )

                    elif typ == "Recomb":
                        result = (item for item in results if item['type'] == "Recomb").next()
                        result['rate'].append( float(rate) )
                        result['ess'].append( float(ess) )

                    elif typ == "Migr":
                        result = (item for item in results
                                  if (item['type'] == "Migr" and
                                      item['from_pop'] == int(from_pop) and
                                      item['to_pop'] == int(to_pop) and
                                      item['epoch'] == int(epoch))).next()
                        result['rate'].append( float(rate) )
                        result['ess'].append( float(ess) )
                    elif typ == "LogL":
                        result = (item for item in results if item['type'] == "LogL").next()
                        result['rate'].append( float(rate) )

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
                id                       = Column(Integer, primary_key=True)
                date                     = Column(String)
                name                     = Column(String)
                simulate_command         = Column(String)
                inference_command        = Column(String)
                np                       = Column(Integer)
                lag                      = Column(Float)
                num_samples              = Column(Integer)
                sequence_length          = Column(Integer)
                chunks                   = Column(Integer)
                missing_leaves           = Column(String)
                recombination_rate       = Column(Float)
                mutation_rate            = Column(Float)
                ancestral_aware          = Column(Boolean)
                phased                   = Column(Boolean)
                infer_recombination      = Column(Boolean)
                pattern                  = Column(String)
                initial_Ne_values        = Column(String)
                initial_migr_values      = Column(String)
                bias_heights             = Column(String)
                bias_strengths           = Column(String)
                guided_recomb_alpha      = Column(Float)
                guided_recomb_beta       = Column(Float)
                dataseed                 = Column(Integer)
                infseed                  = Column(Integer)
                smcsmc_runtime           = Column(Float)
                scrm_version             = Column(String)
                smcsmc_version           = Column(String)
                int_parameter            = Column(Integer)   # generic optional parameter
                str_parameter            = Column(String)    # generic optional parameter
                # add values for analysing dividing data into chunks. This would require changes to Result too.
                # add functionality for different mutation/recombination rates for simulation and inference
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

            # commit tables
            Base.metadata.create_all(engine, checkfirst=True)

        # add data for this experiment
        Session.configure(bind=engine)
        session = Session()
        name = self.name if self.name != None else self.prefix.split('/')[-1]
        this_exp = Experiment( name                     = name,
                               date                     = time.strftime("%d/%m/%Y")+" "+time.strftime("%X"),
                               simulate_command         = self.pop.simulate_command,
                               inference_command        = self.inference_command,
                               np                       = self.np,
                               lag                      = self.lag,
                               num_samples              = self.pop.num_samples,
                               sequence_length          = int(self.pop.sequence_length),
                               chunks                   = int(self.chunks),
                               missing_leaves           = str(self.missing_leaves),
                               recombination_rate       = self.pop.recombination_rate,
                               mutation_rate            = self.pop.mutation_rate,
                               ancestral_aware          = self.ancestral_aware,
                               phased                   = self.phased,
                               infer_recombination      = self.infer_recombination,
                               pattern                  = self.popt+" -tmax "+self.tmax if self.popt else self.popt,
                               initial_Ne_values        = ' '.join(map(str,self.smcsmc_initial_pop_sizes)) if self.smcsmc_initial_pop_sizes else str(self.pop.population_sizes),
                               initial_migr_values      = ' '.join(map(str,self.smcsmc_initial_migr_rates)) if self.smcsmc_initial_migr_rates else str(self.pop.migration_rates),
                               bias_heights             = ' '.join(map(str,self.bias_heights)) if self.bias_heights else self.bias_heights,
                               bias_strengths           = ' '.join(map(str,self.bias_strengths)) if self.bias_strengths else self.bias_strengths,
                               guided_recomb_alpha      = self.guided_recomb_alpha,
                               guided_recomb_beta       = self.guided_recomb_beta,
                               dataseed                 = self.pop.seed[0],
                               infseed                  = self.seed[0],
                               smcsmc_runtime           = self.smcsmc_runtime,
                               scrm_version             = self.scrm_version,
                               smcsmc_version           = self.smcsmc_version,
                               int_parameter            = self.int_parameter,
                               str_parameter            = self.str_parameter )
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
                                     iter = int(elts[0]), epoch = int(elts[1]), start = float(elts[2]),
                                     end = float(elts[3]), type = elts[4], frm = int(elts[5]),
                                     to = int(elts[6]), opp = float(elts[7]), count = float(elts[8]),
                                     rate = float(elts[9]), ne = float(elts[10]), ess = float(elts[11])))
        session.commit()
        session.close()




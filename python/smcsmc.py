#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import logging
import gzip
import subprocess
import time
from collections import defaultdict
import populationmodels
import processrecombination


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


#
# Todo: recomb.gz is emitted from position 1 rather than startpos; is it expected to start from 1?
# If so, leave for now?
#
# Check why we don't use process and use .recomb.gz files
#

class Smcsmc:
    def __init__(self, opts):
        self.theta = None
        self.rho = None
        self.startpos = 1
        self.length = None
        self.segfile = None
        self.outprefix = None
        self.logfile = None
        self.emiters = 0
        self.chunks = 1
        self.infer_recomb = True
        self.do_m_step = True
        self.alpha = 0.0            # posterior mix-in; 0 means use prior
        self.beta = 4               # smoothness parameter; see processrecombination.py
        self.maxNE = 1e99
        self._processes = []
        self.threads = True
        
        self.smcsmcpath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../smcsmc'))
        
        # define command line options, in 3 classes:
        self.classes = ["Population-model related options (fixed)",
                        "Population-model related options (inferred; initial values)",
                        "Inference-related options",
                        "Options related to the EM process"]

        self.options = [
            (0, '-t', 'th',    'Set mutation rate (expected number of mutations) for the locus'),
            (0, '-nsam', 'n',  'Set the total number of samples to n'),
            (0, '-I', 'n s1..sn','Use an n-population island model with si individuals sampled'),
            (0, '-eI','t s1..sn','Sample s1..sn indiviuals from their populations at time t*4N0'),
            (0, '-ej', 't i j','Speciation event at t*4N0; creates population j from population i (not implemented)'),
            (0, '-N0', 'N',    'Set (unscaled) default population size to N (10000; scales the output)'),
            
            (1, '-r', 'r l',   'Set per locus recombination rate and locus length'),
            (1, '-eN', 't n',  'Set the present day size of all populations to n*N0'),
            (1, '-en', 't i n','Change the size of population i to n*N0 at time t*4N0'),
            (1, '-eM', 't m',    'Change the symmetric backward migration rate to m/(npop-1) at time t*4N0'),
            (1, '-em', 't i j m', 'Change the backward migration rate, from population i to population j, to m/(npop-1) at time t*4N0'),
            (1, '-ema', 't s11 s12 ...', 'Set backward migration rate matrix at time t*4N0'),
            
            (2, '-seg', 'f',   'Input .seg file'),
            (2, '-startpos', 'x', 'First locus to process (1)'),
            (2, '-o', 'f',     'Output prefix'),
            (2, '-Np', 'n',    'Number of particles'),
            (2, '-seed', 's',  'Random number seed'),
            (2, '-ancestral_aware', '', 'Assume that haplotype 0 is ancestral'),
            (2, '-bias_heights', 't0..tn', 'Set recombination bias times to h0..hn * 4N0'),
            (2, '-bias_strengths','s1..sn','Set recombination bias strenghts'),
            (2, '-calibrate_lag','s','Explanation required!'),
            
            (3, '-EM', 'n',     'Number of EM iterations-1 ({})'.format(self.emiters)),
            (3, '-cap', 'n',    'Set (unscaled) upper bound on effective population size'),
            (3, '-chunks','n',  'Number of chunks computed in parallel ({})'.format(self.chunks)),
            (3, '-no_infer_recomb', '', 'Do not infer recombination rate'),
            (3, '-no_m_step','','Do not update parameters (but do infer recombination guide)'),
            (3, '-alpha', 't',  'Fraction of posterior recombination to mix in to recombination guide ({})'.format(self.alpha)),
            (3, '-smcsmcpath', 'f','Path to smcsmc executable ({})'.format(self.smcsmcpath)),
            (3, '-nothreads','','Calculate chunks sequentially rather than in parallel threads'),
            (3, '-log', 'f',    'Log file (stdout)'),
            (3, '-@', 'f',      'File with line-, space-, or tab-separated options'),
            (3, '-v', '',       'Show version of underlying smcsmc and scrm implementations'),  
            (3, "-help", '',    'This help')]

        self.opts = opts
        
        
    def print_help_and_exit(self):
        if len(self.opts) == 1 and self.opts[0] == "-v":
            subprocess.check_call(" ".join([self.smcsmcpath] + self.opts), shell=True)
            sys.exit(0)
        if set(self.opts).isdisjoint( set(['-help','--help','-?']) ) and len(self.opts)>0:
            return
        print ("SMC2 - Sequential Monte Carlo Sequentially Markovian Coalescent - Demographic Inference with Particle Filters")
        print ("       Donna Henderson, Sha (Joe) Zhu and Gerton Lunter\n")
        for idx, clas in enumerate(self.classes):
            print ("\n"+clas)
            for _, option, args, hlp in [opt for opt in self.options if opt[0] == idx]:
                print (" {:24} {}".format( option + " " + args, hlp ))
        sys.exit(0)
          

    def load_option_file(self):
        for idx,opt in enumerate(self.opts):
            if opt == "-@":
                fname = self.opts[idx+1]
                fopts = ' '.join([line
                                  for line in open(fname,'r').read().strip().split('\n')
                                  if not line.startswith('#')]).split()
                self.opts = self.opts[:idx] + fopts + self.opts[idx+2:]

                  
    def parse_opts(self):
        opts = self.opts
        idx, unparsed_opts, parsed_opts = 0, [], []
        while idx < len(opts):
            prev_idx = idx
            if opts[idx] == '-t':
                self.theta = float(opts[idx+1])
                idx += 2
            elif opts[idx] == '-r':
                self.rho = float(opts[idx+1])
                self.length = int(float(opts[idx+2]))
                idx += 3
            elif opts[idx] == '-seg':
                self.segfile = opts[idx+1]
                unparsed_opts += opts[idx:idx+2]   # also pass to smc^2
                idx += 2
            elif opts[idx] == '-o':
                self.outprefix = opts[idx+1]
                idx += 2
            elif opts[idx] == '-EM':
                self.emiters = int(opts[idx+1])
                idx += 2
            elif opts[idx] == '-chunks':
                self.chunks = int(opts[idx+1])
                idx += 2
            elif opts[idx] in ['-maxNE','-cap']:
                self.maxNE = float(opts[idx+1])
                idx += 2
            elif opts[idx] == '-startpos':
                self.startpos = int(opts[idx+1])
                idx += 2
            elif opts[idx] == '-no_infer_recomb':
                self.infer_recomb = False
                idx += 1
            elif opts[idx] == '-no_m_step':
                self.do_m_step = False
                idx += 1
            elif opts[idx] == '-alpha':
                self.alpha = float(opts[idx+1])
                idx += 2
            elif opts[idx] == '-smcsmcpath':
                self.smcsmcpath = opts[idx+1]
                idx += 2
            elif opts[idx] == '-log':
                self.logfile = opts[idx+1]
                idx += 2
            elif opts[idx] == '-nothreads':
                self.threads = False
                idx += 1
            elif opts[idx] in ["-sr","-M","-m","-es","-eps","-n","-g","-eg","-G","-eG","-T","-O","-L","-oSFS","-SC","-init"]:
                raise ValueError("SCRM option '{}' is illegal in SMCSMC; for options -m -M -n use alternatives -em -eM -en".format(opts[idx]))
            else:
                unparsed_opts.append( opts[idx] )
                idx += 1
                prev_idx = idx
            parsed_opts += opts[prev_idx : idx]
        self.opts = unparsed_opts
        if self.logfile == None:
            logger.addHandler( logging.StreamHandler() )
        else:
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            fh = logging.FileHandler( self.logfile )
            fh.setLevel( logging.DEBUG )
            fh.setFormatter(formatter)
            logger.addHandler( fh )
        logger.info("Finished parsing EM-related options")
        logger.info("Parsed these options: '{}'".format(' '.join(parsed_opts)))
        logger.info("Options passed to populationmodels: '{}'".format(' '.join(unparsed_opts)))

        if self.theta == None:     raise ValueError("Option -t required")
        if self.rho == None:       raise ValueError("Option -r required")
        if self.length == None:    raise ValueError("Option -r required")
        if self.segfile == None:   raise ValueError("Option -seg required")
        if self.outprefix == None: raise ValueError("Option -o required")
        if self.emiters < 0:       raise ValueError("-EM set to {}; must be >=0".format(self.emiters))
        if self.length < 1:        raise ValueError("Locus length ({}) must be >=1".format(self.length))
        if self.chunks < 1:        raise ValueError("Number of chunks ({}) must be >=1".format(self.chunks))
        if not os.path.exists(self.smcsmcpath):
            raise ValueError("Can't find executable "+self.smcsmcpath+"; check -smcsmcpath")
        logger.info("Running inference over {} bp in {} chunks and {} EM iterations".format(self.length, self.chunks, self.emiters+1))

        
    def prepare_folder(self, iteration):
        self.template = os.path.abspath( self.outprefix + "/emiter{}/chunk{}" )
        self.out_template = os.path.abspath( self.template + ".out" )
        self.log_template = os.path.abspath( self.template + ".log" )
        self.stdout_template = os.path.abspath( self.template + ".stdout" )
        self.stderr_template = os.path.abspath( self.template + ".stderr" )
        self.recomb_template = os.path.abspath( self.template + ".recomb.gz" )
        self.recomb_guide_template = os.path.abspath( self.template + ".recomb_guide.gz" )
        dirname = os.path.dirname( self.template.format(iteration, "dummy") )
        logger.info("Preparing directory {} for EM iteration {}".format( dirname, iteration ) )
        if not os.path.isdir( dirname ):
            try:
                os.makedirs( os.path.abspath( dirname ) )
            except:
                raise ValueError("Cannot find or create directory '{}'".format( os.path.abspath( dirname ) ) )
                                 
            
    def parse_outfile(self, outfname, data=None):
        of = open(outfname,'r')
        of.readline()
        if data == None: data = defaultdict(float)
        iters = set()
        for line in of:
            elts = line.strip().split()
            itr, epoch, frm, to = map(int, [elts[0],elts[1],elts[5],elts[6]])
            start, end, opp, count, rate, ne, ess = map(float, [elts[2],elts[3],elts[7],elts[8],elts[9],elts[10],elts[11]])
            typ = elts[4]
            iters.add(itr)
            if len(iters) > 1: raise ValueError("Found multiple iterations in .out file {}; expected only one".format(outfname))
            key = (typ,epoch,frm,to)
            data[(key,"Opp")] += opp
            data[(key,"Count")] += count
            data[(key,"Wt")] += max(0.0, (1.0/ess - 1e-10)) * opp
            data[(key,"Start")] = start
            data[(key,"End")] = end
        return data


    def write_outfile(self, outfname, data, iteration):
        of = open(outfname,'w')
        of.write("  Iter  Epoch       Start         End   Type   From     To         Opp       Count        Rate          Ne         ESS\n")
        for key in sorted(data.keys()):
            if key[1] == "Count":
                start, end, opp, count, wt = map(lambda label:data[(key[0],label)], ["Start","End","Opp","Count","Wt"])
                typ,epoch,frm,to = key[0]
                if typ == "LogL": opp, wt = 1.0, 1.0    # just add log likelihoods, don't average them
                rate = count / (opp + 1e-30)
                ne = (opp + 1e-10) / (2.0 * count + 1e-30) if typ == "Coal" else 0.0
                ess = 1.0 / (wt / (opp + 1e-30))
                of.write("{:6d} {:>6d} {:11.5g} {:11.5g} {:>6s}  {:>5d}  {:>5d} {:11.5g} {:11.5g} {:11.5g} {:11.5g} {:11.5g}\n".format(
                    iteration, epoch, start, end, typ, frm, to, opp, count, rate, ne, ess))
        of.close()


    def wait_for_outfile(self, iteration, chunk):
        sleeptime = 30
        while True:
            if os.path.exists( self.out_template.format(iteration,chunk) ):
                break
            if sleeptime == 30: logger.info("Waiting for chunk {} to appear...".format(chunk))
            time.sleep( sleeptime )
            sleeptime *= 1.02
        

    def merge_outfiles(self):
        of = open( os.path.abspath( self.outprefix + ".out" ), 'w' )
        of.write("  Iter  Epoch       Start         End   Type   From     To         Opp       Count        Rate          Ne         ESS\n")
        for em_iter in range(0, smcsmc.emiters+1):
            ifile = open( self.out_template.format(em_iter, "final"), 'r')
            ifile.readline()
            for line in ifile:
                of.write(line)
        of.close()
        
    
    def m_step(self, data):
        if not self.do_m_step:
            logger.info("Skipping M step -- re-using initial parameters")
            return
        # set population sizes
        for epoch in range(len(self.pop.change_points)):
            for pop in range(self.pop.num_populations):
                key = ("Coal",epoch,pop,-1)
                rate = data[ (key,"Count") ] / (data[ (key,"Opp") ] + 1e-30)
                popsize = min( self.maxNE / self.pop.N0, 1.0 / (2.0 * rate * self.pop.N0) )
                self.pop.population_sizes[ epoch ][ pop ] = popsize
                logger.info("Setting population size (epoch {} population {}) to {}".format(epoch, pop, popsize))
        # set migration rates
        for epoch in range(len(self.pop.change_points)):
            for frompop in range(self.pop.num_populations):
                for topop in range(self.pop.num_populations):
                    if frompop != topop:
                        key = ("Migr",epoch,frompop,topop)
                        rate = data[ (key,"Count") ] / data[ (key,"Opp") ]
                        # Note: pop.migration_rates are in units of 4Ne, while pop.mutation_rate and pop.recombination_rate are in
                        #       true units (events per nt per generation).  Reporting the 
                        self.pop.migration_rates[ epoch ][ frompop ][ topop ] = rate * 4 * self.pop.N0
                        logger.info("Setting migration rate (epoch {} population {}->{}) to {}  ({})".
                                    format(epoch, frompop, topop, rate, rate*4*self.pop.N0))
        # set recombination rate
        if self.infer_recomb:
            key = ("Recomb",-1,-1,-1)
            rate = data[ (key,"Count") ] / data[ (key,"Opp") ]
            self.pop.recombination_rate = rate
            logger.info("Setting recombination rate to {}".format(rate))

            
    def e_step(self, iteration, chunk):
        start = self.startpos + (self.length * chunk) // self.chunks
        end   = self.startpos + (self.length * (chunk+1)) // self.chunks
        self.pop.startpos = start
        self.pop.sequence_length = end - start
        self.pop.mutation_rate = self.theta / (4 * self.pop.N0 * self.length)
        self.pop.recombination_rate = self.rho / (4 * self.pop.N0 * self.length)
        command = [self.smcsmcpath, self.pop.core_command_line()]
        if not self.infer_recomb:
            command.append( "-xr 1-{}".format( len(self.pop.change_points) ) )
        if self.recombination_guide != None:
            command.append( "-guide {}".format( self.recombination_guide ) )
        command += self.smcsmc_opts
        command += ["-nsam", str(len(self.pop.sample_populations)),
                    "-startpos", str(self.pop.startpos),
                    "-EM 0",
                    "-log",
                    "-o",self.template.format(iteration,chunk),
                    ">",self.stdout_template.format(iteration,chunk),
                    "2>",self.stderr_template.format(iteration,chunk)]
        # final sanity check
        if (os.path.exists( self.stdout_template.format(iteration,chunk) ) or
            os.path.exists( self.log_template.format(iteration,chunk) ) or
            os.path.exists( self.stderr_template.format(iteration,chunk) ) ):
            raise ValueError("Found {} (or .log, or .stderr) -- it appears an inference process is already running!".
                             format(self.stdout_template.format(iteration,chunk)))
        # launch
        command = " ".join(command)
        logger.info("Running: "+command)
        if self.threads:
            p = subprocess.Popen(command, shell=True)
            self._processes.append(p)    # add to list of processes, so that we can wait for them to be done
        else:
            subprocess.call(command, shell=True)
            
                                 
    def do_iteration(self, iteration):
        # check if already done (combined .out file exists); if so finish
        self.prepare_folder( iteration )
        if os.path.exists( self.out_template.format(iteration,"final") ):
            logger.info("Results for iteration {} already exist -- continuing".format(iteration+1))
            if iteration == self.emiters:
                logger.info("All iterations were already present -- results not recomputed")
            return
        
        # build default population class.  Just need to set sequence_length, and we're good to go
        self.pop = populationmodels.Population()
        self.pop.mutation_rate = self.theta / (4 * self.pop.N0)
        self.pop.recombination_rate = self.rho / (4 * self.pop.N0)
        self.smcsmc_opts = self.pop.parse_command_line( " ".join(self.opts) )
        logger.info("Options passed to smcsmc: '{}'".format(' '.join(self.smcsmc_opts)))

        # parse previous .out file, and use previous .recomb.gz file
        self.recombination_guide = None
        if iteration > 0:
            self.m_step( self.parse_outfile( self.out_template.format(iteration-1,"final" ) ) )
            if self.alpha > 0.0:
                # process the .recomb files to recomb-guideing files
                logger.info("Processing recombination guide files...")
                for chunk in range(self.chunks):
                    lr = processrecombination.LocalRecombination( self.recomb_template.format(iteration-1,chunk) )
                    lr.smooth( self.alpha, self.beta )
                    self.recombination_guide = self.recomb_guide_template.format(iteration, chunk)
                    fout = gzip.open( self.recombination_guide, 'w' )
                    lr.write_data( fout )
                    fout.close()

        # launch the chunks (either return immediately, or actually process)
        # also check if already in progress, and raise error if so
        for chunk in range(self.chunks):
            self.e_step( iteration, chunk )

        # wait for child processes to finish
        while len(self._processes)>0:
            logger.info("Waiting for child process {} to finish...".format(len(self._processes)))
            returncode = self._processes[-1].wait()
            if returncode != 0:
                raise ValueError("smcsmc executable for iteration {}, chunk {}, returned non-zero exit code".format(
                    iteration+1,len(self._processes)-1))
            self._processes = self._processes[:-1]

        # wait for all .out files to appear (used if e_step submits jobs to sge)
        for chunk in range(self.chunks):
            self.wait_for_outfile(iteration, chunk)
        logger.info("Found all chunks; calculating sufficient statistics")

        # combine the .out files
        data = None
        for chunk in range(self.chunks):
            data = self.parse_outfile( self.out_template.format(iteration,chunk), data )
        self.write_outfile( self.out_template.format(iteration,"final"), data, iteration )

        
#
# main
#

smcsmc = Smcsmc( sys.argv[1:] )
smcsmc.print_help_and_exit()
smcsmc.load_option_file()
smcsmc.parse_opts()
for em_iter in range(0, smcsmc.emiters+1):
    smcsmc.do_iteration(em_iter)
smcsmc.merge_outfiles()

#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import logging
import gzip
import subprocess
import time
import math
from collections import defaultdict

import populationmodels
import processrecombination
import execute

#
# pattern + start/end time (generations)
#


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
        self.mutation_rate = None  # per nt per gen
        self.recombination_rate = None
        self.theta = None          # 4 N0 mutation_rate L
        self.rho = None            # 4 N0 recombination_rate L
        self.N0 = None
        self.trange = None
        self.length = None
        self.startpos = 1
        self.segfile = None
        self.watterson_estimator = None
        self.outprefix = None
        self.logfile = None
        self.patt_args = None
        self.segfiles = []
        self.emiters = 0
        self.chunks = 1
        self.infer_recomb = True
        self.do_m_step = True
        self.alpha = 0.0            # posterior mix-in; 0 means use prior; <0 means remove .recomb.gz files
        self.beta = 4               # smoothness parameter; see processrecombination.py
        self.maxNE = 1e99
        self._processes = []
        self.threads = True
        self.cluster = False
        self.cconfig = None
        self.maxgap = 200000
        self.minseg = 500000
        self.samples_in_segfile = None
        
        self.smcsmcpath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../smcsmc'))
        
        # define command line options, in 3 classes:
        self.classes = ["Population-model related options (fixed)",
                        "Population-model related options (inferred; initial values)",
                        "Inference-related options",
                        "Options related to the EM process"]

        self.options = [
            (0, '-nsam', 'n',   '[*] Set the total number of samples to n'),
            (0, '-N0', 'N',     '[+] Set (unscaled) initial population size to N'),
            (0, '-mu', 's',     '[+] Set mutation rate (per nucleotide per generation)'),
            (0, '-t', 'th',     '[+] Set mutation rate (expected number of mutations) for the locus (=4 N0 mu L )'),
            (0, '-length', 'L', '[+] Set locus length (nucleotides)'),
            (0, '-I', 'n s1..sn','Use an n-population island model with si individuals sampled'),
            (0, '-ej', 't i j', 'Speciation event at t*4N0; creates population i from population j in forward direction'),
            (0, '-eI','t s1..sn','Sample s1..sn indiviuals from their populations at time t*4N0 generations'),

            (1, '-rho', 'rh',  '[+] Set recombination rate (per nucleotide per generation)'),
            (1, '-r', 'r L',   '[+] Set initial per locus recombination rate (=4 N0 L rho) and locus length (L)'),
            (1, '-eN', 't n',  'Change the size of all populations to n*N0 at time t*4N0'),
            (1, '-en', 't i n','Change the size of population i to n*N0 at time t*4N0'),
            (1, '-eM', 't m',    'Change the symmetric backward migration rate to m/(npop-1) at time t*4N0'),
            (1, '-em', 't i j m', 'Change the backward migration rate, from population i to population j, to m/(npop-1) at time t*4N0'),
            (1, '-ema', 't s11 s12 ...', 'Set backward migration rate matrix at time t*4N0'),
            
            (2, '-o', 'f',     '[*] Output prefix'),
            (2, '-seg', 'f',   '[+] Input .seg file'),
            (2, '-segs', 'f1 f2 ...', '[+] Input .seg files (will be merged into a single .seg file'),
            (2, '-maxgap', 'n', 'Split .seg files over gaps larger than maxgap (500 kb)'),
            (2, '-minseg', 'n', 'After splitting ignore segments shorter than minseg (500 kb)'),
            (2, '-startpos', 'x', 'First locus to process (1)'),
            (2, '-P', 's e p', 'Divide time interval [s,e] (in generations) equally on log scale, using pattern p'),
            (2, '-Np', 'n',    'Number of particles'),
            (2, '-seed', 's',  'Random number seed'),
            (2, '-ancestral_aware', '', 'Assume that haplotype 0 is ancestral'),
            (2, '-bias_heights', 't0..tn', 'Set recombination bias times to h0..hn * 4N0'),
            (2, '-bias_strengths','s1..sn','Set recombination bias strenghts'),
            (2, '-calibrate_lag','s','Explanation required!'),
            (2, '-apf', 'b',    'Auxiliary particle filter: none (0), singletons (1), cherries (2)  (0)'),
            
            (3, '-EM', 'n',     'Number of EM iterations-1 ({})'.format(self.emiters)),
            (3, '-cap', 'n',    'Set (unscaled) upper bound on effective population size'),
            (3, '-chunks','n',  'Number of chunks computed in parallel ({})'.format(self.chunks)),
            (3, '-no_infer_recomb', '', 'Do not infer recombination rate'),
            (3, '-no_m_step','','Do not update parameters (but do infer recombination guide)'),
            (3, '-alpha', 't',  'Fraction of posterior recombination to mix in to recombination guide ({}); negative removes files'.format(self.alpha)),
            (3, '-c', '',       'Use qsub to submit job(s) to cluster (overrides use of threads)'),
            (3, '-C', 'opts',   'Qsub config parameter(s) e.g. "-P project.prj -q long.q"; overrides ./qsub.conf'),
            (3, '-nothreads','','Calculate chunks sequentially rather than in parallel threads'),
            (3, '-smcsmcpath', 'f','Path to smcsmc executable ({})'.format(self.smcsmcpath)),
            (3, '-log', 'f',    'Log file (stdout)'),
            (3, '-@', 'f',      'File with line-, space-, or tab-separated options'),
            (3, '-v', '',       'Show version of underlying smcsmc and scrm implementations'),  
            (3, "-help", '',    'This help')]

        self.opts = opts
        
        
    def print_help_and_exit(self):
        if len(self.opts) == 1 and self.opts[0] == "-v":
            subprocess.check_call(" ".join([self.smcsmcpath] + self.opts), shell=True)
            sys.exit(0)
        if set(self.opts).isdisjoint( set(['-help','--help','-h', '-?']) ) and len(self.opts)>0:
            return
        print ("SMC2 - Sequential Monte Carlo Sequentially Markovian Coalescent - demographic inference with particle filters")
        print ("       Donna Henderson, Sha (Joe) Zhu and Gerton Lunter\n")
        for idx, clas in enumerate(self.classes):
            print ("\n"+clas)
            for _, option, args, hlp in [opt for opt in self.options if opt[0] == idx]:
                print (" {:24} {}".format( option + " " + args, hlp ))
        print ("\n([*], required; [+], optional but one of group is required)")
                
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
            elif opts[idx] == '-mu':
                self.mutation_rate = float(opts[idx+1])
                idx += 2
            elif opts[idx] == '-r':
                self.rho = float(opts[idx+1])
                self.length = int(float(opts[idx+2]))
                idx += 3
            elif opts[idx] == '-rho':
                self.recombination_rate = float(opts[idx+1])
                idx += 2
            elif opts[idx] == '-length':
                self.length = int(float(opts[idx+1]))
                idx += 2
            elif opts[idx] == '-N0':
                self.N0 = float(opts[idx+1])
                unparsed_opts += opts[idx:idx+2]   # also pass to smc^2
                idx += 2
            elif opts[idx] == '-seg':
                self.segfile = opts[idx+1]
                unparsed_opts += opts[idx:idx+2]   # also pass to smc^2
                idx += 2
            elif opts[idx] == '-segs':
                idx, self.segfile = self.parse_segfiles( idx, opts )
                unparsed_opts += ['-seg',self.segfile]
            elif opts[idx] == '-maxgap':
                self.maxgap = int(float(opts[idx+1]))
                idx += 2
            elif opts[idx] == '-minseg':
                self.minseg = int(float(opts[idx+1]))
                idx += 2
            elif opts[idx] == '-o':
                self.outprefix = opts[idx+1]
                idx += 2
            elif opts[idx] == '-EM':
                self.emiters = int(opts[idx+1])
                idx += 2
            elif opts[idx] == '-P':
                self.patt_args = opts[idx+1:idx+4]
                idx += 4
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
            elif opts[idx] == '-c':
                self.cluster = True
                idx += 1
            elif opts[idx] == '-C':
                self.cconfig = opts[idx+1]
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
        self.parsed_opts = parsed_opts
        if self.cluster:
            if self.cconfig is None:
                try:
                    self.cconfig = open("./qsub.conf").read().replace('\n',' ')
                except:
                    self.cconfig = ""
            execute.qsub_config = self.cconfig   # assign to global variable within module 'execute'


    def set_pattern(self):
        if self.patt_args is None: return
        if len(self.patt_args) != 3: raise ValueError("-P should get 3 arguments (start and end generation; and pattern)")
        start, end = map(float, self.patt_args[:2])
        patt = self.patt_args[2]
        # extract all timed options to smcsmc, and separate them from all others (including -I)
        opts = self.opts
        opt_indices = [i for i,o in enumerate(opts + ["-"]) if o.startswith('-')]
        time_opts, remain_opts = [], []
        for ii, idx in enumerate(opt_indices[:-1]):
            if opts[idx] in ('-eI','-ej','-eM','-ema','-em','-eN','-en'):
                time_opts.append( (float(opts[idx+1]), opts[idx : opt_indices[ii+1]]))
            else:
                remain_opts += opts[idx : opt_indices[ii+1] ]
        # parse patt and generate times
        mask = [1]
        try:
            for pat in patt.split('+'):
                a, b = map(int, pat.split('*'))
                masklet = [1] + ([0]*(b-1))
                mask += masklet * a
        except:
            raise ValueError("Problem parsing pattern '{}'".format(patt))
        times = [0] + [start * math.exp( math.log(end/start) * i / (len(mask)-2.0)) / (4 * self.N0)
                       for i in range(len(mask))]
        # generate -eN commands
        new_time_opts = []
        for idx, time in enumerate(times[:-1]):
            if mask[idx] == 1:
                # add -eN command setting the population size at this epoch to the last set population size, or 1.0 if none was set
                last_en_opt = sorted( [(-1.0, ['-eN', '0.0', '1.0'])] + [ to for to in time_opts if to[0] <= time and to[1][0] == '-eN'] )[-1]
                new_time_opts.append( (time, ['-eN {} {}'.format(time, last_en_opt[1][2])]) )
        # set times of other commands to the appropriate epoch time
        # remove -eN commands line opts.  Note: -en commands are left in, and might introduce new epochs; best not to use -en and -P simultaneously
        for to in time_opts:
            if to[1][0] == '-eN': continue
            newtime = sorted( [t for t in times if t <= to[0]] )[-1]
            new_time_opts.append( (newtime, [to[1][0], str(newtime)] + to[1][2:]) )
        new_time_opts.sort()
        self.opts = remain_opts + [' '.join(opt) for t,opt in new_time_opts]
        logger.info("Pattern mask: {}".format(' '.join(map(str,mask))))
        logger.info("Population structure options: {}".format( ' '.join( [' '.join(opt) for t,opt in new_time_opts] ) ) )


    def set_environment(self):
        # workaround for qrsh -V not passing through LD_LIBRARY_PATH
        if "LLP" in os.environ:
            os.environ["LD_LIBRARY_PATH"] = os.environ["LLP"]
            

    def open_file(self, file, mode):
        if file.endswith('.gz'):
            return gzip.open(file,mode)
        return open(file,mode)


    def read_segline(self, line):
        elts = line.strip().split('\t')
        if len(elts) == 3:
            return elts
        elif len(elts) == 6:
            return [elts[0],elts[1],elts[5]]
        raise ValueError("Error parsing a .seg line (expected 3 or 6 fields): '{}'".format(' '.join(elts)))


    def define_chunks(self):
        last = self.startpos
        chunkstarts = [last]
        weighted_length = 1e-10
        segregating_sites = 0
        harmonic = None
        for line in self.open_file(self.segfile, 'r'):
            elts = self.read_segline(line)
            start = int(elts[0])
            size = int(elts[1])
            alleles = len(elts[2])
            missing = elts[2].count('.')
            if start + size < self.startpos:
                continue
            if self.length != None and start > self.startpos + self.length:
                last = start
                break
            if missing == alleles:
                continue
            if harmonic is None:
                harmonic = [ sum(1.0/(i) for i in range(1,k+1)) for k in range(alleles)]
            weighted_length += size * harmonic[ alleles - missing - 1 ]
            segregating_sites += (elts[2].count('0') > 0 and elts[2].count('1') > 0)
            if start - last >= self.maxgap:
                chunkstarts.append( last )
                chunkstarts.append( start )
            last = start + size
        chunkstarts.append( last )
        if self.length is None:
            self.length = last - self.startpos
        else:
            self.minseg = min(self.minseg, self.length)
            if self.length > last - self.startpos:
                if self.length > last - self.startpos + 1000000:
                    raise ValueError("Asked to process {} bp, but .seg file(s) provided only cover positions {}-{} ({} bp)".format(
                        self.length, self.startpos, last, last - self.startpos))
                else:
                    logger.info("Warning: asked to process {} bp, but .seg file(s) provided only cover positions {}-{} ({} bp)".format(
                        self.length, self.startpos, last, last - self.startpos))
        smallchunks = [i for i in range(len(chunkstarts)//2) if chunkstarts[2*i+1]-chunkstarts[2*i] < self.minseg]
        smallchunklen = sum( chunkstarts[2*i+1]-chunkstarts[2*i] for i in smallchunks )
        for i in reversed( smallchunks ):
            chunkstarts = chunkstarts[:2*i] + chunkstarts[2*(i+1):]
        chunkcounts = [1] * (len(chunkstarts) // 2)
        self.watterson_estimator = segregating_sites / weighted_length
        logger.info("Found {} segments in data separated by gaps larger than {} bp (including chromosome boundaries)".format(len(chunkcounts), self.maxgap))
        logger.info("Skipped {} small segments covering {} bp".format(len(smallchunks), smallchunklen))
        if len(chunkcounts) > self.chunks:
            logger.info("Note: using more chunks ({}) than asked for ({}); to reduce, decrease the number of chromosomes and/or increase -maxgap".format(
                    len(chunkcounts), self.chunks))
        if len(chunkcounts) == 0:
            raise ValueError("No segments left - nothing to do...")
        # iteratively increase the number of pieces for the chunk with the current largest piece size
        while sum(chunkcounts) < self.chunks:
            chunksizes = [ (chunkstarts[2*i+1]-chunkstarts[2*i])/c for i,c in enumerate(chunkcounts) ]
            largest = [ i for i,c in enumerate(chunksizes) if c == max(chunksizes) ][0]
            chunkcounts[ largest ] += 1
        # calculate start/end positions of chunks
        self.chunkstarts = []
        self.chunkends = []
        for i,c in enumerate(chunkcounts):
            for ii in range(c):
                chunksize = chunkstarts[2*i+1] - chunkstarts[2*i]
                self.chunkstarts.append( chunkstarts[2*i] + (chunksize * ii) // c )
                self.chunkends.append( chunkstarts[2*i] + (chunksize * (ii+1)) // c )
                logger.info(" Chunk {}: start {} length {}".format(len(self.chunkstarts)-1,
                                                                   self.chunkstarts[-1],
                                                                   self.chunkends[-1]-self.chunkstarts[-1]))


    def validate(self):
        if self.segfile == None:   raise ValueError("Option -seg or -segs required")
        if self.outprefix == None: raise ValueError("Option -o required")
        if not os.path.exists(self.smcsmcpath):
            raise ValueError("Can't find executable "+self.smcsmcpath+"; check -smcsmcpath")
        smcsmc.prepare_folder(0, log=False)

        if self.logfile == None:
            self.logfile = "{}/result.log".format( self.outprefix )
            #logger.addHandler( logging.StreamHandler() )
        formatter = logging.Formatter('%(asctime)s - %(message)s')
        fh = logging.FileHandler( self.logfile )
        fh.setLevel( logging.DEBUG )
        fh.setFormatter(formatter)
        logger.addHandler( fh )
        logger.info("Parsed these options: '{}'".format(' '.join(self.parsed_opts)))
        logger.info("Options passed to populationmodels: '{}'".format(' '.join(self.opts)))



    def validate_parameters(self):
        # infer mutation rate from -t and N0
        if self.mutation_rate == None:
            if self.theta != None and self.N0 != None:
                logger.info("Setting mutation rate (per nt per generation) from -t and N0")
                self.mutation_rate = self.theta / (4 * self.N0 * self.length)

        # infer recombination rate from -r and N0
        if self.recombination_rate == None:
            if self.rho != None and self.N0 != None:
                logger.info("Setting recombination rate (per nt per generation) from -t and N0")
                self.recombination_rate = self.rho / (4 * self.N0 * self.length)

        # estimate N0 (for initial population size, and for scaling time) from watterson's estimate of theta
        if self.N0 == None:
            if self.mutation_rate != None:
                logger.info("Setting N0 from mutation rate and Watterson's estimate of theta")
                self.N0 = self.watterson_estimator / (4 * self.mutation_rate)

        # set -t and -r parameters, if not provided
        if self.theta == None:
            if self.N0 != None:
                logger.info("Setting -t parameter from N0, mutation rate, and length")
                self.theta = 4 * self.N0 * self.mutation_rate * self.length
        if self.rho == None:
            if self.N0 != None:
                logger.info("Setting -r parameter from N0, recombination rate, and length")
                self.rho = 4 * self.N0 * self.recombination_rate * self.length

        if self.N0 == None: raise ValueError("N0 required -- use -N0, or (implicitly) -mu")
        if self.mutation_rate == None: raise ValueError("Mutation rate required -- use -mu, or (implicitly) -t / -N0")
        if self.recombination_rate == None: raise ValueError("Recombination rate required -- use -rho, or (implicitly) -r / -N0")
        if self.theta == None: raise ValueError("-t value required -- use -t or (implicitly) -mu / -N0")
        if self.rho == None: raise ValueError("-r value required -- use -r or (implicitly) -rho / -N0")
        if self.length == None: raise ValueError("Option -r required")
        
        # ensure parameters are consistent
        if abs( 4*self.N0*self.mutation_rate*self.length / self.theta - 1.0 ) > 1e-5:
            raise ValueError("Values for -t and -N0, -mu, -length not consistent")
        if abs( 4*self.N0*self.recombination_rate*self.length / self.rho - 1.0 ) > 1e-5:
            raise ValueError("Values for -r and -N0, -rho, -length not consistent")

        if self.emiters < 0:       raise ValueError("-EM set to {}; must be >=0".format(self.emiters))
        if self.length < 1:        raise ValueError("Locus length ({}) must be >=1".format(self.length))
        if self.chunks < 1:        raise ValueError("Number of chunks ({}) must be >=1".format(self.chunks))
        logger.info("Parameter value:-")
        logger.info("  mutation rate:        {}".format(self.mutation_rate))
        logger.info("  recombination rate:   {}".format(self.recombination_rate))
        logger.info("  N0:                   {}".format(self.N0))
        logger.info("  theta (4 N0 mu):      {}".format(4 * self.N0 * self.mutation_rate))
        logger.info("  Watterson's estimate: {}".format(self.watterson_estimator))
        logger.info("  locus size:           {}".format(self.length))
        logger.info("Running inference over {} bp in {} chunks and {} EM iterations".format(self.length, len(self.chunkstarts), self.emiters+1))
        logger.info("Unique chunk sizes (ordered): {}".format( ' '.join( map(str, sorted(list(set((self.chunkends[i]-start for i, start in enumerate(self.chunkstarts)))))))))
    

    def parse_segfiles(self, idx0, files):
        if self.segfile != None: raise ValueError("You can only provide the --seg or --segs option once!")
        segfiles = []
        for idx in range(idx0 + 1, len(files) + 1):
            if idx == len(files): break
            fname = files[idx]
            if fname.startswith('-'): break
            segfiles.append( fname )
        self.segfiles = segfiles
        self.mergedseg_name = "{}/merged.seg".format(self.outprefix)
        self.mergedmap_name = "{}/merged.map".format(self.outprefix)
        return idx, self.mergedseg_name


    def process_segfiles(self):
        if len(self.segfiles) == 0: return
        fout_name = self.mergedseg_name
        fmap_name = self.mergedmap_name
        fout = open(fout_name,'w')
        chrom_map = []
        start = 1
        for idx, fname in enumerate(self.segfiles):
            logger.info("Merging .seg files -- processing {}...".format(fname))
            name = fname[:-3] if fname.endswith('.gz') else fname
            name = name[:-4] if name.endswith('.seg') else name
            name = name.replace('/','.')
            chrom_map.append( (name, start) )
            for line in self.open_file(fname,'r'):
                elts = self.read_segline(line)
                pos = int(elts[0])
                fout.write("{}\t{}\t{}\n".format( start+pos-1, elts[1], elts[2] ) )
            abs_pos = start + pos + int(elts[1]) - 1
            new_abs_pos = abs_pos + self.maxgap
            new_abs_pos = 1000000 * ((new_abs_pos // 1000000) + 1) + 1
            gap = new_abs_pos - abs_pos
            if idx < len(self.segfiles)-1:
                ## add gap
                fout.write("{}\t{}\t{}\n".format( abs_pos, gap, "."*len(elts[2]) ) )
                start = new_abs_pos
        fout.close()
        fout = self.open_file(fmap_name,'w')
        for name, start in chrom_map:
            fout.write("{}\t{}\n".format(name,start))
        fout.close()

        
    def prepare_folder(self, iteration, log=True):
        self.template = os.path.abspath( self.outprefix + "/emiter{}/chunk{}" )
        self.out_template = os.path.abspath( self.template + ".out" )
        self.log_template = os.path.abspath( self.template + ".log" )
        self.stdout_template = os.path.abspath( self.template + ".stdout" )
        self.stderr_template = os.path.abspath( self.template + ".stderr" )
        self.recomb_template = os.path.abspath( self.template + ".recomb.gz" )
        self.recomb_guide_template = os.path.abspath( self.template + ".recomb_guide.gz" )
        dirname = os.path.dirname( self.template.format(iteration, "dummy") )
        if log: logger.info("Preparing directory {} for EM iteration {}".format( dirname, iteration ) )
        if not os.path.isdir( dirname ):
            try:
                os.makedirs( os.path.abspath( dirname ) )
            except:
                raise ValueError("Cannot find or create directory '{}'".format( os.path.abspath( dirname ) ) )
                                 
            
    def parse_outfile(self, outfname, data=None, clump=None):
        of = open(outfname,'r')
        if data == None: data = defaultdict(float)
        header = { idx:elt for idx, elt in enumerate( of.readline().strip().split() ) }
        iters = set()
        for line in of:
            elts = { header[idx] : elt for idx, elt in enumerate(line.strip().split()) }
            for name, val in list(elts.iteritems()):  # Leave "Type" as String
                if name in ["Iter","Epoch","From","To","Clump"]:
                    elts[name] = int(val)
                elif name in ["Start","End","Opp","Count","Rate","Ne","ESS","InitRate"]:
                    elts[name] = float(val)
            iters.add( elts["Iter"] )
            if len(iters) > 1: raise ValueError("Found multiple iterations in .out file {}; expected only one".format(outfname))
            keys = [(elts["Type"],elts["Epoch"],elts["From"],elts["To"], -1)]
            if "Clump" in elts:
                clump = elts["Clump"]
            if clump != None:
                keys.append((elts["Type"],elts["Epoch"],elts["From"],elts["To"], clump))
            for key in keys:
                data[(key,"Opp")] += elts["Opp"]
                data[(key,"Count")] += elts["Count"]
                data[(key,"Wt")] += max(0.0, (1.0/elts["ESS"] - 1e-10)) * elts["Opp"]
                data[(key,"Start")] = elts["Start"]
                data[(key,"End")] = elts["End"]
                if "InitRate" in elts:
                    data[(key,"InitRate")] = elts["InitRate"]
        return data


    def write_outfile(self, outfname, data, iteration):
        of = open(outfname,'w')
        of.write("  Iter  Epoch       Start         End   Type   From     To            Opp          Count           Rate             Ne         ESS  Clump\n")
        for key in sorted(data.keys(), key = lambda elt : (elt[0][-1]>=0,elt) ):
            if key[1] == "Count":
                start, end, opp, count, wt = map(lambda label:data[(key[0],label)], ["Start","End","Opp","Count","Wt"])
                typ,epoch,frm,to,clump = key[0]
                if typ == "LogL": opp, wt = 1.0, 1.0    # just add log likelihoods, don't average them
                rate = count / (opp + 1e-30)
                ne = (opp + 1e-10) / (2.0 * count + 1e-30) if typ == "Coal" else 0.0
                ess = 1.0 / (wt / (opp + 1e-30))
                of.write("{:6d} {:>6d} {:11.5g} {:11.5g} {:>6s}  {:>5d}  {:>5d} {:14.8g} {:14.8g} {:14.8g} {:14.8g} {:11.5g} {:>6d}\n".format(
                    iteration, epoch, start, end, typ, frm, to, opp, count, rate, ne, ess, clump))
        of.close()


    def have_outfile(self, iteration, chunk):
        if not os.path.exists( self.out_template.format(iteration,chunk) ):
            return False
        try:
            data = self.parse_outfile( self.out_template.format(iteration,chunk) )
            for key in data:
                if key[0][0] == "LogL": 
                    return True
        except:
            pass
        return False


    def wait_for_outfile(self, iteration, chunk):
        sleeptime = 30
        while True:
            if self.have_outfile( iteration, chunk ):
                logger.info("Found chunk {}".format(chunk))
                break
            if sleeptime == 30: logger.info("Waiting for chunk {} to appear...".format(chunk))
            time.sleep( sleeptime )
            sleeptime *= 1.02


    def merge_outfiles(self):
        of = open( os.path.abspath( self.outprefix + "/result.out" ), 'w' )
        write_header = True
        for em_iter in range(smcsmc.emiters, -1, -1):
            ifile = open( self.out_template.format(em_iter, "final"), 'r')
            header = ifile.readline()
            if write_header:
                of.write(header)
                write_header = False
            for line in ifile:
                elts = line.split()
                if elts[-1] == "-1":  # only pass through aggregate data
                    of.write(line)
        of.close()
        
    
    def m_step(self, data):
        if not self.do_m_step:
            logger.info("Skipping M step -- re-using initial parameters")
            return
        # set population sizes
        for epoch in range(len(self.pop.change_points)):
            for pop in range(self.pop.num_populations):
                key = ("Coal",epoch,pop,-1,-1)
                rate = data[ (key,"Count") ] / (data[ (key,"Opp") ] + 1e-30)
                popsize = min( self.maxNE / self.pop.N0, 1.0 / (2.0 * rate * self.pop.N0) )
                self.pop.population_sizes[ epoch ][ pop ] = popsize
                logger.info("Setting population size (epoch {} population {}) to {}".format(epoch, pop, popsize))
        # set migration rates
        for epoch in range(len(self.pop.change_points)):
            for frompop in range(self.pop.num_populations):
                for topop in range(self.pop.num_populations):
                    if frompop != topop:
                        key = ("Migr",epoch,frompop,topop,-1)
                        rate = data[ (key,"Count") ] / data[ (key,"Opp") ]
                        # Note: pop.migration_rates are in units of 4Ne, while pop.mutation_rate and pop.recombination_rate are in
                        #       true units (events per nt per generation).  Reporting the 
                        self.pop.migration_rates[ epoch ][ frompop ][ topop ] = rate * 4 * self.pop.N0
                        logger.info("Setting migration rate (epoch {} population {}->{}) to {}  ({})".
                                    format(epoch, frompop, topop, rate, rate*4*self.pop.N0))
        # set recombination rate
        if self.infer_recomb:
            key = ("Recomb",-1,-1,-1,-1)
            rate = data[ (key,"Count") ] / data[ (key,"Opp") ]
            self.pop.recombination_rate = rate
            logger.info("Setting recombination rate to {}".format(rate))

            
    def e_step(self, iteration, chunk):
        start = self.chunkstarts[chunk]
        end   = self.chunkends[chunk]
        self.pop.startpos = start
        self.pop.sequence_length = end - start
        self.pop.mutation_rate = self.theta / (4 * self.pop.N0 * self.length)
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
        # launch.  Note, no check is made to see if another process is already writing to 
        # the output files (e.g. self.stdout_template.format(iteration,chunk)).  This is
        # hard to do on a cluster system, so it's up to the user to make sure they don't
        # launch processes writing to the same file in parallel!
        # However do check if a complete .out file is already present, and do not overwrite if so
        if self.have_outfile(iteration, chunk):
            logger.info("Results for iteration {} chunk {} already present -- continuing".format(iteration+1,chunk))
            return                       # data is complete
        command = " ".join(command)
        logger.info("Running: "+command)
        if self.threads and not self.cluster:
            p = subprocess.Popen(command, shell=True)
            self._processes.append(p)    # add to list of processes, so that we can wait for them to be done
        else:
            execute.submit(command, shell=True)
            
                                 
    def do_iteration(self, iteration):
        # check if already done (combined .out file exists); if so finish
        self.prepare_folder( iteration )
        if self.have_outfile( iteration, "final" ):
            logger.info("Results for iteration {} already exist -- continuing".format(iteration+1))
            if iteration == self.emiters:
                logger.info("All iterations were already present -- results not recomputed")
            return
        
        # build default population class.  Just need to set sequence_length, and we're good to go
        self.pop = populationmodels.Population()
        self.pop.recombination_rate = self.recombination_rate
        self.pop.mutation_rate = self.mutation_rate
        self.smcsmc_opts = self.pop.parse_command_line( " ".join(self.opts) )
        logger.info("Options passed to smcsmc: '{}'".format(' '.join(self.smcsmc_opts)))

        # parse previous .out file, and use previous .recomb.gz file
        self.recombination_guide = None
        if iteration > 0:
            self.m_step( self.parse_outfile( self.out_template.format(iteration-1,"final" ) ) )
            if self.alpha > 0.0:
                # process the .recomb files to recomb-guideing files
                logger.info("Processing recombination guide files...")
                for chunk in range(len(self.chunkstarts)):
                    lr = processrecombination.LocalRecombination( self.recomb_template.format(iteration-1,chunk) )
                    lr.smooth( self.alpha, self.beta )
                    self.recombination_guide = self.recomb_guide_template.format(iteration, chunk)
                    fout = gzip.open( self.recombination_guide, 'w' )
                    lr.write_data( fout )
                    fout.close()
            elif self.alpha < 0.0:
                # remove the (large) .recomb files
                for chunk in range(len(self.chunkstarts)):
                    try:
                        os.remove( self.recomb_template.format(iteration-1,chunk) )
                    except:
                        pass

        # launch the chunks (either return immediately, or actually process)
        # also check if already in progress, and raise error if so
        for chunk in range(len(self.chunkstarts)):
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
        for chunk in range(len(self.chunkstarts)):
            self.wait_for_outfile(iteration, chunk)
        logger.info("Found all chunks; calculating sufficient statistics")

        # combine the .out files
        data = None
        for chunk in range(len(self.chunkstarts)):
            data = self.parse_outfile( self.out_template.format(iteration,chunk), data, chunk )
        self.write_outfile( self.out_template.format(iteration,"final"), data, iteration )

        
#
# main
#

smcsmc = Smcsmc( sys.argv[1:] )
smcsmc.print_help_and_exit()
smcsmc.load_option_file()
smcsmc.parse_opts()
smcsmc.validate()
smcsmc.process_segfiles()
smcsmc.set_environment()
smcsmc.define_chunks()
smcsmc.validate_parameters()
smcsmc.set_pattern()
for em_iter in range(0, smcsmc.emiters+1):
    smcsmc.do_iteration(em_iter)
smcsmc.merge_outfiles()

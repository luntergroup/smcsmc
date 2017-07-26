#!/usr/bin/env python
from __future__ import print_function

import sys
import os
import tempfile
import subprocess
import itertools

from optparse import OptionParser

"""
Population structures
"""


defaults = {'N0':                 10000,
            'mutation_rate':      2.5e-8,
            'recombination_rate': 1e-8,
            'sequence_length':    1e6,
            'num_samples':        2,
            'scrmpath': 'scrm'}

class Population:

    def __init__(self,
                 N0                   = defaults['N0'],
                 mutation_rate        = defaults['mutation_rate'],
                 recombination_rate   = defaults['recombination_rate'],
                 sequence_length      = defaults['sequence_length'],
                 num_samples          = defaults['num_samples'],
                 change_points        = [0, .01, 0.06, 0.2, 1, 2],
                 population_sizes     = [[1], [0.1], [1], [0.5], [1], [2]],
                 num_populations      = 1,
                 migration_rates      = None,
                 sample_populations   = None,
                 sample_times         = None,
                 migration_commands   = None,
                 split_command        = "",
                 seed                 = (1,),
                 filename             = None,
                 scrmpath             = defaults['scrmpath'] ):

        self.N0 = N0
        self.mutation_rate = mutation_rate
        self.recombination_rate = recombination_rate
        self.startpos = 1
        self.sequence_length = sequence_length
        self.num_samples = num_samples
        self.change_points = change_points              # list of times, in 4Ne units
        self.population_sizes = population_sizes        # list, or list of lists (per population)
        self.migration_rates = migration_rates          # None, or list of lists of lists
        self.num_populations = num_populations
        self.sample_populations = sample_populations    # encoded as sample_populations_ in model.h
        self.sample_times = sample_times                # encoded as sample_times_ in scrm/model.h
        self.migration_commands = migration_commands    # to implement speciation events
        self.split_command = split_command
        self.seed = seed
        self.filename = filename
        self.scrmpath = scrmpath

    def parse_command_line(self, cmdline):
        """ Parses options relating to population model: -nsam -em -eM -en -eN -eI -I -ej
            Unrecognized options are returned unchanged
            Note: need to set mutation_rate, recombination_rate, sequence_length separately! 
            TODO: add -ej (split) """
        self.num_samples = -1
        self.change_points = []
        self.population_sizes = []
        self.migration_rates = []
        self.sample_populations = []
        self.sample_times = []
        self.num_populations = 1
        opts, unparsed_opts = cmdline.split(), []
        idx = 0
        while idx < len(opts):
            if opts[idx] == "-N0":
                self.N0 = int(opts[idx+1])
                idx += 2
            elif opts[idx] == "-nsam":
                self.num_samples = int(opts[idx+1])
                idx += 2
            elif opts[idx] == "-I":
                self.num_populations = int(opts[idx+1])
                if len(self.change_points) != 0: raise ValueError("Options -I must go before -eN / -en / -eM / -em / -ema")
                self._parse_time( 0, ["0.0"] )
                idx = self._parse_sample( 0, idx+2, opts )
            elif opts[idx] == "-eI":
                idx = self._parse_time( idx+1, opts )
                idx = self._parse_sample( len(self.change_points)-1, idx, opts )
            elif opts[idx] == "-ej":
                time = float( opts[ idx+1 ] )
                source, sink = int(opts[idx+2]), int(opts[idx+3])
                self.split_command = "-ej {} {} {}".format( time, source, sink )
                idx += 4
            elif opts[idx] == "-eM":
                idx = self._parse_time( idx+1, opts )
                rate = float( opts[idx] ) / (self.num_populations - 1)
                idx += 1
                for i in range(self.num_populations):
                    for j in range(self.num_populations):
                        if i!=j:
                            self.migration_rates[-1][i][j] = rate
            elif opts[idx] == "-ema":
                idx = self._parse_time( idx+1, opts )
                for i in range(self.num_populations):
                    for j in range(self.num_populations):
                        rate = float( opts[idx] )
                        idx += 1
                        self.migration_rates[-1][i][j] = rate
            elif opts[idx] == "-em":
                idx = self._parse_time( idx+1, opts )
                i, j, rate = int(opts[idx]), int(opts[idx+1]), float( opts[idx+2] )
                idx += 3
                self.migration_rates[-1][i-1][j-1] = rate
            elif opts[idx] == "-eN":
                idx = self._parse_time( idx+1, opts )
                for i in range(self.num_populations):
                    self.population_sizes[-1][i] = float(opts[idx])
                idx += 1
            elif opts[idx] == "-en":
                idx = self._parse_time( idx+1, opts )
                i = int(opts[idx])
                self.population_sizes[-1][i-1] = float(opts[idx+1])
                idx += 2
            else:
                unparsed_opts.append( opts[idx] )
                idx += 1
        if self.sample_populations == []: self.sample_populations = None
        if self.sample_times == []: self.sample_times = None
        self._finalize_and_validate()
        return unparsed_opts
                
    def _add_migration_rate(self, time_idx, i, j, rate):
        self.migration_rates[time_idx][i][j] = rate        

    def _parse_time(self, opt_idx, opts ):
        if self.change_points == []:
            self.change_points = [float(opts[opt_idx])]
            self.population_sizes = [ [0.0] * self.num_populations ]
            self.migration_rates = [[[0.0]*self.num_populations for _ in range(self.num_populations)]]
        elif self.change_points[-1] != float(opts[opt_idx]):
            self.change_points.append( float(opts[opt_idx]) )
            self.population_sizes.append( self.population_sizes[-1][:] )
            self.migration_rates.append( [vec[:] for vec in self.migration_rates[-1]] )
        return opt_idx + 1

    def _parse_sample(self, time_idx, opt_idx, opts ):
        for idx in range(self.num_populations):
            num_sampled = int(opts[ opt_idx + idx ])
            self.sample_times += [ self.change_points[time_idx] ] * num_sampled
            self.sample_populations += [ 1 + idx ] * num_sampled
        return opt_idx + self.num_populations

    def _finalize_and_validate(self):

        assert type(self.change_points) == list
        assert type(self.population_sizes) == list

        # set defaults where necessary
        if self.sample_populations is None:
            self.sample_populations = [1] * self.num_samples
        if self.sample_times is None:
            self.sample_times = [0] * self.num_samples
        if self.migration_rates is None:
            self.migration_rates = [None] * len(self.change_points)
        if self.migration_commands is None:
            self.migration_commands = [None] * len(self.change_points)
        for idx, migr_rates in enumerate(self.migration_rates):
            if migr_rates is None:
                migr_rates = [ [0] * self.num_populations
                               for dummy in range(self.num_populations) ]
                self.migration_rates[idx] = migr_rates
            assert len(migr_rates) == self.num_populations

        # check sanity of samples, their times, and their population
        assert len(self.sample_populations) == self.num_samples
        assert len(self.sample_times) == self.num_samples
        assert min( self.sample_times ) == 0.0
        populations = set( self.sample_populations )
        for i in range(1, self.num_populations + 1):
            assert i in populations
        for i in range(self.num_samples - 1):
            assert self.sample_times[i] <= self.sample_times[i+1]
            if self.sample_times[i] == self.sample_times[i+1]:
                assert self.sample_populations[i] <= self.sample_populations[i+1]

        # check change_points, population_sizes and migration rates
        # also implement simplified interface when populations have identical sizes
        for idx, change_point in enumerate(self.change_points):
            if idx == 0:
                assert change_point == 0.0
            else:
                assert change_point > self.change_points[idx-1]
        if len(self.population_sizes) != len(self.change_points):
            raise ValueError("Expected {} population sizes; got {}".format(
                len(self.change_points),self.population_sizes))
        if len(self.migration_rates) != len(self.change_points):
            raise ValueError("Expected {} migration rate matrices; got {}".format(
                len(self.migration_rates),self.population_sizes))
        # allow one-population models be expressed as a simple list of population sizes
        if type(self.population_sizes[0]) != list:
            self.population_sizes = [ [size] * self.num_populations
                                      for size in self.population_sizes ]
        for idx, migr_rates in enumerate(self.migration_rates):
            for popidx, migr_vec in enumerate(migr_rates):
                assert len(migr_vec) == self.num_populations
                assert migr_vec[popidx] == 0
                assert min(migr_vec) == 0
        if float( self.split_command.split()[1] ) not in [ float(t) for t in self.change_points]:
            raise ValueError("Time used for -ej ({}) not a valid change point! (change points: {})".format(self.split_command.split()[1],self.change_points))

        # require the population size at time 0 to be set explicitly
        assert( self.change_points[0] == 0.0 )


    def _create_sample_command_line_options(self):

        self._finalize_and_validate()
        expression = ""
        if self.num_populations == 1 and max( self.sample_times ) == 0.0:
            return expression
        times = sorted( list( set( self.sample_times ) ) )
        for time in times:
            popsizes = [ len( [ None
                                for idx, sample_time in enumerate(self.sample_times)
                                if sample_time == time and self.sample_populations[idx] == pop ] )
                         for pop in range( 1, self.num_populations + 1 ) ]
            if time == 0.0:
                expression = "-I {} {}".format(self.num_populations,
                                               " ".join(map(str,popsizes)))
            else:
                expression = expression + " -eI {} {}".format(time,
                                                              " ".join(map(str,popsizes)))
        return expression


    def _create_popsize_migration_command_line_options(self,
                                                      inference_popsizes = None,
                                                      inference_migrationrates = None,
                                                      inference_changepoints = None,
                                                      inference_migrationcommands = None):
        self._finalize_and_validate()
        my_popsizes =          self.population_sizes   if inference_popsizes is None          else inference_popsizes
        my_migrationrates =    self.migration_rates    if inference_migrationrates is None    else inference_migrationrates
        my_changepoints =      self.change_points      if inference_changepoints is None      else inference_changepoints
        my_migrationcommands = self.migration_commands if inference_migrationcommands is None else inference_migrationcommands

        expression = []
        for i in range(len(my_changepoints)):
            time = my_changepoints[i]
            popsizes = my_popsizes[i]
            if len(set(popsizes)) == 1:
                expression.append("-eN {} {}".format(time, popsizes[0]))
            else:
                for idx, popsize in enumerate(popsizes):
                    expression.append("-en {} {} {}".format(time, idx+1, popsize))
            migration_matrix = my_migrationrates[i]
            all_rates = set( migration_matrix[j][k]
                             for j,k in itertools.product( range(self.num_populations),
                                                           range(self.num_populations) )
                             if j != k )
            if len(all_rates) == 1:
                rate = all_rates.pop()
                if rate != 0:
                    total_symmetric_rate = rate * (self.num_populations - 1)
                    expression.append("-eM {} {}".format(time, total_symmetric_rate))
            else:
                if len(all_rates) > 1:
                    expression.append("-ema {}".format(time))
                    # note that the migration rates are listed in the order
                    #   m_00 m_01 m_02 ... m_10 m_11 m_12 ...
                    # where m_ij is the rate FROM population i TO population j BACKWARDS in time,
                    # which seems to contradict the help text of scrm.
                    for j in range(self.num_populations):
                        for k in range(self.num_populations):
                            expression.append("{}".format(migration_matrix[j][k]))
            migration_command = my_migrationcommands[i]
            if migration_command is not None:
                expression.append(migration_command)
        return " ".join(expression)


    def core_command_line( self,
                           inference_popsizes = None,
                           inference_migrationrates = None,
                           inference_changepoints = None,
                           inference_migrationcommands = None):
        """ builds core command line common to simulation and inference;
            without the initial "num_samples num_loci" parameters, and without
            output-related commands """

        self.mutations      = 4 * self.N0 * self.mutation_rate * self.sequence_length
        self.recombinations = 4 * self.N0 * self.recombination_rate * self.sequence_length

        command = "-N0 {N0} -t {muts} -r {recs} {seqlen} {sample} {split} {popmigr}".format(
            N0 = self.N0,
            muts = self.mutations,
            recs = self.recombinations,
            seqlen = self.sequence_length,
            sample = self._create_sample_command_line_options(),
            split = self.split_command,
            popmigr = self._create_popsize_migration_command_line_options(
                inference_popsizes,
                inference_migrationrates,
                inference_changepoints,
                inference_migrationcommands
            )
        )
        return command


    def simulate(self, missing_leaves = [], phased = True, debug = False):

        # make file name if required
        if self.filename == None: final_filename = command.replace(' ','_') + ".seg"
        else:                     final_filename = self.filename

        temp_file = tempfile.NamedTemporaryFile( dir = os.path.dirname( os.path.abspath( final_filename ) ) )
        filename = temp_file.name + ".seg"
        temp_file.close()
        scrmfilename = filename + ".scrm"

        # -T        print Newick trees;
        # -L        print TMRCA and local tree length for each segment;
        # -p 10     use 10 digits precision;
        # -l 300000 limit exact LD window to 300 kb
        command = "{cmd} {nsam} {nloci} {core} {seed} -T -L -p 10 -l 300000 > {fn}".format(
            cmd = self.scrmpath,
            nsam = self.num_samples,
            nloci = 1,
            core = self.core_command_line(),
            seed = "-seed {}".format(' '.join(map(str,self.seed))),
            fn = scrmfilename)
        self.simulate_command = command

        # ensure no race conditions can occur when running many simulations in parallel
        if os.path.exists( final_filename ): return

        # execute
        if debug: print (self.simulate_command)
        returnvalue = subprocess.check_call(command, shell=True)
        if returnvalue > 0: raise ValueError("Problem executing " + command)

        # generate the .seg and .seg.recomb files
        self.convert_scrm_to_seg( scrmfilename, filename, missing_leaves, phased )
        # don't generate .seg.recomb files; these are only useful for plotting to check that inference works
        #self.convert_scrm_to_recomb( scrmfilename, filename + ".recomb")

        # done; remove scrm output file
        os.unlink( scrmfilename )

        # finally, move to true destination (if possible)
        self.move_into_place( filename, final_filename )
        #self.move_into_place( filename + ".recomb", final_filename + ".recomb" )


    def move_into_place(self, filename, final_filename ):
        if os.path.exists( final_filename ):
            os.unlink( filename )                     # appeared while we were simulating -- remove
        else:
            try:
                os.rename( filename, final_filename ) # move into destination
            except:
                os.unlink( filename )                 # assume failure means file has appeared -- remove


    def convert_scrm_to_seg(self, infilename, outfilename, missing_leaves, phased ):
        # convert scrm output file to .seg file for input to smcsmc

        scrmfile = open(infilename, "r")
        data = None
        positions = None
        for line in scrmfile:
            if line.startswith('positions:'):
                positions = list(map(float,line.strip().split()[1:]))
                data = []
            elif type(data) == type([]):
                if positions == None:
                    raise ValueError("Unexpected output from scrm -- no 'positions' line found")
                dataline = line.strip()
                # make data go missing...
                if len(data) in missing_leaves:
                    dataline = '.' * len(dataline)
                data.append(dataline)
                if len(data[-1]) != len(positions):
                    raise ValueError("Unexpected number of mutations: got {}, expected {}".format(
                        len(data[-1]),len(positions)))
        scrmfile.close()

        if data == None or len(data) != self.num_samples:
            raise ValueError("Unexpected data from scrm")

        outfile = open(outfilename,'w')
        row = "{pos}\t{distance}\tT\tF\t1\t{genotype}\n"

        positions = list(map( lambda realpos : int(realpos * self.sequence_length + 0.5),
                              positions ))
        positions = [1] + positions

        if phased:
            for idx in range(len(positions) - 1):
                outfile.write( row.format( pos = positions[idx],
                                           distance = positions[idx+1] - positions[idx],
                                           genotype = ''.join( [sequence[idx]
                                                                for sequence in data] ) ) )
        else:
            for idx in range(len(positions) - 1):
                alleles_list = [sequence[idx] for sequence in data]
                for haplotype_idx in range( len(alleles_list)/2 ):
                    if alleles_list[haplotype_idx*2] != alleles_list[haplotype_idx*2+1]:
                        alleles_list[haplotype_idx*2  ] = "/"
                        alleles_list[haplotype_idx*2+1] = "/"
                alleles = ''.join( alleles_list )
                outfile.write( row.format( pos = positions[idx],
                                           distance = positions[idx+1] - positions[idx],
                                           genotype = alleles ) )

        outfile.write( row.format( pos = positions[-1],
                                   distance = self.sequence_length - positions[-1],
                                   genotype = "." * len(data) ) )

        outfile.close()


    def convert_scrm_to_recomb(self, infilename, outfilename, segsize = 100, fourne = 40000 ):
        # convert scrm output file to .recomb file

        scrmfile = open(infilename, "r")
        pos = 0
        recombinations = []
        treelengths = []
        for line in scrmfile:
            if line.startswith('['):
                # newick line; parse segment length and add to recombination positions
                segment = int(line[1:].split(']')[0])
                recombinations.append( segment )
            elif line.startswith('time'):
                # tmrca + total tree length line
                treelengths.append( float(line[:-1].split()[2]) * fourne )

        if len(recombinations) != len(treelengths):
            raise ValueError("Unexpected number of tree lengths; got {}, expected {}".format(
                len(treelengths),
                len(recombinations)))
        scrmfile.close()

        # now convert to segments.  this is done by a series of generators, that each generate a
        # sequence of (pos, segment_size, value) or (pos, segment_size, value1, valud2) tuples.
        # The first two generators use the list just read in to do this, the others
        # consume data from an upstream generator.  In this way it is relatively straightforward to
        # do a fairly complex transformation of the data: convert into segments, discretize to
        # 100 bp blocks, and merge blocks with the same value together (and add a header).
        # It's possibly not the simplest way to do this, but it nicely separates the various steps
        # of the transformation.  This style is also memory-efficient -- no intermediate data
        # stores are necessary; trivial advantage in this case, but nice for very large data sets.
        class Conversion:
            def __init__(self, lengths, heights, segsize):
                self.lengths = lengths
                self.heights = heights
                self.segsize = segsize
            def segment_gen_opp(self):
                pos = 0
                for idx, length in enumerate(self.lengths):
                    yield pos, length, self.heights[idx]
                    pos += length
            def segment_gen_rec(self):
                pos = 0
                for idx, length in enumerate(self.lengths):
                    yield pos, length-1, 0
                    pos += length-1
                    yield pos, 1, 1
                    pos += 1
            def segment_discretize(self, gen):
                pos = 0
                integral = 0
                try:
                    while True:
                        newpos, length, height = next(gen)
                        # generate complete segments until current overlaps endpoint of
                        # [newpos, newpos+length)
                        while pos + self.segsize < newpos + length:
                            # add height over appropriate subinterval
                            if newpos < pos + self.segsize:
                                integral += height * (pos + self.segsize - max(pos, newpos))
                            # report current interval
                            yield pos, self.segsize, integral / float(self.segsize)
                            # next interval
                            pos += self.segsize
                            integral = 0
                        # cannot generate further segments; add last bit of current, get new one
                        integral += max(0, min( pos + self.segsize, newpos + length)
                                           - max( pos, newpos )) * height
                except StopIteration:
                    # no further segments -- yield last one, and end
                    yield pos, self.segsize, integral / float(self.segsize)
            def merge(self, gen1, gen2):
                for pos, length, height in gen1:
                    try:
                        pos2, length2, height2 = next(gen2)
                    except StopIteration:
                        pos2, length2, height2 = pos1, length1, 0
                    assert pos == pos2 and length == length2
                    yield pos, length, height, height2
            def collapse(self, gen):
                pos, length, height1, height2 = "locus", "size", "opp_per_nt", "recomb"
                for newpos, newlength, newheight1, newheight2 in gen:
                    if newheight1 == height1 and newheight2 == height2:
                        length += newlength
                    else:
                        yield pos, length, height1, height2
                        pos, length, height1, height2 = newpos, newlength, newheight1, newheight2
                yield pos, length, height1, height2

        conversion = Conversion(recombinations, treelengths, segsize)

        opportunity =   conversion.segment_discretize( conversion.segment_gen_opp() )
        recombination = conversion.segment_discretize( conversion.segment_gen_rec() )
        combined =      conversion.merge( opportunity, recombination )
        collapsed =     conversion.collapse( combined )

        # output list
        outfile = open( outfilename, 'w' )
        for pos, length, opp, rec in collapsed:
            if type(opp) == float: opp = "{:11.2f}".format(opp)
            outfile.write( "{pos}\t{length}\t{opp}\t{rec}\n".format(
                pos=pos,
                length=length,
                opp=opp,
                rec=rec))
        outfile.close()



#
# some models
#

class Pop2(Population):

    def __init__(self, **kwargs):
        Population.__init__(self, **kwargs)


class Pop4(Population):

    def __init__(self, **kwargs):

        pop4_defaults = {'change_points':    [0, .01, 0.06, 0.2,  1, 2],
                         'population_sizes': [1, .1,  1,    0.5,  1, 2],
                         'num_samples': 4}

        # set kwargs to Pop4-default values, but allow caller to override
        for key, value in pop4_defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        Population.__init__(self, **kwargs)



class PopSingleConst(Population):

    def __init__(self, **kwargs):

        defaults = {'change_points':    [0, .5, 1.0],
                    'population_sizes': [1, 1,  1],
                    'num_samples': 4}

        # set kwargs to default values, but allow caller to override
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        Population.__init__(self, **kwargs)



class PopSingleExpand(Population):

    def __init__(self, **kwargs):

        defaults = {'change_points':    [0, .02],
                    'population_sizes': [2, 1],
                    'num_samples': 4}

        # set kwargs to default values, but allow caller to override
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        Population.__init__(self, **kwargs)


class PopSingleShrink(Population):

    def __init__(self, **kwargs):

        defaults = {'change_points':    [0, .02],
                    'population_sizes': [.5, 1],
                    'num_samples': 4}

        # set kwargs to default values, but allow caller to override
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        Population.__init__(self, **kwargs)


class PopSingleBottleneck(Population):

    def __init__(self, **kwargs):

        defaults = {'change_points':    [0, .02, .04],
                    'population_sizes': [1, .5, 1],
                    'num_samples': 4}

        # set kwargs to default values, but allow caller to override
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        Population.__init__(self, **kwargs)



class TwoPopUniDirMigr(Population):

    def __init__(self, **kwargs):

        defaults = {'change_points':    [0, .1, .5],
                    'population_sizes': [[1,1], [1,1], [1,1]],
                    'num_populations':  2,
                    'migration_rates':  [ [[0,0],[1,0]],   # -em 0   2 1 1
                                          [[0,0],[1,0]],   # -em 0.1 2 1 1
                                          [[0,0],[1,0]] ], # -em .5  2 1 1
                    'num_samples':      8}

        # set kwargs to default values, but allow caller to override
        for key, value in defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        Population.__init__(self, **kwargs)


class TwoPopBiDirMigr(Population):

    def __init__(self, **kwargs):

        pop2migr_defaults = {'change_points':    [0, .1, .5],
                             'population_sizes': [[1,1], [1,1], [1,1]],
                             'migration_rates':  [ [[0,0.5],[0.5,0]],
                                                   [[0,0.5],[0.5,0]],
                                                   [[0,0.5],[0.5,0]] ],
                             'num_samples': 8,
                             'num_populations': 2}

        # set kwargs to default values, but allow caller to override
        for key, value in pop2migr_defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        Population.__init__(self, **kwargs)


class TwoPopSplitNoMigr(Population):

    def __init__(self, **kwargs):

        pop2migr_defaults = {'change_points':      [0, .1, .5],
                             'population_sizes':   [[1,1], [1,1], [1,1]],
                             'num_populations':    2,
                             'migration_commands': [None,None,"-ej 0.5 2 1"],
                             'num_samples': 8}

        # set kwargs to default values, but allow caller to override
        for key, value in pop2migr_defaults.items():
            if key not in kwargs:
                kwargs[key] = value

        Population.__init__(self, **kwargs)




#
# main code for standalone execution
#
if __name__ == "__main__":
    models = {"Population": Population}
    labels = str(tuple(models.keys())).replace("'","")
    parser = OptionParser()
    parser.add_option("-o", "--output", dest="outfile",
                      help="write data to FILE.seg", metavar="FILE.seg")
    parser.add_option("-m", "--model", dest="model",
                      help="choose model {}".format(labels), default = "Population")
    parser.add_option("-n", "--numsamples", type="int", dest="num_samples",
                      help="number of samples to simulate")
    parser.add_option("-l", "--seqlen", type="int", dest="seqlen",
                      help="sequence length", default = 1e6)
    parser.add_option("-s", "--seed", type="int", dest="seed",
                      help="random number seed", default = 1)

    (options, args) = parser.parse_args()

    if len(args) != 0:
        raise ValueError("Unexpected arguments:" + str(args))

    if options.num_samples == None:
        print ("Error: -n, --numsamples is required.  Use --help for help")
        sys.exit(1)

    popclass = models[ options.model ]
    pop = popclass( filename = options.outfile,
                    num_samples = options.num_samples,
                    sequence_length = options.seqlen,
                    seed = (options.seed, 1, 1))
    pop.simulate()

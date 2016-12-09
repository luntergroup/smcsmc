#!/usr/bin/env python
from __future__ import print_function

import sys
import subprocess
import os

from optparse import OptionParser

"""
Population structures
"""


defaults = {'N0': 10000,
            'mutation_rate': 2.5e-8,
            'recombination_rate': 5e-9,
            'sequence_length': 1e6,
            'num_samples': 2,
            'scrmpath': 'scrm'}

class Population:

    def __init__(self,
                 N0                   = defaults['N0'],
                 mutation_rate        = defaults['mutation_rate'],
                 recombination_rate   = defaults['recombination_rate'],
                 sequence_length      = defaults['sequence_length'],
                 num_samples          = defaults['num_samples'],
                 change_points        = [.01, 0.06, 0.2, 1, 2],
                 population_sizes     = [0.1, 1, 0.5, 1, 2],
                 migration_commands   = [None, None, None, None, None],
                 additional_commands  = "",
                 seed                 = (3647837471,),
                 filename             = None,
                 scrmpath             = defaults['scrmpath'] ):

        self.N0 = N0
        self.mutation_rate = mutation_rate
        self.recombination_rate = recombination_rate
        self.sequence_length = sequence_length
        self.num_samples = num_samples
        self.change_points = change_points
        self.population_sizes = population_sizes
        self.migration_commands = migration_commands
        self.additional_commands = additional_commands
        self.seed = seed
        self.filename = filename
        self.scrmpath = scrmpath

        if len(population_sizes) != len(change_points):
            raise ValueError("Expected {} population sizes; got {}".format(len(change_points),population_sizes))
        if len(migration_commands) != len(change_points):
            raise ValueError("Expected {} migration commands; got {}".format(len(change_points),migration_commands))


    def simulate(self, missing_leaves = []):

        num_loci = 1
        self.mutations      = 4 * self.N0 * self.mutation_rate * self.sequence_length
        self.recombinations = 4 * self.N0 * self.recombination_rate * self.sequence_length

        command = "{nsam} {nloci} -t {muts} -r {recs} {seqlen} {add}".format(nsam = self.num_samples,
                                                                             nloci = num_loci,
                                                                             muts = self.mutations,
                                                                             recs = self.recombinations,
                                                                             seqlen = self.sequence_length,
                                                                             add = self.additional_commands)

        # add population structure and migration
        for i in range(len(self.change_points)):
            time = self.change_points[i]
            popsize = self.population_sizes[i]
            migration_command = self.migration_commands[i]
            if popsize != None:
                command += " -eN {time} {popsize}".format(time=time, popsize=popsize)
            if migration_command != None:
                command += " " + migration_command

        # make file name if required
        if self.filename == None:
            filename = command.replace(' ','_') + ".seg"
        else:
            filename = self.filename

        # complete the command
        command = "{cmd} {args} -seed {seed}".format(cmd = self.scrmpath,
                                                     args = command,
                                                     seed = ' '.join(map(str,self.seed)))

        # print Newick trees;
        # print TMRCA and local tree length for each segment;
        # use 10 digits precision;
        # limit exact LD window to 300 kb
        # output
        scrmfilename = filename + ".scrm"
        command += " -T -L -p 10 -l 300000 > " + scrmfilename
        print(command)
        # execute
        returnvalue = subprocess.check_call(command, shell=True)

        if returnvalue > 0:
            raise ValueError("Problem executing " + command)

        self.convert_scrm_to_seg( scrmfilename, filename, missing_leaves )

        self.convert_scrm_to_recomb( scrmfilename, filename + ".recomb")

        # done; remove scrm output file
        os.unlink( scrmfilename )


    def convert_scrm_to_seg(self, infilename, outfilename, missing_leaves ):
        # convert scrm output file to .seg file

        scrmfile = open(infilename, "r")
        data = None
        positions = None
        for line in scrmfile:
            if line.startswith('positions:'):
                positions = list(map(float,line.strip().split()[1:]))
                data = []
            elif type(data) == type([]):
                if positions == None:
                    raise ValueError("Unexpected output from scrm -- not 'positions' line found")
                dataline = line.strip()
                # make data go missing...
                if len(data) in missing_leaves:
                    dataline = '.' * len(dataline)
                data.append(dataline)
                if len(data[-1]) != len(positions):
                    raise ValueError("Unexpected number of mutations: got {}, expected {}".format(len(data[-1]),len(positions)))
        scrmfile.close()

        if data == None or len(data) != self.num_samples:
            raise ValueError("Unexpected data from scrm")

        outfile = open(outfilename,'w')
        row = "{pos}\t{distance}\tT\tF\t1\t{genotype}\n"

        positions = list(map( lambda realpos : int(realpos * self.sequence_length + 0.5), positions ))
        positions.append( int(self.sequence_length) )

        outfile.write( row.format( pos=1, distance = positions[0]-1, genotype = "." * len(data) ) )
        for idx in range(len(positions) - 1):
            outfile.write( row.format( pos=positions[idx],
                                       distance = positions[idx+1] - positions[idx],
                                       genotype = ''.join( [sequence[idx] for sequence in data] ) ) )
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
            raise ValueError("Unexpected number of tree lengths; got {}, expected {}".format(len(treelengths),len(recombinations)))
        scrmfile.close()

        # now convert to segments.  this is done by a series of generators, that each generate a sequence of (pos, segment_size, value)
        # or (pos, segment_size, value1, valud2) tuples.  The first two generators use the list just read in to do this, the others
        # consume data from an upstream generator.  In this way it is relatively straightforward to do a fairly complex transformation
        # of the data: convert into segments, discretize to 100 bp blocks, and merge blocks with the same value together (and add a header).
        # it's possibly not the simplest way to do this, but it nicely separates the various steps of the transformation.  This style is
        # also very memory-efficient -- no intermediate data stores are necessary; a trivial advantage in this case, but nice for very large
        # data sets.
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
                        # generate complete segments until current overlaps endpoint of [newpos, newpos+length)
                        while pos + self.segsize < newpos + length:
                            # add height over appropriate subinterval
                            if newpos < pos + self.segsize:
                                integral += height * (pos + self.segsize - max(pos, newpos))
                            # report current interval
                            yield pos, self.segsize, integral / float(self.segsize)
                            # next interval
                            pos += self.segsize
                            integral = 0
                        # cannot generate further segments -- add last bit of current, and get new one
                        integral += max(0, min( pos + self.segsize, newpos + length) - max( pos, newpos )) * height
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
            outfile.write( "{pos}\t{length}\t{opp}\t{rec}\n".format(pos=pos, length=length, opp=opp, rec=rec))
        outfile.close()



#
# some models
#
class Pop1(Population):

    def __init__(self, **kwargs):
        Population.__init__(self, **kwargs)


class Pop1Migr(Population):

    def __init__(self,
                 change_points        = [.01, 0.06, 0.2, 0.6,  1, 2],
                 population_sizes     = [0.1, 1,    0.5, None, 1, 2],
                 migration_commands   = [None, None, None, "-ej 0.6 2 1", None, None],
                 **kwargs):

        nsam = kwargs.get('num_samples', defaults['num_samples'])

        Population.__init__(self,
                            change_points = change_points,
                            population_sizes = population_sizes,
                            migration_commands = migration_commands,
                            additional_commands = "-I 2 {} {} 0.0".format( int(nsam/2), int(nsam/2) ),
                            **kwargs)

class Pop4(Population):

    def __init__(self,
                 change_points        = [.01, 0.06, 0.2,  1, 2],
                 population_sizes     = [0.1, 1,    0.5,  1, 2],
                 migration_commands   = [None, None, None, None, None],
                 num_samples          = 4,
                 **kwargs):

        nsam = kwargs.get('num_samples', defaults['num_samples'])

        Population.__init__(self,
                            change_points = change_points,
                            population_sizes = population_sizes,
                            migration_commands = migration_commands,
                            num_samples = num_samples,
                            **kwargs)


class PopSingleConst(Population):

    def __init__(self,
                 change_points        = [.5, 1.0],
                 population_sizes     = [1, 1],
                 migration_commands   = [None, None],
                 num_samples          = 2,
                 **kwargs):

        nsam = kwargs.get('num_samples', defaults['num_samples'])

        Population.__init__(self,
                            change_points = change_points,
                            population_sizes = population_sizes,
                            migration_commands = migration_commands,
                            num_samples = num_samples,
                            **kwargs)
        self.target_min = [0,     9500, 9500]
        self.target_max = [1e+10, 11500, 11500]


class PopSingleExpand(Population):

    def __init__(self,
                 change_points        = [0, .5, 1.0],
                 population_sizes     = [5, 5, 1],
                 migration_commands   = [None, None, None],
                 num_samples          = 2,
                 **kwargs):

        nsam = kwargs.get('num_samples', defaults['num_samples'])

        Population.__init__(self,
                            change_points = change_points,
                            population_sizes = population_sizes,
                            migration_commands = migration_commands,
                            num_samples = num_samples,
                            **kwargs)
        self.target_min = [0,     49000, 9500]
        self.target_max = [1e+10, 51000, 11500]


class TwoPopUniDirMigr(Population):

    def __init__(self,
                 change_points        = [0, .5, 1.0],
                 population_sizes     = [5, 5, 1],
                 migration_commands   = ["-em 0 2 1 1", "-em 0.5 2 1 1", "-em 1.0 2 1 1"],
                 num_samples          = 8,
                 **kwargs):

        nsam = kwargs.get('num_samples', defaults['num_samples'])

        Population.__init__(self,
                            change_points = change_points,
                            population_sizes = population_sizes,
                            migration_commands = migration_commands,
                            num_samples = num_samples,
                            additional_commands = "-I 2 {} {}".format( int(num_samples/2), int(num_samples/2) ),
                            **kwargs)
        self.target_min = [[0,     49000, 9500],
                           [0,     49000, 9500]]
        self.target_max = [[1e+10, 51000, 11500],
                           [1e+10, 51000, 11500]]
        self.migr_target_min = [[0, 0, 0],
                                [0, 0, 0]]
        self.migr_target_max = [[3e-5, 3e-5, 3e-5],
                                [3e-5, 3e-5, 3e-5]]


class PopSingleShrink(Population):

    def __init__(self,
                 change_points        = [0, .5, 1.0],
                 population_sizes     = [5,  5, 10],
                 migration_commands   = [None, None],
                 num_samples          = 2,
                 **kwargs):

        nsam = kwargs.get('num_samples', defaults['num_samples'])

        Population.__init__(self,
                            change_points = change_points,
                            population_sizes = population_sizes,
                            migration_commands = migration_commands,
                            num_samples = num_samples,
                            **kwargs)
        self.target_min = [0,     49000, 95000]
        self.target_max = [1e+10, 51000, 115000]

#
# main code for standalone execution
#
if __name__ == "__main__":
    models = {"Pop1": Sim1,
              "Pop1Migr": Sim1Migr}
    labels = str(tuple(models.keys())).replace("'","")
    parser = OptionParser()
    parser.add_option("-o", "--output", dest="outfile", help="write data to FILE.seg", metavar="FILE.seg")
    parser.add_option("-m", "--model", dest="model", help="choose model {}".format(labels), default = "Pop1")
    parser.add_option("-n", "--numsamples", type="int", dest="num_samples", help="number of samples to simulate")
    parser.add_option("-l", "--seqlen", type="int", dest="seqlen", help="sequence length", default = 1e6)
    parser.add_option("-s", "--seed", type="int", dest="seed", help="random number seed", default = 1)

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

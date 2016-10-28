#!/usr/bin/env python

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
            'num_samples': 2}

class Population:

    def __init__(self,
                 N0                   = defaults['N0'],
                 mutation_rate        = defaults['mutation_rate'],
                 recombination_rate   = defaults['recombination_rate'],
                 sequence_length      = 1e6,
                 num_samples          = defaults['num_samples'],
                 change_points        = [.01, 0.06, 0.2, 1, 2],
                 population_sizes     = [0.1, 1, 0.5, 1, 2],
                 migration_commands   = [None, None, None, None, None],
                 additional_commands  = "",
                 seed                 = (1,1,1),
                 filename             = None ):

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


        if len(population_sizes) != len(change_points):
            raise ValueError("Expected {} population sizes; got {}".format(len(change_points),population_sizes))
        if len(migration_commands) != len(change_points):
            raise ValueError("Expected {} migration commands; got {}".format(len(change_points),migration_commands))


    def simulate(self):

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
            filename = command.replace(' ','_') + ".scrmdata"
        else:
            filename = self.filename

        # complete the command
        command = "{cmd} {args} -seed {seed[0]} {seed[1]} {seed[2]}".format(cmd = "scrm",
                                                                         args = command,
                                                                         seed = self.seed)

        # print Newick trees; use 10 digits precision; limit exact LD window to 300 kb
        command += " -T -p 10 -l 300000 > " + filename + ".tmp"

        # execute
        returnvalue = subprocess.check_call(command, shell=True)

        if returnvalue > 0:
            raise ValueError("Problem executing " + command)

        scrmfile = open(filename + ".tmp", "r")
        data = None
        positions = None
        for line in scrmfile:
            if line.startswith('positions:'):
                positions = map(float,line.strip().split()[1:])
                data = []
            elif type(data) == type([]):
                data.append(line.strip())
                if len(data[-1]) != len(positions):
                    raise ValueError("Unexpected number of mutations: got {}, expected {}".format(len(data[-1]),len(positions)))
        scrmfile.close()

        if data == None or len(data) != self.num_samples:
            raise ValueError("Unexpected data from scrm")

        outfile = open(filename,'w')
        row = "{pos}\t{distance}\tT\tF\t1\t{genotype}\n"

        positions = map( lambda realpos : int(realpos * self.sequence_length + 0.5), positions )
        positions.append( int(self.sequence_length) )
        
        outfile.write( row.format( pos=1, distance = positions[0]-1, genotype = ".." ) )
        for idx in range(len(positions) - 1):
            outfile.write( row.format( pos=positions[idx],
                                       distance = positions[idx+1] - positions[idx],
                                       genotype = ''.join( [sequence[idx] for sequence in data] ) ) )
        outfile.close()

        os.unlink( filename + ".tmp" )


#
# some models
#
class Sim1(Population):

    def __init__(self, **kwargs):
        Population.__init__(self, **kwargs)


class Sim1Migr(Population):

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




if __name__ == "__main__":
    models = {"Sim1": Sim1,
              "Sim1Migr": Sim1Migr}
    parser = OptionParser()
    parser.add_option("-o", "--output", dest="outfile", help="write data to FILE.seg", metavar="FILE.seg")
    parser.add_option("-m", "--model", dest="model", help="choose model (Sim1, Sim1Migr)", default = "Sim1")
    parser.add_option("-n", "--numsamples", type="int", dest="num_samples", help="number of samples to simulate")
    parser.add_option("-l", "--seqlen", type="int", dest="seqlen", help="sequence length", default = 1e6)
    parser.add_option("-s", "--seed", type="int", dest="seed", help="random number seed", default = 1)

    (options, args) = parser.parse_args()

    if len(args) != 0:
        raise ValueError("Unexpected arguments:" + str(args))

    if options.num_samples == None:
        print "Error: -n, --numsamples is required.  Use --help for help"
        sys.exit(1)

    popclass = models[ options.model ]
    pop = popclass( filename = options.outfile,
                    num_samples = options.num_samples,
                    sequence_length = options.seqlen,
                    seed = (options.seed, 1, 1))
    pop.simulate()

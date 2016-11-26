from __future__ import print_function

import unittest
import os
import subprocess

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
        self.tmax = 4
        self.seed = (3647837471,)
        self.success = False

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
        popt = "-p 1*3+15*4+1"
        npopt = "-Np {np}".format(np=self.np)
        emopt = "-EM {em}".format(em=self.em)
        tmaxopt = "-tmax {tmax}".format(tmax=self.tmax)
        seedopt = "-seed {seed}".format(seed=' '.join(map(str,self.seed)))
        segopt = "-seg {}".format( self.segfile )
        return "../smcsmc {nsam} {t} {r} {np} {em} {p} {tmax} {seed} {seg}".format(
            nsam = nsamopt,
            t = topt,
            r = ropt,
            np = npopt,
            em = emopt,
            p = popt,
            tmax = tmaxopt,
            seed = seedopt,
            seg = segopt)

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
                   'Recomb':[]}
        for line in open(self.outfile):
            elts = line.strip().split()
            if elts[0] == "Iter": continue
            it, epoch, start, end, typ, ne = elts[0],elts[1],elts[2],elts[3],elts[4],elts[10]
            if int(epoch) == -1:
                # recombination rate is assigned to epoch -1
                epoch = 0
            if len(results[typ]) <= int(epoch):
                results[typ].append( {'start': float(start),
                                      'end': float(end),
                                      'estimates': []} )
            assert( len(results[typ][int(epoch)]['estimates']) == int(it) )
            results[typ][int(epoch)]['estimates'].append( float(ne) )
        return results


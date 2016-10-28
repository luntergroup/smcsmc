import unittest
import os
import subprocess

import simulate


class TestGeneric(unittest.TestCase):

    def setUp(self):
        self.simulate()
        self.copy_opts()

    def simulate(self):
        self.filename = None
        self.outfile = None
        self.pop = simulate.Sim1( filename = self.filename )
        self.pop.simulate()

    def copy_opts(self):
        self.nsamopt = "-nsam {}".format(self.pop.num_samples)
        self.topt = "-t {}".format(self.pop.mutations)
        self.ropt = "-r {} {}".format(self.pop.recombinations,self.pop.sequence_length)
        self.npopt = "-Np 1000"
        self.emopt = "-EM 0"
        self.popt = "-p 1*3+15*4+1"
        self.tmaxopt = "-tmax 4"
        self.seedopt = "-seed 1 1 1"
        self.segopt = "-seg {}".format( self.filename )

    def command(self):
        return "../smcsmc {nsam} {t} {r} {np} {em} {p} {tmax} {seed} {seg}".format(
            nsam = self.nsamopt,
            t = self.topt,
            r = self.ropt,
            np = self.npopt,
            em = self.emopt,
            p = self.popt,
            tmax = self.tmaxopt,
            seed = self.seedopt,
            seg = self.segopt)

    def infer(self, case = 0):
        self.outfile = self.filename + ".out" + str(case)
        cmd = self.command() + " > " + self.outfile
        returnvalue = subprocess.check_call(cmd, shell = True)
        self.assertTrue( returnvalue == 0 )
    
    def tearDown(self):
        #if filename != None: os.unlink( filename )
        #if outfile != None:  os.unlink( outfile )
        pass
        



class TestTwoSample(TestGeneric):

    def simulate(self):
        self.filename = "2sample.seg"
        self.seqlen = 1e7
        self.pop = simulate.Sim1( filename = self.filename, sequence_length = self.seqlen )
        self.pop.simulate()
        
    def copy_opts(self):
        TestGeneric.copy_opts(self)

    def test_inference(self):
        self.infer( case = 1 )
        
        



if __name__ == "__main__":
    unittest.main()

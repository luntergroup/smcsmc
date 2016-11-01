import unittest
import os
import subprocess

import populationmodels


class TestGeneric(unittest.TestCase):

    # called every time an instance of TestGeneric is made -- set inference defaults
    def setUp(self):
        self.np = 1000
        self.em = 0
        self.tmax = 4
        self.seed = [1,1,1]
        self.success = False

    # called every time an instance of TestGeneric is destroyed -- remove output file
    def tearDown(self):
        if self.success and 'caseprefix' in self.__dict__:
            if self.caseprefix != None:
                for suffix in ['Resample','.outXXX','.logXXX','.stdout','.stderr']:
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
        seedopt = "-seed {seed[0]} {seed[1]} {seed[2]}".format(seed=self.seed)
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



class TestTwoSample(TestGeneric):

    # Called by unittest.main when TestTwoSample is created.
    # Responsible for per-class creation of simulated data.
    # It is implemented as a class method to avoid simulating the same data multiple times
    @classmethod
    def setUpClass(cls):
        cls.prefix = "2sample"
        cls.segfile = cls.prefix + ".seg"
        cls.seqlen = 1e7
        cls.scrmpath = "../scrm"
        cls.pop = populationmodels.Pop1( filename = cls.segfile, sequence_length = cls.seqlen, scrmpath=cls.scrmpath )
        cls.success = True   # set ourselves up for success

        print ("simulating for",cls.prefix,"...")
        cls.pop.simulate()

    # remove simulated data
    @classmethod
    def tearDownClass(cls):
        if cls.segfile != None and cls.success:
            os.unlink( cls.segfile )
            cls.filename = None

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
            
    # actual test
    def test_inference(self):
        self.em = 4
        self.np = 1000
        self.seed = [1,1,1]
        self.infer( case = 1 )

        # Ne target ranges for the varous epochs
        target_min = [0,     1000, 700,  3000, 9000, 9000, 8000,7000, 5500, 5000, 4500, 4500, 4500, 4500, 5500, 10000, 5000]
        target_max = [1e+10, 2000, 1100, 5000, 12000,11000,9500,8000, 6500, 6000, 5500, 5500, 5500, 6500, 7000, 15000, 8000]
        results = self.readResults()

        for epoch, emresults in enumerate(results['Coal']):
            # get the last EM estimate
            result = emresults['estimates'][-1]

            # see if it is within range
            out_of_range = 0
            print ("Checking epoch",epoch,result,target_min[epoch],target_max[epoch])
            #self.assertTrue( result >= target_min[epoch] )
            #self.assertTrue( result <= target_max[epoch] )
            if result < target_min[epoch] or result > target_max[epoch]:
                print ("Out of range!")
                out_of_range += 1

        self.assertTrue( out_of_range < 3 )
        self.success = True


if __name__ == "__main__":
    unittest.main()

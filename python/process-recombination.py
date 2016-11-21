from __future__ import print_function

import sys

class LocalRecombination:
    
    def __init__(self, infile):
        self.read_data(infile)
        self.size   = self.data[-1][1] + self.data[-1][2]
        self.rate   = self.calculate_rate()
        self.leaves = len(self.data[0][4:])

    def calculate_rate(self):
        opportunity, recombinations = 0, 0
        for data in self._forward_data_generator( encode=False ):
            opportunity += data[1] * data[2]
            recombinations += data[1] * sum(data[3:])
        return recombinations / opportunity

    def read_data(self, infile, iter=0):
        self.data = []
        gcd = 0
        for line in open(infile,'r').readlines():
            if line.startswith('iter'): continue
            eltsS = line.strip().split()
            elts = list(map(int, eltsS[:3])) + list(map(float, eltsS[3:]))
            if elts[0] < iter: continue
            if elts[0] > iter: break
            if gcd == 0:
                gcd = elts[1]
            else:
                a, b = elts[1], gcd
                while b > 0:
                    a, b = b, a % b
                gcd = a
            self.data.append( elts )
        self.step = gcd

    def smooth(self, window, alpha):
        """ Smooths the posterior recombination using windiw size 'window',
            and assign a fraction alpha to the posterior, 1-alpha to the uniform prior """

        self.window = window
        self.smoothed_data = []
        for leaf in range(self.leaves):
            self.smoothed_data.append( self._smooth_leaf( leaf, alpha ) )

    def write_data(self, outfile):
        headerline = "locus\tsize\trecomb_rate" + ''.join(["\t{}".format(leaf+1) for leaf in range(self.leaves)]) + "\n"
        outfile.write(headerline)
        start, curvalues = None, None
        for idx in range(len(self.smoothed_data[0])):
            values = [ self.smoothed_data[leaf][idx] for leaf in range(self.leaves) ]
            if values != curvalues:
                if start != None:
                    self._write_line(outfile, curvalues, start, idx)
                start, curvalues = idx, values
        self._write_line(outfile, curvalues, start, idx+1)
            
    def _write_line(self, outfile, curvalues, start, idx):
        rate = sum(curvalues)
        curvalues = [ v / rate for v in curvalues ]
        dataline = "{}\t{}\t{:12.6e}".format(start * self.step,
                                             (idx - start) * self.step,
                                             rate)
        dataline += ''.join(["\t{:5.3f}".format(v) for v in curvalues]) + "\n"
        outfile.write(dataline)
            
    def _encode(self, value):
        if value < 0.1 * self.rate / self.leaves:
            return 0
        elif value < self.rate / self.leaves:
            return 1
        else:
            return 2
        
    def _forward_data_generator(self, encode=False):
        """ convert into low - medium - high codes """
        for elts in self.data:
            locus, size, opp = elts[1:4]
            if encode:
                ratecodes = [ self._encode(count / opp) for count in elts[4:] ]
            else:
                ratecodes = elts[4:]
            for loc in range(locus, locus+size, self.step):
                yield [loc, 100, opp] + ratecodes
    
    def _smooth_leaf(self, leaf, alpha):
        """ smooth using an oil-bleed algorithm """
        leafidx = leaf + 3
        shoulderbins = int(self.window / self.step)
        bins = int(self.size / self.step)
        data = [ d[leafidx] for d in self._forward_data_generator( encode=True ) ]

        # bleed out the high recombination bins by the required window size + 1
        for i in range( shoulderbins ):
            self._bleed( data, 2 )
        # bleed out the low recombination bins by 1
        self._bleed( data, 0 )

        # find streaks of low or medium, and or high rec. rates, and calc avg rates
        startidx = None
        curstate = None
        opp, recombs = 0, 0
        rates = []
        for idx, datum in enumerate(self._forward_data_generator( encode=False )):
            newstate = data[idx] == 2
            if newstate != curstate:
                # assign avg rate to previous streak
                if startidx != None:
                    newrate = alpha * (recombs / opp) + (1.0-alpha)*self.rate
                    for i in range(startidx, idx):
                        rates.append( newrate )
                startidx, curstate, opp, recombs = idx, newstate, 0, 0
            opp += datum[1] * datum[2]
            recombs += datum[1] * datum[leafidx]
        newrate = alpha * (recombs / opp) + (1.0-alpha)*self.rate
        for i in range(startidx, idx+1):
            rates.append( newrate )
        return rates

    def _bleed(self, data, state):
        start = None
        for idx, s in enumerate(data):
            if s == state:
                if start == None:
                    # start of streak
                    start = idx
            else:
                if start != None:
                    # end of streak -- implement bleed
                    if start > 0:
                        data[start - 1] = state
                    start = None
                    data[idx] = state
        if start != None and start > 0:
            data[start - 1] = state


if __name__ == "__main__":
    
    filename = sys.argv[1]

    lr = LocalRecombination( filename )

    lr.smooth( 200, 0.8 )

    lr.write_data( sys.stdout )


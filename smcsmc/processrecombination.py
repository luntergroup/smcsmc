from __future__ import print_function

import sys
import gzip
import math
import heapq
import bisect
import itertools


# seems to work nicely; for beta=4 segment count is about half number of recombinations
# next: sub-segment the leaves; write averages

# NOTE: will not work properly with the newfangled output files that include time-weighted counts!


class LocalRecombination:
    def __init__(self, infile):
        self.data = self._read_data(infile)
        self.size = self.data[-1][1] + self.data[-1][2]
        self.rate = self._calculate_rate()
        self.leaves = len(self.data[0][4:])

    def _calculate_rate(self):
        opportunity, recombinations = 0, 0
        for data in self._unmerge_data_generator():
            opportunity += data[1] * data[2]
            recombinations += data[1] * sum(data[3:])
        return recombinations / opportunity

    def _cusum(self, leaf=None):
        cusum = []
        curate = 0.0
        curate2 = 0.0
        for data in self._unmerge_data_generator():
            if leaf == None:
                datum = sum(data[3:]) / data[2] - self.rate
            else:
                datum = data[3 + leaf] / data[2] - self.rate / self.leaves
            curate += datum
            cusum.append(curate)
        return cusum

    def _unmerge_data_generator(self):
        # generates opportunity and per-leaf rates from input data, but in uniformly-sized windows
        for elts in self.data:
            locus, size, opp = elts[1:4]  # per-nt opportunity
            counts = elts[4:]  # per-nt counts
            for loc in range(locus, locus + size, self.step):
                yield [loc, self.step, opp] + counts

    def _smooth_column(self, B, leaf=None):
        bidx, previdx, acc = 0, 0, 0
        for idx, data in enumerate(self._unmerge_data_generator()):
            if bidx < len(B) and idx == B[bidx]:
                # found breakpoint; report last segment
                for sidx in range(previdx, idx):
                    yield acc / (idx - previdx)
                # prepare for next segment
                previdx, acc, bidx = idx, 0, bidx + 1
            opp = data[2]
            if leaf == None:
                acc += sum(data[3:]) / opp
            else:
                acc += data[3 + leaf] / opp
        for sidx in range(previdx, idx + 1):
            yield acc / (idx + 1 - previdx)

    def _read_data(self, infile, iter=0):
        data = []
        gcd = 0
        if infile.upper().endswith(".GZ"):
            rfile = gzip.open(infile, "r")
        else:
            rfile = open(infile, "r")
        for line in rfile:
            if line.startswith("iter"):
                continue
            eltsS = line.strip().split()
            elts = list(map(int, eltsS[:3])) + list(map(float, eltsS[3:]))
            if elts[0] < iter:
                continue
            if elts[0] > iter:
                break
            if gcd == 0:
                curpos = elts[0]
                gcd = elts[2]
            else:
                if curpos != elts[1]:
                    raise ValueError(
                        "Found gaps or overlaps in input file, line '{}'".format(
                            line.strip()
                        )
                    )
                a, b = elts[2], gcd
                while b > 0:
                    a, b = b, a % b
                gcd = a
            curpos += elts[2]
            data.append(elts)
        self.step = gcd
        rfile.close()
        return data

    # requires self.smoothed_data
    def write_data(self, outfile):
        headerline = (
            "locus\tsize\trecomb_rate"
            + "".join(["\t{}".format(leaf + 1) for leaf in range(self.leaves)])
            + "\n"
        )
        outfile.write(headerline)
        start, curvalues = None, None
        for idx in range(len(self.smoothed_data[0])):
            values = [self.smoothed_data[leaf][idx] for leaf in range(self.leaves)]
            if values != curvalues:
                if start != None:
                    self._write_line(outfile, curvalues, start, idx)
                start, curvalues = idx, values
        self._write_line(outfile, curvalues, start, idx + 1)

    def _write_line(self, outfile, curvalues, start, idx):
        rate = sum(curvalues)
        curvalues = [v / (rate + 1e-30) for v in curvalues]
        dataline = "{}\t{}\t{:9.3e}".format(
            start * self.step, (idx - start) * self.step, rate
        )
        dataline += "".join(["\t{:5.3f}".format(v) for v in curvalues]) + "\n"
        outfile.write(dataline)

    # minimize |X^b_{s,e}|.  Note, using half-open segment convention for [s,e)
    def argmax_xbse(self, s, e, cusum):
        n = float(e - s)
        sumleft = 0
        previous_cusum = 0 if s == 0 else cusum[s - 1]
        sumright = cusum[e - 1] - previous_cusum
        maxsum, maxb = -1, -1
        # loop to find change point (first point of right-hand segment)
        for b in range(s + 1, e):
            # update sumleft and sumright
            here = cusum[b - 1] - previous_cusum
            previous_cusum = cusum[b - 1]
            sumleft += here
            sumright -= here
            f1 = math.sqrt((e - b) / (n * (b - s)))
            f2 = math.sqrt((b - s) / (n * (e - b)))
            xbse = f1 * sumleft - f2 * sumright
            # we can compute the square of the difference, without the square roots, but hey
            if abs(xbse) > maxsum:
                maxsum = abs(xbse)
                maxb = b
        return maxsum, maxb  # segments are [s,b) and [b,e)

    # smooth using Wild Binary Segmentation (WBS) algorithm -- stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf
    def _wbs(self, cusum, beta, B=None):
        # change points
        if B == None:
            B = []
        # build test segments - not randomly though, that seems silly
        testsegs = []
        for l in [
            2,
            3,
            4,
            6,
            9,
            13,
            20,
            30,
            40,
            60,
            90,
            130,
            200,
            300,
            400,
            600,
            900,
            1300,
            2000,
        ]:
            for s in range(0, len(cusum), l / 2):
                if s + l < len(cusum):
                    testsegs.append((s, s + l))
        # add test segments from B
        for s, e in itertools.izip([0] + B, B + [len(cusum)]):
            testsegs.append((s, e))
        # build data structure
        F = []
        for s, e in testsegs:
            value, b = self.argmax_xbse(s, e, cusum)
            F.append((-value, b, s, e))
        heapq.heapify(F)
        # main loop of WBS algorithm
        while len(F) > 0:
            value, bk, s, e = heapq.heappop(F)
            if -value < beta * self.rate:
                break
            # check that [s,e) does not contain any change points that have already been identified
            idx1 = bisect.bisect_right(B, s)
            idx2 = bisect.bisect_left(B, e)
            if idx1 != idx2:
                continue
            # found good breakpoint
            B.append(bk)
            B.sort()
            # print("X added ",bk,"; ",len(B)," brkpts; ",len(F)," segs remaining; zeta=",-value," rate=",self.rate)
        return B

    def smooth(self, alpha, beta):
        assert alpha >= 0 and alpha <= 1
        assert beta > 0
        # calculate change points for overall recombination rate
        cusum = self._cusum()
        B = self._wbs(cusum, beta)
        data = [self._smooth_column(B)]
        # iteratively calculate denser super-grid for the per-leaf rates
        Bprime = B[:]
        for leaf in range(self.leaves):
            cusum = self._cusum(leaf)
            Bprime = self._wbs(cusum, beta, Bprime)
        data += [self._smooth_column(Bprime, leaf) for leaf in range(self.leaves)]
        # store
        self.smoothed_data = [[] for _ in range(self.leaves)]
        for allcols in itertools.izip(*data):
            smoothrate = allcols[0]
            relrates = [r / (sum(allcols[1:]) + 1e-30) for r in allcols[1:]]
            smoothedrates = [
                alpha * (r * smoothrate) + (1 - alpha) * self.rate / self.leaves
                for r in relrates
            ]
            for idx, r in enumerate(smoothedrates):
                self.smoothed_data[idx].append(r)
        return


if __name__ == "__main__":
    # just a quick example of usage:
    try:
        filename = sys.argv[1]
        alpha = float(sys.argv[2])  # mixin fraction
        beta = float(sys.argv[3])  # min ratio of zeta and recombination rate
    except:
        print("Usage: {} <file.recomb.gz> alpha beta".format(sys.argv[0]))
        print("\n  alpha (0-1) : proportion of posterior to mix in; 0 = use prior")
        print(
            "  beta (>0)   : smoothing factor; 1 is ok; use larger for more smoothing"
        )
        sys.exit(1)
    lr = LocalRecombination(filename)
    lr.smooth(alpha, beta)
    lr.write_data(sys.stdout)

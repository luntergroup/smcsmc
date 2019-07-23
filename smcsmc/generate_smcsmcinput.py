#!/usr/bin/env python3

#
# Built on Stefan Schiffel's generate-multihetsep.py, part of msmc-tools (commit 12758d9)
#

import sys
import gzip
import string
import copy
import argparse
import io
import collections
import os
import logging
import pdb

class MaskIterator:
    def __init__(self, filename, negative=False):
        if filename[-3:] == ".gz":
            self.file = io.TextIOWrapper(gzip.open(filename, "r"))
        else:
            self.file = open(filename, "r") #io.TextIOWrapper(open(filename, "r"))
        self.eof = False
        self.lastPos = 1
        self.negative = negative
        self.readLine()
  
    def readLine(self):
        try:
            line = next(self.file)
            fields = line.strip().split()
            if len(fields) == 2:
                self.start = int(fields[0])
                self.end = int(fields[1])
            else:
                self.start = int(fields[1]) + 1
                self.end = int(fields[2])
        except StopIteration:
            self.eof = True
  
    def getVal(self, pos):
        assert pos >= self.lastPos
        self.lastPos = pos
        while pos > self.end and not self.eof:
            self.readLine()
        if pos >= self.start and pos <= self.end:
            return True if not self.negative else False
        else:
            return False if not self.negative else True

class MergedMask:
    def __init__(self, mask_iterators):
        self.maskIterators = mask_iterators
  
    def getVal(self, pos):
        if len(self.maskIterators) == 0:
            return True
        else:
            return all((m.getVal(pos) for m in self.maskIterators))

class VcfIterator:
    def __init__(self, filename):
        self.file = io.TextIOWrapper(gzip.open(filename, "r"))
    
    def __iter__(self):
        return self
    
    def __next__(self):
        line = next(self.file)
        while line[0] == "#":
            line = next(self.file)
        fields = line.strip().split()
        chrom = fields[0]
        pos = int(fields[1])
        alleles = [fields[3]]
        for alt_a in fields[4].split(","):
            alleles.append(alt_a)
        geno = fields[9][:3]
        phased = geno[1] == "|"
        return (chrom, pos, tuple(alleles), (int(geno[0]), int(geno[2])), phased)

class OrderedAlleles:
    def __init__(self):
        self.ordered_alleles = []
    
    def addGenotype(self, a1, a2, phasing, ref):
        self.ref = ref
        if len(self.ordered_alleles) == 0:
            self.ordered_alleles = [[a1, a2]]
            if not phasing and a1 != a2:
                self.ordered_alleles.append([a2, a1])
        else:
            new = []
            for o in self.ordered_alleles:
                new.append(o + [a1, a2])
                if not phasing and a1 != a2:
                    new.append(o + [a2, a1])
            self.ordered_alleles = new
  
    def phase(self, trio):
        child, father, mother = trio
        new = [] 
        for o in self.ordered_alleles:
            child_1, child_2 = o[2 * child : 2 * (child + 1)]
            pat = o[2 * father]
            mat = o[2 * mother]
            if child_1 == pat and child_2 == mat or child_2 == pat and child_1 == mat:
                new.append(o)
        if len(new) > 0:
            self.ordered_alleles = new
        self.ordered_alleles = unique(self.ordered_alleles)
    
    def getPrint(self, trios):
        child_indices = []
        for i in range(len(self.ordered_alleles[0])):
            for child, father, mother in trios:
                if i == 2 * child or i == 2 * child + 1:
                    child_indices.append(i)
        print_indices = [i for i in range(len(self.ordered_alleles[0])) if i not in child_indices]
        stripped_alleles = unique([[o[i] for i in print_indices] for o in self.ordered_alleles])
        # convert to smcsmc form
        smcsmc_allele = [ "10//"[ 2*(len(set(o[i] for o in stripped_alleles))-1) + (stripped_alleles[0][i] == self.ref) ]
                          for i in range(len(stripped_alleles[0])) ]
        return ''.join(smcsmc_allele)

def unique(list_of_lists):
    return list(set([tuple(l) for l in list_of_lists]))

class JoinedVcfIterator:
    def __init__(self, filenames, trios):
        self.vcfIterators = [VcfIterator(f) for f in filenames]
        self.current_lines = [next(v) for v in self.vcfIterators]
        self.trios = trios
    
    def __iter__(self):
        return self
    
    def __next__(self):
        minIndices = self.getMinIndices()
        chrom = self.current_lines[minIndices[0]][0]
        pos = self.current_lines[minIndices[0]][1]
        ref = self.current_lines[minIndices[0]][2][0]
      
        ordered_alleles = OrderedAlleles()
        
        for i, l in enumerate(self.current_lines):
            if i not in minIndices:
                ordered_alleles.addGenotype(ref, ref, True, ref)
            else:
                alleles, geno, phased = l[2:5]
                ordered_alleles.addGenotype(alleles[geno[0]], alleles[geno[1]], phased, ref)
                try:
                    self.current_lines[i] = next(self.vcfIterators[i])
                except StopIteration:
                    self.current_lines[i] = None
        for trio in self.trios:
            ordered_alleles.phase(trio)
        return (chrom, pos, ordered_alleles.getPrint(self.trios))
    
    def getMinIndices(self):
        activeLines = [(i, l) for i, l in enumerate(self.current_lines) if l]
        if len(activeLines) == 0:
            raise StopIteration
        if len(activeLines) == 1:
            return [activeLines[0][0]]
        else:
            minIndices = [activeLines[0][0]]
            minPos = activeLines[0][1][1]
            for a in activeLines[1:]:
                if a[1][1] == minPos:
                    minIndices.append(a[0])
                if a[1][1] < minPos:
                    minPos = a[1][1]
                    minIndices = [a[0]]
            return minIndices
    
 
def run_multihetsep(files, output, mask = None, negative_mask = None, trio = None, minsize = 1000):
    """
    This is the right one

    .. todo: Docuemnt this
    """

    fout = gzip.GzipFile(output, 'w')
    
    trios = []
    if trio is not None:
        trios = [tuple(map(int, t.split(","))) for t in trio]
 
    nrIndividuals = len(files)
    nrHaplotypes = 2 * (nrIndividuals - len(trios))
 
    joinedVcfIterator = JoinedVcfIterator(files, trios)
    
    maskIterators = []
    if mask is not None:
        for f in mask:
            if len(maskIterators) < nrIndividuals:
                logging.debug("adding mask for individual {}: {}\n".format(len(maskIterators)+1,f))
            else:
                logging.debug("adding mask: {}\n".format(f))
            maskIterators.append(MaskIterator(f))
    if len(maskIterators) < nrIndividuals and len(maskIterators) > 0:
        raise ValueError("Must have at least {} masks\n".format(nrIndividuals))
    if negative_mask is not None:
        for nm in negative_mask:
            sys.stderr.write("adding negative mask: {}\n".format(nm))
            maskIterators.append(MaskIterator(nm, True))

    mergedMask = MergedMask(maskIterators)
    #pdb.set_trace()
    def is_segregating(alleles):
        orders = alleles.split(",")
        for o in orders:
            o = [a for a in o if a != '.']
            if len(o)>0:
                for a in o[1:]:
                    if a != o[0] or a == '/':
                        return True
        return False

    def pattern(mi, mm, pos):
        if len(mi) == 0:
            return tuple([True for i in range(nrIndividuals)])
        else:
            if mm.getVal(pos):
                return tuple([m.getVal(pos) for m in mi[:nrIndividuals]])
            else:
                return (False,) * nrIndividuals

    def dup(lst):
        return [ lst[i//2] for i in range(2*len(lst)) ]

    #pdb.set_trace()
    pos = 1
    last_pos = 1
    nr_called_by_pattern = collections.defaultdict(int)
    for chrom, snp_pos, alleles in joinedVcfIterator:
        # sys.stderr.write("{}\t{}\t{}\n".format(chrom, snp_pos, alleles))
        while pos < snp_pos:
            pos += 1
            patt = pattern( maskIterators, mergedMask, pos )
            nr_called_by_pattern[patt] += 1
            if pos % 1000000 == 0:
                logging.debug("processing pos {}".format(pos), file=sys.stderr) 
        if any( patt ): 
            calledAllele = ''.join([['.',a][m] for a,m in zip(alleles,dup(patt))])  # replace masked alleles with '.'
            if is_segregating(calledAllele):
                removePatterns = [p for p,c in nr_called_by_pattern.items()         # find patterns that we do NOT want to output
                                  if c < minsize and any(p) and p != patt]          #   (but treat as missing instead)
                if len(removePatterns) == 1:
                    removePatterns = []                                             # no point removing just one pattern, which will be replaced by the empty pattern
                removePatterns = unique( [(False,)*nrIndividuals] + removePatterns )
                numEmpty = sum( nr_called_by_pattern[p] for p in removePatterns )
                for p in removePatterns:
                    del nr_called_by_pattern[p]
                if numEmpty > 0:
                    nr_called_by_pattern[ (False,)*nrIndividuals ] = numEmpty
                c = chrom 
                for p,count in nr_called_by_pattern.items():                        # output the patterns, including the empty record, but excluding the one for the mutation
                    if p != patt:
                        out = str(last_pos) + '\t' + str(count) + '\t' + ''.join(['.0'[m] for m in dup(p)]) + '\n'
                        fout.write(out.encode('ascii'))
                        #pdb.set_trace()
                        #print(last_pos, count, ''.join(['.0'[m] for m in dup(p)]), sep="\t")
                        last_pos += count
                # output mutation
                count = nr_called_by_pattern[ patt ]
                #print(last_pos, count, calledAllele, sep="\t")
                out = str(last_pos) + '\t' + str(count) + '\t' + str(calledAllele) + '\n'
                fout.write(out.encode('ascii'))
                last_pos += count
                assert last_pos == snp_pos
                nr_called_by_pattern = collections.defaultdict(int)            
      
      
def split_vcfs(input, vcfdir, key, chroms = range(1,23)):
    """
    Splits samples into VCFs by themselves.

    :param list of tuples input: A list of (vcf, sample_name) pairings.
    :param str vcfdir: Path to the temporary folder where you wish to store the individual VCFs.
    :param str key: Prefix for the files.
    :param list chroms: Chromosomes to process."""
    for chrom in chroms:
        for vcf, sample in input:
            fname =  "{}/tmp{}.{}.chr{}.vcf.gz".format(vcfdir, key, sample, chrom)

            try:
                try_open = gzip.GzipFile( fname, 'r')
                print("Found:\t\t", fname, "\tso not doing anything...")
                have_files = True
            except:
                have_files = False

            if not have_files:
                if not os.path.exists(vcfdir):
                    os.makedirs(vcfdir)

                fout = gzip.GzipFile (fname, 'w')
                fin = gzip.GzipFile( vcf.format(chrom), 'r')
                print("Reading:\t", vcf.format(chrom))

                cols = []

                for line in fin:
                    #try:
                     #   elts = line.strip().split('\t')
                    #except TypeError:
                    elts = line.decode().strip().split('\t')
                    if line.startswith(b'#CHROM'):
                            col = [i for i, e in enumerate(elts) if e == sample]
                            if len(col) == 0:
                                raise ValueError("Could not find individual {}".format(sample))
                            cols.append(col[0])
                    if line.startswith(b'#') and not line.startswith(b'#CHROM'):
                            fout.write(line)
                    else: 
                            # filter out hom ref calls, and indel calls 
                            if (not elts[cols[0]].startswith("0|0")) and '.' not in elts[cols[0]][:3] and len(elts[3]) == 1 and len(elts[4]) == 1:
                                 
                                fout.write('\t'.join(elts[:9] + [elts[cols[0]]] + ["\n"]).encode('ascii'))
                    





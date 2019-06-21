import tskit
import os
from .model import *
import pdb
import numpy as np
import pandas as pd
import glob
import gzip
import io

def prune_tree_sequence(tree_sequence_path, num_samples):
    """
    Subsamples from a Tree Sequence.

    Parameters
    ----------
    tree_sequence_path : str
        File path to tree sequence file produced by `ts.dump`.
    num_samples : int
        The number of samples to sample.

    Returns
    -------
    ts : tree_sequence
        A pruned tree sequence.

    See Also
    --------
    ts_to_seg
    
    This function is identical to the one in PopSim analysis repository.
    """
    ts = tskit.load(tree_sequence_path)
    if ts.num_samples > num_samples:
        subset = np.random.choice(ts.samples(), num_samples, replace=False)
        ts = ts.simplify(subset)
    return ts

def ts_to_seg(path, n = None):
    """
    Converts a tree sequence into a seg file for use by :code:`smcsmc.run_smcsmcs()`. This is especially
    useful if you are simulating data from :code:`msprime` and would like to directly 
    use it in :code:`smcsmc`. For details of how to do this, please see the tutorial on simulation using :code:`msprime`.

    Provide the path to the tree sequence, and the suffix will be replaced by :code:`.seg`. This code is adapted from PopSim.

    :param str path: Full file path to the tree sequence created by :code:`ts.dump`.
    :param list n: If more than one sample of haplotypes is being analysed simulateously, provide it here as a list. Otherwise, simply provide the number of haplotypes as a single-element list. 
    """

    if n is None:
        ts = tskit.load(pathe)
        dirr = os.path.dirname(path)
        filen = os.path.basename(path)
        sep = filen.split(".")
        output = os.path.join(dirr,".".join(sep) + ".seg")
        fi = open(output, "w")
        prev = 1
        cur = 0
        for var in ts.variants():
            cur = int(var.site.position)
            if cur > prev:
                geno = ''.join(map(str,var.genotypes))
                fi.write(f"{prev}\t{cur-prev}\t{geno}\n")
            prev = cur
        fi.close()
    else: 
        for sample_size in n:
            ts = smcsmc.utils.prune_tree_sequence(path, sample_size)
            dirr = os.path.dirname(path)
            filen = os.path.basename(path)
            sep = filen.split(".")
            chrom = sep[0]
            sep.insert(0,str(sample_size))
            output = os.path.join(dirr,".".join(sep) + ".seg")
            fi = open(output, "w")
            prev = 1
            cur = 0
            for var in ts.variants():
                cur = int(var.site.position)
                if cur > prev:
                    geno = ''.join(map(str,var.genotypes))
                    fi.write(f"{prev}\t{cur-prev}\t{geno}\n")
                prev = cur
            fi.close()
    return None

def dict_to_args(smcsmc_params):
    '''
    Converts a dictionary of arguments into smcsmc compatible input.

    Parameters
    ----------
    smcsmc_params : dict
        A dictionary of arguments. See documentation for details.

    Returns
    -------
    args : list of arguments
        A list of arguments suitable for processing by `smcsmc.Smcsmc`.

    See Also
    --------
    run_smcsmc : Use these arguments to run `smcsmc`.
    '''
    args = []
    [args.extend(['-' + k, v]) for k, v in smcsmc_params.items()]

    # Remove the empty flags to expose the booleans
    args = [arg for arg in args if arg != ''] 

    # Remove any entries turned off by None flags.
    idx = [i for i,x in enumerate(args) if x == "None" or x is None]
    idx = idx + [i - 1 for i in idx]

    # By deleting in revrese, avoid moving idxs
    idx.sort(reverse=True)
    for i in idx:
        del args[i]
 
    args = [a for b in args for a in b.split()]

    return(args)

def run_smcsmc(args):
    """
    Run SMCSMC.

    Parameters
    ----------
    args : dict
        A dictionary of arguments. See the documentation for a complete list of options. To see how they are parsed, look at `smcsmc.utils.dict_to_args`.

    Returns
    -------
    None

    See Also
    --------
    dict_to_args

    """
    args = dict_to_args(args)
    run = Smcsmc(args)
    run.print_help_and_exit()
    run.load_option_file()
    run.parse_opts()
    run.validate()
    run.process_segfiles()
    run.set_environment()
    run.define_chunks()
    run.validate_parameters()
    run.set_pattern()
    for em_iter in range(0, run.emiters+1):
        run.do_iteration(em_iter)
    run.merge_outfiles()


class Expectations:
    def __init__(self, height, bin):
        log_height = np.log(height)
        boundaries = np.linspace(0, log_height, bin)

        start = boundaries[:-1]
        end = boundaries[1:]

        self.df = pd.DataFrame( {
            "Start" : np.exp(start),
            "End" : np.exp(end),
            "Events": 0
            } )

    def count_events(self, path_to_tree):
        ts = tskit.load(path_to_tree)
        out = []
        for tree in ts.trees():
            for noden in tree.nodes():
                out.append(ts.node(noden).time)
        return out

    def update_table(self, counts, gen_time=29):
        for i in counts:
            if i > 1:
                self.df.Events[max(self.df[self.df.Start <= i*gen_time].index)] += 1


    def update(self, paths):
        for path in glob.glob(paths):
            counts = self.count_events(path)
            self.update_table(counts)


class VcfIter:
    def __init__(self, filename):
        """
        Initiate a VCF class for reading raw genotype data.

        :param str filename: Path to the (optionally GZIPPED) VCF file. 
                            We support VCF versions above 4.0
        :return iterator: An iterator over the lines of the VCF in a structured way."""
        self.file = io.TextIOWrapper(gzip.open(filename, mode = 'rb'))

    def __iter__(self):
        return self

    def __next__(self):
        """
        Next item in the VCF, return is a tuple where the elements, in order, are:

        0. chrom (str)
        1. pos (int)
        2. alleles (list of str)
        3. geno 1 (int) 
        4. geno 2 (int)
        5. phased (bool)

        """
        line = next(self.file) # Grab the next line of the VCF. Is this safe?
        while line[0] == "#":
            line = next(self.file)
        items = line.strip().split()

        chrom =         items[0]
        pos =           int(items[1]) 
        
        alleles =       [items[3]]
        for alt_a in items[4].split(","): 
            alleles.append(alt_a) 

        # There should always be at least one allele
        assert ( len(alleles) > 0 )

        geno =          items[9][:3] 
        phased = geno[1] == "|"

        return (chrom, pos, tuple(alleles), (int(geno[0]), int(geno[2])), phased)


class MaskIter:
    def __init__(self, filename):
        if filename[-3:] == ".gz":
            self.file = io.TextIOWrapper(gzip.open(filename, "r"))
        else:
            self.file = open(filename, "r")
        self.eof = False
        self.lastPos = 1
 
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
            return True 
        else:
            return False

class MergedMask:
    """Strictly stolen from Stephen Schiffels"""
    def __init__(self, mask_iterators):
        self.maskIterators = mask_iterators
  
    def getVal(self, pos):
        return all((m.getVal(pos) for m in self.maskIterators))

class MergedVcf:
    def __init__(self, vcf_iterators):
        self.VcfIterators = vcf_iterators

    def __next__(self):
        states = [next(vcf) for vcf in self.VcfIterators]




def vcf_to_seg(vcfs, masks, minsize=1000):
    """
    Takes given samples from given VCFs and creates seg files for 
    input to :code:`smcsmc.run_smcsmc()`. 

    Input VCFs do not have to phased, but you should be aware that 
    if they are not, this will decrease the effectiveness of the lookahead likelihood.

    Provide a list of VCFs to merge together and masks for each individual in bed format.

    :param list vcfs: List of paths to the VCFs which you wish to pull samples from.
    :param str mask: Directory of similarly named mask files. Currently not implemented.
   
    :return: Nothing.
    """ 
    nInd = len(vcfs)

    # Need to progress along the sequence 
    # Filter out hom ref and indels for simplicity

    vcf_iters = MergedVcf([VcfIter(file) for file in vcfs])
    mask_iters = MergedMask([MaskIter(file) for file in masks])

    pos = 1
    last_pos = 1

    for chrom, snp_pos, alelles in vcf_iters:
        while pos < snp_pos:
            pos += 1
            



    
    
    
    


import tskit
import os
from .model import *
import pdb
import numpy as np
import pandas as pd
import glob


def prune_tree_sequence(tree_sequence_path, num_samples):
    """
    take in a tree sequence, and a number of samples
    less than the number of samples in the tree, 
    then simplify the tree sequence on that subset.

    Identical to the one in msmc. 
    """
    ts = tskit.load(tree_sequence_path)
    if ts.num_samples > num_samples:
        subset = np.random.choice(ts.samples(), num_samples, replace=False)
        ts = ts.simplify(subset)
    return ts


def ts_to_seg(path, n = None):
    """
    take one .trees file and write out 
    path.seg which acts as a single input to smcsmc. 

    Based off of the msmc multihep version.
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
            ts = prune_tree_sequence(path, sample_size)
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

    Requires a list of arguements akin to optparse, basically as if
    they had been entered on the command line. I will get around to 
    making this cleaner, but I want to get it working first.
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


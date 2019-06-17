import tskit
import os
from .model import *
import pdb
import numpy as np
import pandas as pd
import glob
import gzip

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

def vcf_to_seg(vcfs, samples, key, vcfdir="tmp", mask = None, chroms = range(1,23)):
    """
    Takes given samples from given VCFs and creates seg files for 
    input to :code:`smcsmc.run_smcsmc()`. 

    .. warning:: 
        This function is not fully implemented yet and will not work as expected (or at all!)  

    Input VCFs do not have to phased, but you should be aware that 
    if they are not, this will decrease the effectiveness of the lookahead likelihood
    calculation. 

    Provide a list of VCFs to merge together and a list of sample names 
    from each to include. These must match up identically. 

    Additionally provide mask files, and an output directory to place
    the VCFs. By default all somatic chromosomes are anlaysed, but 
    you can change this by providing an iterator to the chroms
    argument.

    :param list vcfs: List of paths to the VCFs which you wish to pull samples from.
    :param list samples: List of sample names to extract from the VCFs. The order must be the same.
    :param str key: Name of your output seg files.
    :param str vcfdir: Serves as a staging ground for intermediate seg files.  
    :param str mask: Directory of similarly named mask files. Currently not implemented.
    :param iterable chroms: An iterator of the chromosomes you wish to convert. 

    :return: Nothing.
    """ 
    raise NotImplementedError
    assert (len(vcfs)==len(samples)), "Your lists of samples and VCFs must be the same length." 
    #for chrom in chroms:
    #    for vcf, sample in zip(vcfs, samples):
    #        fname =  f"{vcfdir}/tmp{key}.{sample}.chr{chrom}.vcf.gz"

    #        try:
    #            try_open = gzip.GzipFile( fname, 'r')
    #            have_files = True
    #            print(f"You have already created {fname}")
    #        except:
    #            have_files = False

    #        if not have_files:
    #            if not os.path.exists(vcfdir):
    #                os.makedirs(args.vcfdir)

    #            fout = gzip.GzipFile (fname, 'w')
    #            fin = gzip.GzipFile( vcf.format(chrom), 'r')
    #            print "Reading:\t", vcf.format(chrom)

    #            cols = []

    #            for line in fin:
    #                elts = line.strip().split('\t')
    #                if line.startswith('#CHROM'):
    #                        col = [i for i, e in enumerate(elts) if e == sample]
    #                        if len(col) == 0:
    #                            raise ValueError("Could not find individual {}".format(sample))
    #                        cols.append(col[0])
    #                if line.startswith('#') and not line.startswith('#CHROM'):
    #                        fout.write(line)
    #                else:
    #                        # filter out hom ref calls, and indel calls
    #                        if (not elts[cols[0]].startswith("0|0")) and len(elts[3]) == 1 and len(elts[4]) == 1:
    #                            fout.write('\t'.join(elts[:9] + [elts[cols[0]]] + ["\n"]))




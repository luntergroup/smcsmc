import tskit
import os
from .model import *
import pdb
import numpy as np
import pandas as pd
import glob
from smcsmc.trees2tskit import trees2tskit
from collections import namedtuple
from tqdm import tqdm


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


def ts_to_seg(path, n, mask_file=None):
    """
    Converts a tree sequence into a seg file for use by :code:`smcsmc.run_smcsmcs()`. This is especially
    useful if you are simulating data from :code:`msprime` and would like to directly
    use it in :code:`smcsmc`. For details of how to do this, please see the tutorial on simulation using :code:`msprime`.

    Provide the path to the tree sequence, and the suffix will be replaced by :code:`.seg`. This code is adapted from PopSim.

    :param str path: Full file path to the tree sequence created by :code:`ts.dump`.
    :param list n: If more than one sample of haplotypes is being analysed simulateously, provide it here as a list. Otherwise, simply provide the number of haplotypes as a single-element list.
    :param str mask_file: Path to a bed file that you wish to mask out.
    """

    def mask(mask_dict, chr):
        for d in mask_dict[chr]:
            yield d

    def get_next_m(iter):
        try:
            m = next(miter)
        except StopIteration:
            m = pd.Interval(left=0, right=1)

        return m

    if mask_file:
        mask_table = pd.read_csv(mask_file, sep="\t", header=None)
        mask_dict = {}
        for key in mask_table[0].unique():
            sub = mask_table[mask_table[0] == key]
            mask_dict[key] = pd.IntervalIndex.from_arrays(sub[1], sub[2])

    for sample_size in n:
        ts = smcsmc.utils.prune_tree_sequence(path, sample_size)
        dirr = os.path.dirname(path)
        filen = os.path.basename(path)
        sep = filen.split(".")
        chrom = sep[0]
        sep.insert(0, str(sample_size))
        output = os.path.join(dirr, ".".join(sep) + ".seg")
        fi = open(output, "w")
        prev = 1
        cur = 0
        last_masked = 0

        if mask_file:
            miter = mask(mask_dict, chrom)

            m = get_next_m(miter)
            last_m = pd.Interval(left=0, right=1)

            for var in ts.variants():
                cur = int(var.site.position)

                if m.left > prev:
                    if m.right < cur:
                        geno = "." * sample_size
                        # Write the missing segment
                        right = m.left - prev
                        fi.write(f"{m.left}\t{m.left-prev}\t{geno}")
                        # And fill in the missing bit with the previous
                        # genotype
                        prev_geno = "".join(map(str, var.genotypes))
                        right = m.right - m.left
                        fi.write(f"{m.right}\t{right}\t{prev_geno}\n")

                        pdb.set_trace()
                        last_m = m
                        m = get_next_m(miter)

                    elif m.right > cur and m.left < cur:
                        geno = "." * sample_size
                        right = m.right - m.left
                        fi.write(f"{m.left}\t{right}\t{geno}\n")
                        cur = m.right
                        last_m = m

                        m = get_next_m(miter)

                if cur > prev and cur > last_m.right:

                    if last_m.right < cur and prev > last_m.left and last_m.left > 0:
                        geno = "".join(map(str, var.genotypes))
                        fi.write(f"{last_m.right}\t{cur-last_m.right}\t{geno}\n")
                        last_m = pd.Interval(left=0, right=1)

                    geno = "".join(map(str, var.genotypes))
                    fi.write(f"{prev}\t{cur-prev}\t{geno}\n")

                prev = cur
                prev_var = var

            fi.close()

        else:
            for var in ts.variants():
                cur = int(var.site.position)
                if cur > prev:
                    geno = "".join(map(str, var.genotypes))
                    fi.write(f"{prev}\t{cur-prev}\t{geno}\n")
                prev = cur
            fi.close()


def dict_to_args(smcsmc_params):
    """
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
    """
    args = []
    [args.extend(["-" + k, v]) for k, v in smcsmc_params.items()]

    # Remove the empty flags to expose the booleans
    args = [arg for arg in args if arg != ""]

    # Remove any entries turned off by None flags.
    idx = [i for i, x in enumerate(args) if x == "None" or x is None]
    idx = idx + [i - 1 for i in idx]

    # By deleting in revrese, avoid moving idxs
    idx.sort(reverse=True)
    for i in idx:
        del args[i]

    args = [a for b in args for a in b.split()]

    return args


def run_smcsmc(args):
    """
    The main entry point to :code:`smcsmc`, this function runs the whole analysis portion from start to finish. The one single argument is a dictionary of arguments as detailed on the :ref:`args`.

    .. tip::

        It's a really good idea to run :code:`smcsmc` in a :code:`tmux` session on the login node of your cluster if you are doing large analyses. The main loop takes very little memory, and spawns off cluster jobs if it is configured to do so (:ref:`cluster`).

    :param dict args: A dictionary of arguments as per :ref:`args`.

    An entirely equivalent entry point is the command line interface called :code:`smc2`. See :ref:`sec_cli` for examples.

    :code:`smcsmc` requires **segment files** as input. Below we detail three main ways to create them.

    - From VCF files
    - From :code:`msprime` sufficient tree sequences
    - From :code:`SCRM` simulations

    If you are using a different data type and would like help converting it to *seg* format, please let us know. For more details on converting files and interpreting the output of the algorithm, please see :ref:`getting_started`.
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
    for em_iter in range(0, run.emiters + 1):
        run.do_iteration(em_iter)
    run.merge_outfiles()


class Expectations:
    def __init__(self, height, bin):
        log_height = np.log(height)
        boundaries = np.linspace(0, log_height, bin)

        start = boundaries[:-1]
        end = boundaries[1:]

        self.df = pd.DataFrame(
            {"Start": np.exp(start), "End": np.exp(end), "Events": 0}
        )

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
                self.df.Events[max(self.df[self.df.Start <= i * gen_time].index)] += 1

    def update(self, paths):
        for path in glob.glob(paths):
            counts = self.count_events(path)
            self.update_table(counts)


def vcf_to_seg(
    input, output, masks=None, tmpdir=".tmp", key="tmp", chroms=range(1, 23)
):
    """
    This function converts data from a VCF to the segment type required for :code:`smcsmc`. You can also (optionally) include masks
    for the VCF (for example, from 1000 genomes) to indicate callable regions. This function first creates a number of
    intermediary VCF files in :code:`tmpdir`, identifiable by their :code:`key` and sample IDs. The function, by default, will loop over all
    chromosomes and create seperate :code:`seg.gz` files for each of them, however you can specify only a subset of chromosomes by providing a list
    to :code:`chroms`.

    It is preferable for VCFs to be phased, but it is by no means necessary. Phasing helps to improve the effectiveness of the lookahead likelihood.

    .. warning::
        This function does not take a lot of memory, but it can run for quite a while with genome-scale VCFs.

    The format of the input argument is as follows:

    .. code-block:: python

        input = [
            ("path_to_vcf", "sample_ID_1"),
            ("path_to_vcf", "sample_ID_2)
        ]

    For each individual that you would like to include in the segment file, you must specify its
    identifier and the path to its VCF. This means, in some cases, that you will be repeating the VCF paths. That's okay. Give the pairing of the VCF path and Sample ID (in that order) as a tuple, and give the list of tuples as the input argument. This is the same order in which individual's genotypes will appear in the seg file.

    **Chromosomes**

    *If you have one VCF with many chromsoomes*, specify the chromosomes that you want to use in the anlaysis in the :code:`chrom` option, as mentioned above. *If you have many VCFs for each chromosome*, you will have to run this function separately for each, specifying the correct chromsome each time.

    **Masks**

    Masks are bed files which indicate the callable regions from a VCF file, and may be included. If you do include masks, please make sure to include at least as many masks as individuals.

    **Example**

    If you wanted to convert two individuals from two different VCFs, whilst specifying masks, for chromosomes 5,6, and 7:

    .. code-block:: python

        smcsmc.vcf_to_seg(
            input = [
                ("/path/to/vcf1", "name_of_individual_1"),
                ("/path/to/vcf2", "name_of_individual_2")],
            output = "chrs567.seg.gz",
            masks = [
                "/path/to/mask1.bed.gz",
                "/path/to/mask2.bed.gz"],
            key = "chrs567",
            chroms = [5,6,7]
        )

    Which would create :code:`chr567.seg.gz` in the current working directory.

    :param list input: List of tuples, where the first element is the path to the VCF file and the second is the individual to be included.
    :param str output: Path to the output segment file. Gzipped if the suffix indicates so.
    :param list masks: List of masks for each of the individuals given in :code:`input`. These are **positive masks**. Masks are given as bed files, optionally gzipped.
    :param str tempdir: Directory to write intermediary vcf files. This is *not cleaned* after runs, so make sure you know where it is! This is to preserve the files for any further conversions.
    :param str key: Unique identifier of this run.
    :param list chroms: Either a list or range of chromosomes codes to use in this particular run.


    .. todo::

        It would be good to have a "cleanup = True" option to get rid of the intermediary files. This would be highly inefficient for rerunning but maybe worth it for some people.

    :return: Nothing.
    """
    smcsmc.generate_smcsmcinput.split_vcfs(input, tmpdir, key, chroms)

    for chr in chroms:
        # for n_sample in range(len(input)):
        names = [t[1] for t in input]
        file = [f"{tmpdir}/tmp{key}.{name}.chr{chr}.vcf.gz" for name in names]
        try:
            mask = [mask.format(chr) for mask in masks]
        except TypeError:
            mask = None
        smcsmc.generate_smcsmcinput.run_multihetsep(file, output, mask)


def convert_position(pos, map):
    # I'm not sure why this is happening
    if pos == 0.0:
        pos = 1.1
        # print("A weird position...")

    idx = map[map[1] < pos][1].idxmax()
    new_pos = int(round(pos - map[1][idx]))
    chr = map[0][idx].split(".")[-1][3:]
    return (new_pos, chr)


def find_segments(
    path,
    frm,
    to,
    time_range,
    hap,
    g=29,
    suffix=".trees.gz",
    pos_key=None,
    d=False,
    nomap=False,
    print_summary=False,
):
    """
    Find segments which have migrated from one population to another within a specified time range.

    :param str path: Path to the output folder which contains the trees of interest.
    :param int frm: Source population (backwards in time).
    :param int to: Sink population (backwards in time).
    :param tuple time_range: Tuple with (start, end) in years.
    :param int g: Years to a generation.
    :param str suffix: String to search for trees files.
    :param bool nomap: If there is no map."""
    Segment = namedtuple("Segment", "chr left right orientation")
    pattern = path + "*" + suffix
    files = glob.glob(pattern)

    if not nomap:
        if pos_key is None:
            key = pd.read_csv(path + "../merged.map", header=None, sep="\t")
        else:
            key = pd.read_csv(pos_key, header=None, sep="\t")
    else:
        key = pd.DataFrame({0: ["file", "file2"], 1: [1, np.Inf]})

    # pdb.set_trace()
    # subset = [Segment(convert_position(seg.left, key)[1],convert_position(seg.left, key)[0], convert_position(seg.right, key)[0], '+', seg.time, seg.descendants) for file in files for seg in trees2tskit(file, hap=hap,d = d).migrationlist if seg.source == frm and seg.dest == to and seg.time > time_range[0] and seg.time < time_range[1]]
    # segments = [segment for subsegments in segments for segment in subsegments]
    # subset = [Segment(convert_position(seg.left, key)[1],convert_position(seg.left, key)[0], convert_position(seg.right, key)[0], '+', seg.time)  for seg in segments if seg.source == frm and seg.dest == to and seg.time > time_range[0] and seg.time < time_range[1]]

    subset = []

    h = False
    if h:
        files = files[:2]

    for file in tqdm(files):
        for seg in trees2tskit(file, hap=hap, d=d).migrationlist:
            if (
                seg.source == frm
                and seg.dest == to
                and seg.time > time_range[0]
                and seg.time < time_range[1]
            ):
                subset = subset + [
                    Segment(
                        convert_position(seg.left, key)[1],
                        convert_position(seg.left, key)[0],
                        convert_position(seg.right, key)[0],
                        "+",
                    )
                ]

    df = pd.DataFrame(subset)

    if print_summary:
        print("Summary:")
        print("\tMean length: " + str(np.mean(df["right"] - df["left"])))
        print("\tTotal: " + str(np.sum(df["right"] - df["left"]) / 3e9))
        print("\tN: " + str(len(df)))
        print("\tTime: " + str(1 / (1e-8 * np.mean(df["right"] - df["left"]))))

    return df


def bed_to_marey(df, out):
    left = pd.DataFrame({"map": "Chromosome " + df["chr"], "phys": df["left"]})
    right = pd.DataFrame({"map": "Chromosome " + df["chr"], "phys": df["right"]})
    all = pd.concat([left, right])
    all.to_csv(out, sep=" ", index=False)


def smcsmc_available():
    try:
        return True
    except ImportError:
        return False

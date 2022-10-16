import argparse
import gzip
import os


def split_vcfs(input, vcfdir, key, chroms=range(1, 23)):
    """
    Splits samples into VCFs by themselves.

    :param list of tuples input: A list of (vcf, sample_name) pairings.
    :param str vcfdir: Path to the temporary folder where you wish to store the individual VCFs.
    :param str key: Prefix for the files.
    :param list chroms: Chromosomes to process."""
    for chrom in chroms:
        for vcf, sample in input:
            fname = "{}/tmp{}.{}.chr{}.vcf.gz".format(vcfdir, key, sample, chrom)

            try:
                try_open = gzip.GzipFile(fname, "r")
                print("Found:\t\t", fname, "\tso not doing anything...")
                have_files = True
            except:
                have_files = False

            if not have_files:
                if not os.path.exists(vcfdir):
                    os.makedirs(vcfdir)

                fout = gzip.GzipFile(fname, "w")
                fin = gzip.GzipFile(vcf.format(chrom), "r")
                print("Reading:\t", vcf.format(chrom))

                cols = []

                for line in fin:
                    elts = line.strip().split("\t")
                    if line.startswith("#CHROM"):
                        col = [i for i, e in enumerate(elts) if e == sample]
                        if len(col) == 0:
                            raise ValueError(
                                "Could not find individual {}".format(sample)
                            )
                        cols.append(col[0])
                    if line.startswith("#") and not line.startswith("#CHROM"):
                        fout.write(line)
                    else:
                        # filter out hom ref calls, and indel calls
                        if (
                            (not elts[cols[0]].startswith("0|0"))
                            and len(elts[3]) == 1
                            and len(elts[4]) == 1
                        ):
                            fout.write("\t".join(elts[:9] + [elts[cols[0]]] + ["\n"]))

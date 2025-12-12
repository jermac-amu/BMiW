import os
import argparse
import pathlib
import errno
import numpy as np
import re

from Bio import GenBank as GB

parser = argparse.ArgumentParser(prog="BaktaAnalyser",
                                    description="A tool for retrieving some simple statistics about results from annotation tool Bakta.")

parser.add_argument("-i",
                    "--i",
                    "--in",
                    type=pathlib.Path,
                    required=True,
                    help="Name of a .gbff file from Bakta output.",
                    dest="input")


# Verify that input files exist and can be read
#   args - An argparse.Namespace with arguments for this program.
def verify_files(args):
    inabsp = args.input.resolve()

    if not inabsp.exists():
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), str(inabsp))
    elif not inabsp.is_file():
        raise ValueError(f"File \"{str(inabsp)}\" is not a regular file")
    elif not os.access(inabsp, os.R_OK):
        raise PermissionError(f"Permission denied to read regular file \"{str(inabsp)}\"")


# Read a gbff file and retrieve basic information about its contents
#   file - A io.TextIOWrapper object exposing contents of a gbff file (usually output of open()).
#
#   returns: A numpy.ndarray of tuples with following contents:
#       gene start position (int),
#       gene end position (int),
#       flag informing whether the function of gene is known (bool),
def stats_gbff(file):

    gene_dat = np.empty(100, dtype=np.dtype("i8,i8,b1"))
    curr_gene = [-1,-1,False]
    num_genes = 0
    for record in GB.parse(file):
        for feat in record.features:
            if "gene" in feat.key:
                if num_genes > 0:
                    if num_genes > gene_dat.size: # We ran out of space in array, resize
                        gene_dat = np.resize(gene_dat, gene_dat.size*2)

                    gene_dat[[num_genes-1]] = tuple(curr_gene)
                    curr_gene = [-1,-1,False]

                num_genes += 1

                m = re.search(r"[0-9]+\.\.>?[0-9]+", feat.location)
                loc = m.group(0).replace(">", "")

                curr_gene[0], curr_gene[1] = [int(n) for n in loc.split("..")]

            elif "CDS" in feat.key:
                for qual in feat.qualifiers:
                    if "product" in qual.key:
                        if re.search(r"hypothetical\sprotein", qual.value) == None:
                            curr_gene[2] = True
                            break

    # Resize array to its final size
    gene_dat = np.resize(gene_dat,num_genes)

    gene_dat[[num_genes-1]] = tuple(curr_gene)

    return gene_dat


if __name__ == "__main__":
    args = parser.parse_args()

    verify_files(args)

    f = open(args.input)
    gene_dat = stats_gbff(f)
    f.close()

    num = len(gene_dat)
    known = sum([1 for s,e,k in gene_dat if k == True])
    unknown = num - known

    gene_lengths = sorted([e-s+1 for s,e,k in gene_dat])
    if num % 2 == 0:
        median = ( gene_lengths[num//2-1] + gene_lengths[num//2] ) / 2
    else:
        median = gene_lengths[num//2]

    print(f"Number of genes with known function: {known}")
    print(f"Number of genes with unknown function: {unknown}")
    print(f"Median of gene length: {median}")

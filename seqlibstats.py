import os
import argparse
import pathlib
import warnings
import errno
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(prog="SeqLibStats",
                                    description="A tool for retrieving some simple statistics of DNA sequencing results.")

parser.add_argument("-i",
                    "--i",
                    "--in",
                    action="extend",
                    nargs='+',
                    type=pathlib.Path,
                    required=True,
                    help="Names of one or more input fastq files.",
                    dest="input")

parser.add_argument("-H",
                    "--H",
                    "--histogram",
                    nargs="?",
                    const=pathlib.Path("./read_length_histogram.png"),
                    type=pathlib.Path,
                    required=False,
                    help="Generate a histogram of read lengths. Optionally, this argument can be followed by a file name of the histogram.",
                    dest="hist")

# Verify that input/output files exist and can be read/written
#   args - An argparse.Namespace with arguments for this program.
def verify_files(args):
    for filep in args.input: # Verify that input files exist and can be read

        absp = filep.resolve()
        if not absp.exists():
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), str(absp))
        elif not absp.is_file():
            raise ValueError(f"File \"{str(absp)}\" is not a regular file")
        elif not os.access(absp, os.R_OK):
            raise PermissionError(f"Permission denied to read regular file \"{str(absp)}\"")

    if args.hist is not None: # If histogram output file will be needed and exists, verify that it can be written to
        if args.hist.exists() and not os.access(args.hist, os.W_OK):
            raise PermissionError(f"Permission denied to write to regular file \"{str(args.hist.resolve())}\"")


# Read a FASTQ file and retrieve basic information about its contents
#   file - A io.TextIOWrapper object exposing contents of a FASTQ file (usually output of open()).
#   filename - (optional) Name of the file being read, used for warning/error messages.
#
#   returns: A tuple of:
#       number of reads (int),
#       array of read lengths (numpy.ndarray),
#       array of read GC contents (numpy.ndarray)
def stats_fastq(file, filename=""):

    read_lengths=np.empty(100, dtype=np.uint32) # array of read lengths; 32 bits is far more than enough for even largest singular reads
    gc_percents=np.empty(100, dtype=np.float32) # array of GC content percents

    num_reads = 0        # number of reads in file
    len_current = 0      # length of current read
    len_qual_current = 0 # length of quality string of current read
    gc_current = 0       # number of G and C in current read
    idx_current = 0      # index of line with current read's header
    location = None      # which FASTQ field is currently being read
    for idx, line in enumerate(file):
        line = line.strip().upper() # strip trailing whitespaces, uppercase to simplify recognition of G and C
        if len(line) == 0: continue # empty line

        elif (location is None or location == "quality") and line[0] == '@' and len_qual_current >= len_current: # beginning of next read (current line is header)
            if len_qual_current > len_current:
                warnings.warn(f"In file{' '+filename}: In line {idx_current}: Following read has sequence and quality strings of differing lengths. This is not valid in FASTQ format and might lead to inaccurate results.")

            idx_current = idx
            num_reads += 1
            if num_reads > read_lengths.size: # We ran out of space in arrays, resize
                read_lengths = np.resize(read_lengths, read_lengths.size*2)
                gc_percents = np.resize(gc_percents, gc_percents.size*2)

            if location is not None: # Only save previous read's stats if there was a previous read
                read_lengths[num_reads-2] = len_current
                gc_percents[num_reads-2] = gc_current/len_current

            len_current = 0
            len_qual_current = 0
            gc_current = 0

            location = "sequence" # we are about to read sequence

        elif location == "sequence":
            if line[0] == '+': # we are about to read quality
                location = "quality"
            else:              # still reading sequence
                len_current += len(line)
                gc_current += line.count('G') + line.count('C')

        else:
            len_qual_current += len(line)

    if len_qual_current > len_current:
        warnings.warn(f"In file{' '+filename}: In line {idx_current}: Following read has sequence and quality strings of differing lengths. This is not valid in FASTQ format and might lead to inaccurate results.")

    # Resize arrays to their final size
    read_lengths = np.resize(read_lengths,num_reads)
    gc_percents = np.resize(gc_percents, num_reads)

    read_lengths[num_reads-1] = len_current
    gc_percents[num_reads-1] = gc_current/len_current

    return (num_reads, read_lengths, gc_percents)


# Draw a histogram of read length distribution.
#   read_lengths - numpy.ndarray of read lengths
def read_length_histogram(read_lengths):
    plt.figure(figsize=(24,13.5), dpi=80)

    bins = min(100, (max(read_lengths) - min(read_lengths))//2)
    plt.hist(read_lengths, bins=bins, align="right", density=True)

    step = (max(read_lengths) - min(read_lengths)) // 25
    plt.xticks(np.arange(min(read_lengths), max(read_lengths)+1, step))

    plt.xlabel("Read length (bp)")
    plt.ylabel("Density")

    plt.suptitle("Histogram of read length distribution")


if __name__ == "__main__":
    args = parser.parse_args()

    verify_files(args)

    num_reads = 0
    read_lengths = None
    gc_percents = None

    for filep in args.input:
        filep = filep.resolve()
        print(f"Reading {filep}...")

        file = filep.open(mode='r')

        n, lths, gc = stats_fastq(file, filep)

        num_reads += n

        if read_lengths is not None:
            read_lengths = np.concatenate((read_lengths, lths), dtype=np.uint32)
        else:
            read_lengths = lths.copy()
        del lths

        if gc_percents is not None:
            gc_percents = np.concatenate((gc_percents, gc), dtype=np.float32)
        else:
            gc_percents = gc.copy()
        del gc

    read_lengths.sort(kind="quicksort")

    print(f"\nTotal number of reads: {num_reads:,}")
    print(f"Total number of nucleotides: {read_lengths.sum(dtype=np.uint64):,}")
    print(f"Mean read length: {read_lengths.mean(dtype=np.float32):,.2f}")
    if num_reads % 2 == 0:
        median = (read_lengths[num_reads//2-1] + read_lengths[num_reads//2])/2
    else:
        median = read_lengths[num_reads//2]

    print(f"Median read length: {median:,}")
    print(f"Mean GC pair content: {gc_percents.mean(dtype=np.float32)*100:.2f}%")

    if args.hist is not None:
        read_length_histogram(read_lengths)
        plt.savefig(args.hist)

        print(f"\nHistogram saved as \"{args.hist}\".")

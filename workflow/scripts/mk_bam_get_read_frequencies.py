# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Author : Katsantoni Maria
# Company: Mihaela Zavolan group, Biozentrum, Basel
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#   This script is part of the RCRUNCH pipeline. RCRUNCH pipeline detects
#   significantly enriched (in reads) peaks, which correspond to binding sites
#   of an RBP (CLIP experiment analysis).
#
#   This script, calculates the frequency of each read (or pair of reads),
#   using a bam file as input.
# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
import sys
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import pysam


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():
    """ Main function """

    __doc__ = "Obtain sliding windows of read coverage from bam file."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)
    # ----------------------------------------------------------
    parser.add_argument(
        "--bamfile",
        dest="bamfile",
        help="Reads alignment",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--paired",
        dest="paired",
        help="Reads alignment",
        required=True)

    parser.add_argument(
        "--outfile",
        dest="outfile",
        help="Output filename",
        required=True)

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # get the arguments
    # -------------------------------------------------------------------------

    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    bamfile = options.bamfile
    paired = int(options.paired)
    outfile = options.outfile
    coverage_profile(bamfile, paired, outfile)
    return


def coverage_profile(bamfile, paired, outfile):
    '''Calculate frequencies of reads'''
    read_frequencies = {}
    sam = pysam.AlignmentFile(bamfile, "rb")
    if paired == 1:
        for read in sam.fetch():
            if read.opt("NH") == 1:
                continue
            if read.qname in read_frequencies:
                continue
            else:
                read_frequencies[read.qname] = 1 / read.opt("NH")
    elif paired == 2:
        for read in sam.fetch():
            if read.opt("NH") == 1:
                continue
            if read.is_proper_pair:
                if read.is_read1:
                    if read.qname in read_frequencies:
                        continue
                    else:
                        read_frequencies[read.qname] = 1 / read.opt("NH")

    sam.close()
    read_frequencies_df = pd.Series(
        read_frequencies, name='Frequency')
    read_frequencies_df.to_csv(
        outfile,
        index=True,
        header=True,
        sep='\t')
    return


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)

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
#   This script, discards the reads mapping to the - strand in the case of the
#   transcriptome. Works for single-end and paired-end.
# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
import sys
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
        help="aligned reads (bam format)",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--outfile",
        dest="outfile",
        help="output filename",
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
    outfile = options.outfile
    discard_duplicates(bamfile, outfile)
    return


def discard_duplicates(bamfile, outfile):
    '''Calculate frequencies of reads'''
    sam = pysam.AlignmentFile(bamfile, "rb")
    filtered = pysam.AlignmentFile(outfile, "wb", template=sam)
    for read in sam.fetch():
        if read.is_duplicate:
            continue
        filtered.write(read)
    sam.close()
    filtered.close()
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

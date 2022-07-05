# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import sys

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------


def main():
    """ Main function """

    __doc__ = "Merge sliding windows of read coverage from bam files."

    parser = ArgumentParser(description=__doc__,
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--transcriptome",
        help="number of reads in rolling windows of transcriptome",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--genome",
        help="number of reads in rolling windows of genome",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--out",
        dest="out",
        help="Output filename",
        required=True,
        metavar="FILE")

    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    genome_df = pd.read_csv(
        options.genome,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')

    transcriptome_df = pd.read_csv(
        options.transcriptome,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')

    final_table = pd.concat(
        [genome_df, transcriptome_df],
        axis=0,
        join='outer',
        ignore_index=False)

    final_table.to_csv(options.out, sep='\t', index=False, header=True)

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

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Author : Katsantoni Maria
# Company: Mihaela Zavolan group, Biozentrum, Basel
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#   This script is part of the RCRUNCH pipeline. RCRUNCH pipeline detects
#   significantly enriched (in reads) peaks, which correspond to binding sites
#   of an RBP (CLIP experiment analysis).
#   This script creates the files as described in documentation called:
#   1. (training/test)set: contains a random half of the total peak sequences
#   2. (training/test)_bg:  contains 10 randomisations of the nucleotides
        # for each of the sequences in the test_set
#   3. (training/test)_pool: contains the test_set + the test_bg sequences
#   Background sequences are created by randomisation of the nucleotides of
#   each peak sequence. (How does that help in case of low comlexity regions?)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import os


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():
    """ Main function """

    __doc__ = "Divide the peaks into \
    training and test set and create Phylogibbs input files"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--phylogenetic_peaks",
        dest="phylogenetic_peaks",
        help="multiple alignments of the peak sequences",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--out_folder",
        dest="out_folder",
        help="output folder",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--prefix",
        dest="prefix",
        help="prefix name of the files",
        required=False,
        default=10)

    parser.add_argument(
        "--random_seq_per_peak_num",
        dest="random_seq_per_peak_num",
        help="number of random sequences to be created for \
        each peak sequence",
        required=False,
        default=10)
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

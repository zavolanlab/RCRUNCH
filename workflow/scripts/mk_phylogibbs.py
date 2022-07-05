# ----------------------------------------------------------------------
# RCRUNCH - Processing workflow of CLIP data
#         - Detection of binding peaks and motifs
# Author : Katsantoni Maria (modify and adapted for RCRUNCH)
# Company: Mihaela Zavolan group, Biozentrum, Basel
#
# Script: Run Phylogibbs

import sys
import numpy as np
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import subprocess


def main():
    """Run Phylogibbs to detect de novo motifs"""

    __doc__ = "Run Phylogibbs to detect de novo motifs of given size"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--peaks",
        help="fasta file with peak sequences",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--motif_length",
        help="motif size predicted by phylogibbs",
        required=True)

    parser.add_argument(
        "--phylogenetic_tree",
        help="phylogenetic tree of the organism at hand",
        required=False)

    parser.add_argument(
        "--outfile",
        help="output filename for the motif",
        required=True)

    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    initial_path = os.getcwd()
    peaks_path = os.path.join(initial_path, options.peaks)
    outfile = os.path.join(initial_path, options.outfile)
    if options.phylogenetic_tree:
        phylogibbs__phyl_process(
            peak_path=peaks_path,
            motif_length=options.motif_length,
            outfile=outfile,
            phylogenetic_tree=options.phylogenetic_tree)
    else:
        phylogibbs_process(
            peak_path=peaks_path,
            motif_length=options.motif_length,
            outfile=outfile)
    return


def phylogibbs_process(peak_path=None, motif_length=None, outfile=None):
    '''
     call phylogibbs to predict motifs of specific length in set of peaks
    '''
    peak_file = open(peak_path, 'r')
    total_peaks = 0
    for line in peak_file.readlines():
        if '>' in line:
            total_peaks += 1
    num_peak_successes = np.floor(total_peaks * 0.7)
    phylogibbs = [
        'phylogibbs',
        '-r',
        '-D', '0',
        '-m', str(motif_length),
        '-z', '2',
        '-y', str(int(num_peak_successes)),
        '-N', '1',
        '-f', peak_path,
        '-q',
        '-o', outfile]
    subprocess.call(phylogibbs)
    if not os.path.exists(outfile):
        out = open(outfile, 'w')
        out.close()
    return


def phylogibbs__phyl_process(peak_path=None, motif_length=None,
                             outfile=None, phylogenetic_tree=None):
    '''
     call phylogenetic phylogibbs to predict motifs
     of specific length in set of peaks
    '''
    peak_file = open(peak_path, 'r')
    total_peaks = 0
    for line in peak_file.readlines():
        if '>>' in line:
            total_peaks += 1
    num_peak_successes = np.floor(total_peaks * 0.7)
    phylogibbs = [
        'phylogibbs',
        '-r',
        '-D', '1',
        '-L', '"' + phylogenetic_tree + '"',
        '-m', str(motif_length),
        '-z', '2'
        '-y', str(num_peak_successes),
        '-N', '1',
        '-q',
        '-f', peak_path,
        '-o', outfile]
    subprocess.call(phylogibbs)
    if not os.path.exists(outfile):
        out = open(outfile, 'w')
        out.close()
    return


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)

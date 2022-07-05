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
#   This script merges the list of known motifs (e.g from AtTRACT database)
#   with the de novo detected motifs (by Phylogibbs) and the refined de novo
#   motifs (Motevo).
# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
import sys
import pandas as pd
import re
import io
import os
import numpy as np
from argparse import ArgumentParser, RawTextHelpFormatter


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():
    """ Main function """

    __doc__ = "Merge known and de novo motifs, TRANSFAC \
    (motevo compatible format)"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--denovo_motifs",
        dest="denovo_motifs",
        help="List of de novo motifs",
        required=True)

    parser.add_argument(
        "--known_motifs",
        dest="known_motifs",
        help="file with motifs in trasfac format",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--outfile",
        dest="outfile",
        help="output file with known and de novo \
        motifs in transfac format",
        required=True,
        metavar="FILE")

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

    # fka;fk
    output_file = open(options.outfile, 'w')
    with open(options.known_motifs, 'r') as input_file:
        for line in input_file.readlines():
            output_file.write(line)

    # for new_motif in new_motifs:
    #     output

def read_pwm(path):
    pwm = re.compile(r"^(\/\/\nNA.*\n([^\/]+))", re.MULTILINE)
    with open(path, 'r') as myfile:
        filein = str(myfile.read())
    motifs = []
    for match in pwm.finditer(filein):
        found = match.group(2)
        motif = pd.read_csv(
            io.StringIO(found), sep='\s+', comment='#', engine='python')

        # Drop non-informative columns from beginning and end
        # of the PWM
        information = motif['inf'].values
        index = motif.index.values
        non_inf = []
        counter = 0
        while information[counter] <= 0.25:
            non_inf.append(index[counter])
            counter += 1
        counter = -1
        while information[counter] <= 0.25:
            non_inf.append(index[counter])
            counter -= 1
        motif.drop(non_inf, inplace=True)
        motif.drop(['cons', 'inf'], axis=1, inplace=True)
        #######################################################

        motif['PO'] = np.arange(len(motif))
        motif.set_index('PO', drop=True, inplace=True)
        motif = motif.astype('float32')
        a = pd.Series(motif.sum(axis='columns'))
        a[a == 0] = 1
        motif = motif.divide(a, axis='rows')
        motif = motif * 100
        motif = motif.astype('int')
        motifs.append(motif)
    return motifs


# ___________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)


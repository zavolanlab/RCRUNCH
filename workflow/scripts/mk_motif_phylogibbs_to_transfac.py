# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# import needed (external) modules
# ----------------------------------------------------------------------------------------------------------------------
import sys
import pandas as pd
import re
import io
import os
from argparse import ArgumentParser, RawTextHelpFormatter

# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# Main function
# ----------------------------------------------------------------------------------------------------------------------


def main():
    """ Convert the two motifs detected by phylogibbs into motevo
    compatible motif format. These will be refined in the refinement
    mode of Motevo."""

    __doc__ = "Phylogibbs motif format to Motevo compatible format."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--phylogibbs_motifs",
        dest="phylogibbs_motifs",
        help="phylogibbs motif formated file",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--prefix",
        dest="prefix",
        help="prefix name for each of the output files",
        required=False,
        default="")

    parser.add_argument(
        "--out1",
        dest="out1",
        help="output file for the first motif",
        required=True)

    parser.add_argument(
        "--out2",
        dest="out2",
        help="output file for second motif",
        required=True)

    # __________________________________________________________________________________________________________________
    # ------------------------------------------------------------------------------------------------------------------
    # get the arguments
    # ------------------------------------------------------------------------------------------------------------------
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # -------------------------------------------------------

    motifs = read_pwm(options.phylogibbs_motifs)

    # Setups the template that the client requested with 3-5
    # rows of information
    # Followed by 0 blank rows and the dataframe
    template1 = """//\nNA """ + str(
        options.prefix) + """_1\n{}//"""

    with open(options.out1, 'w') as fp:
        fp.write(template1.format(motifs[0].to_csv(
            sep='\t',
            index=True, header=True)))

    template2 = """//\nNA """ + str(
        options.prefix) + """_2\n{}//"""
    with open(options.out2, 'w') as fp:
        fp.write(template2.format(motifs[-1].to_csv(
            sep='\t',
            index=True, header=True)))


def read_pwm(path):
    pwm = re.compile(r"^(\/\/\nNA.*\n([^\/]+))", re.MULTILINE)
    with open(path, 'r') as myfile:
        filein = str(myfile.read())
    motifs = []
    for match in pwm.finditer(filein):
        found = match.group(2)
        motif = pd.read_csv(
            io.StringIO(found), sep='\s+', comment='#', engine='python')
        cons = motif['cons'].values
        inf = motif['inf'].values
        motif.drop(['cons', 'inf'], axis=1, inplace=True)
        motif.set_index('PO', drop=True, inplace=True)
        motif = motif.astype('float32')
        a = pd.Series(motif.sum(axis='columns'))
        a[a == 0] = 1
        motif = motif.divide(a, axis='rows')
        motif = motif * 100
        motif = motif.round(3)
        motif['cons'] = cons
        motif['inf'] = inf
        motif.index = motif.index.map("{:02}".format)
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

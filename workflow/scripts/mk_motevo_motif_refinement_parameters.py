# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# import needed (external) modules
# ----------------------------------------------------------------------------------------------------------------------
import sys
from argparse import ArgumentParser, RawTextHelpFormatter

# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# Main function
# ----------------------------------------------------------------------------------------------------------------------

def main():
    """ Find enriched regions in comparison to a background sample.
        The enrichment is the outcome of significant binding of the RBP
        used as target for the eCLIP experiment."""

    __doc__ = "Obtain sliding windows of read coverage from bam file."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--phylogenetic_tree",
        dest="phylogenetic_tree",
        help="phylogenetic tree of the organism at hand",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--out",
        dest="out",
        help="output file containing the parameters needed \
        for running motevo in refinement mode",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--genome_tag",
        dest="genome_tag",
        help="genome name (eg:hg38)",
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

    phylogenetic_tree = open(options.phylogenetic_tree, 'r')
    phylogenetic_tree = str(phylogenetic_tree.read())
    # Create the parameters file for running Motevo in
    # motif refinement mode as explained in the manual of motevo
    parameters = open(options.out, 'w')
    parameters.write(
        '''refspecies ''' + options.genome_tag + '''
TREE ''' + phylogenetic_tree + '''

Mode WMREF

minposteriorWM 0.5
wmdiff 0.0001

markovorderBG 1
bgprior 0.99
bg A 0.25
bg T 0.25
bg G 0.25
bg C 0.25

restrictparses 0
minposterior 0.01 ''')

    parameters.close()



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

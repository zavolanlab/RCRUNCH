# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# import needed (external) modules
# ----------------------------------------------------------------------------------------------------------------------
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import os


# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# Main function
# ----------------------------------------------------------------------------------------------------------------------

def main():
    """ Trim uninformative columns (Information lower
    than 0.25 from de novo predicted motifs"""

    __doc__ = "Trim uninormative columns from predicted motifs."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--inmotif",
        dest="inmotif",
        help="motif to be trimmed",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--outmotif",
        dest="outmotif",
        help="trimmed motif",
        required=True,
        metavar="FILE")

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
    # read only the table info. ignore the other lines
    df = pd.read_csv(
        options.inmotif,
        header=1,
        sep='\t',
        index_col=0,
        comment='//',
        engine='python')
    # get the name included in the non-table line
    # use it again when i output the table
    myfile = open(options.inmotif, 'r')
    name = myfile.readlines()[1].strip()
    myfile.close()

    remove = []
    # get the index of the uninformative rows
    # at the beggining of the motif
    df.dropna(inplace=True)
    for index, row in df.iterrows():
        if float(row['inf']) <= 0.5:
            remove.append(index)
        else:
            break
    # get the index of the uninformative rows
    # at the end of the motif
    for index, row in df[::-1].iterrows():
        if float(row['inf']) <= 0.5:
            remove.append(index)
        else:
            break
    # drop the rows that were uninformative
    df.drop(remove, inplace=True)
    # reindex the table
    df.index = (np.arange(1, len(df) + 1))
    df.index.name = "P0"
    # round the decimals to a max of 3
    df[['A', 'C', 'G', 'T', 'inf']] = df[['A', 'C', 'G', 'T', 'inf']].round(3)
    # make the index values have a 0 in front if they have one digit
    df.index = df.index.map("{:02}".format)

    # output the file with the preceding and ending lines that
    # are not part of the table
    template = """//\n""" + str(name) + """\n{}//"""
    with open(options.outmotif, 'w') as fp:
        fp.write(template.format(df.to_csv(
            sep='\t', index=True, header=True)))


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

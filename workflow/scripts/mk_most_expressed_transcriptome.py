# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
import sys
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import re


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():
    """ Main function """

    __doc__ = "Obtain the set of most expressed transcripts of \
                each gene using output of Salmon."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--quantification_file",
        dest="quantification_file",
        help="Salmon quantification file",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--transcriptome",
        dest="transcriptome",
        help="Fasta file of the transcriptome",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--ensembl_csv",
        dest="ensembl_csv",
        help="Ensembl table of characteristics",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--ensembl_gtf",
        dest="ensembl_gtf",
        help="Ensembl table of characteristics",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--out_gtf",
        dest="out_gtf",
        help="modified gtf of customised transcriptome",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--out_fasta",
        dest="out_fasta",
        help="modified fasta of customised transcriptome",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--transcript_info",
        dest="transcript_info",
        help="customised transcriptome information",
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

    gene_info = pd.read_csv(
        options.ensembl_csv,
        header=0,
        sep='\t',
        index_col=0,
        comment='#',
        engine='python')
    gene_info = gene_info[['transcript_id', 'gene_id']]
    gene_info.dropna(axis=0, how='any', inplace=True)
    gene_info.drop_duplicates(subset=['transcript_id'], inplace=True)
    gene_info.set_index('transcript_id', drop=True, inplace=True)
    # find the most expressed transcript based on the quantification
    expression_table = pd.read_csv(
        options.quantification_file,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')
    expression_table['Name'] = expression_table['Name'].apply(
        lambda x: x.split('.')[0])
    expression_table.set_index('Name', drop=True, inplace=True)
    expression_table = expression_table.merge(
        gene_info,
        how='inner',
        left_index=True,
        right_index=True)
    expression_table['count_max'] = expression_table.groupby(
        ['gene_id'])['NumReads'].transform(max)
    most_expressed = expression_table[
        expression_table['NumReads'] == expression_table['count_max']]
    most_expressed['count_max'].fillna(0)
    most_expressed = most_expressed[most_expressed['NumReads'] >= 1]
    most_expressed.drop('count_max', axis=1, inplace=True)

    most_expressed_info = pd.DataFrame(most_expressed['Length'])
    # modified_names = [i.split('.')[0] for i in most_expressed_info.index.values]
    # most_expressed_info['name'] = modified_names
    # most_expressed_info.set_index('name', inplace=True, drop=True)

    most_expressed_info.to_csv(
        options.transcript_info,
        sep='\t',
        index=True,
        header=False)

    true_trs = modify_gtf(
        options.ensembl_gtf,
        most_expressed.index.values,
        options.out_gtf)

    modify_fasta(
        options.transcriptome,
        true_trs,
        options.out_fasta)
    return
# _________________________________________________________________________
# -------------------------------------------------------------------------
# Functions
# -------------------------------------------------------------------------


def modify_gtf(ensembl_file, highest_tr, outpath):
    true_transcripts = []
    outfile = open(outpath, 'w')
    with open(ensembl_file) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('#'):
            outfile.write("%s" % line)
        transcript = re.match('.*transcript_id\s\"(\w+)\".*', line)
        if transcript is None:
            continue
        if transcript.group(1) in highest_tr:
            outfile.write("%s" % line)
            true_transcripts.append(transcript.group(1))
    outfile.close()
    return true_transcripts


def modify_fasta(transcr, highest_tr, outpath):
    outfile = open(outpath, 'w')
    with open(transcr) as f:
        contents = f.read()
    lines = re.split('\n\>', contents)
    for line in lines:
        transcript = re.match('^\>?(\w+)\W', line)
        if transcript is None:
            continue
        if transcript.group(1) in highest_tr:
            if line.startswith('>'):
                name = str(line)
                # names = name.split(' ')
                # name1 = name.split('.')[0]
                # names[0] = name1
                # names = ' '.join(names)
                outfile.write("%s\n" % name)
            else:
                outfile.write(">%s\n" % line)
    outfile.close()
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

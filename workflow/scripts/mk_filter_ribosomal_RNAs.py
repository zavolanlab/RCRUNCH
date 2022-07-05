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
#   This script, discards the reads mapping to the ribosomal_rnas.
# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import pysam
from collections import defaultdict
import pandas as pd


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
        "--rRNAs",
        dest="rRNAs",
        help="rRNA mappings \
        to be excluded",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--rRNA_out",
        dest="rRNA_out",
        help="output filename",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--annotation",
        dest="annotation",
        help="ensembl annotation table",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--paired",
        dest="paired",
        help="paired",
        required=True)

    parser.add_argument(
        "--sense",
        dest="sense",
        help="sense",
        required=True)

    parser.add_argument(
        "--flag",
        dest="flag",
        help="flag to indicate successful run",
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

    bamfile = options.bamfile
    rRNA_out = options.rRNA_out
    annotation_path = options.annotation
    paired = options.paired
    sense = options.sense

    rRNAs_df = pd.read_csv(
        options.rRNAs,
        header=None,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')
    rRNAs = rRNAs_df[0].tolist()
    sys.stdout.write('Start\n')
    sys.stdout.flush()

    annotation = load_annotation(annotation_path)
    rRNA_df = get_ensembl_rRNAs(annotation)

    discard_rRNAs(
        bamfile=bamfile,
        rRNAs=rRNAs,
        ensembl_df=rRNA_df,
        rRNA_out=rRNA_out,
        paired=paired,
        sense=sense)

    myfile = open(options.flag, 'w')
    myfile.close()

    sys.stdout.write('Done!\n')
    sys.stdout.flush()
    return


def load_annotation(ensembl):
    df = pd.read_csv(
        ensembl,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        usecols=['seqname', 'feature', 'transcript_id', 'gene_biotype', 'strand', 'start', 'end'],
        engine='python')
    df[['start', 'end']] = df[['start', 'end']].astype('int')
    return df


def get_ensembl_rRNAs(df):
    df = df[(df['gene_biotype'] == 'rRNA') & (df['feature'] == 'exon')].copy(deep=True)
    return df


def discard_rRNAs(bamfile=None, rRNAs=None, ensembl_df=None,
                  rRNA_out=None, paired=None, sense=None):
    '''Find reads mapping to rRNAs'''
    paired = int(paired)
    sense = int(sense)
    bam = pysam.AlignmentFile(bamfile, "rb")
    rRNA_outfile = open(rRNA_out, 'w')
    reads = set()
    for rRNA in rRNAs:
        for read in bam.fetch(rRNA):
            reads.add(read.query_name)
    if paired == 1:
        if sense == 1:
            for index, row in ensembl_df.iterrows():
                chromosome = str(row['seqname'])
                start = int(row['start'])
                end = int(row['end'])
                strand = row['strand']
                try:
                    bam.fetch(chromosome, 1, 1)
                except:
                    continue
                for read in bam.fetch(chromosome, start, end):
                    if strand == '+':
                        if not read.is_reverse:
                            reads.add(read.query_name)
                    if strand == '-':
                        if read.is_reverse:
                            reads.add(read.query_name)
        elif sense == 2:
            for index, row in ensembl_df.iterrows():
                chromosome = str(row['seqname'])
                start = int(row['start'])
                end = int(row['end'])
                strand = row['strand']
                try:
                    bam.fetch(chromosome, 1, 1)
                except:
                    continue
                for read in bam.fetch(chromosome, start, end):
                    if strand == '+':
                        if read.is_reverse:
                            reads.add(read.query_name)
                    if strand == '-':
                        if not read.is_reverse:
                            reads.add(read.query_name)
    elif paired == 2:
        if sense == 1:
            for index, row in ensembl_df.iterrows():
                chromosome = str(row['seqname'])
                start = int(row['start'])
                end = int(row['end'])
                strand = row['strand']
                try:
                    bam.fetch(chromosome, 1, 1)
                except:
                    continue
                for read in bam.fetch(chromosome, start, end):
                    if strand == '+':
                        if (read.is_read1) and (not read.is_reverse):
                            reads.add(read.query_name)
                        elif (read.is_read2) and (read.is_reverse):
                            reads.add(read.query_name)
                    if strand == '-':
                        if (read.is_read1) and (read.is_reverse):
                            reads.add(read.query_name)
                        elif (read.is_read2) and (not read.is_reverse):
                            reads.add(read.query_name)
        elif sense == 2:
            for index, row in ensembl_df.iterrows():
                chromosome = str(row['seqname'])
                start = int(row['start'])
                end = int(row['end'])
                strand = row['strand']
                try:
                    bam.fetch(chromosome, 1, 1)
                except:
                    continue
                for read in bam.fetch(chromosome, start, end):
                    if strand == '+':
                        if (read.is_read1) and (read.is_reverse):
                            reads.add(read.query_name)
                        elif (read.is_read2) and (not read.is_reverse):
                            reads.add(read.query_name)
                    if strand == '-':
                        if (read.is_read1) and (not read.is_reverse):
                            reads.add(read.query_name)
                        elif (read.is_read2) and (read.is_reverse):
                            reads.add(read.query_name)
    names = '\n'.join(list(reads))
    rRNA_outfile.write(names)
    sys.stdout.write('Parsed once the bam\n')
    sys.stdout.flush()

    bam.close()
    rRNA_outfile.close()
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

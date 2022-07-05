# ----------------------------------------------------------------------
# RCRUNCH - Processing workflow of CLIP data
#         - Detection of binding peaks and motifs
# Author : Katsantoni Maria
# Company: Mihaela Zavolan group, Biozentrum, Basel
#
# Script: Remove non coding RNA categories

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import pysam
import pandas as pd
import numpy as np
from io import StringIO
from csv import writer
import HTSeq
from pandas.io.common import EmptyDataError


def main():
    """ Main function """

    __doc__ = "Obtain sliding windows of read coverage from bam file."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "--bamfile",
        help="aligned reads (bam format)",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--RNA_central",
        help="RNA central non coding RNAs",
        nargs='?',
        const='',
        required=False)

    parser.add_argument(
        "--outfile",
        help="output filename",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--paired",
        help="paired",
        required=True)

    parser.add_argument(
        "--sense",
        help="sense",
        required=True)

    parser.add_argument(
        "--ncRNAs",
        nargs='*',
        help="ncRNA biotypes to be excluded",
        required=False,
        default='')

    parser.add_argument(
        "--flag",
        help="flag to indicate successful run",
        required=True,
        metavar="FILE")

    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    bamfile = options.bamfile
    outfile = options.outfile
    rnacentral_path = options.RNA_central
    paired = options.paired
    sense = options.sense

    sys.stdout.write('Start\n')
    sys.stdout.flush()

    if (not rnacentral_path) | (len(options.ncRNAs) == 0):
        ncRNAs = open(outfile, 'w')
        ncRNAs.write('No ncRNAs\n')
        ncRNAs.close()
    else:
        ncRNA_rnacentral = rna_central_annotation(rnacentral_path, options.ncRNAs)
        ncRNA_rnacentral_nr = merge_overlapping_rrnas(ncRNA_rnacentral)
        discard_rRNAs(
            bamfile=bamfile,
            rna_central_df=ncRNA_rnacentral_nr,
            outfile=outfile,
            paired=paired,
            sense=sense)
    myfile = open(options.flag, 'w')
    myfile.close()
    sys.stdout.write('Done!\n')
    sys.stdout.flush()
    return


def rna_central_annotation(rna_central, ncRNAs):
    output = StringIO()
    csv_writer = writer(output)
    gff = HTSeq.GFF_Reader(rna_central)
    for gff_line in gff:
        if gff_line.type == 'transcript':
            if gff_line.attr['type'] in ncRNAs:
                line = pd.Series([
                    gff_line.iv.chrom,
                    gff_line.type,
                    gff_line.iv.start,
                    gff_line.iv.end,
                    gff_line.iv.strand])
                csv_writer.writerow(line)
    output.seek(0)
    try:
        df = pd.read_csv(output)
        df.columns = ['seqname', 'feature', 'start', 'end', 'strand']
    except EmptyDataError:
        df = pd.DataFrame()
    
    # df = pd.read_csv(
    #     rna_central,
    #     header=None,
    #     sep='\t',
    #     index_col=None,
    #     comment='#',
    #     engine='python')
    # df.columns = ['seqname', 'source', 'feature', 'start', 'end',
    #               'score', 'strand', 'frame', 'attribute']
    # df = df[(df['feature'].str.contains('transcript')) & (
    #     (df['attribute'].str.contains('type=rRNA')) |
    #     (df['attribute'].str.contains('type=snRNA')) |
    #     (df['attribute'].str.contains('type=tRNA')))].copy(deep=True)
    # df = df[['seqname', 'feature', 'start', 'end', 'strand']].copy(deep=True)
    # print('df', df)
    return df


def merge_overlapping_rrnas(df):
    if df.empty:
        return df
    df.sort_values(
        by=['strand', 'seqname', 'start', 'end'],
        ascending=[True, True, True, True],
        inplace=True)
    output = StringIO()
    csv_writer = writer(output)
    chromosome_bef = ''
    start_bef = 0
    end_bef = 0
    strand_bef = ''
    chr_range_bef = set()
    for index, row in df.iterrows():
        chromosome = str(row['seqname'])
        start = int(row['start'])
        end = int(row['end'])
        strand = row['strand']
        chr_range = set(list(np.arange(start, end + 1)))
        if ((chromosome == chromosome_bef) & (strand == strand_bef)):
            if chr_range.isdisjoint(chr_range_bef) is False:
                start_bef = np.min([start, start_bef, end, end_bef])
                end_bef = np.max([start, start_bef, end, end_bef])
                chr_range_bef = set(list(np.arange(start_bef, end_bef + 1)))
                continue
        if chromosome_bef:
            rrna = pd.Series([chromosome_bef, start_bef, end_bef, strand_bef])
            csv_writer.writerow(rrna)
        chromosome_bef = chromosome
        strand_bef = strand
        start_bef = start
        end_bef = end
        chr_range_bef = chr_range
    output.seek(0)
    output_df = pd.read_csv(output)
    output_df.columns = ['seqname', 'start', 'end', 'strand']
    print(output_df)
    return output_df


def discard_rRNAs(bamfile=None, rna_central_df=None,
                  outfile=None, paired=None, sense=None):
    '''
        Find reads mapping to non coding RNAs of specific categories
    '''
    paired = int(paired)
    sense = int(sense)
    bam = pysam.AlignmentFile(bamfile, "rb")
    ncRNAs = open(outfile, 'w')
    reads = set()
    if rna_central_df.empty:
        ncRNAs.write('No ncRNAs\n')
        ncRNAs.close()
        return

    if (paired == 1) & (sense == 1):
        conditional = '((strand == "+") & (not read.is_reverse)) | ' + \
            '((strand == "-") & read.is_reverse)'

    elif (paired == 1) & (sense == 2):
        conditional = '((strand == "+") & read.is_reverse) | ' + \
            '((strand == "-") & (not read.is_reverse))'

    elif (paired == 2) & (sense == 1):
        conditional = '((strand == "+") & (((not read.is_reverse) & ' + \
            'read.is_read1) |(read.is_reverse & read.is_read2))) | ' + \
            '((strand == "-") & ((read.is_reverse & read.is_read1) | ' + \
            '((not read.is_reverse) & read.is_read2)))'

    elif (paired == 2) & (sense == 2):
        conditional = '((strand == "+") & ((read.is_reverse & ' + \
            'read.is_read1) |((not read.is_reverse) & read.is_read2))) | ' + \
            '((strand == "-") & (((not read.is_reverse) & read.is_read1) | ' + \
            '(read.is_reverse & read.is_read2)))'

    for index, row in rna_central_df.iterrows():
        chromosome = str(row['seqname'])
        start = int(row['start'])
        end = int(row['end'])
        strand = row['strand']
        try:
            bam.fetch(chromosome, 1, 1)
        except:
            continue
        for read in bam.fetch(chromosome, start, end):
            if eval(conditional):
                reads.add(read.query_name)

    names = '\n'.join(list(reads))
    ncRNAs.write(names)
    sys.stdout.write('Parsed once the bam\n')
    sys.stdout.flush()
    bam.close()
    ncRNAs.close()
    return


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)

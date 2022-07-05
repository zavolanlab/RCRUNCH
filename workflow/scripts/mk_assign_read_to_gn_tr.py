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
#   In this script the same reads are mapped to the genome or transcriptome.
#   Using the AS scores found in STAR mappings, the read is kept only in the
#   trasncriptome /genome depending on the highest AS score. If the read has
#   the same score on the genome/transcriptome then it is preferentially
#   assigned to the transcriptome.
# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import numpy as np
import pandas as pd
import pysam
from io import StringIO
from csv import writer
import multiprocessing as mp
import os
import logging
import math
logger = mp.log_to_stderr(logging.DEBUG)


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
        "--annotation",
        dest="annotation",
        help="csv file of annotation",
        metavar="FILE",
        required=True)

    parser.add_argument(
        "--tr_bam",
        dest="tr_bam",
        help="reads aligned to transcriptome",
        metavar="FILE",
        required=True)

    parser.add_argument(
        "--gn_bam",
        dest="gn_bam",
        help="reads aligned to genome",
        metavar="FILE",
        required=True)

    parser.add_argument(
        "--out_folder",
        dest="out_folder",
        help="output folder for temp small bam files",
        required=True)

    parser.add_argument(
        "--filtered_tr_bam",
        dest="filtered_tr_bam",
        help="reads that align best to transcriptome",
        required=True)

    parser.add_argument(
        "--filtered_gn_bam",
        dest="filtered_gn_bam",
        help="reads that align best to genome",
        required=True)

    parser.add_argument(
        "--paired",
        dest="paired",
        help="paired end sequencing or single-end",
        required=True)

    parser.add_argument(
        "--sense",
        dest="sense",
        help="1 if mate1 is the sense, otherwise 2",
        required=True)

    parser.add_argument(
        "--chromosome_info",
        dest="chromosome_info",
        help="File containing chromosome names\
        and lengths provided during STAR indexing, required",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--threads",
        dest="threads",
        help="number of cores used",
        required=True)
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

    annotation = options.annotation
    tr_bam = options.tr_bam
    gn_bam = options.gn_bam
    filtered_gn_bam = options.filtered_gn_bam
    filtered_tr_bam = options.filtered_tr_bam
    paired = options.paired
    sense = options.sense
    out_folder = options.out_folder
    threads = options.threads

    # chromosome_info = pd.read_csv(
    #     options.chromosome_info,
    #     header=None,
    #     sep='\t',
    #     index_col=None,
    #     comment='#',
    #     engine='python')
    # chromosome_info.columns = ['chromosome', 'chromosome_length']
    # chromosome_info.set_index('chromosome', inplace=True, drop=True)

    annotation_info = dictionary_of_tables(annotation)

    names = []
    arguments = []
    for chromosome in annotation_info.keys():
        #     # if chromosome not in chromosome_info.index:
        #     #     continue
        #     chr_length = chromosome_info.at[chromosome, 'chromosome_length']
        # #     chr_middle = math.floor(chr_length / 2)
        #     splits = [[0, chr_middle], [chr_middle, chr_length]]
        # for i in np.arange(2):
        chromosome = str(chromosome)
        # start = splits[i][0]
        # end = splits[i][-1]
        name = os.path.join(out_folder, chromosome + '.bam')
        each_arguments = {}
        each_arguments['gn_bam'] = gn_bam
        each_arguments['tr_bam'] = tr_bam
        each_arguments['annotation'] = annotation_info[chromosome]
        each_arguments['chromosome'] = chromosome
        each_arguments['gn_out'] = name
        each_arguments['paired'] = paired
        each_arguments['sense'] = sense
        # each_arguments['start'] = start
        # each_arguments['end'] = end

        names.append(name)
        arguments.append(each_arguments)

    with mp.Pool(processes=len(annotation_info.keys())) as pool:
        gn_reads = []
        for gn_reads_per_chr in pool.imap_unordered(remove_genomic_reads, arguments):
            gn_reads.append(gn_reads_per_chr)

    merge_bamfiles(filtered_gn_bam, names, threads)
    total_gn_reads = pd.concat(gn_reads)

    tr_arguments = {}
    tr_arguments['tr_bam'] = tr_bam
    tr_arguments['genomic_reads'] = total_gn_reads
    tr_arguments['tr_out'] = filtered_tr_bam
    tr_arguments['paired'] = paired
    tr_arguments['sense'] = sense
    remove_transcriptomic_reads(tr_arguments)
    return


def merge_bamfiles(out_file_name, input_file_names, threads=1):
    threads = int(threads)
    if threads > 1:
        args = ["-@", str(threads)]
    else:
        args = []

    args.append(out_file_name)
    args.extend(input_file_names)
    pysam.merge(*args)


def load_annotation(path):
    '''read annotation csv file'''
    df = pd.read_csv(
        path,
        header=0,
        sep='\t',
        index_col=0,
        comment='#',
        usecols=['seqname', 'feature', 'transcript_id',
                 'strand', 'start', 'end'],
        engine='python')
    df[['start', 'end']] = df[['start', 'end']].astype('int')
    return df


def dictionary_of_tables(path):
    '''create a dictionary of tables,
    where each key is the chromosome and
    value is the annotation csv chromosome
    specific part'''
    annotation = load_annotation(path)
    annotation_dict = {}
    for i in annotation.groupby(['seqname']):
        name = str(i[0])
        annotation_dict[name] = i[1]
    return annotation_dict


def coordinate_list(x):
    '''function that returns a list of positions,
    based on start and end values'''
    if x[1] > x[2]:
        return np.arange(x[2], x[1] + 1)
    else:
        return np.arange(x[1], x[2] + 1)


def group_by_tr_id(df):
    '''merge all the features (exons, 3utrs etc) per transcript
    and create union of the coordinates '''
    df = df.groupby('transcript_id').agg(
        strand=('strand', 'first'),
        positions=("positions", lambda x: (list(set.union(*map(set, x))))))
    return df


def convert_tr_to_gn_coordinates(annotation_df):
    '''convert transcriptomic coordinates to genomic ones
    Return list of positions where the read maps to'''
    annotation_df = annotation_df[(annotation_df['feature'] != 'gene') & (
        annotation_df['feature'] != 'transcript')]
    new_df = annotation_df.copy(deep=True)
    new_df['positions'] = new_df.apply(coordinate_list, axis=1)
    new_df.drop(['feature', 'start', 'end'], axis=1, inplace=True)
    new_df = group_by_tr_id(new_df)
    return new_df


def get_transcriptomic_reads(tr_bam, annotation_info, paired, sense):
    """ get all the transcriptomic reads that map
    to specific chromosome (the annotation_info is a
    subtable of annotation containing one chromosome)"""

    output = StringIO()
    csv_writer = writer(output)
    annotation = convert_tr_to_gn_coordinates(annotation_info)
    annotation_index = set(annotation.index.values)

    bam = pysam.AlignmentFile(tr_bam, "rb")
    paired = int(paired)
    sense = int(sense)

    if (paired == 1) & (sense == 1):
        for read in bam.fetch():
            name = read.reference_name.split('.')[0]
            if name in annotation_index:
                if not read.is_reverse:
                    gn_strand = annotation.at[name, 'strand']
                    gn_coord = sorted(annotation.at[name, 'positions'])
                    tr_pos = list(read.get_reference_positions())
                    tr_pos = tr_pos[math.floor(len(tr_pos) / 2)]
                    gn_pos = list(read.get_reference_positions())
                    gn_pos = gn_pos[math.floor(len(gn_pos) / 2)]
                    if gn_strand == '-':
                        gn_pos = gn_coord[-gn_pos]
                    else:
                        gn_pos = gn_coord[gn_pos]
                    csv_writer.writerow(
                        [read.query_name, read.get_tag('AS'), gn_strand,
                            gn_pos,
                            read.reference_name, tr_pos, True])

    elif (paired == 1) & (sense == 2):
        for read in bam.fetch():
            name = read.reference_name.split('.')[0]
            if name in annotation_index:
                if read.is_reverse:
                    gn_strand = annotation.at[name, 'strand']
                    gn_coord = sorted(annotation.at[name, 'positions'])
                    tr_pos = list(read.get_reference_positions())
                    tr_pos = tr_pos[math.floor(len(tr_pos) / 2)]
                    gn_pos = list(read.get_reference_positions())
                    gn_pos = gn_pos[math.floor(len(gn_pos) / 2)]
                    if gn_strand == '-':
                        gn_pos = gn_coord[-gn_pos]
                    else:
                        gn_pos = gn_coord[gn_pos]
                    csv_writer.writerow(
                        [read.query_name, read.get_tag('AS'), gn_strand,
                            gn_pos,
                            read.reference_name, tr_pos, True])

    elif (paired == 2) & (sense == 1):
        for read in bam.fetch():
            name = read.reference_name.split('.')[0]
            if name in annotation_index:
                if ((read.is_read1) & (read.is_reverse)):
                    continue
                elif (read.is_read2) & (not read.is_reverse):
                    continue
                elif not read.is_proper_pair:
                    continue
                else:
                    gn_strand = annotation.at[name, 'strand']
                    gn_coord = sorted(annotation.at[name, 'positions'])
                    tr_pos = list(read.get_reference_positions())
                    tr_pos = tr_pos[math.floor(len(tr_pos) / 2)]
                    gn_pos = list(read.get_reference_positions())
                    gn_pos = gn_pos[math.floor(len(gn_pos) / 2)]
                    if gn_strand == '-':
                        gn_pos = gn_coord[-gn_pos]
                    else:
                        gn_pos = gn_coord[gn_pos]
                    csv_writer.writerow(
                        [read.query_name, read.get_tag('AS'), gn_strand,
                            gn_pos,
                            read.reference_name, tr_pos, read.is_read1])

    elif (paired == 2) & (sense == 2):
        for read in bam.fetch():
            name = read.reference_name.split('.')[0]
            if name in annotation_index:
                if (read.is_read1) & (not read.is_reverse):
                    continue
                elif (read.is_read2) & (read.is_reverse):
                    continue
                elif not read.is_proper_pair:
                    continue
                else:
                    gn_strand = annotation.at[name, 'strand']
                    gn_coord = sorted(annotation.at[name, 'positions'])
                    tr_pos = list(read.get_reference_positions())
                    tr_pos = tr_pos[math.floor(len(tr_pos) / 2)]
                    gn_pos = list(read.get_reference_positions())
                    gn_pos = gn_pos[math.floor(len(gn_pos) / 2)]
                    if gn_strand == '-':
                        gn_pos = gn_coord[-gn_pos]
                    else:
                        gn_pos = gn_coord[gn_pos]
                    csv_writer.writerow(
                        [read.query_name, read.get_tag('AS'), gn_strand,
                            gn_pos,
                            read.reference_name, tr_pos, read.is_read1])
                    # if gn_strand == '-':
                    #     gn_coord.sort(reverse=True)
                    # gn_pos = [gn_coord[i] for i in read.get_reference_positions()]
                    # csv_writer.writerow(
                    #     [read.query_name, read.get_tag('AS'), gn_strand,
                    #         gn_pos[math.floor(len(gn_pos) / 2)],
                    #         read.reference_name, tr_pos[math.floor(len(tr_pos) / 2)]])
    bam.close()
    final = pd.DataFrame()
    output.seek(0)
    first_char = output.read(1)
    if first_char:
        output.seek(0)
        final = pd.read_csv(output)
        final.columns = [
            'read_name', 'as_score', 'gn_strand', 'gn_middle', 'tr_name', 'tr_middle', 'read1']
        final.set_index('read_name', inplace=True, drop=True)
    return final


def remove_genomic_reads(arguments):
    """
        Keep a table with those reads that map with a higher AS score
        to genome. Create new genomic bamfile where reads with highest
        AS mapping to transcriptome are excluded.
    """
    gn_bam = arguments['gn_bam']
    tr_bam = arguments['tr_bam']
    annotation = arguments['annotation']
    chromosome = str(arguments['chromosome'])
    gn_out = arguments['gn_out']
    paired = int(arguments['paired'])
    sense = int(arguments['sense'])

    sys.stdout.write(
        chromosome + 'start_remove_genomic_reads\n')
    sys.stdout.flush()

    transcriptomic_reads = get_transcriptomic_reads(
        tr_bam, annotation, paired, sense)
    output = StringIO()
    csv_writer = writer(output)
    chromosome = str(chromosome)

    if transcriptomic_reads.empty:
        tr_reads = set()
    else:
        tr_reads = set(transcriptomic_reads.index.values)

    bam = pysam.AlignmentFile(gn_bam, "rb")
    filtered_genomic = pysam.AlignmentFile(gn_out, "wb", template=bam)
    try:
        bam.fetch(chromosome, 1, 2)
    except:
        bam.close()
        final = pd.DataFrame()
        filtered_genomic.close()
        return final

    for read in bam.fetch(chromosome):
        if paired == 2:
            if not read.is_proper_pair:
                continue
        if read.query_name in tr_reads:
            overlap = 0
            df_part = transcriptomic_reads.loc[[read.query_name], :]
            for row in df_part.itertuples():
                if (row.gn_middle in set(
                        read.get_reference_positions())) & (row.read1 == read.is_read1):
                    overlap += 1
                    if read.get_tag('AS') - row.as_score > 3:
                        csv_writer.writerow(pd.Series(row))
                        filtered_genomic.write(read)
                    break
            # no overlap write in gene
            if overlap == 0:
                filtered_genomic.write(read)
        else:
            filtered_genomic.write(read)

    bam.close()
    filtered_genomic.close()

    final = pd.DataFrame()
    output.seek(0)
    first_char = output.read(1)
    if first_char:
        output.seek(0)
        final = pd.read_csv(output)
        final.columns = [
            'read_name', 'as_score', 'gn_strand', 'gn_middle', 'tr_name', 'tr_middle', 'read1']
        final.set_index('read_name', inplace=True, drop=True)
    sys.stdout.write(chromosome + 'finish_remove_genomic_reads \n')
    sys.stdout.flush()

    return final


def remove_transcriptomic_reads(arguments):
    """
        Remove reads from transcriptomic bam that map to genome with
        higher AS score
    """
    tr_bam = arguments['tr_bam']
    gn_reads = arguments['genomic_reads']
    tr_out = arguments['tr_out']
    paired = int(arguments['paired'])
    sense = int(arguments['sense'])
    if gn_reads.empty:
        gn_reads_set = set()
    else:
        gn_reads_set = set(gn_reads.index.values)
    bam = pysam.AlignmentFile(tr_bam, "rb")
    filtered_transcriptomic = pysam.AlignmentFile(tr_out, "wb", template=bam)

    if (paired == 1) & (sense == 1):
        for read in bam.fetch():
            if read.is_reverse:
                continue
            elif read.query_name in gn_reads_set:
                success = 0
                df_part = gn_reads.loc[[read.query_name], :]
                df_part = df_part[(df_part['tr_name'] == read.reference_name)]
                for row in df_part.itertuples():
                    if row.tr_middle in set(read.get_reference_positions()):
                        success += 1
                if success == 0:
                    filtered_transcriptomic.write(read)
            else:
                filtered_transcriptomic.write(read)
    elif (paired == 1) & (sense == 2):
        for read in bam.fetch():
            if not read.is_reverse:
                continue
            elif read.query_name in gn_reads_set:
                success = 0
                df_part = gn_reads.loc[[read.query_name], :]
                df_part = df_part[(df_part['tr_name'] == read.reference_name)]
                for row in df_part.itertuples():
                    if row.tr_middle in set(read.get_reference_positions()):
                        success += 1
                if success == 0:
                    filtered_transcriptomic.write(read)
            else:
                filtered_transcriptomic.write(read)
    elif (paired == 2) & (sense == 1):
        for read in bam.fetch():
            if not read.is_proper_pair:
                continue
            elif (read.is_read1) & (read.is_reverse):
                continue
            elif (read.is_read2) & (not read.is_reverse):
                continue
            elif read.query_name in gn_reads_set:
                success = 0
                df_part = gn_reads.loc[[read.query_name], :]
                df_part = df_part[(df_part['tr_name'] == read.reference_name)]
                for row in df_part.itertuples():
                    if row.tr_middle in set(read.get_reference_positions()):
                        success += 1
                if success == 0:
                    filtered_transcriptomic.write(read)
            else:
                filtered_transcriptomic.write(read)
    elif (paired == 2) & (sense == 2):
        for read in bam.fetch():
            if not read.is_proper_pair:
                continue
            elif (read.is_read1) & (not read.is_reverse):
                continue
            elif (read.is_read2) & (read.is_reverse):
                continue
            elif read.query_name in gn_reads_set:
                success = 0
                df_part = gn_reads.loc[[read.query_name], :]
                df_part = df_part[(df_part['tr_name'] == read.reference_name)]
                for row in df_part.itertuples():
                    if row.tr_middle in set(read.get_reference_positions()):
                        success += 1
                if success == 0:
                    filtered_transcriptomic.write(read)
            else:
                filtered_transcriptomic.write(read)
    bam.close()
    filtered_transcriptomic.close()
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

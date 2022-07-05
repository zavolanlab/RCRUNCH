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
#   This script utilises the significant peaks detected in the previous step
#   and, based on the mi, sigma, produces a fasta file containing the sequences
#   of the peaks.
#
#   The first 1000 peaks (if there are less than 1000, then all) peaks are
#   converted into sequences for motif prediction.
# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import random
from pyfasta import Fasta
import pybedtools
from Bio import SeqIO

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():
    """ Main function """

    __doc__ = "Create motevo input files of ALL peaks"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--peaks_file",
        dest="peaks_file",
        help="File with the significant RBP binding peaks.",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--all_regions",
        dest="all_regions",
        help="File with the all regions. Used to sample background",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--genome_fasta",
        dest="genome_fasta",
        help="Fasta file of the genome to get \
        sequences of the genomic peaks.",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--transcriptome_fasta",
        dest="transcriptome_fasta",
        help="Fasta file of the transcriptome to get \
            sequences of the transcriptomic peaks.",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--chromosomes_t",
        dest="chromosomes_t",
        help="File containing chromosome names and \
            lengths provided during STAR indexing",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--chromosomes_g",
        dest="chromosomes_g",
        help="File containing chromosome names and \
            lengths provided during STAR indexing",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--peak_size",
        dest="peak_size",
        help="Size of the peak",
        required=False)

    parser.add_argument(
        "--tempdir",
        dest="tempdir",
        help="temporary directory for \
        pybed calculations",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--out_folder",
        dest="out_folder",
        help="output folder",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--crosslink_type",
        dest="crosslink_type",
        help="specify the whole peak, or have a crosslink site",
        choices=['peak_center', 'crosslink'],
        default='crosslink',
        required=False)

    parser.add_argument(
        "--genome_tag",
        dest="genome_tag",
        help="genome name (eg:hg38)",
        required=True)

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

    pybedtools.helpers.set_tempdir(options.tempdir)
    # ----------------------------------------------------------
    # obtain the peaks
    # if there are more than 1000
    # keep only the first 1000 sorted by zscore
    peaks = pd.read_csv(
        options.peaks_file, header=0,
        sep='\t', index_col=None,
        comment='#', engine='c',
        dtype={"chromosome": str, "start": int, "end": int,
               "strand": str, "mi": float, "ro": float,
               "sigma": float, "crosslink": float, "crosslink_weight": float,
               "z_score": float, "nj": float, "mj": float})
    peaks.sort_values(by=['z_score'], ascending=False, inplace=True)

    # if there are more than 1000 peaks, keep the first 1000

    if len(peaks) > 1000:
        peaks = peaks.iloc[:1000]
    peaks.sort_values(by=['z_score'], axis=0, ascending=False,
                      inplace=True, kind='quicksort')
    peaks.reset_index(drop=True, inplace=True)

    all_regions = pd.read_csv(
        options.all_regions, header=0, sep='\t',
        index_col=0, comment='#', engine='c',
        dtype={"chromosome": str, "start": int, "end": int,
               "name": str, "type": str, "Reads_f": float,
               "Reads_b": float, "strand": str, "startbg": int, "endbg": int,
               "z_score": float, "P_nimi": float, "p_bgi": float,
               "cumsum_p_bgi": float, "T": int, "FDR": float})
    all_regions.sort_values(by=['z_score'], ascending=True, inplace=True)
    all_regions = all_regions[all_regions['z_score'] <= 0]
    all_regions.index = np.arange(len(all_regions))

    # ----------------------------------------------------------
    # based on the mi and sigma values get the peak coordinates
    new_peaks = peaks[['chromosome', 'strand', 'z_score', 'sigma']].copy(deep=True)
    if options.crosslink_type == 'peak_center':
        if options.peak_size:
            peak_size = int(int(options.peak_size) / 2)
            new_peaks['start'] = peaks['start'] + \
                peaks['mi'] - peak_size
            new_peaks['end'] = peaks['start'] + \
                peaks['mi'] + peak_size
        else:
            new_peaks['start'] = peaks['start'] + \
                peaks['mi'] - peaks['sigma']
            new_peaks['end'] = peaks['start'] + \
                peaks['mi'] + peaks['sigma']

    elif options.crosslink_type == 'crosslink':
        if options.peak_size:
            peak_size = int(int(options.peak_size) / 2)
            new_peaks['start'] = peaks['start'] + \
                peaks['crosslink'] - peak_size
            new_peaks['end'] = peaks['start'] + \
                peaks['crosslink'] + peak_size
        else:
            new_peaks['start'] = peaks['start'] + \
                peaks['crosslink'] - peaks['sigma']
            new_peaks['end'] = peaks['start'] + \
                peaks['crosslink'] + peaks['sigma']

    if options.chromosomes_t:
        transcriptome = pd.read_csv(
            options.chromosomes_t,
            header=None,
            sep='\t',
            index_col=None,
            comment='#',
            engine='python')

        transcriptome.columns = ['chromosome', 'chromosome_length']
        transcriptome_names = list(transcriptome['chromosome'].values)
        transcriptome.set_index('chromosome', inplace=True, drop=True)

        ft = Fasta(
            options.transcriptome_fasta,
            key_fn=lambda key: key.split()[0])
    else:
        transcriptome_names = []

    genome = pd.read_csv(
        options.chromosomes_g,
        header=None,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python',
        dtype={0: str, 1: int})

    genome.columns = ['chromosome', 'chromosome_length']
    genome.set_index('chromosome', inplace=True, drop=True)

    fg = Fasta(
        options.genome_fasta,
        key_fn=lambda key: key.split()[0])

    for index, row in new_peaks.iterrows():
        new_peaks.loc[index, 'start'] = get_start(row['start'])
        info = fg
        if row['chromosome'] in transcriptome_names:
            new_peaks.loc[index, 'end'] = get_end(
                row['end'], row['chromosome'], transcriptome)
        else:
            new_peaks.loc[index, 'end'] = get_end(
                row['end'], row['chromosome'], genome)

    new_peaks['end'] = new_peaks['end'].astype(int)
    new_peaks['start'] = new_peaks['start'].astype(int)
    new_peaks.index = np.arange(len(new_peaks))
    new_peaks = new_peaks[['chromosome', 'start', 'end',
                           'sigma', 'z_score', 'strand']]
    new_peaks.sort_values(by=['chromosome', 'start', 'end'],
                          axis=0, ascending=[True, True, True],
                          inplace=True, kind='quicksort')
    new_peaks.to_csv(
        os.path.join(options.tempdir, 'rcrunch.bed'),
        header=False, index=False, sep='\t')
    new_peaks = merge_overlapping_peaks(
        os.path.join(options.tempdir, 'rcrunch.bed'),
        'rcrunch')
    new_peaks.reset_index(drop=True, inplace=True)
    pybedtools.cleanup()

    # ----------------------------------------------------------
    # randomise rows in the table
    new_peaks = new_peaks.loc[
        random.sample(list(new_peaks.index), len(new_peaks))]
    new_peaks.index = np.arange(len(new_peaks))

    # ----------------------------------------------------------
    #  Subsets construction
    test_set = open(os.path.join(
        options.out_folder,
        'test.fasta'), 'w')

    test_set_bg = open(os.path.join(
        options.out_folder,
        'test_bg.fasta'), 'w')

    # ----------------------------------------------------------
    # convert peak coordinates to fasta squences
    # (genomic & transcriptomic)

    # ----------------------------------------------------------
    # use pyfasta to find the sequences that correspond to the
    # coordinates of the peaks
    # 1. obtain the genomic sequences
    peak_counter = 0
    for index, row in new_peaks.iterrows():
        chromosome = str(row['chromosome'])
        start = get_start(row['start'])
        strand = row['strand']
        score = row['z_score']
        info = fg
        if chromosome in transcriptome_names:
            end = get_end(row['end'], chromosome, transcriptome)
            strand = '+'
            info = ft
        else:
            end = get_end(row['end'], chromosome, genome)
        size = end - start
        sequence = info.sequence({
            'chr': chromosome,
            'start': start,
            'stop': end,
            'strand': strand}, one_based=False)
        line = f'>{chromosome}_{start}_{end}_{strand}_{score}\n{sequence}\n'
        test_set.write(line)
        for j in np.arange(1, 11):
            for x in np.arange(1, 11):
                bg_row = all_regions.sample(n=1).squeeze()
                bg_chromosome = str(bg_row['chromosome'])
                bg_start = get_start(bg_row['start'])
                if bg_row['end'] - bg_start > 10 + size:
                    break
            try:
                bg_start = random.choice(np.arange(
                    bg_start, bg_row['end'] - size - 1))
            except:
                sys.stderr.write(f'***\n{pd.DataFrame(row).T}\n{pd.DataFrame(peaks.iloc[index]).T}\n{size}\n{pd.DataFrame(bg_row).T}')
            bg_end = bg_start + size
            bg_strand = bg_row['strand']
            bg_info = fg
            bg_score = bg_row['z_score']
            if bg_chromosome in transcriptome_names:
                bg_end = get_end(
                    bg_end, bg_chromosome, transcriptome)
                bg_strand = '+'
                bg_info = ft
            else:
                bg_end = get_end(bg_end, bg_chromosome, genome)

            bg_sequence = bg_info.sequence({
                'chr': bg_chromosome,
                'start': bg_start,
                'stop': bg_end,
                'strand': bg_strand}, one_based=False)
            bg_line = f'>{bg_chromosome}_{bg_start}_{bg_end}_{bg_strand}_{bg_score}_random_sequence{str(peak_counter)}\n{bg_sequence}\n'
            test_set_bg.write(bg_line)
            peak_counter += 1

    test_set.close()
    test_set_bg.close()


def merge_overlapping_peaks(peak_df, tool=''):
    '''merge overlapping peaks, keep the highest
    score observed and the outer boundaries'''
    overlap = pybedtools.BedTool(peak_df)
    if tool == 'rcrunch':
        overlap = overlap.merge(
            s=True, c="4,5", o="first,max").to_dataframe(
            names=['chromosome', 'start', 'end', 'strand',
                   'sigma', 'z_score'])
    return overlap


def get_start(start):
    if start <= 0:
        start = 0
    return int(start)


def get_end(end, chromosome, info):
    if end >= int(
            info.loc[chromosome, 'chromosome_length']):
        end = int(
            info.loc[chromosome, 'chromosome_length'] - 1)
    return int(end)


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interruptions
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)

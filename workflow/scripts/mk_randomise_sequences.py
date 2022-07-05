
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Author : Katsantoni Maria
# Company: Mihaela Zavolan group, Biozentrum, Basel
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
#   This script is part of the RCRUNCH pipeline. RCRUNCH pipeline detects
#   significantly enriched (in reads) peaks, which correspond to binding sites
#   of an RBP (CLIP experiment analysis).
#   This script creates the files as described in documentation called:
#   1. (training/test)set: contains a random half of the total peak sequences
#   2. (training/test)_bg:  contains 10 randomisations of the nucleotides
#      for each of the sequences in the test_set
#   3. (training/test)_pool: contains the test_set + the test_bg sequences
#   Background sequences are created by randomisation of the nucleotides of
#   each peak sequence. (How does that help in case of low comlexity regions?)

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
from ushuffle import shuffle, Shuffler
from Bio import SeqIO


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():
    """ Main function """

    __doc__ = "Divide the peaks into \
    training and test set and create Phylogibbs input files"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--phylogenetic_peaks",
        dest="phylogenetic_peaks",
        help="multiple alignments of the peak sequences",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--phylogenetic_out_folder",
        dest="phylogenetic_out_folder",
        help="output folder for the phylogenetic random sequences",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--peaks",
        dest="peaks",
        help="peak sequences",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--peaks_bg",
        dest="peaks_bg",
        help="peak sequences for the bg",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--out_folder",
        dest="out_folder",
        help="output folder",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--prefix",
        dest="prefix",
        help="prefix name of the files",
        required=False)

    parser.add_argument(
        "--species",
        dest="species",
        help="species name for the motevo alignments",
        required=False)

    parser.add_argument(
        "--random_seq_per_peak_num",
        dest="random_seq_per_peak_num",
        help="number of random sequences to be created for \
        each peak sequence",
        required=False,
        default=10)
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

   

    # ----------------------------------------------------------
    if options.peaks:
        peaks_set = open(os.path.join(
            options.out_folder,
            options.prefix + '.fasta'), 'w')
        peaks_bg = open(os.path.join(
            options.out_folder,
            options.prefix + '_bg.fasta'), 'w')
        peaks_pool = open(os.path.join(
            options.out_folder,
            options.prefix + '_pool.fasta'), 'w')

        motevo_peaks_set = open(os.path.join(
            options.out_folder,
            options.prefix + '_motevo.fasta'), 'w')
        motevo_peaks_bg = open(os.path.join(
            options.out_folder,
            options.prefix + '_bg_motevo.fasta'), 'w')
        motevo_peaks_pool = open(os.path.join(
            options.out_folder,
            options.prefix + '_pool_motevo.fasta'), 'w')

        table_peak = read_peaks(options.peaks, options.species)

        fasta = get_phylogibbs_fasta(table_peak)
        peaks_set.write(fasta)
        peaks_pool.write(fasta)
        fasta = get_motevo_fasta(table_peak, '>>')
        motevo_peaks_set.write(fasta)
        motevo_peaks_pool.write(fasta)

        if options.peaks_bg:
            table_peak = read_peaks(options.peaks_bg, options.species)
            new_peaks = table_peak
            fasta = get_phylogibbs_fasta(
                new_peaks,
                '>')
            peaks_bg.write(fasta)
            peaks_pool.write(fasta)
            fasta = get_motevo_fasta(
                new_peaks,
                '>>')
            motevo_peaks_bg.write(fasta)
            motevo_peaks_pool.write(fasta)

        peaks_set.close()
        peaks_bg.close()
        peaks_pool.close()
        motevo_peaks_set.close()
        motevo_peaks_bg.close()
        motevo_peaks_pool.close()
        sys.stdout.write('Done!')

    # obtain fasta sequences from the phylogenetic multiple alignments
    elif options.phylogenetic_peaks:
        randomisation_number = int(options.random_seq_per_peak_num)
        phylogenetic_peaks_set = open(os.path.join(
            options.phylogenetic_out_folder,
            options.prefix + '.fasta'), 'w')
        phylogenetic_peaks_bg = open(os.path.join(
            options.phylogenetic_out_folder,
            options.prefix + '_bg.fasta'), 'w')
        phylogenetic_peaks_pool = open(os.path.join(
            options.phylogenetic_out_folder,
            options.prefix + '_pool.fasta'), 'w')

        motevo_phyl_peaks_set = open(os.path.join(
            options.phylogenetic_out_folder,
            options.prefix + '_motevo.fasta'), 'w')
        motevo_phyl_peaks_bg = open(os.path.join(
            options.phylogenetic_out_folder,
            options.prefix + '_bg_motevo.fasta'), 'w')
        motevo_phyl_peaks_pool = open(os.path.join(
            options.phylogenetic_out_folder,
            options.prefix + '_pool_motevo.fasta'), 'w')

        with open(options.phylogenetic_peaks) as peaks:
            # for each of the peaks create a table with one column being the
            # title (starting with '>') and one the sequence
            for peak in myreadlines(peaks, '\n\\\\\n'):
                if not peak:
                    continue
                table_peak = pd.DataFrame(columns=['title', 'sequence'])
                counter = 0
                peak_info = peak.split('\n')

                for line in peak_info:
                    line = line.strip('\n')
                    if (line.startswith('>>')) or (line.startswith('>')):
                        counter += 1
                        table_peak.at[counter, 'title'] = line
                        continue
                    elif line:
                        table_peak.at[counter, 'sequence'] = line
                    else:
                        continue

                fasta_ma = get_peak_fasta(table_peak)
                phylogenetic_peaks_set.write(fasta_ma + '\n')
                phylogenetic_peaks_pool.write(fasta_ma + '\n')
                motevo_phyl_peaks_set.write(fasta_ma + '\n')
                motevo_phyl_peaks_pool.write(fasta_ma + '\n')

                nucleotide_table = table_peak['sequence'].apply(
                    lambda x: pd.Series(list(x)))
                for i in np.arange(1, randomisation_number + 1):
                    table1 = randomise_sequences(nucleotide_table)
                    table1['title'] = table_peak['title']
                    table1 = table1[['title', 'sequence']]
                    fasta_ma = get_peak_fasta(
                        table1, ('_random_sequence' + str(i)))
                    phylogenetic_peaks_bg.write(fasta_ma + '\n')
                    phylogenetic_peaks_pool.write(fasta_ma + '\n')
                    motevo_phyl_peaks_bg.write(fasta_ma + '\n')
                    motevo_phyl_peaks_pool.write(fasta_ma + '\n')

        phylogenetic_peaks_set.close()
        phylogenetic_peaks_bg.close()
        phylogenetic_peaks_pool.close()
        motevo_phyl_peaks_set.close()
        motevo_phyl_peaks_bg.close()
        motevo_phyl_peaks_pool.close()
        sys.stdout.write('Done phylogenetic!')
    return


def read_peaks(infile, species):
    with open(infile, 'r') as peaks:
        peak_number = 0
        table_peak = pd.DataFrame(columns=['title', 'sequence'])
        for peak in SeqIO.parse(peaks, 'fasta'):
            title = str(peak.id)
            title = title.replace('>', '')
            title = species + '_' + title
            seq = str(peak.seq)
            peak_number += 1
            table_peak.at[peak_number, 'title'] = title
            table_peak.at[peak_number, 'sequence'] = seq
    return table_peak


def get_motevo_fasta(sequence_table, prefix='>', random_index=''):
    fasta = ''
    for index, row in sequence_table.iterrows():
        title = prefix + row['title'] + random_index
        sequence = row['sequence'].replace('n', 'N')
        fasta += '\n'.join([title, sequence]) + '\n\n'
    return fasta


def get_phylogibbs_fasta(sequence_table, prefix='>', random_index=''):
    '''Modify format to be compatible for PHYLOGIBBS use'''
    fasta = ''
    for index, row in sequence_table.iterrows():
        title = prefix + row['title'] + random_index
        sequence = row['sequence'].replace('N', 'X')
        sequence = sequence.replace('n', 'X')
        fasta += '\n'.join([title, sequence]) + '\n\n'
    return fasta


def randomise_shuffle_kmers(sequence_table, kmer=2):
    for index, row in sequence_table.iterrows():
        sequence = row['sequence'].encode('utf-8')
        shuffler = Shuffler(sequence, 2)
        sequence_new = shuffler.shuffle()
        sequence_table.loc[index, 'sequence'] = sequence_new.decode('utf-8')
    return sequence_table


# def randomise_all_seqs(sequences, length):
#     new_peak_seq = []
#     new_seq = random.sample(sequences, len(sequences))
#     x = 0
#     for length in lengths:
#         subseq = ''.join(new_seq[x: x + length])
#         subseq = subseq.replace('N', 'X')
#         new_peak_seq.append(subseq)
#         x += length
#     return new_peak_seq


def get_peak_fasta(sequence_table, random_index=''):
    fasta = ''
    for index, row in sequence_table.iterrows():
        title = row['title'] + random_index
        sequence = row['sequence'].replace('N', 'X')
        fasta += '\n'.join([title, sequence]) + '\n'
    return fasta


def randomise_sequences(sequence_table, nucl=1):
    sequence_table = sequence_table.T
    gaps = ''
    gaps_counter = 0
    dinucl = ''
    dinucl_counter = 0
    seq = []
    for index, row in sequence_table.iterrows():
        if row[1] == '-':
            gaps += str(int(index)) + '_'
            gaps_counter += 1
            if dinucl_counter > 0:
                dinucl = dinucl.strip('_')
                seq.append(dinucl)
                dinucl = ''
                dinucl_counter = 0
        else:
            if gaps_counter > 0:
                gaps = gaps.strip('_')
                seq.append(gaps)
                gaps = ''
                gaps_counter = 0
            else:
                dinucl += str(int(index)) + '_'
                dinucl_counter += 1

                if dinucl_counter == nucl:
                    dinucl = dinucl.strip('_')
                    seq.append(dinucl)
                    dinucl = ''
                    dinucl_counter = 0
    if (dinucl_counter > 0) & (dinucl_counter < nucl):
        dinucl = dinucl.strip('_')
        seq.append(dinucl)
    if gaps_counter > 0:
        gaps = gaps.strip('_')
        seq.append(gaps)

    random_sequence = np.random.choice(
        seq,
        len(seq),
        replace=False)

    flattened = []
    for i in random_sequence:
        values = i.split('_')
        for y in values:
            flattened.append(int(y))
    print(flattened)
    sequence_table = sequence_table.loc[flattened, :]
    sequence_table = sequence_table.T
    sequence_table_new = pd.DataFrame()
    for index, row in sequence_table.iterrows():
        sequence_table_new.at[index, 'sequence'] = row.str.cat(sep='')

    return sequence_table_new


def myreadlines(file, delimiter='\n', bufsize=4096):
    buf = ''
    while True:
        newbuf = file.read(bufsize)
        if not newbuf:
            yield buf
            return
        buf += newbuf
        lines = buf.split(delimiter)
        for line in lines[:-1]:
            yield line
        buf = lines[-1]


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

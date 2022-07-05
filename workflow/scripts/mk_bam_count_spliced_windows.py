# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import pysam
import multiprocessing as mp
import os
import logging
from io import StringIO
from csv import writer
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
        "--bam_foreground",
        dest="bam_foreground",
        help="Foreground reads aligned to genome",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--bam_foreground_frequencies",
        dest="bam_foreground_frequencies",
        help="Read frequencies in foreground",
        required=False,
        default=0)

    parser.add_argument(
        "--sense_f",
        dest="sense_f",
        help="Foreground sense mate",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--bam_background",
        dest="bam_background",
        help="Background reads aligned to genome",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--bam_background_frequencies",
        dest="bam_background_frequencies",
        help="Read frequencies in background",
        required=False,
        default=0)

    parser.add_argument(
        "--sense_b",
        dest="sense_b",
        help="Background sense mate",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--paired_f",
        dest="paired_f",
        help="Foreground pe or se",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--paired_b",
        dest="paired_b",
        help="Background pe or se",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--chromosomes",
        dest="chromosomes",
        help="File containing chromosome names\
        and lengths provided during STAR indexing, required",
        required=True,
        metavar="FILE")

    # ----------------------------------------------------------
    parser.add_argument(
        "--out",
        dest="out",
        help="Output folder path",
        required=True)

    parser.add_argument(
        "--prefix",
        dest="prefix",
        help="prefix",
        required=True)

    parser.add_argument(
        "--window_f",
        dest="window_f",
        help="Length of the sliding windows",
        required=False,
        default="300")

    parser.add_argument(
        "--window_b",
        dest="window_b",
        help="Length of the sliding windows",
        required=False,
        default="300")

    parser.add_argument(
        "--data_type",
        dest="data_type",
        help="genome or transcriptome",
        required=True)

    parser.add_argument(
        "--step_size",
        dest="step_size",
        default=150,
        help="Step size for the generation of sliding windows.\
        If different from '--window-size', overlapping or \
        gapped sliding windows will be generated.")

    parser.add_argument(
        "--cutoff",
        dest="cutoff",
        help="Cutoff of coverage for window \
        to be included",
        required=False,
        default="1")

    parser.add_argument(
        "--threads",
        dest="threads",
        help="Number of threads for multiprocessing",
        required=False,
        default="1")

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

    # -------------------------------------------------------------------------
    #           Cast options
    # -------------------------------------------------------------------------
    bam_foreground = options.bam_foreground
    bam_background = options.bam_background
    chromosomes = options.chromosomes
    out = options.out
    prefix = options.prefix
    window_f = int(float(options.window_f))
    step_size = int(float(options.step_size))
    window_b = int(float(options.window_b))
    sense_f = int(float(options.sense_f))
    sense_b = int(float(options.sense_b))
    paired_f = int(float(options.paired_f))
    paired_b = int(float(options.paired_b))
    data_type = options.data_type
    threads = int(float(options.threads))
    cutoff = int(float(options.cutoff))

    global fg_frequencies
    if options.bam_foreground_frequencies:
        fg_frequencies = pd.read_csv(
            options.bam_foreground_frequencies,
            header=0,
            sep='\t',
            index_col=0,
            comment='#',
            engine='python')
        fg_frequencies = fg_frequencies['Frequency'].to_dict()
    else:
        fg_frequencies = {}

    global bg_frequencies
    if options.bam_background_frequencies:
        bg_frequencies = pd.read_csv(
            options.bam_background_frequencies,
            header=0,
            sep='\t',
            index_col=0,
            comment='#',
            engine='python')
        bg_frequencies = bg_frequencies['Frequency'].to_dict()
    else:
        bg_frequencies = {}

    processes = 4 * threads
    number_of_windows = 20000

    # obtain the chromosome lengths and names
    # used for creating windows
    chromosome_info = pd.read_csv(
        chromosomes,
        header=None,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')

    chromosome_info.columns = ['chromosome', 'chromosome_length']
    # get a table of all the possible windows
    # format should be :
    # chr   start  end chr:start-end:strand  coverage  strand
    # 1       0       300     1:0-300:+       0       +
    # 1       150     450     1:150-450:+     0       +

    # -------------------------------------------------------------------------
    #           GENOMIC APPROACH
    # -------------------------------------------------------------------------

    # -----------------------------------------------------------
    # Foreground and background windows calculated at the same time. Decisions
    # on spliced windows based both on foreground and background.

    # format should be :
    # chr  start  end  chr:st1-end1;stn-endn:strand coverage_fg coverage_bg strand type
    # 1     0     300     1:0-150;200-350:+              0            0      +     spliced
    # 1     150   450     1:150-450:+                    0            0      +     genomic



    sliding_windows_parameters_fg = []
    sliding_windows_parameters_bg = []



    for row in chromosome_info.itertuples():
        sliding_windows_parameters_fg.append([
            row.chromosome,
            row.chromosome_length,
            window_f,
            step_size])

        sliding_windows_parameters_bg.append([
            row.chromosome,
            row.chromosome_length,
            window_b,
            step_size])
    try:
        pool = mp.Pool(processes=processes)
        sliding_windows_fg = pool.map(
            sliding_windows,
            sliding_windows_parameters_fg)
    finally:
        pool.close()
        pool.join()

    try:
        pool = mp.Pool(processes=processes)
        sliding_windows_bg = pool.map(
            sliding_windows,
            sliding_windows_parameters_bg)
    finally:
        pool.close()
        pool.join()

    parameters = []
    count = 0

    windows_fg = split_tables(sliding_windows_fg, number_of_windows)
    sliding_windows_fg = []
    for each_window in windows_fg:
        parameters.append([
            bam_foreground,
            sense_f,
            out,
            prefix,
            each_window,
            'fg',
            paired_f,
            count,
            cutoff])
        count += 1

    windows_fg = []

    windows_bg = split_tables(sliding_windows_bg, number_of_windows)
    sliding_windows_bg = []
    for each_window in windows_bg:
        parameters.append([
            bam_background,
            sense_b,
            out,
            prefix,
            each_window,
            'bg',
            paired_b,
            count,
            cutoff])
        count += 1
    windows_bg = []

    try:
        pool = mp.Pool(processes=processes)
        pool.map(coverage_profile_gn, parameters)

    finally:
        pool.close()
        pool.join()



def sliding_windows(params):
    chromosome_name = params[0]
    chromosome_length = params[1]
    window_length = int(params[2] / 2)
    step_size = params[3]

    if (chromosome_length <= step_size):
        window = pd.DataFrame()
        window[0] = [chromosome_name]
        window[1] = [0]
        window[2] = [chromosome_length]
        return window

    else:
        windows_center = np.arange(step_size, chromosome_length, step_size)
        window = pd.DataFrame()
        window['start'] = windows_center - window_length
        window['end'] = windows_center + window_length
        window['name'] = chromosome_name
        window = window[['name', 'start', 'end']]
        window.loc[window.start < 0, 'start'] = 0
        window.loc[window.end > chromosome_length, 'end'] = \
            chromosome_length
        window.columns = [0, 1, 2]
    return window


def split_tables(dfs, length):
    dfs_new = []
    for window in dfs:
        if len(window) > length:
            number_of_slices = np.floor(len(window) / length)
            windows = np.array_split(window, number_of_slices)
            dfs_new.extend(windows)
        else:
            dfs_new.append(window)
    return dfs_new


def estimate_frequency(name, fragments, cutoff):
    unique_fragments = list(fragments)
    multimapper_number = 0
    total_number = len(unique_fragments)
    unique_fragments_number = len(unique_fragments)
    if unique_fragments_number > cutoff:
        if name == 'fg':
            if fg_frequencies:
                for i in unique_fragments:
                    if i in fg_frequencies.keys():
                        unique_fragments_number = unique_fragments_number \
                            - 1 + fg_frequencies[i]
                        multimapper_number += 1
        if name == 'bg':
            if bg_frequencies:
                for i in unique_fragments:
                    if i in bg_frequencies.keys():
                        unique_fragments_number = unique_fragments_number \
                            - 1 + bg_frequencies[i]
                        multimapper_number += 1
    return unique_fragments_number, multimapper_number, total_number


def coverage_profile_gn(params):
    # bamfile, chromosome, start, end, strand, sense_mate):
    '''Build coverage profile of a given enriched region'''
    # open the alignment file
    sam = pysam.AlignmentFile(params[0], "rb")
    # get the parameter info
    sense_mate = params[1]
    outfile = params[2]
    prefix_name = params[3]
    windows = params[4]
    name = params[5]
    paired = params[6]
    counter = params[7]
    cutoff = params[8]

    output = StringIO()
    csv_writer = writer(output)
    success = 0

    for row in windows.itertuples(index=False):
        chromosome = row[0]
        start = int(row[1])
        end = int(row[2])
        segment = list(np.arange(start, end))

        reads_codes_plus = set()
        reads_codes_minus = set()

        for read in sam.fetch(chromosome, start, end + 1):
            successful_block = 0
            for block in read.get_blocks():
                if any(
                    i in segment for i in list(
                        np.arange(block[0], block[-1]))):
                    successful_block += 1
                    break
            if successful_block > 0:
                if paired == 1:
                    if sense_mate == 1:
                        if read.is_reverse:
                            reads_codes_minus.add(read.qname)
                        else:
                            reads_codes_plus.add(read.qname)
                    elif sense_mate == 2:
                        if read.is_reverse:
                            reads_codes_plus.add(read.qname)
                        else:
                            reads_codes_minus.add(read.qname)

                elif paired == 2:
                    if read.is_proper_pair:
                        if sense_mate == 1:
                            if read.is_read1:
                                if read.is_reverse:
                                    reads_codes_minus.add(read.qname)
                                else:
                                    reads_codes_plus.add(read.qname)
                            elif read.is_read2:
                                if read.is_reverse:
                                    reads_codes_plus.add(read.qname)
                                else:
                                    reads_codes_minus.add(read.qname)
                        elif sense_mate == 2:
                            if read.is_read1:
                                if read.is_reverse:
                                    reads_codes_plus.add(read.qname)
                                else:
                                    reads_codes_minus.add(read.qname)
                            elif read.is_read2:
                                if read.is_reverse:
                                    reads_codes_minus.add(read.qname)
                                else:
                                    reads_codes_plus.add(read.qname)
        if reads_codes_plus:
            unique_fragments_pl, multimappers_pl, total_pl = estimate_frequency(
                name,
                reads_codes_plus,
                cutoff)

            if unique_fragments_pl > cutoff:
                success += 1
                df_pl = pd.Series([
                    chromosome,
                    start,
                    end,
                    str(chromosome) + ':' + str(start) + '-' + str(end) + ':+',
                    unique_fragments_pl,
                    '+',
                    multimappers_pl,
                    total_pl])
                csv_writer.writerow(df_pl)

        if reads_codes_minus:
            unique_fragments_mn, multimappers_mn, total_mn = estimate_frequency(
                name,
                reads_codes_minus,
                cutoff)

            if unique_fragments_mn > cutoff:
                success += 1
                df_mn = pd.Series([
                    chromosome,
                    start,
                    end,
                    str(chromosome) + ':' + str(start) + '-' + str(end) + ':-',
                    unique_fragments_mn,
                    '-',
                    multimappers_mn,
                    total_mn])
                csv_writer.writerow(df_mn)
    sam.close()
    del windows
    if success == 0:
        return 'Done'

    output.seek(0)
    df_output = pd.read_csv(output)

    pd.DataFrame(df_output).to_csv(
        os.path.join(
            outfile,
            str(prefix_name +
                '_' + chromosome +
                '_' + str(counter) +
                '.' + name + '.bed')),
        index=False,
        header=False,
        sep='\t')
    return 'Done!'


def coverage_profile_tr(params):
    # bamfile, chromosome, start, end, strand, sense_mate):
    '''Build coverage profile of a given enriched region'''
    # open the alignment file
    sam = pysam.AlignmentFile(params[0], "rb")
    # get the parameter info
    sense_mate = params[1]
    outfile = params[2]
    prefix_name = params[3]
    windows = params[4]
    name = params[5]
    paired = int(params[6])
    counter = params[7]
    cutoff = int(params[8])

    output = StringIO()
    csv_writer = writer(output)
    success = 0

    for row in windows.itertuples(index=False):
        chromosome = row[0]
        start = int(row[1])
        end = int(row[2])
        segment = list(np.arange(start, end))

        reads_codes_plus = set()

        for read in sam.fetch(chromosome, start, end):
            successful_block = 0
            for block in read.get_blocks():
                if any(
                    i in segment for i in list(
                        np.arange(block[0], block[-1]))):
                    successful_block += 1
                    break

            if successful_block > 0:
                if paired == 1:
                    if sense_mate == 1:
                        if not read.is_reverse:
                            reads_codes_plus.add(read.qname)

                    elif sense_mate == 2:
                        if read.is_reverse:
                            reads_codes_plus.add(read.qname)

                elif paired == 2:
                    if read.is_proper_pair:
                        if sense_mate == 1:
                            if read.is_read1:
                                if not read.is_reverse:
                                    reads_codes_plus.add(read.qname)
                            elif read.is_read2:
                                if read.is_reverse:
                                    reads_codes_plus.add(read.qname)
                        elif sense_mate == 2:
                            if read.is_read1:
                                if read.is_reverse:
                                    reads_codes_plus.add(read.qname)
                            elif read.is_read2:
                                if not read.is_reverse:
                                    reads_codes_plus.add(read.qname)
        if reads_codes_plus:
            unique_fragments_pl, multimappers_pl, total_pl = estimate_frequency(
                name,
                reads_codes_plus,
                cutoff)
            if unique_fragments_pl > cutoff:
                success += 1
                df_pl = pd.Series([
                    chromosome,
                    start,
                    end,
                    str(chromosome) + ':' + str(start) + '-' + str(end) + ':+',
                    unique_fragments_pl,
                    '+',
                    multimappers_pl,
                    total_pl])
                csv_writer.writerow(df_pl)
    sam.close()
    del windows
    if success == 0:
        return 'Done'

    output.seek(0)
    df_output = pd.read_csv(output)

    pd.DataFrame(df_output).to_csv(
        os.path.join(
            outfile,
            str(prefix_name +
                '_' + chromosome +
                '_' + str(counter) +
                '.' + name + '.bed')),
        index=False,
        header=False,
        sep='\t')

    return 'Done!'


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

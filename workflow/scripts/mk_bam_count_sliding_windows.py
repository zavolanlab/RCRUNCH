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
        help="prefix",
        required=True)

    parser.add_argument(
        "--window_f",
        help="Length of the sliding windows",
        required=False,
        default="300")

    parser.add_argument(
        "--window_b",
        help="Length of the sliding windows",
        required=False,
        default="300")

    parser.add_argument(
        "--data_type",
        help="genome or transcriptome",
        required=True)

    parser.add_argument(
        "--step_size",
        default=150,
        help="Step size for the generation of sliding windows.\
        If different from '--window-size', overlapping or \
        gapped sliding windows will be generated.")

    parser.add_argument(
        "--cutoff",
        help="Cutoff of coverage for window \
        to be included",
        required=False,
        default="1")

    parser.add_argument(
        "--threads",
        help="Number of threads for multiprocessing",
        required=False,
        default="1")

    parser.add_argument(
        "--background",
        help="Type of background. If foreground is used for background \
             estimation use the shifted_fg",
        required=False,
        default="standard",
        choices=['standard', 'shifted_fg'])

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
    background_type = options.background

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

    processes = threads
    number_of_windows = 1000

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

    sliding_windows_parameters = []

    for row in chromosome_info.itertuples():
        sliding_windows_parameters.append([
            row.chromosome,
            row.chromosome_length,
            window_f,
            window_b,
            step_size,
            background_type])
    try:
        pool = mp.Pool(processes=processes)
        windows = pool.map(
            sliding_windows,
            sliding_windows_parameters)
    finally:
        pool.close()
        pool.join()
    sys.stdout.write(f'Obtain sliding windows: {str(len(sliding_windows_parameters))}\n')
    sys.stdout.write(f'Obtain sliding windows done!\n')
    sys.stdout.flush()

    parameters = []
    windows = split_tables(windows, number_of_windows)
    for each_window in windows:
        parameters.append([
            bam_foreground,
            bam_background,
            sense_f,
            sense_b,
            each_window,
            paired_f,
            paired_b,
            cutoff,
            data_type,
            background_type])
    sys.stdout.write(f'Number of subsets to be multiprocessed: {str(len(parameters))}\n')
    sys.stdout.flush()
    windows = []

    with mp.Pool(processes=processes) as pool:
        with open(os.path.join(out, str(prefix) + '.bed'), 'w') as output:
            writer_final = writer(output, delimiter='\t')
            writer_final.writerow([
                'chromosome',
                'start',
                'end',
                'name',
                'type',
                'Reads_f',
                'Reads_b',
                'strand',
                'startbg',
                'endbg'])
            for rows in pool.imap_unordered(coverage_profile, parameters):
                writer_final.writerows(rows)
    return


def sliding_windows(params):
    chr_name = str(params[0])
    chr_len = int(params[1])
    fg_len = int(params[2] / 2)
    bg_len = int(params[3] / 2)
    step_size = int(params[4])
    background_type = params[5]
    if (chr_len <= step_size):
        window = pd.DataFrame()
        window['name'] = [chr_name]
        window['start'] = [0]
        window['end'] = [chr_len]
        window['startbg'] = [0]
        window['endbg'] = [chr_len]
        return window
    else:
        windows_center = np.arange(
            fg_len,
            chr_len + np.floor(step_size / 2),
            step_size)
        window = pd.DataFrame()
        window['name'] = [chr_name] * len(windows_center)
        window['start'] = windows_center - fg_len
        window['end'] = windows_center + fg_len
        if background_type == 'standard':
            window['startbg'] = windows_center - bg_len
            window['endbg'] = windows_center + bg_len
        elif background_type == 'shifted_fg':
            window['startbg'] = windows_center - fg_len - bg_len
            window['endbg'] = windows_center + fg_len + bg_len
        # windows close to boundaries --special cases
        window.loc[window.start < 0, 'end'] = 2 * fg_len
        window.loc[window.end > chr_len, 'start'] = chr_len - 2 * fg_len
        window.loc[window.start < 0, 'start'] = 0
        window.loc[window.end > chr_len, 'end'] = chr_len
        if background_type == 'standard':
            window.loc[window.startbg < 0, 'endbg'] = 2 * bg_len
            window.loc[window.endbg > chr_len, 'startbg'] = \
                chr_len - (2 * bg_len)

        elif background_type == 'shifted_fg':
            window.loc[window.startbg < 0, 'endbg'] = 2 * fg_len + 2 * bg_len
            window.loc[window.endbg > chr_len, 'startbg'] = \
                chr_len - 2 * fg_len - 2 * bg_len
        window.loc[window.startbg < 0, 'startbg'] = 0
        window.loc[window.endbg > chr_len, 'endbg'] = chr_len
    window.columns = [0, 1, 2, 3, 4]
    return window


def split_tables(dfs, length):
    dfs_new = []
    for window in dfs:
        window = window.sample(frac=1).reset_index(drop=True)
        if len(window) > length:
            number_of_slices = np.floor(len(window) / length)
            windows = np.array_split(window, number_of_slices)
            dfs_new.extend(windows)
        else:
            dfs_new.append(window)
    return dfs_new


def estimate_frequency(name, unique_fragments, cutoff):
    multimapper_number = 0
    total_number = len(unique_fragments)
    unique_fragments_number = total_number
    if unique_fragments_number >= cutoff:
        if name == 'fg':
            if fg_frequencies:
                for i in unique_fragments:
                    try:
                        unique_fragments_number = unique_fragments_number \
                            - 1 + fg_frequencies[i]
                        multimapper_number += 1
                    except KeyError:
                        continue
        elif name == 'bg':
            if bg_frequencies:
                for i in unique_fragments:
                    try:
                        unique_fragments_number = unique_fragments_number \
                            - 1 + bg_frequencies[i]
                        multimapper_number += 1
                    except KeyError:
                        continue
    return unique_fragments_number, multimapper_number, total_number


def coverage_profile(params):
    '''Build coverage profile of a given enriched region'''
    samfg = pysam.AlignmentFile(params[0], "rb")
    sambg = pysam.AlignmentFile(params[1], "rb")
    # get the parameter info
    sense_fg = params[2]
    sense_bg = params[3]
    windows = params[4]
    pairedfg = params[5]
    pairedbg = params[6]
    cutoff = params[7]
    approach = params[8]
    background_type = params[9]
    df = []
    for row in windows.itertuples(index=False):
        chromosome = str(row[0])
        start = int(row[1])
        end = int(row[2])
        startbg = int(row[3])
        endbg = int(row[4])
        window_options = {}
        window_options['fg'] = [start, end, samfg, sense_fg, pairedfg]
        window_options['bg'] = [startbg, endbg, sambg, sense_bg, pairedbg]
        coverages = {}
        for key in window_options:
            each_start = window_options[key][0]
            each_end = window_options[key][1]
            sam = window_options[key][2]
            sense_mate = window_options[key][3]
            paired = window_options[key][4]
            segment = set(list(np.arange(each_start, each_end)))
            reads_codes_plus = set()
            reads_codes_minus = set()
            if paired == 1:
                if sense_mate == 1:
                    for read in sam.fetch(chromosome, each_start, each_end + 1):
                        if read.is_reverse:
                            if read.reference_end in segment:
                                reads_codes_minus.add(read.qname)
                        else:
                            if read.reference_start in segment:
                                reads_codes_plus.add(read.qname)            
                elif sense_mate == 2:
                    for read in sam.fetch(chromosome, start, end + 1):
                        if read.is_reverse:
                            if read.reference_end in segment:
                                reads_codes_plus.add(read.qname)
                        else:
                            if read.reference_start in segment:
                                reads_codes_minus.add(read.qname)
            elif paired == 2:
                if sense_mate == 1:
                    for read in sam.fetch(chromosome, each_start, each_end + 1):
                        if read.is_read1:
                            if read.is_proper_pair:
                                if not read.is_reverse:
                                    if read.reference_start in segment:
                                        reads_codes_plus.add(read.qname)
                                else:
                                    if read.reference_end in segment:
                                        reads_codes_minus.add(read.qname)
                elif sense_mate == 2:
                    for read in sam.fetch(chromosome, each_start, each_end + 1):
                        if read.is_read2:
                            if read.is_proper_pair:
                                if read.is_reverse:
                                    if read.reference_end in segment:
                                        reads_codes_minus.add(read.qname)
                                else:
                                    if read.reference_start in segment:
                                        reads_codes_plus.add(read.qname)
            unique_fr_pl = 0
            multimappers_pl = 0
            total_pl = 0
            unique_fr_mn = 0
            multimappers_mn = 0
            total_mn = 0
            if reads_codes_plus:
                unique_fr_pl, multimappers_pl, total_pl = estimate_frequency(
                    key,
                    reads_codes_plus,
                    cutoff)
            if reads_codes_minus:
                unique_fr_mn, multimappers_mn, total_mn = estimate_frequency(
                    key,
                    reads_codes_minus,
                    cutoff)
            coverages[key] = [
                unique_fr_pl, multimappers_pl, total_pl,
                unique_fr_mn, multimappers_mn, total_mn]

        if background_type == 'standard':
            if (coverages['fg'][0] >= cutoff):
                df.append([
                    chromosome,
                    start,
                    end,
                    f'{str(chromosome)}:{str(start)}-{str(end)}:{str(startbg)}-{str(endbg)}:+',
                    approach,
                    coverages['fg'][0],
                    coverages['bg'][0],
                    '+',
                    startbg,
                    endbg])
            if approach == 'genome':
                if (coverages['fg'][3] >= cutoff):
                    df.append([
                        chromosome,
                        start,
                        end,
                        f'{str(chromosome)}:{str(start)}-{str(end)}:{str(startbg)}-{str(endbg)}:-',
                        approach,
                        coverages['fg'][3],
                        coverages['bg'][3],
                        '-',
                        startbg,
                        endbg])

        elif background_type == 'shifted_fg':
            if (coverages['fg'][0] >= cutoff):
                df.append([
                    chromosome,
                    start,
                    end,
                    f'{str(chromosome)}:{str(start)}-{str(end)}:{str(startbg)}-{str(endbg)}:+',
                    approach,
                    coverages['fg'][0],
                    coverages['bg'][0] - coverages['fg'][0],
                    '+',
                    startbg,
                    endbg])
            if approach == 'genome':
                if (coverages['fg'][3] >= cutoff):
                    df.append([
                        chromosome,
                        start,
                        end,
                        f'{str(chromosome)}:{str(start)}-{str(end)}:{str(startbg)}-{str(endbg)}:-',
                        approach,
                        coverages['fg'][3],
                        coverages['bg'][3] - coverages['fg'][3],
                        '-',
                        startbg,
                        endbg])
    samfg.close()
    sambg.close()
    del windows
    return df


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

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
from collections import defaultdict
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
        "--bed_foreground",
        dest="bed_foreground",
        help="Foreground reads aligned to genome",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--foreground_frequencies",
        dest="foreground_frequencies",
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
        "--bed_background",
        dest="bed_background",
        help="Background reads aligned to genome",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--background_frequencies",
        dest="background_frequencies",
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
    bed_foreground = options.bed_foreground
    bed_background = options.bed_background
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
    if options.foreground_frequencies:
        fg_frequencies = pd.read_csv(
            options.foreground_frequencies,
            header=0,
            sep='\t',
            index_col=0,
            comment='#',
            engine='python')
        fg_frequencies = fg_frequencies['Frequency'].to_dict()
    else:
        fg_frequencies = {}

    global bg_frequencies
    if options.background_frequencies:
        bg_frequencies = pd.read_csv(
            options.background_frequencies,
            header=0,
            sep='\t',
            index_col=0,
            comment='#',
            engine='python')
        bg_frequencies = bg_frequencies['Frequency'].to_dict()
    else:
        bg_frequencies = {}

    processes = int(options.threads)
    number_of_windows = 50000

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
    sys.stdout.write(
        f'Obtain sliding windows: {str(len(sliding_windows_parameters))}\n')
    sys.stdout.write('Obtain sliding windows done!\n')
    sys.stdout.flush()
    countsfg = get_total_coverage(
        bed_foreground, chromosome_info, sense_f, paired_f, 'fg')
    countsbg = get_total_coverage(
        bed_background, chromosome_info, sense_b, paired_b, 'bg')
    parameters = []
    windows = split_tables(windows, number_of_windows)
    sys.stdout.write(f'Windows:{str(len(windows))}\n')
    sys.stdout.flush()
    for each_window in windows:
        parameters.append([
            countsfg,
            countsbg,
            each_window,
            cutoff,
            data_type,
            background_type])
    sys.stdout.write(
        f'Number of subsets to be multiprocessed: {str(len(parameters))}\n')
    sys.stdout.flush()
    windows = []

    pool = mp.Pool(processes=processes)
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
        # try:
        for rows in pool.imap_unordered(coverage_profile, parameters):
            writer_final.writerows(rows)
        # except:
        #     print("Outer exception caught!")
    pool.close()
    pool.join()
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


def get_total_coverage(bedfile, chromosome_info, sense, paired, sample_type):
    counts = defaultdict(dict)
    strands = ['+', '-']
    for chromosome in chromosome_info['chromosome'].values:
        chromosome = str(chromosome)
        counts[chromosome] = defaultdict(dict)
        for strand in strands:
            counts[chromosome][strand] = defaultdict(int)
    # iterate through the bedfiles and assign reads to appropriate chromosome
    with open(bedfile, "r") as infile:
        if sample_type == 'fg':
            if paired == 1:
                if sense == 1:
                    for line in infile:
                        fields = line.split()
                        try:
                            if(fields[5] == '+'):
                                try:
                                    fragment_frequency = fg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][fields[5]][int(fields[1])] += fragment_frequency
                            else:
                                try:
                                    fragment_frequency = fg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][fields[5]][int(fields[2])] += fragment_frequency
                        except KeyError:
                            print(line)
                elif sense == 2:
                    for line in infile:
                        fields = line.split()
                        try:
                            if(fields[5] == '+'):
                                each_strand = '-'
                                try:
                                    fragment_frequency = fg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][each_strand][int(fields[1])] += fragment_frequency
                            else:
                                each_strand = '+'
                                try:
                                    fragment_frequency = fg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][each_strand][int(fields[2])] += fragment_frequency
                        except KeyError:
                            print(line)
            elif paired == 2:
                if sense == 1:
                    for line in infile:
                        fields = line.split()
                        if fields[3].endswith('/2'):
                            continue
                        try:
                            if(fields[5] == '+'):
                                try:
                                    fragment_frequency = fg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][fields[5]][int(fields[1])] += fragment_frequency
                            else:
                                try:
                                    fragment_frequency = fg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][fields[5]][int(fields[2])] += fragment_frequency
                        except KeyError:
                            print(line)
                elif sense == 2:
                    for line in infile:
                        fields = line.split()
                        if fields[3].endswith('/1'):
                            continue
                        try:
                            if(fields[5] == '+'):
                                try:
                                    fragment_frequency = fg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][fields[5]][int(fields[1])] += fragment_frequency
                            else:
                                try:
                                    fragment_frequency = fg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][fields[5]][int(fields[2])] += fragment_frequency
                        except KeyError:
                            print(line)
        elif sample_type == 'bg':
            if paired == 1:
                if sense == 1:
                    for line in infile:
                        fields = line.split()
                        try:
                            if(fields[5] == '+'):
                                try:
                                    fragment_frequency = bg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][fields[5]][int(fields[1])] += fragment_frequency
                            else:
                                try:
                                    fragment_frequency = bg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][fields[5]][int(fields[2])] += fragment_frequency
                        except KeyError:
                            print(line)
                elif sense == 2:
                    for line in infile:
                        fields = line.split()
                        try:
                            if(fields[5] == '+'):
                                each_strand = '-'
                                try:
                                    fragment_frequency = bg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][each_strand][int(fields[1])] += fragment_frequency
                            else:
                                each_strand = '+'
                                try:
                                    fragment_frequency = bg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][each_strand][int(fields[2])] += fragment_frequency
                        except KeyError:
                            print(line)
            elif paired == 2:
                if sense == 1:
                    for line in infile:
                        fields = line.split()
                        if fields[3].endswith('/2'):
                            continue
                        try:
                            if(fields[5] == '+'):
                                try:
                                    fragment_frequency = bg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][fields[5]][int(fields[1])] += fragment_frequency
                            else:
                                try:
                                    fragment_frequency = bg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][fields[5]][int(fields[2])] += fragment_frequency
                        except KeyError:
                            print(line)
                elif sense == 2:
                    for line in infile:
                        fields = line.split()
                        if fields[3].endswith('/1'):
                            continue
                        try:
                            if(fields[5] == '+'):
                                try:
                                    fragment_frequency = bg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][fields[5]][int(fields[1])] += fragment_frequency
                            else:
                                try:
                                    fragment_frequency = bg_frequencies[
                                        fields[3].split('/')[0]]
                                except KeyError:
                                    fragment_frequency = 1
                                counts[fields[0]][fields[5]][int(fields[2])] += fragment_frequency
                        except KeyError:
                            print(line)
    return counts


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


def coverage_profile(params):
    '''Build coverage profile of a given enriched region'''
    countsfg = params[0]
    countsbg = params[1]
    # get the parameter info
    windows = params[2]
    cutoff = params[3]
    approach = params[4]
    background_type = params[5]
    df = []
    if approach == 'genome':
        for row in windows.itertuples(index=False):
            chromosome = row[0]
            start = int(row[1])
            end = int(row[2])
            startbg = int(row[3])
            endbg = int(row[4])
            wincountplus = 0
            for i in range(start, end + 1):
                wincountplus += countsfg[chromosome]['+'][i]
            if wincountplus >= cutoff:
                wincountplusbg = 0
                for i in range(startbg, endbg + 1):
                    wincountplusbg += countsbg[chromosome]['+'][i]
                if background_type == 'standard':
                    df.append([
                        chromosome,
                        start,
                        end,
                        f'{str(chromosome)}:{str(start)}-{str(end)}:{str(startbg)}-{str(endbg)}:+',
                        approach,
                        wincountplus,
                        wincountplusbg,
                        '+',
                        startbg,
                        endbg])
                elif background_type == 'shifted_fg':
                    df.append([
                        chromosome,
                        start,
                        end,
                        f'{str(chromosome)}:{str(start)}-{str(end)}:{str(startbg)}-{str(endbg)}:+',
                        approach,
                        wincountplus,
                        wincountplusbg - wincountplus,
                        '+',
                        startbg,
                        endbg])
            wincountminus = 0
            for i in range(start, end + 1):
                wincountminus += countsfg[chromosome]['-'][i]
            if wincountminus >= cutoff:
                wincountminusbg = 0
                for i in range(startbg, endbg + 1):
                    wincountminusbg += countsbg[chromosome]['-'][i]
                if background_type == 'standard':
                    df.append([
                        chromosome,
                        start,
                        end,
                        f'{str(chromosome)}:{str(start)}-{str(end)}:{str(startbg)}-{str(endbg)}:-',
                        approach,
                        wincountminus,
                        wincountminusbg,
                        '-',
                        startbg,
                        endbg])
                elif background_type == 'shifted_fg':
                    df.append([
                        chromosome,
                        start,
                        end,
                        f'{str(chromosome)}:{str(start)}-{str(end)}:{str(startbg)}-{str(endbg)}:-',
                        approach,
                        wincountminus,
                        wincountminusbg - wincountminus,
                        '-',
                        startbg,
                        endbg])
    else:
        for row in windows.itertuples(index=False):
            chromosome = row[0]
            start = row[1]
            end = row[2]
            startbg = row[3]
            endbg = row[4]
            wincountplus = 0
            for i in range(start, end + 1):
                wincountplus += countsfg[chromosome]['+'][i]
            if wincountplus >= cutoff:
                wincountplusbg = 0
                for i in range(startbg, endbg + 1):
                    wincountplusbg += countsbg[chromosome]['+'][i]
                if background_type == 'standard':
                    df.append([
                        chromosome,
                        start,
                        end,
                        f'{str(chromosome)}:{str(start)}-{str(end)}:{str(startbg)}-{str(endbg)}:+',
                        approach,
                        wincountplus,
                        wincountplusbg,
                        '+',
                        startbg,
                        endbg])
                elif background_type == 'shifted_fg':
                    df.append([
                        chromosome,
                        start,
                        end,
                        f'{str(chromosome)}:{str(start)}-{str(end)}:{str(startbg)}-{str(endbg)}:+',
                        approach,
                        wincountplus,
                        wincountplusbg - wincountplus,
                        '+',
                        startbg,
                        endbg])
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

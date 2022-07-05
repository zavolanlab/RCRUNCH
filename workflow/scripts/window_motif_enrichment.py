# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Author : Katsantoni Maria
# Company: Mihaela Zavolan group, Biozentrum, Basel
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import os
from pyfasta import Fasta
import multiprocessing as mp
import logging
import subprocess
import glob
from subprocess import Popen, PIPE, CalledProcessError
from contextlib import closing
from time import sleep
logger = mp.log_to_stderr(logging.DEBUG)

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():
    """ Main function """

    __doc__ = "get motif enrichment"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    # parser.add_argument(
    #     "--genome_fasta",
    #     dest="genome_fasta",
    #     help="Fasta file of the genome to get \
    #     sequences of the genomic peaks.",
    #     required=True,
    #     metavar="FILE")

    # parser.add_argument(
    #     "--transcriptome_fasta",
    #     dest="transcriptome_fasta",
    #     help="Fasta file of the transcriptome to get \
    #         sequences of the transcriptomic peaks.",
    #     required=False,
    #     metavar="FILE")

    # parser.add_argument(
    #     "--chromosomes_t",
    #     dest="chromosomes_t",
    #     help="File containing chromosome names and \
    #         lengths provided during STAR indexing",
    #     required=False,
    #     metavar="FILE")

    # parser.add_argument(
    #     "--chromosomes_g",
    #     dest="chromosomes_g",
    #     help="File containing chromosome names and \
    #         lengths provided during STAR indexing",
    #     required=True,
    #     metavar="FILE")

    # parser.add_argument(
    #     "--genome_tag",
    #     dest="genome_tag",
    #     help="genome name (eg:hg38)",
    #     required=True)
    parser.add_argument(
        "--sequences",
        dest="sequences",
        help="folder containing the fasta sequences",
        nargs='+',
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--outpath",
        dest="outpath",
        help="output file containing the parameters needed for \
              running motevo in refinement mode",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--wm",
        dest="wm",
        help="motif paths",
        required=True,
        metavar="FILE")
    
    # parser.add_argument(
    #     "--window_length",
    #     dest="window_length",
    #     help="size of subsequences",
    #     required=False,
    #     default="400")
     
    # parser.add_argument(
    #     "--step_size",
    #     dest="step_size",
    #     help="size of overlapping windows",
    #     required=False,
    #     default="300")

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

    # if options.chromosomes_t:
    #     transcriptome = pd.read_csv(
    #         options.chromosomes_t,
    #         header=None,
    #         sep='\t',
    #         index_col=None,
    #         comment='#',
    #         engine='python')

    #     transcriptome.columns = ['chromosome', 'chromosome_length']
    #     transcriptome_names = list(transcriptome['chromosome'].values)
    #     # transcriptome.set_index('chromosome', inplace=True, drop=True)

    #     ft = Fasta(
    #         options.transcriptome_fasta,
    #         key_fn=lambda key: key.split()[0])
    # else:
    #     transcriptome_names = []

    # genome = pd.read_csv(
    #     options.chromosomes_g,
    #     header=None,
    #     sep='\t',
    #     index_col=None,
    #     comment='#',
    #     engine='python')

    # genome.columns = ['chromosome', 'chromosome_length']
    # # genome.set_index('chromosome', inplace=True, drop=True)

    # fg = Fasta(
    #     options.genome_fasta,
    #     key_fn=lambda key: key.split()[0])
    # if not os.path.exists(options.outpath):
    #     os.makedirs(options.outpath)
    # temporary_path = options.outpath.replace(os.path.basename(options.outpath), 'temp')
    # temp_seq_path = os.path.join(temporary_path, 'sequences')
    # if not os.path.exists(temporary_path):
    #     os.makedirs(temporary_path)
    # if not os.path.exists(temp_seq_path):
    #     os.makedirs(temp_seq_path)
    # # break the windows into subsets since these are more than 1 million
    # # and motevo works optimally with small input sequences (less than 75000)
    # ##################################################################
    # number_of_windows = 20000
    # window_length = int(options.window_length) 
    # step_size = int(options.step_size) 
    # sliding_windows_parameters = []
    # for row in genome.itertuples():
    #     # if str(row.chromosome) != '21':
    #     #     continue
    #     sliding_windows_parameters.append([
    #         row.chromosome,
    #         row.chromosome_length,
    #         window_length,
    #         step_size,
    #         'genome'])
    # if options.chromosomes_t:
    #     for row in transcriptome.itertuples():
    #         sliding_windows_parameters.append([
    #             row.chromosome,
    #             row.chromosome_length,
    #             window_length,
    #             step_size,
    #             'transcriptome'])

    # try:
    #     processes = 30  #* threads
    #     pool = mp.Pool(processes=processes)
    #     new_peaks_all = pool.map(
    #         sliding_windows,
    #         sliding_windows_parameters)
    # finally:
    #     pool.close()
    #     pool.join()
    parameters = []
    for sequence in options.sequences:
    
    # counter = 0
    # new_peaks_all = split_tables(new_peaks_all, number_of_windows)
    # for new_peaks in new_peaks_all:
    # ##################################################################
    # ##################################################################
    # # parameters = []
    # # for new_peaks in pd.read_csv(
    # #         options.all_regions,
    # #         header=0,
    # #         sep='\t',
    # #         index_col=0,
    # #         comment='#',
    # #         engine='python',
    # #         chunksize=5000):
    #     counter += 1
    #     with open(os.path.join(temp_seq_path, 'test_' + str(counter) + '.fasta'), 'w') as test_set:
    #         for index, row in new_peaks.iterrows():
    #             chromosome = str(row['chromosome'])
    #             start = int(row['start'])
    #             end = int(row['end'])
    #             if chromosome in transcriptome_names:
    #                 strand = '+'
    #                 info = ft
    #             else:
    #                 strand = row['strand']
    #                 info = fg
    #             sequence = info.sequence({
    #                 'chr': chromosome,
    #                 'start': start,
    #                 'stop': end,
    #                 'strand': strand}, one_based=False)
    #             sequence = sequence.replace('n', 'N')

    #             line = f'>>{options.genome_tag}_{chromosome}_{start}_{end}_{strand}\n{sequence}\n\n'
    #             test_set.write(line)
       
        sequence = os.path.abspath(sequence)

        # ------------------------------------------------------------------
        # this part is due to Motevo which runs in the local folder
        # and does not take input path as an argument
        initial_path = os.getcwd()
        # os.chdir(options.outpath)
        # outpath = os.getcwd()
        # -------------------------------------------------------------------
        
        # wm = os.path.join(initial_path, options.wm)
        parameters.append(
            [wm, abs_path, options.genome_tag,
                os.path.abspath(options.outpath),
                os.path.abspath(temporary_path), counter])
    try:
        pool = mp.Pool()
        results = pool.imap_unordered(motevo_process, parameters)
    finally:
        pool.close()
        pool.join()

    wm = os.path.basename(options.wm)
    outfile_name = os.path.join(options.outpath, f'{wm}_posterior_sites.bed')
    merge_posteriors = pd.DataFrame()
    for posterior in glob.glob(os.path.join(temporary_path, wm, '*', 'posterior_sites')):
        try:
            each_posterior = pd.read_csv(posterior, sep='\s+', header=None, index_col=None)
        except pd.io.common.EmptyDataError:
            sys.stdout.write(f'{posterior} was empty')
            continue
        each_posterior = each_posterior.drop_duplicates()
        each_posterior.columns = ['position', 'strand', 'posterior', 'motif', 'window']
        each_posterior[['genome', 'chromosome', 'windowstart', 'windowend', 'strand1']] = each_posterior['window'].str.split(pat="_", expand=True)
        each_posterior[['start', 'end']] = each_posterior['position'].str.split(pat="-", expand=True)
        each_posterior[['start', 'end', 'windowstart', 'windowend']] = each_posterior[['start', 'end', 'windowstart', 'windowend']].astype('int')
        each_posterior['start'] = each_posterior['windowstart'] + each_posterior['start']
        each_posterior['end'] = each_posterior['windowstart'] + each_posterior['end']
        each_posterior.drop(['strand1', 'window', 'position','windowstart', 'windowend'], axis=1, inplace=True)
        each_posterior = each_posterior[['chromosome','start', 'end', 'motif', 'posterior', 'strand', 'genome']]
        if merge_posteriors.empty:
            merge_posteriors = each_posterior
        else:
            merge_posteriors = pd.concat([merge_posteriors, each_posterior])
    merge_posteriors = merge_posteriors.drop_duplicates()
    merge_posteriors.to_csv(outfile_name, index=False, header=False, sep='\t')
    success = open(os.path.join(initial_path, options.outpath, 'done'), 'w')
    success.close()
    return

############################################
############################################
def sliding_windows(params):
    chr_name = str(params[0])
    chr_len = int(params[1])
    window_length = int(params[2] / 2)
    step_size = int(params[3])
    annotation = params[4]
    if (chr_len <= step_size):
        window = pd.DataFrame()
        window['chromosome'] = [chr_name]
        window['start'] = [0]
        window['end'] = [chr_len]
        return window
    else:
        windows_center = np.arange(step_size, chr_len, step_size)
        window = pd.DataFrame()
        window['chromosome'] = [chr_name] * len(windows_center)
        window['start'] = windows_center - window_length
        window['end'] = windows_center + window_length
        # windows close to boundaries --special cases
        window.loc[window.start < 0, 'end'] = 2 * window_length
        window.loc[window.end > chr_len, 'start'] = chr_len - 2 * window_length
        window.loc[window.start < 0, 'start'] = 0
        window.loc[window.end > chr_len, 'end'] = chr_len
        window['strand'] = '+'
    if annotation== 'transcriptome':
        window1 = window.copy(deep=True)
        window1['strand'] = '-'
        window = pd.concat([window, window1])
    return window


def split_tables(dfs, length):
    dfs_new = []
    for window in dfs:
        # randomise rows, so that windows with high coverage are dispersed
        window = window.sample(frac=1).reset_index(drop=True)
        if len(window) > length:
            number_of_slices = np.floor(len(window) / length)
            windows = np.array_split(window, number_of_slices)
            dfs_new.extend(windows)
        else:
            dfs_new.append(window)
    return dfs_new
############################################
############################################

def motevo_process(params):
    '''prior and posterior calls to Motevo and beta optimisation'''
    wm = params[0]
    test = params[1]
    genome_tag = params[2]
    outpath = params[3]
    temppath = params[4]
    counter = str(params[5])   
    wmname = wm.split('/')[-1] 
    sys.stdout.write(f'start {wmname}\n')
    sys.stdout.flush()
    # specify the wm path
    each_path = os.path.join(temppath, wmname, counter)
    if not os.path.exists(each_path):
        os.makedirs(each_path)
    os.chdir(each_path)
    # sys.stdout.write(os.getcwd())
    # sys.stdout.flush()
    # Calculate the posteriors for each of the motifs
    bgprior = 0.99
    parameters_posterior = 'posterior_params'
    # make the motevo parameter file for the posteriors
    motevo_parameters_posterior(bgprior, parameters_posterior, genome_tag)
    posterior_motevo = ['motevo', test, parameters_posterior, wm]

    # call Motevo
    with open('posterior_report', "wb") as outfile:
        try:
            posteriorcall = Popen(posterior_motevo, stdout=outfile, stderr=PIPE)
            posteriorcall.wait()
            stdout, stderr = posteriorcall.communicate()
            return 'done'
        except CalledProcessError as err:
            print(f"Error ocurred: posterior {wmname} {err}")
            return None


def motevo_parameters_posterior(bgprior, parampath, genome_tag,
                                ufewmprior=0, ufewmfile=None,
                                ufewmlen=8, ufeprint=0, markovorderbg=0,
                                bga=0.25, bgt=0.25, bgc=0.25, bgg=0.25,
                                restrictparses=0, printsiteals=0,
                                minposterior=0.01):
    tree = str('TREE (' + genome_tag + ':1);')
    with open(parampath, 'w') as paramfile:
        paramfile.write(
            'refspecies ' + genome_tag + '\n' +
            tree + '\n' +
            'Mode TFBS\n' +
            'EMprior 0\n' +
            'markovorderBG ' + str(markovorderbg) + '\n' +
            'bgprior ' + str(bgprior) + '\n' +
            'bg A ' + str(bga) + '\n' +
            'bg T ' + str(bgt) + '\n' +
            'bg G ' + str(bgg) + '\n' +
            'bg C ' + str(bgc) + '\n' +
            'restrictparses ' + str(restrictparses) + '\n' +
            'singlestrand 1' + '\n' +
            'sitefile posterior_sites\n' +
            'priorfile posteriors\n' +
            'printsiteals ' + str(printsiteals) + '\n' +
            'minposterior ' + str(minposterior) + '\n')
    return


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

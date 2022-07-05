# -----------------------------------------------------------------------------
# Author : Katsantoni Maria
# Company: Mihaela Zavolan, Biozentrum, Basel
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# This script is part of the RCRUNCH pipeline, which is used for analysing CLIP
# data. In this step any regions deemed as significantly enriched in binding
# sites are further processed by this script, which fits peaks in those regions
# and assesses whether these new fitted peaks correspond to significant binding
# sites.
#
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Import needed (external) modules
# -----------------------------------------------------------------------------
import sys
import numpy as np
import pandas as pd
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from argparse import ArgumentParser, RawTextHelpFormatter
import pysam
import multiprocessing as mp
from io import StringIO
from csv import writer
import logging
import time
matplotlib.use('Agg')
logger = mp.log_to_stderr(logging.DEBUG)
sns.set_style(style='white')


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():
    """ Main function """

    __doc__ = "Fit individual binding sites in regions enriched\
     in binding sites"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    # ENRICHED REGIONS STEP INFO
    parser.add_argument(
        "--enriched_regions",
        dest="enriched_regions",
        help="table of the detected enriched regions",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--parameters",
        dest="parameters",
        help="parameters file obtained during enrichment detection",
        required=True,
        metavar="FILE")
    # ---------------------------------------------------

    # INFO BASED ON MAPPING(OR GUESSED)
    parser.add_argument(
        "--fragment_size",
        dest="fragment_size",
        help="fragment size",
        required=False,
        default=100,
        metavar="Number")
    # ---------------------------------------------------

    parser.add_argument(
        "--window_size",
        dest="window_size",
        help="length of the sliding windows",
        required=True,
        metavar="Number")
    # ---------------------------------------------------

    # GENOME RELATED VARIABLES
    # genome chrom names
    parser.add_argument(
        "--chromosomes_g",
        dest="chromosomes_g",
        help="File containing chrom names and lengths \
        provided during STAR indexing, required",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--bam_foreground_g",
        dest="bam_foreground_g",
        help="foreground reads aligned to genome",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--bam_background_g",
        dest="bam_background_g",
        help="background reads aligned to genome",
        required=True,
        metavar="FILE")
    # ---------------------------------------------------

    # TRANSCRIPTOME RELATED VARIABLES
    # transcript names
    parser.add_argument(
        "--chromosomes_t",
        dest="chromosomes_t",
        help="File containing transcript names and lengths \
        provided during STAR indexing, required",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--bam_foreground_t",
        dest="bam_foreground_t",
        help="foreground reads aligned to transcriptome",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--bam_background_t",
        dest="bam_background_t",
        help="background reads aligned to transcriptome",
        required=False,
        metavar="FILE")
    # ---------------------------------------------------

    parser.add_argument(
        "--bam_fg_fq",
        dest="bam_foreground_frequencies",
        help="Read frequencies in foreground",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--bam_bg_fq",
        dest="bam_background_frequencies",
        help="Read frequencies in background",
        required=False,
        metavar="FILE")

    # LIBRARY TYPE AND SENSE MATE
    parser.add_argument(
        "--paired_f",
        dest="paired_f",
        choices=['1', '2'],
        help="whether data is paired in the foreground,\
        1=single-end, 2=paired-end ",
        required=True,
        metavar="Number")

    parser.add_argument(
        "--paired_b",
        dest="paired_b",
        choices=['1', '2'],
        help="whether data is paired in the background,\
        1=single-end, 2=paired-end ",
        required=True,
        metavar="Number")

    parser.add_argument(
        "--sense_f",
        dest="sense_f",
        choices=['1', '2'],
        help=" sense mate of the foreground",
        required=True,
        metavar="Number")

    parser.add_argument(
        "--sense_b",
        dest="sense_b",
        choices=['1', '2'],
        help=" sense mate of the background",
        required=True,
        metavar="Number")

    # ---------------------------------------------------

    # VARIABLES NOT USED ALWAYS
    parser.add_argument(
        "--log_likelihood_cutoff",
        dest="log_likelihood_cutoff",
        help="log likelihood threshold value to stop the EM",
        required=False,
        default="0.00001")

    parser.add_argument(
        "--length_cutoff",
        dest="length_cutoff",
        help="peak length cutoff",
        required=False,
        metavar="Number",
        default=2)

    parser.add_argument(
        "--background",
        help="Type of background. If foreground is used for background \
             estimation use the shifted_fg",
        required=False,
        default="standard",
        choices=['standard', 'shifted_fg'])

    parser.add_argument(
        "--verbose",
        action='store_true',
        help="Plot for each region and calculations ")

    parser.add_argument(
        "--threads",
        help="Plot for each region and calculations ",
        required=False,
        default=12)

    # ---------------------------------------------------
    # OUTPUT
    parser.add_argument(
        "--out",
        dest="out",
        help="output folder path",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--out_folder",
        dest="out_folder",
        help="destination folder for the plots",
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

    parameters = pd.read_csv(
        options.parameters,
        header=0,
        sep='\t',
        index_col=0,
        comment='#',
        engine='python')

    enriched_region_windows = pd.read_csv(
        options.enriched_regions,
        header=0,
        sep='\t',
        index_col=0,
        comment='#',
        engine='python')
 

    # enriched_region_windows = enriched_region_windows.head(1000)
    # 22      49515450        49537050 
    # enriched_region_windows = enriched_region_windows[(enriched_region_windows['chromosome']=='20') & (enriched_region_windows['start']>36000000)] # & (enriched_region_windows['start']<37050190)]
    # enriched_region_windows = enriched_region_windows[(enriched_region_windows['chromosome']=='1') & (enriched_region_windows['start']==30931350)] # & (enriched_region_windows['start']<37050190)]

    print(enriched_region_windows)
    # if len(enriched_region_windows) == 1:
    #    enriched_region_windows = enriched_region_windows.to_frame().T

    # --------------------------------------------------------------------------
    # enriched_region_windows = \
    #     enriched_region_windows.sample(frac=1).reset_index(drop=True)

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

    chrom = pd.read_csv(
        options.chromosomes_g,
        header=None,
        sep='\t',
        index_col=0,
        comment='#',
        engine='python')

    chrom_names = list(chrom.index.values.astype(str))
    chrom.columns = ['length']

    gene_list = enriched_region_windows[
        enriched_region_windows['chromosome'].astype(str).isin(
            chrom_names)]
    gene_list['size'] = gene_list['end'] - gene_list['start']
    gene_list.sort_values(by='size', ascending=False, inplace=True)
    # gene_list.sort_values(by='size', ascending=True, inplace=True)
    # gene_list = gene_list.iloc[1002:1003]

    paired_f = int(options.paired_f)
    sense_f = int(options.sense_f)
    # split very big windows into smaller ones
    count = 0
    split_params = []
    split_windows = []
    intact_windows = []
    for index, row in gene_list.iterrows():
        count += 1
        if row['end'] - row['start'] > 800:
            if str(row['chromosome']) in chrom_names:
                split_params.append(
                    [row.copy(deep=True), options.bam_foreground_g,
                     paired_f, sense_f, options.window_size, 400])
        else:
            intact_windows.append(row.copy(deep=True))
    if intact_windows:
        gene_list1 = pd.concat(intact_windows, axis=1).T
    else:
        gene_list1 = pd.DataFrame()
    time1 = time.time()
    if split_big_regions:
        processes = int(options.threads)
        pool = mp.Pool(processes=processes)
        for x in pool.imap_unordered(
                split_big_regions, split_params):
            split_windows += x
        if split_windows:
            gene_list2 = pd.concat(split_windows, axis=1).T
        else:
            gene_list2 = pd.DataFrame()
    time2 = time.time() - time1
    print(time2)
    gene_list = pd.concat([gene_list1, gene_list2])
    gene_list1 = []
    gene_list2 = []

    # gene_list = gene_list[(gene_list['start'] == 36050195) ] #& (gene_list['end'] <= 49529000)]


    total_parameters = []
    for index, row in gene_list.iterrows():
        total_parameters.append([
            row,
            # each row looks like:
            # 17      X       154562250       154562550       +
            # number     chrom  start_position  end_position    strand
            options.bam_foreground_g,
            options.bam_background_g,
            parameters,
            # parameters file looks like:
            # nr sigma ro mi z_score fdr_cutoff log likelihood total_b total_f
            # 0 0.28 0.99  -0.19 3.94  0.1 -3288185.48  13894575.0  25780502.0
            options.out_folder,
            options.fragment_size,
            options.log_likelihood_cutoff,
            options.window_size,
            options.sense_f,
            options.sense_b,
            options.length_cutoff,
            options.paired_f,
            options.paired_b,
            options.background,
            options.verbose])

    # --------------------------------------------------------------------------
    if options.chromosomes_t:
        transcriptome = pd.read_csv(
            options.chromosomes_t,
            header=None,
            sep='\t',
            index_col=0,
            comment='#',
            engine='python')

        transcriptome_names = list(transcriptome.index.values)
        transcriptome.columns = ['length']

        transcript_list = enriched_region_windows[
            enriched_region_windows['chromosome'].astype(str).isin(
                transcriptome_names)]

        for index, row in transcript_list.iterrows():
            total_parameters.append([
                row,
                # each row looks like:
                # 33      ENST00000638183.1       2850    3150    +
                # number     transcript  start_position  end_position    strand
                options.bam_foreground_t,
                options.bam_background_t,
                parameters,
                options.out_folder,
                options.fragment_size,
                options.log_likelihood_cutoff,
                options.window_size,
                options.sense_f,
                options.sense_b,
                options.length_cutoff,
                options.paired_f,
                options.paired_b,
                options.background,
                options.verbose])
    
    processes = int(options.threads)
    pool = mp.Pool(processes=processes)
    df_total = list(pool.imap_unordered(
        second_EM_define_peaks,
        total_parameters))

    if len(df_total) > 0:
        df_final = pd.concat(
            df_total,
            axis=0,
            join='outer',
            ignore_index=True)
        df_final.sort_values(by=['z_score'], ascending=False, inplace=True)
    else:
        df_final = pd.DataFrame()
    df_final.to_csv(
        os.path.join(options.out),
        sep='\t',
        index=False,
        header=True)
    del enriched_region_windows
    time3 = time.time() - time1 - time2
    print(time3)
    return


def split_big_regions(params):
    [row, bam, paired, sense, window_size, step] = params
    window_size = int(window_size)
    rows = []
    if (paired == 1) & (sense == 1):
        coverage_f, fragments_f = coverage_paired_1_sense_1(
            bam, str(row['chromosome']),
            row['start'], row['end'], row['strand'], 'fg')
    elif (paired == 1) & (sense == 2):
        coverage_f, fragments_f = coverage_paired_1_sense_2(
            bam, str(row['chromosome']),
            row['start'], row['end'], row['strand'], 'fg')
    elif (paired == 2) & (sense == 1):
        coverage_f, fragments_f = coverage_paired_2_sense_1(
            bam, str(row['chromosome']),
            row['start'], row['end'], row['strand'], 'fg')
    elif (paired == 2) & (sense == 2):
        coverage_f = coverage_paired_2_sense_2_initial(
            bam,  str(row['chromosome']),
            row['start'], row['end'], row['strand'], 'fg', window_size)
    new_end = row['start']
    j = 0
    for i in np.arange(step, row['end'] - row['start'] - window_size, step):
        row1 = row.copy(deep=True)
        if coverage_f['coverage'].values[j + window_size: i].min() >= 10:
            continue
        new_start = new_end
        new_end = coverage_f.iloc[j + window_size: i].sort_values(
            by=['coverage'], ascending=True).index[0]
        j = new_end
        new_end = coverage_f.loc[new_end].coordinates
        row1.loc['start'] = new_start
        row1.loc['end'] = new_end
        row1.loc['startbg'] = new_start
        row1.loc['endbg'] = new_end
        rows.append(row1)
    row1 = row.copy(deep=True)
    new_start = new_end
    new_end = row['end']
    row1.loc['start'] = new_start
    row1.loc['end'] = new_end
    row1.loc['startbg'] = new_start
    row1.loc['endbg'] = new_end
    rows.append(row1)
    return rows


def second_EM_define_peaks(params):
    '''EM procedure to fit gaussians to the enriched coverage profiles'''
    # The parameters consist of:
    # 1. the enriched region
    enriched_region = params[0]
    # 1a. chrom where the enriched region is found
    chrom = str(enriched_region['chromosome'])
    # 1b. coordinates of start of enriched region
    region_start = int(enriched_region['start'])
    # 1c. coordinates of end of enriched region
    region_end = int(enriched_region['end'])
    # 1d. strand of enriched region
    strand = enriched_region['strand']
    # 1e. coordinates of start of bg enriched region
    region_start_b = int(enriched_region['startbg'])
    # 1f. coordinates of end of bg enriched region
    region_end_b = int(enriched_region['endbg'])
    del enriched_region
    # 2. foreground bam file (used for coverage profiles, using pysam)
    bam_f = params[1]
    # 3. the foreground bam file (used for coverage profiles, using pysam)
    bam_b = params[2]
    # 4. parameters of the previous EM step
    parameters = params[3]
    # 5. folder to output the figures of the coverage profiles
    plotfolder = params[4]
    # 6. fragment size
    fragment_size = float(params[5])
    # 7. cutoff for the log_l convergence
    # (affects the number of cycles during the EM)
    log_likelihood_cutoff = float(params[6])
    # 9. window size (refers to the foreground window size)
    window_size = params[7]
    # 10. sense mate of the foreground sample
    sense_f = int(params[8])
    # 11. sense mate of the background sample
    sense_b = int(params[9])
    # 12. length cutoff for the peak width
    length_cutoff = int(params[10])
    # 9a. windows that have an overlap with the enriched region
    # do this to decrease the number of rows because the windows file
    # contains too many rows
    paired_f = int(params[11])
    paired_b = int(params[12])
    background_type = params[13]
    verbose = params[14]
    # counter for plots at different stages of the fitting)
    # define the coverage profile of the enriched region
    # foreground
    coverage_f = pd.DataFrame()
    if (paired_f == 1) & (sense_f == 1):
        coverage_f, fragments_f = coverage_paired_1_sense_1(
            bam_f, chrom, region_start, region_end, strand, 'fg')
    elif (paired_f == 1) & (sense_f == 2):
        coverage_f, fragments_f = coverage_paired_1_sense_2(
            bam_f, chrom, region_start, region_end, strand, 'fg')
    elif (paired_f == 2) & (sense_f == 1):
        coverage_f, fragments_f = coverage_paired_2_sense_1(
            bam_f, chrom, region_start, region_end, strand, 'fg')
    elif (paired_f == 2) & (sense_f == 2):
        coverage_f, fragments_f = coverage_paired_2_sense_2(
            bam_f, chrom, region_start, region_end, strand, 'fg')
    else:
        pass
    # background
    coverage_b = pd.DataFrame()
    if (paired_b == 1) & (sense_b == 1):
        coverage_b, fragments_b, codes_b = coverage_paired_1_sense_1(
            bam_b, chrom, region_start_b, region_end_b, strand, 'bg', True)
    elif (paired_b == 1) & (sense_b == 2):
        coverage_b, fragments_b, codes_b = coverage_paired_1_sense_2(
            bam_b, chrom, region_start_b, region_end_b, strand, 'bg', True)
    elif (paired_b == 2) & (sense_b == 1):
        coverage_b, fragments_b, codes_b = coverage_paired_2_sense_1(
            bam_b, chrom, region_start_b, region_end_b, strand, 'bg', True)
    elif (paired_b == 2) & (sense_b == 2):
        coverage_b, fragments_b, codes_b = coverage_paired_2_sense_2(
            bam_b, chrom, region_start_b, region_end_b, strand, 'bg', True)
    else:
        pass
    coverage_b.columns = [
        'position_bg', 'coverage_bg', 'multimappers_bg',
        'total_bg', 'coordinates', 'crosslink_bg']
    # low_coverage = [ ]
    # counter = 0
    # while counter == 0:
    #     for index, row in coverage_f.iterrows():
    #         if row['coverage'] < 2:
    #             low_coverage.append(index)
    #         else:
    #             counter+=1
    # counter = 0
    # while counter == 0:
    #     for index in reversed(coverage_f.index):
    #         if coverage_f.loc[index, 'coverage'] < 2:
    #             low_coverage.append(index)
    #         else:
    #             counter+=1
    # coverage_f = coverage_f[~coverage_f.index.isin(low_coverage)]
    # coverage_f.reset_index(drop=True, inplace=True)

    enrichment_per_position = coverage_f.merge(
        coverage_b, how='outer', right_on='coordinates', left_on='coordinates')
    enrichment_per_position.fillna(0, inplace=True)

    if background_type == 'shifted_fg':
        fragments_b = fragments_b - fragments_f
        if fragments_b < 0:
            fragments_b = 0
    if verbose is True:
        enrichment_per_position.to_csv(
            os.path.join(
                plotfolder,
                f'{str(chrom)}_{str(region_start)}_{strand}_coverage.csv'),
            sep='\t',
            index=True,
            header=True)

    initials = []
    best_iteration = 0
    log_max = None
    iterations = pd.DataFrame()
    iteration_number = (2 * np.floor(
        len(coverage_f) / (float(fragment_size))))

    for count in np.arange(iteration_number):
        # 1. Gaussian initialisation
        gaussians = initialise_gaussians(coverage_f, fragment_size)
        
        if gaussians.empty:
            gaussians = pd.DataFrame(columns=[
                'chromosome', 'start', 'end', 'strand',
                'mi', 'ro', 'sigma', 'crosslink',
                'crosslink_weight', 'z_score', 'nj', 'mj'])
            return gaussians
        gaussians['iterations'] = count
        initials.append(gaussians.copy(deep=True))
        # 2. Expectation maximisation to optimise mi, sigma, ro
        gaussians = expectation_maximisation(
            gaussians, coverage_f, log_likelihood_cutoff)
        gaussians.dropna(axis=0, how='any', inplace=True)
        if gaussians.empty:
            continue
        gaussians['chromosome'] = chrom
        gaussians['start'] = region_start
        gaussians['end'] = region_end
        gaussians['strand'] = strand
        if log_max is None:
            log_max = gaussians['log_l'].values[0]
            iterations = gaussians
        else:
            if gaussians['log_l'].values[0] >= log_max:
                iterations = gaussians
                best_iteration = count
                log_max = gaussians['log_l'].values[0]
    # ------------------------------------------------------------------------
    # 3. Choose optimal iteration
    # ------------------------------------------------------------------------
    gaussians = iterations
    if gaussians.empty:
        gaussians = pd.DataFrame(columns=[
            'chromosome', 'start', 'end', 'strand',
            'mi', 'ro', 'sigma', 'crosslink',
            'crosslink_weight', 'z_score', 'nj', 'mj'])
        return gaussians
    # fitted_curves(
    #     gaussians, enrichment_per_position, chrom,
    #     region_start, 0, window_size, plotfolder, plot=True)

    for i in np.arange(len(gaussians)):
        gaussian_number = len(gaussians)
        fitted_curves(
            gaussians, enrichment_per_position, chrom,
            region_start, i, window_size, plotfolder, plot=False)
        gaussians = remove_redundant_peaks(gaussians)
        if (len(gaussians) == gaussian_number) | (len(gaussians) == 0):
            break
        gaussians = expectation_maximisation(
            gaussians, coverage_f, log_likelihood_cutoff)
        gaussians.dropna(axis=0, how='any', inplace=True)
        if gaussians.empty:
            break

    if gaussians.empty:
        gaussians = pd.DataFrame(columns=[
            'chromosome', 'start', 'end', 'strand',
            'mi', 'ro', 'sigma', 'crosslink',
            'crosslink_weight', 'z_score', 'nj', 'mj'])
        return gaussians
    initial_df = pd.concat(initials)

    # gaussians = merge_overlapping_peaks_part(gaussians)
    gaussians['chromosome'] = chrom
    gaussians['start'] = region_start
    gaussians['end'] = region_end
    gaussians['strand'] = strand
    initial_df = initial_df[initial_df['iterations'] == best_iteration]

    if verbose is True:
        # a) plot the initial random gaussians
        # counter += 1
        # fitted_curves(
        #     initial_df, enrichment_per_position, chrom,
        #     region_start, counter, window_size, plotfolder, plot=True)

        # b) plot the fitted curves versus the coverage profile
        counter = 0
        counter += 1
        fitted_curves(
            gaussians, enrichment_per_position, chrom,
            region_start, counter, window_size, plotfolder, plot=True)

    # ------------------------------------------------------------------------
    # 4. Renormalise coverage profiles, recalculate z-scores
    # ------------------------------------------------------------------------
    # make new z-scores comparable to those of the enriched regions step
    # calculate new z-score for each of the peaks
    # keep only the peaks that have a new z-score
    # above the cutoff of the enriched region detection step
    gaussians = renormalise_number_of_reads(
        gaussians, enrichment_per_position, codes_b,
        fragments_f, fragments_b)
    gaussians = z_scores_peaks(gaussians, parameters, length_cutoff)
    gaussians = gaussians[gaussians['z_score'] >= float(parameters['z_score'])]
    if gaussians.empty:
        gaussians = pd.DataFrame(columns=[
            'chromosome', 'start', 'end', 'strand',
            'mi', 'ro', 'sigma', 'crosslink',
            'crosslink_weight', 'z_score', 'nj', 'mj'])
        return gaussians

    gaussians = max_crosslink(gaussians, coverage_f, fragment_size)
    

    # gaussians = gaussians[gaussians['z_score'] >= float(parameters['z score'])]
    # If none of the peaks in an enriched region has a z-score
    # above the cutoff, the enriched region is ignored
    gaussians['chromosome'] = chrom
    gaussians['start'] = region_start
    gaussians['end'] = region_end
    gaussians['strand'] = strand
    gaussians = gaussians[[
        'chromosome', 'start', 'end', 'strand',
        'mi', 'ro', 'sigma', 'crosslink', 'crosslink_weight', 'height',
        'z_score', 'nj', 'mj']]
    return gaussians


def expectation_maximisation(gaussians, coverage_f, cutoff=0.0001):
    log_ls = []
    for em_cycle in np.arange(300):
        gaussians, log_l = EM_step_define_peaks(
            coverage_f,
            gaussians)
        gaussians['log_l'] = log_l
        if em_cycle < 10:
            log_ls.append(log_l)
            continue
        elif (abs((log_l - log_ls[-1]) / log_ls[-1]) < cutoff):
            break
        log_ls.append(log_l)
    if em_cycle > 90:
        print(em_cycle)
    return gaussians


def estimate_frequency(name, unique_fragments):
    unique_fragments_number = len(unique_fragments)
    if name == 'fg':
        if fg_frequencies:
            for i in unique_fragments:
                try:
                    unique_fragments_number = unique_fragments_number \
                        - 1 + fg_frequencies[i]
                except KeyError:
                    continue
    elif name == 'bg':
        if bg_frequencies:
            for i in unique_fragments:
                try:
                    unique_fragments_number = unique_fragments_number \
                        - 1 + bg_frequencies[i]
                except KeyError:
                    continue
    return unique_fragments_number


def weight_multimappers(name, unique_fragments):
    multimapper_number = 0
    total_number = len(unique_fragments)
    unique_fragments_number = len(unique_fragments)
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


def log_likelihood(enriched_region, gaussians):
    '''log_l calculations without fitting'''
    gaussians = gaussians.copy(deep=True)
    sigma_squared = gaussians['sigma'].pow(2)
    a = gaussians['ro'] / np.sqrt(2 * np.pi * sigma_squared)
    b = pd.concat(
        [enriched_region['position']] * len(gaussians),
        axis=1,
        ignore_index=True).subtract(gaussians['mi'], axis='columns')
    sigma_part = b.pow(2)
    division_part = 2 * (sigma_squared)
    b = np.exp(-sigma_part.divide(division_part, axis='columns'))
    aj_bj = b.multiply(a, axis='columns')
    Sum_j = aj_bj.sum(axis='columns')
    Sum_Ci_ro = (1 - gaussians['ro'].sum()) / len(enriched_region)
    if Sum_Ci_ro < 0:
        Sum_Ci_ro = 0
    Sum_j = Sum_j + Sum_Ci_ro
    log_likelihood = np.sum(
        np.log(
            Sum_j.replace(0, np.nan)
            ) * enriched_region['coverage']
        )
    return log_likelihood


def max_crosslink(gaussians, enriched_regions, fragment_size=80):
    fragment_size = fragment_size / 2
    for index, row in gaussians.iterrows():
        start = row['mi'] - (2 * row['sigma'])
        end = row['mi'] + (2 * row['sigma'])
        if start <= 0:
            start = 0
        if end >= len(enriched_regions):
            end = len(enriched_regions)
        if row['strand'] == '+':
            chunk_enriched_regions = enriched_regions.loc[
                start: end].sort_values(
                    by=['crosslink', 'position'],
                    ascending=[False, True])
        else:
            chunk_enriched_regions = enriched_regions.loc[
                start: end].sort_values(
                    by=['crosslink', 'position'],
                    ascending=[False, False])
        gaussians.loc[index, 'crosslink'] = \
            chunk_enriched_regions['position'].values[0]
        gaussians.loc[index, 'crosslink_weight'] = \
            chunk_enriched_regions['crosslink'].values[0]
    return gaussians


def coverage_paired_2_sense_2_initial(
        bamfile, chrom, start, end, strand, name, window_size, cutoff=10):
    '''Build coverage profile of a given enriched region'''
    window_size = int(window_size)
    sam = pysam.AlignmentFile(bamfile, "rb")
    counter = 0
    output = StringIO()
    csv_writer = writer(output)
    if start == 0:
        start = 1
    csv_writer.writerow(
        ['position', 'coverage', 'coordinates'])
    if strand == '+':
        for i in np.arange(start, end + 1):
            counter += 1
            if ((counter < window_size) | (counter > (end - start - window_size)) | (i % 5 != 0) ):
                csv_writer.writerow([counter, cutoff, i])
                continue
            reads_codes = set()
            for read in sam.fetch(chrom, i, i + 1):
                if len(reads_codes) >= cutoff:
                    break
                if (read.is_read2):
                    if (not read.is_reverse):
                        if (read.is_proper_pair):
                            if i in list(read.get_reference_positions()):
                                reads_codes.add(read.qname)
            csv_writer.writerow([counter, len(reads_codes), i])
    elif strand == '-':
        for i in np.arange(start, end + 1):
            counter += 1
            if ((counter < window_size) | (counter > (end - start - window_size)) | (i % 5 != 0) ):
                csv_writer.writerow([counter, cutoff, i])
                continue
            reads_codes = set()
            for read in sam.fetch(chrom, i-1, i):
                if len(reads_codes) >= cutoff:
                    break
                if (read.is_read2):
                    if (read.is_reverse):
                        if (read.is_proper_pair):
                            if i in list(read.get_reference_positions()):
                                reads_codes.add(read.qname)
            csv_writer.writerow([counter, len(reads_codes), i])
    sam.close()
    output.seek(0)
    df_output = pd.read_csv(output)
    return df_output


def coverage_paired_2_sense_2(
        bamfile, chrom, start, end, strand, name, reads=False):
    '''Build coverage profile of a given enriched region'''
    peaks = {}
    crosslinks = {}
    region_counts = set()
    for a in np.arange(start, end+1):
        peaks[a] = set()
        crosslinks[a] = set()
    sam = pysam.AlignmentFile(bamfile, "rb")
    region_counts = set()
    if strand == '+':
        for read in sam.fetch(chrom, start, end + 1):
            if (read.is_read2):
                if (not read.is_reverse):
                    region_counts.add(read.qname)
                    for i in list(read.get_reference_positions()):
                        try:
                            peaks[i].add(read.qname)
                        except KeyError:
                            continue
                    try:
                        crosslinks[read.reference_start].add(read.qname)
                    except KeyError:
                        pass
    elif strand == '-':
        for read in sam.fetch(chrom, start, end + 1):
            if (read.is_read2):
                if (read.is_reverse):
                    region_counts.add(read.qname)
                    for i in list(read.get_reference_positions()):
                        try:
                            peaks[i].add(read.qname)
                        except KeyError:
                            continue
                    try:
                        crosslinks[read.reference_end].add(read.qname)
                    except KeyError:
                        pass
    sam.close()
    fragments = estimate_frequency(name, region_counts)
    rows = []
    for key in peaks:
        if len(peaks[key]) == 0:
            rows.append([0, 0, 0, key, 0])
        else:
            weighted, multim, total = weight_multimappers(
                name, peaks[key]
                )
            if len(crosslinks[key]) == 0:
                rows.append([weighted, multim, total, key, 0])
            else:
                weightedc, multimc, totalc = weight_multimappers(
                    name, crosslinks[key]
                    )
                rows.append([weighted, multim, total, key, weightedc])
    df_output = pd.DataFrame(rows, columns=[
        'coverage', 'multimappers',
        'total', 'coordinates', 'crosslink'])
    df_output['position'] = np.arange(1, end - start + 2)
    df_output.reset_index(drop=True, inplace=True)
    df_output = df_output[[
        'position', 'coverage', 'multimappers',
        'total', 'coordinates', 'crosslink']]
    df_output = df_output.astype(
        {"position": int, "coverage": float, "multimappers": int,
         "total": int, "coordinates": int, "crosslink": float})
    if reads:
        return df_output, fragments, peaks
    else:
        return df_output, fragments


def coverage_paired_2_sense_1(
        bamfile, chrom, start, end, strand, name, reads=False):
    '''Build coverage profile of a given enriched region'''
    sam = pysam.AlignmentFile(bamfile, "rb")
    counter = 0
    output = StringIO()
    csv_writer = writer(output)
    codes = {}
    if start == 0:
        start = 1
    csv_writer.writerow(
        ['position', 'coverage', 'multimappers',
         'total', 'coordinates', 'crosslink'])
    region_counts = set()
    if strand == '+':
        for i in np.arange(start, end + 1):
            counter += 1
            reads_codes = set()
            crosslink = set()
            for read in sam.fetch(chrom, i, i + 1):
                if (read.is_read1):
                    if (not read.is_reverse):
                        if (read.is_proper_pair):
                            if i in list(read.get_reference_positions()):
                                reads_codes.add(read.qname)
                                region_counts.add(read.qname)
                            if read.reference_start == i:
                                crosslink.add(read.qname)
            weighted, multim, total = weight_multimappers(
                name, reads_codes
                )
            # estimate frequency of the crosslink aggregate reads
            # (weight in case of multimappers)
            weighted_c, multim_c, total_c = weight_multimappers(
                name, crosslink
                )
            if reads:
                codes[i] = list(reads_codes)
            csv_writer.writerow(
                [counter, weighted, multim, total, i, weighted_c]
                )
    elif strand == '-':
        for i in np.arange(start, end + 1):
            counter += 1
            reads_codes = set()
            crosslink = set()
            for read in sam.fetch(chrom, i-1, i):
                if (read.is_read1):
                    if (read.is_reverse):
                        if (read.is_proper_pair):
                            if i in list(read.get_reference_positions()):
                                reads_codes.add(read.qname)
                                region_counts.add(read.qname)
                            if read.reference_end == i:
                                crosslink.add(read.qname)
            weighted, multim, total = weight_multimappers(
                name, reads_codes
                )
            weighted_c, multim_c, total_c = weight_multimappers(
                name, crosslink
                )
            if reads:
                codes[i] = list(reads_codes)
            csv_writer.writerow(
                [counter, weighted, multim, total, i, weighted_c]
                )
    sam.close()
    output.seek(0)
    fragments = estimate_frequency(name, region_counts)
    df_output = pd.read_csv(output)
    if reads:
        return df_output, fragments, codes
    else:
        return df_output, fragments


def coverage_paired_1_sense_1(
        bamfile, chrom, start, end, strand, name, reads=False):
    '''Build coverage profile of a given enriched region'''
    sam = pysam.AlignmentFile(bamfile, "rb")
    counter = 0
    output = StringIO()
    csv_writer = writer(output)
    codes = {}
    if start == 0:
        start = 1
    csv_writer.writerow(
        ['position', 'coverage', 'multimappers',
         'total', 'coordinates', 'crosslink'])
    region_counts = set()
    if strand == '+':
        for i in np.arange(start, end + 1):
            counter += 1
            reads_codes = set()
            crosslink = set()
            for read in sam.fetch(chrom, i, i + 1):
                if not read.is_reverse:
                    if i in list(read.get_reference_positions()):
                        reads_codes.add(read.qname)
                        region_counts.add(read.qname)
                    if read.reference_start == i:
                        crosslink.add(read.qname)
            weighted, multim, total = weight_multimappers(
                name, reads_codes
                )
            # estimate frequency of the crosslink aggregate reads
            # (weight in case of multimappers)
            weighted_c, multim_c, total_c = weight_multimappers(
                name, crosslink
                )
            if reads:
                codes[i] = list(reads_codes)
            csv_writer.writerow(
                [counter, weighted, multim, total, i, weighted_c]
                )
    elif strand == '-':
        for i in np.arange(start, end + 1):
            counter += 1
            reads_codes = set()
            crosslink = set()
            for read in sam.fetch(chrom, i-1, i):
                if read.is_reverse:
                    if i in list(read.get_reference_positions()):
                        reads_codes.add(read.qname)
                        region_counts.add(read.qname)
                    if read.reference_end == i:
                        crosslink.add(read.qname)
            weighted, multim, total = weight_multimappers(
                name, reads_codes
                )
            weighted_c, multim_c, total_c = weight_multimappers(
                name, crosslink
                )
            if reads:
                codes[i] = list(reads_codes)
            csv_writer.writerow(
                [counter, weighted, multim, total, i, weighted_c]
                )
    sam.close()
    output.seek(0)
    fragments = estimate_frequency(name, region_counts)
    df_output = pd.read_csv(output)
    if reads:
        return df_output, fragments, codes
    else:
        return df_output, fragments


def coverage_paired_1_sense_2(
        bamfile, chrom, start, end, strand, name, reads=False):
    '''Build coverage profile of a given enriched region'''
    sam = pysam.AlignmentFile(bamfile, "rb")
    counter = 0
    output = StringIO()
    csv_writer = writer(output)
    codes = {}
    if start == 0:
        start = 1
    csv_writer.writerow(
        ['position', 'coverage', 'multimappers',
         'total', 'coordinates', 'crosslink'])
    region_counts = set()
    if strand == '+':
        for i in np.arange(start, end + 1):
            counter += 1
            reads_codes = set()
            crosslink = set()
            for read in sam.fetch(chrom, i, i + 1):
                if read.is_reverse:
                    if i in list(read.get_reference_positions()):
                        reads_codes.add(read.qname)
                        region_counts.add(read.qname)
                    if read.reference_start == i:
                        crosslink.add(read.qname)
            weighted, multim, total = weight_multimappers(
                name, reads_codes
                )
            # estimate frequency of the crosslink aggregate reads
            # (weight in case of multimappers)
            weighted_c, multim_c, total_c = weight_multimappers(
                name, crosslink
                )
            if reads:
                codes[i] = list(reads_codes)
            csv_writer.writerow(
                [counter, weighted, multim, total, i, weighted_c]
                )
    elif strand == '-':
        for i in np.arange(start, end + 1):
            counter += 1
            reads_codes = set()
            crosslink = set()
            for read in sam.fetch(chrom, i-1, i):
                if not read.is_reverse:
                    if i in list(read.get_reference_positions()):
                        reads_codes.add(read.qname)
                        region_counts.add(read.qname)
                    if read.reference_end == i:
                        crosslink.add(read.qname)
            weighted, multim, total = weight_multimappers(
                name, reads_codes
                )
            weighted_c, multim_c, total_c = weight_multimappers(
                name, crosslink
                )
            if reads:
                codes[i] = list(reads_codes)
            csv_writer.writerow(
                [counter, weighted, multim, total, i, weighted_c]
                )
    sam.close()
    output.seek(0)
    fragments = estimate_frequency(name, region_counts)
    df_output = pd.read_csv(output)
    if reads:
        return df_output, fragments, codes
    else:
        return df_output, fragments


def initialise_gaussians(enriched_region, fragment_size=80):
    '''initialise gaussians and distribute them evenly to the window'''
    # number of initial gaussians,
    # divide the length of the enriched region by the fragment size
    covered = enriched_region[enriched_region['coverage'] > 2]
    # initalise with at least 2 peaks
    j_number = np.max([np.floor(len(covered) / float(fragment_size)), 2])
    # gaussians = pd.DataFrame(columns=['sigma', 'mi', 'ro'])
    if len(covered) < (fragment_size / 2):
        sys.stderr.write(
            "\nvery small coverage\n")
        gaussians = pd.DataFrame()
        return gaussians
    all_positions = covered['position'].values
    peaks = []
    for indexj in np.arange(j_number):
        # the total length is equally divided between
        # the initialised gaussians
        gaussian_sigma = len(covered) / (4 * j_number)
        # the center of the gaussian is a random position in the region
        gaussian_mi = float(np.random.choice(all_positions))
        all_positions = list(set(all_positions).difference(
            set(np.arange(
                (gaussian_mi - gaussian_sigma),
                (gaussian_mi + gaussian_sigma) + 1
            )))
        )
        peaks.append([gaussian_sigma, gaussian_mi, 1 / (4 * j_number)])
    gaussians = pd.DataFrame(peaks, columns=["sigma", "mi", "ro"])
    gaussians = gaussians.astype({"sigma": float, "mi": float, "ro": float})
    return gaussians


def remove_redundant_peaks(gaussians1):
    gaussians = gaussians1.copy(deep=True)
    gaussians['start'] = gaussians['mi'] - 2 * gaussians['sigma']
    gaussians['end'] = gaussians['mi'] + 2 * gaussians['sigma']
    gaussians.sort_values(by='sigma', ascending=False, inplace=True)
    gaussians = gaussians.astype({"start": np.int, "end": np.int})
    for index, row in gaussians.iterrows():
        for index1, row1 in gaussians.iterrows():
            if index == index1:
                continue
            overlap = range(
                max(row['start'], row1['start']),
                min(row['end'], row1['end']) + 1)
            if len(overlap) == 0:
                continue
            result1 = float(len(overlap) / (row['end'] - row['start']))
            result2 = float(len(overlap) / (row1['end'] - row1['start']))
            if (result1 >= 0.5) | (result2 >= 0.5):
                if row1['height'] > row['height']:
                    gaussians1.drop([index], inplace=True)
                    return gaussians1
                elif row1['height'] <= row['height']:
                    gaussians1.drop([index1], inplace=True)
                    return gaussians1
    return gaussians1


def EM_step_define_peaks(enriched_region, gaussians):
    ''' Expectation maximisation for the peak fitting'''
    gaussians = gaussians.copy(deep=True)
    sigma_squared = np.array(gaussians['sigma']) ** 2
    a = [gaussians['ro'].values / np.sqrt(
        2 * np.pi * sigma_squared)] * len(enriched_region['position'])
    position = np.array([
        enriched_region['position'].values] * len(gaussians)).T
    coverage = np.array([enriched_region['coverage']] * len(gaussians)).T
    each = [np.array(gaussians['mi'])] * len(enriched_region['position'])
    b = position - each
    sigma_part = b ** 2
    sigma_division = [2 * sigma_squared] * len(enriched_region['position'])
    b = np.exp(- sigma_part / sigma_division)
    aj_bj = a * b
    Sum_j = aj_bj.sum(axis=1)
    Sum_Ci = enriched_region['coverage'].sum()
    Sum_Ci_ro = (1 - np.sum(gaussians['ro'].values)) / len(enriched_region)
    if Sum_Ci_ro < 0:
        Sum_Ci_ro = 0
    Sum_j = np.array([Sum_j + Sum_Ci_ro] * len(gaussians)).T
    PG_ji = aj_bj / Sum_j
    PG_jiC = PG_ji * coverage
    PG_jiC_sum = PG_jiC.sum(axis=0)
    iPG_jiC = PG_jiC * position
    np.seterr(all='ignore')
    gaussians['ro'] = PG_jiC_sum / Sum_Ci
    gaussians['mi'] = iPG_jiC.sum(axis=0) / PG_jiC_sum
    gaussians['sigma'] = np.sqrt(
        (sigma_part * PG_jiC).sum(axis=0) / PG_jiC_sum)
    log_likelihood = np.sum(np.log(
        np.nan_to_num(Sum_j[:, 0])) * enriched_region['coverage'].values)
    return gaussians, log_likelihood


def slicedict(d, start, end):
    new_codes = set()
    {new_codes.update(v) for k, v in d.items() if (k >= start) & (k <= end)}
    new_codes = list(new_codes)
    return new_codes


def renormalise_number_of_reads(
        gaussians, enriched_region, coverage_bg,
        fg_fragm, bg_fragm):
    '''estimate number of unique fragments within the peak'''
    gaussians = gaussians.copy(deep=True)
    gaussians['nj'] = 0
    gaussians['mj'] = 0
    for indexj, rowj in gaussians.iterrows():
        start = int(rowj['start'] + rowj['mi'] - rowj['sigma'])
        end = int(rowj['start'] + rowj['mi'] + rowj['sigma'])
        gaussians.loc[indexj, 'nj'] = fg_fragm * (
            rowj['ro'] + (
                (2 * rowj['sigma'] / len(enriched_region)) *
                float(1.0 - gaussians['ro'].sum())
            )) + 1.5
        fragments_b = slicedict(coverage_bg, start, end)
        fragments_b, multim_b, total_b = weight_multimappers('bg', fragments_b)
        gaussians.loc[indexj, 'mj'] = fragments_b + 1.5
    return gaussians


def z_scores_peaks(gaussians, parameters, length_cutoff):
    '''recalculate z-scores for each of the peaks'''
    gaussians['z_score'] = 0
    M = float(parameters['total_b'])
    N = float(parameters['total_f'])
    sigma = float(parameters['sigma'])
    mi = float(parameters['mi'])
    for indexj, rowj in gaussians.iterrows():
        nj = float(rowj['nj'])
        mj = float(rowj['mj'])
        zscore_j = (
            np.log(nj / N) - np.log(mj / M) - mi
        ) / np.sqrt(
            (2 * (sigma ** 2)) + (1 / nj) + (1 / mj))
        gaussians.loc[indexj, 'z_score'] = zscore_j
    return gaussians


def gaussian_fit(x, a, b, c):
    ''' calculate the gaussian fitted to each peak'''
    y = (a / (c * np.sqrt(2 * np.pi))) * np.exp(
        - ((x - b) ** 2) / (2 * (c ** 2)))
    return y


def fitted_curves(gaussians, enriched_regions, chromosome, start,
                  counter, window_size, plotfolder, plot=True):
    '''draw the coverage profiles along with the gaussians'''
    enriched_regions.dropna(inplace=True)
    window_size = int(window_size)
    enriched_regions['coverage'].fillna(0, inplace=True)
    if plot is True:
        fig, ax = plt.subplots(figsize=(len(enriched_regions) * 0.02, 3))
    for index, row in gaussians.iterrows():
        # width of gaussian
        sigma = row['sigma']
        # mean position of the gaussian
        mi = enriched_regions['position'].min() + row['mi']
        # coverage
        ro = row['ro']
        ro_tot = int(ro * enriched_regions['coverage'].sum())
        if ro_tot < 1:
            ro_tot = 1
        # set sampling space for gaussian (x)
        sampling = np.linspace(
            int(mi - (3 * sigma)),
            int(mi + (3 * sigma)),
            2 * ro_tot)
        # fit the gaussians (y)
        fitted_gaussian = gaussian_fit(
            sampling,
            ro_tot,
            mi,
            sigma)
        # plot the gaussian
        gaussians.loc[index, 'height'] = fitted_gaussian.max()
        if plot is True:
            plt.plot(sampling, fitted_gaussian, 'k', color='red', linewidth=3)
    if plot is True:
        sns.lineplot(
            x=enriched_regions['position'],
            y=enriched_regions['coverage'],
            ax=ax,
            color='black',
            lw=3)
        sns.lineplot(
            x=enriched_regions['position'],
            y=enriched_regions['coverage_bg'],
            ax=ax,
            color='grey',
            lw=3)
        sns.lineplot(
            x=enriched_regions['position'],
            y=enriched_regions['crosslink'],
            ax=ax,
            color='green',
            lw=3)
        ax.xaxis.set_major_locator(mticker.MultipleLocator(100))
        ax.xaxis.set_minor_locator(mticker.MultipleLocator(10))
        plt.xlim(enriched_regions['position'].min(), enriched_regions['position'].max())
        plt.ylim(0, np.amax([
            np.amax(fitted_gaussian),
            np.amax(enriched_regions['coverage'].values),
            np.amax(enriched_regions['coverage'].values),
            10]) + 5)
        fig.savefig(
            os.path.join(
                plotfolder,
                f'{str(chromosome)}_{str(start)}_{str(counter)}_peaks.png'
                )
            )
        plt.close(fig)
        plt.clf()
    return gaussians


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

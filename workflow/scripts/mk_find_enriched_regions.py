# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# import needed (external) modules
# ----------------------------------------------------------------------------------------------------------------------
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import numpy as np
import pandas as pd
from io import StringIO
from csv import writer
import matplotlib
import random
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# Main function
# ----------------------------------------------------------------------------------------------------------------------

def main():
    """ Find enriched regions in comparison to a background sample.
        The enrichment is the outcome of significant binding of the RBP
        used as target for the eCLIP experiment."""

    __doc__ = "Obtain sliding windows of read coverage from bam file."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--windows_table",
        dest="windows_table",
        help="foreground and background reads per window",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--out",
        dest="out",
        help="output table file containing the enriched regions",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--out_parameters",
        dest="out_parameters",
        help="output table file containing the fitted parameters",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--all_regions",
        dest="all_regions",
        help="all regions containing at least n reads",
        required=True)

    parser.add_argument(
        "--out_folder",
        dest="out_folder",
        help="destination folder for the plots ",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--log_likelihood_cutoff",
        dest="log_likelihood_cutoff",
        help=" log likelihood threshold value to stop the EM",
        required=False,
        default="0.01")

    parser.add_argument(
        "--fdr_cutoff",
        dest="fdr_cutoff",
        help=" fdr threshold for choosing z_score cutoff",
        required=False,
        default="1")

    # __________________________________________________________________________________________________________________
    # ------------------------------------------------------------------------------------------------------------------
    # get the arguments
    # ------------------------------------------------------------------------------------------------------------------
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    sys.stdout.write('Starting the detection of enriched windows...\n')

    windows = pd.read_csv(
        options.windows_table,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')
    windows.index = np.arange(len(windows))
    windows['Reads_f'] = windows['Reads_f'].astype('float64')
    windows['Reads_b'] = windows['Reads_b'].astype('float64')
    windows = windows[(windows['Reads_f'] >= 1)]
    windows['Reads_f'] = windows['Reads_f'] + 1
    windows['Reads_b'] = windows['Reads_b'] + 1

    n, N = get_fg_reads(windows)
    m, M = get_bg_reads(windows)

    fdr_cutoff = float(options.fdr_cutoff)
    log_likelihood_cutoff = float(options.log_likelihood_cutoff)
    plot_folder = options.out_folder

    table = pd.DataFrame()
    sys.stdout.write('Step1. Starting the expectation maximisation step...\n')
    sys.stdout.write('sigma, ro, mi, log_likelihood\n')
    counter = 0
    log_likelihood = ''
    sigma = random.uniform(0.1, 0.9)
    ro = random.uniform(0.8, 0.9)
    mi = random.uniform(-3, 3)
    sys.stdout.write(f'{sigma}, {ro}, {mi}\n')
    sys.stdout.flush()
    while True:
        try:
            counter += 1
            if counter > 50:
                sys.stdout.write('No convergence of the EM could be achieved.')
                sys.exit(1)
                # break
            sigma, ro, mi, log_likelihood = EM_step_enriched_regions(
                windows,
                sigma=sigma,
                ro=ro,
                mi=mi)
        except ValueError:
            windows = windows[(windows['Reads_f'] >= 3)]
            sigma = random.uniform(0.1, 0.9)
            ro = random.uniform(0.7, 0.98)
            mi = random.uniform(-3, 3)
            continue
        sys.stdout.write(f'{sigma}, {ro}, {mi}, {log_likelihood}\n')
        sys.stdout.flush()
        table.at[counter, 'sigma'] = sigma
        table.at[counter, 'ro'] = ro
        table.at[counter, 'mi'] = mi
        table.at[counter, 'log_likelihood'] = log_likelihood
        if counter < 20:
            continue
        elif abs(log_likelihood - table.iloc[:-1]['log_likelihood'].max()) < log_likelihood_cutoff:
            break
    sys.stdout.write('Step1. Expectation maximisation step finished.\n')
    sys.stdout.flush()

    table[['mi', 'ro', 'sigma']].plot()
    plt.savefig(
        os.path.join(plot_folder, 'expectation_maximisation.png'),
        transparent=True)
    plt.clf()

    sys.stdout.write('plotted em values\n')
    sys.stdout.flush()

    # calculate the z_score for each of the windows
    windows = z_score(windows, mi=mi, sigma=sigma)

    # calculate the false discovery rate(FDR) and calculate zscore cutoff
    z_score_ctf, windows_total = FDR(
        windows, sigma=sigma, ro=ro, mi=mi, fdr_cutoff=fdr_cutoff)

    windows_total.to_csv(
        os.path.join(options.all_regions),
        sep='\t',
        index=True,
        header=True)

    windows_selected = windows_total[windows_total['z_score'] >= 2]
    # merge the final windows that show enrichment and have overlapping regions
    sys.stdout.write('start merging.....\n')
    sys.stdout.flush()

    enriched_regions = merge_overlapping_regions(windows_selected)
    sys.stdout.write('finished merging.....\n')
    sys.stdout.flush()

    # create a parameters table whith the fitted values
    df_parameters = pd.DataFrame()
    df_parameters.loc[0, 'sigma'] = sigma
    df_parameters.loc[0, 'ro'] = ro
    df_parameters.loc[0, 'mi'] = mi
    df_parameters.loc[0, 'z_score'] = z_score_ctf
    df_parameters.loc[0, 'fdr cutoff'] = options.fdr_cutoff
    df_parameters.loc[0, 'log likelihood'] = log_likelihood
    df_parameters.loc[0, 'total_b'] = M
    df_parameters.loc[0, 'total_f'] = N

    df_parameters.to_csv(
        options.out_parameters, sep='\t', index=True, header=True)
    enriched_regions.to_csv(
        options.out, sep='\t', index=True, header=True)

    sys.stdout.write('Detection of enriched windows finshed.\n')
    sys.stdout.flush()

    # create some plots for these results
    sys.stdout.write('Creating some additional plots.\n')
    sys.stdout.flush()
    try:
        plots(windows, z_score_ctf, plot_folder)
    except:
        sys.stderr.write(
            'Finished. Something went wrong with some of the plots.\n')
    return

# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# Functions
# ----------------------------------------------------------------------------------------------------------------------

def get_fg_reads(windows):
    '''calculate total number of reads for foreground'''
    n = windows['Reads_f'].values
    N = float(n.sum())
    return n, N


def get_bg_reads(windows):
    '''calculate total number of reads for background'''
    m = windows['Reads_b'].values
    M = float(m.sum())
    return m, M


def get_W(windows):
    '''calculate W: range of values for the diff in log read densities'''
    n, N = get_fg_reads(windows)
    m, M = get_bg_reads(windows)
    f_nm = np.log(n) + np.log(M) - np.log(m) - np.log(N)
    W = float(f_nm.max() - f_nm.min())
    return W


def calculate_a(windows, sigma=None):
    n, N = get_fg_reads(windows)
    m, M = get_bg_reads(windows)
    a = (2 * (sigma ** 2)) + (1 / m) + (1 / n)
    return a


def calculate_b(windows, mi=None):
    n, N = get_fg_reads(windows)
    m, M = get_bg_reads(windows)
    b = np.log(n) + np.log(M) - np.log(m) - np.log(N) - mi
    return b


def calculate_c(windows):
    n, N = get_fg_reads(windows)
    m, M = get_bg_reads(windows)
    c = np.log(n) + np.log(M) - np.log(m) - np.log(N)
    return c


def get_Pnimi(windows, sigma=None, mi=None, ro=None):
    n, N = get_fg_reads(windows)
    m, M = get_bg_reads(windows)
    nominator = calculate_b(windows, mi=mi)
    denominator = calculate_a(windows, sigma=sigma)
    P_nimi = np.exp(
        np.log(
            1 / np.sqrt(2 * np.pi * denominator)) + (
            - (nominator ** 2) / (2 * denominator)))
    return P_nimi


def get_p_bgi(windows, sigma=None, mi=None, ro=None):
    P_nimi = get_Pnimi(windows, sigma=sigma, mi=mi, ro=ro)
    W = get_W(windows)
    p_bgi = (ro * P_nimi) / ((ro * P_nimi) + ((1 - ro) * (1 / W)))
    return p_bgi


def get_log_likelihood(windows, sigma=None, mi=None, ro=None):
    P_nimi = get_Pnimi(windows, sigma=sigma, mi=mi, ro=ro)
    W = get_W(windows)
    log_likelihood = sum(np.log((ro * P_nimi) + ((1 - ro) / W)))
    return log_likelihood


def EM_step_enriched_regions(windows, sigma=None, ro=None, mi=None):
    '''Expectation maximisation step to detect enriched
    regions given a background.'''
    p_bgi = get_p_bgi(windows, sigma=sigma, mi=mi, ro=ro)
    a = calculate_a(windows, sigma=sigma)
    c = calculate_c(windows)
    L = float(len(windows))
    mi_new = np.sum(p_bgi * (c / a)) / np.sum(p_bgi / a)
    ro_new = np.sum(p_bgi) / L
    try:
        sigma_new = approximate_sigma(windows, ro, mi)
    except ValueError:
        raise ValueError("derivative of sigma failed!\n")
    log_likelihood = get_log_likelihood(windows, sigma=sigma_new, mi=mi_new, ro=ro_new)
    return sigma_new, ro_new, mi_new, log_likelihood


def approximate_sigma(windows, ro=None, mi=None, s_start=10 ** (-8), cutoff=10 ** (-6)):
    '''Binary search to converge to a sigma value.'''
    s_end = s_start
    derivative_s_start = find_sigma_boundaries(windows, s=s_start, ro=ro, mi=mi)
    derivative_s_end = find_sigma_boundaries(windows, s=s_end, ro=ro, mi=mi)
    count = 0
    while derivative_s_start * derivative_s_end > 0:
        count += 1
        derivative_s_end = find_sigma_boundaries(
            windows, s=s_end, ro=ro, mi=mi)
        derivative_s_start = find_sigma_boundaries(
            windows, s=s_start, ro=ro, mi=mi)
        s_end = s_end * 10
        s_start = s_start / 10
        if count > 30:
            raise ValueError("derivative of sigma has no positive values! recalculate\n")
    # iteratively cut distance in 2
    while abs(s_end - s_start) > cutoff:
        s_new = (s_start + s_end) / 2.0
        derivative_s = find_sigma_boundaries(windows, s=s_new, ro=ro, mi=mi)
        if derivative_s_start * derivative_s < 0:
            s_end = s_new
            derivative_s_end = derivative_s
        elif derivative_s_end * derivative_s < 0:
            s_start = s_new
            derivative_s_start = derivative_s
    sigma_approximated = np.sqrt(s_new / 2)
    return sigma_approximated


def find_sigma_boundaries(windows, s=None, ro=None, mi=None):
    '''Given an s the derivative of s is calculated'''
    sigma = np.sqrt(s / 2)
    p_bgi = get_p_bgi(windows, sigma=sigma, mi=mi, ro=ro)
    b = calculate_b(windows, mi=mi)
    a = calculate_a(windows, sigma=sigma)
    derivative_s = np.sum(p_bgi * (
        ((b ** 2) / (2 * (a ** 2))) - (1.0 / (2 * a))))
    return derivative_s


def z_score(windows, mi=None, sigma=None):
    '''Z-score calculation for each of the windows.'''
    sys.stdout.write(
        'Step2. Calculate the z_score values for each of the windows...\n')
    sys.stdout.flush()
    b = calculate_b(windows, mi=mi)
    a = calculate_a(windows, sigma=sigma)
    windows['z_score'] = b / np.sqrt(a)
    sys.stdout.write('Step2. Z_score values successfully calculated...\n')
    sys.stdout.flush()
    return windows


def FDR(windows, sigma=None, ro=None, mi=None, fdr_cutoff=0.1):
    '''Probability for a given window to be a false positive.'''
    sys.stdout.write(
        'Step3. Calculate the z_score cutoff for given FDR cutoff...\n')
    sys.stdout.flush()

    windows.sort_values(by=['z_score'], ascending=[False], inplace=True)

    windows.index = np.arange(len(windows))
    windows['P_nimi'] = get_Pnimi(windows, sigma=sigma, mi=mi, ro=ro)
    windows['p_bgi'] = get_p_bgi(windows, sigma=sigma, mi=mi, ro=ro)
    windows['cumsum_p_bgi'] = windows['p_bgi'].cumsum()
    windows['T'] = np.arange(1, len(windows) + 1)
    windows['FDR'] = windows['cumsum_p_bgi'] / windows['T']

    sys.stdout.write('Step3. FDR chosen:' + str(fdr_cutoff) + '\n')
    sys.stdout.flush()

    windows_selected = windows[windows['FDR'] <= fdr_cutoff]
    if windows_selected.empty:
        sys.stderr.write(
            'no peaks detected to have high enough z score!\n')
        sys.exit(1)
    else:
        z_score_cutoff = windows_selected['z_score'].min()

    sys.stdout.write('Step3. Z_score cutoff chosen:' + str(z_score_cutoff) + '\n')
    sys.stdout.flush()
    return z_score_cutoff, windows


def merge_overlapping_regions(windows):
    '''merge overlapping regions in the genome that are considered enriched'''
    windows = windows[['chromosome', 'start', 'end', 'strand', 'startbg', 'endbg']].copy(deep=True)
    windows['start'] = windows['start'].astype('int')
    windows['end'] = windows['end'].astype('int')
    windows['startbg'] = windows['startbg'].astype('int')
    windows['endbg'] = windows['endbg'].astype('int')
    windows['chromosome'] = windows['chromosome'].astype('str')
    output = StringIO()
    csv_writer = writer(output)
    csv_writer.writerow(['chromosome', 'start', 'end', 'strand', 'startbg', 'endbg'])
    for strand in ['+', '-']:
        for chromosome in list(set(windows['chromosome'].values)):
            delete = []
            windows1 = windows[(windows['chromosome'] == chromosome) & (windows['strand'] == strand)].copy(deep=True)
            windows1.sort_values(by=['start'], ascending=[True], inplace=True)
            windows1.set_index(np.arange(len(windows1)), inplace=True)
            for index in windows1.index:
                if index in delete:
                    continue
                chrom = windows1.loc[index, 'chromosome']
                strand = windows1.loc[index, 'strand']
                start = windows1.loc[index, 'start']
                end = windows1.loc[index, 'end']
                startbg = windows1.loc[index, 'startbg']
                endbg = windows1.loc[index, 'endbg']
                chr_range = set(list(np.arange(start, end + 1)))
                success = 0
                for index2 in windows1.index:
                    if index2 <= index:
                        continue
                    elif index2 in delete:
                        continue
                    start2 = windows1.loc[index2, 'start']
                    end2 = windows1.loc[index2, 'end']
                    start2bg = windows1.loc[index2, 'startbg']
                    end2bg = windows1.loc[index2, 'endbg']
                    if start2 in chr_range:
                        success += 1
                        start = np.min([start, start2, end, end2])
                        end = np.max([start, start2, end, end2])
                        startbg = np.min([startbg, start2bg, endbg, end2bg])
                        endbg = np.max([startbg, start2bg, endbg, end2bg])
                        chr_range = set(list(np.arange(start, end + 1)))
                        delete.append(index2)
                peak = pd.Series([
                    chrom,
                    start,
                    end,
                    strand,
                    startbg,
                    endbg])
                csv_writer.writerow(peak)
    output.seek(0)
    windows = pd.read_csv(output)
    return windows


def plots(windows, cutoff, plot_folder):
    '''Various plots to visualise the results.'''

    sys.stdout.write('Step4. Creating plots of the results...\n')
    sys.stdout.flush()

    df_plots = windows.copy(deep=True)
    df_plots.index = np.arange(len(df_plots))

    m, std = norm.fit(df_plots['z_score'][(df_plots['z_score'] < 3.5) & (
        df_plots['z_score'] > -3.5)])
    plt.hist(np.array(
        df_plots['z_score']),
        bins=200,
        normed=True,
        histtype='step',
        log=True,
        lw=0.5
    )
    x = np.random.normal(m, std, len(df_plots['z_score'][(df_plots['z_score'] < 3.5) & (
        df_plots['z_score'] > -3.5)]))
    plt.hist(
        x,
        bins=50,
        normed=True,
        histtype='step',
        log=True,
        lw=1,
        color='red')
    plt.xlim([-10, 10])
    plt.title('z score distribution')
    plt.savefig(os.path.join(plot_folder, 'z_score.png'), transparent=True)
    plt.clf()


    plt.hist(
        np.array(df_plots['z_score']),
        bins=100000,
        normed=False,
        histtype='step',
        log=True,
        cumulative=-1,
        lw=3.0)
    plt.title('reverse cumulative distribution of z scores')
    plt.savefig(
        os.path.join(plot_folder, 'reverse_cumulative_z_score.png'),
        transparent=True)
    plt.clf()
    sys.stdout.write('Step4. Successfull plot creation.\n')
    sys.stdout.flush()
    return


# __________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(1)

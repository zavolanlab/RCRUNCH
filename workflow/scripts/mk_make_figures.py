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
#   This script, utilises the results obtained at the various steps of the
#   analysis and produces various figures.
# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
from argparse import ArgumentParser, RawTextHelpFormatter
import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser, RawTextHelpFormatter, FileType
import sys
from scipy.stats import norm
import os
import matplotlib.ticker as plticker
import glob
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D

sns.set_style(style='white')
sns.set_context("poster")
# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():
    """ Main function """

    __doc__ = "Obtain figures from the various steps of the analysis"

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--enriched_regions",
        dest="enriched_regions",
        help="file with the enrichment analysis values \
        of all the windows",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--fdr_cutoff",
        dest="fdr_cutoff",
        help="false discovery rate threshold chosen",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--annotation",
        dest="annotation",
        help="fasta file of the transcriptome to get \
            sequences of the transcriptomic peaks.",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--defined_peaks",
        dest="defined_peaks",
        help="final peaks determined \
        to be true binding events",
        required=False,
        metavar="FILE")

    parser.add_argument(
        "--plot_folder",
        dest="plot_folder",
        help="destination folder for the plots ",
        required=True,
        metavar="FILE")
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
    total_enriched = pd.read_csv(
        options.enriched_regions,
        header=0,
        sep='\t',
        index_col=0,
        comment='#',
        engine='python')

    annotation = pd.read_csv(
        options.enriched_regions,
        header=0,
        sep='\t',
        index_col=0,
        comment='#',
        engine='python')
    annotation[['start', 'end']] = annotation[['start', 'end']].astype('int')

    # #########################################################################
    # Enrichment specific figures
    # _________________________________________________________________________

    sys.stdout.write('z_score distribution in the windows...\n')
    sys.stdout.flush()
    get_windows_scatterplot_colorbar(total_enriched, "z_score", "RdBu_r")
    plt.savefig(
        os.path.join(options.plot_folder, 'scatterplot_zscore.png'),
        transparent=True,
        dpi=300)
    plt.clf()
    sys.stdout.write('z_score distribution in the windows finished 1/10\n')
    sys.stdout.flush()
    # _________________________________________________________________________

    sys.stdout.write('FDR distribution in the windows...\n')
    sys.stdout.flush()
    get_windows_scatterplot_colorbar(total_enriched, "FDR", "RdBu")
    plt.savefig(
        os.path.join(options.plot_folder, 'scatterplot_FDR.png'),
        transparent=True,
        dpi=300)
    plt.clf()
    sys.stdout.write('FDR distribution in the windows finished 2/10\n')
    sys.stdout.flush()
    # _________________________________________________________________________

    sys.stdout.write('Significant windows chosen...\n')
    sys.stdout.flush()
    get_windows_scatterplot(total_enriched, options.fdr_cutoff)
    plt.savefig(
        os.path.join(options.plot_folder, 'scatterplot_significant.png'),
        transparent=True,
        dpi=300)
    plt.clf()
    sys.stdout.write('Significant windows chosen finished 3/10\n')
    sys.stdout.flush()
    # _________________________________________________________________________

    sys.stdout.write('Z-score distribution...\n')
    sys.stdout.flush()
    get_zscore_distribution(total_enriched)
    plt.savefig(
        os.path.join(options.plot_folder, 'z_score_distribution.png'),
        transparent=True,
        dpi=300)
    plt.clf()
    sys.stdout.write('Z-score distribution finished 4/10\n')
    sys.stdout.flush()
    # _________________________________________________________________________

    sys.stdout.write('Z-score FDR scatterplot...\n')
    sys.stdout.flush()
    get_zscore_distribution(total_enriched)
    plt.savefig(
        os.path.join(options.plot_folder, 'z_score_fdr_distribution.png'),
        transparent=True,
        dpi=300)
    plt.clf()
    sys.stdout.write('Z-score FDR scatterplot finished 5/10\n')
    sys.stdout.flush()
    return


def get_windows_scatterplot_colorbar(df, name, clr):
    """Scatterplot color according to column values"""
    windows = pd.read_csv(
        df,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')
    windows = windows.replace(0, np.nan)
    windows.dropna(inplace=True)
    windows['Reads_f'] = np.log(windows['Reads_f'] / windows['Reads_f'].sum())
    windows['Reads_b'] = np.log(windows['Reads_b'] / windows['Reads_b'].sum())
    f, ax = plt.subplots(figsize=(15, 10))
    points = ax.scatter(
        windows['Reads_f'],
        windows['Reads_b'],
        c=windows[name],
        vmin=windows[name].min(),
        vmax=windows[name].max(),
        s=0.7,
        cmap=clr)
    f.colorbar(points)
    # add a diagonal line
    plt.plot(
        [windows['Reads_b'].min(), windows['Reads_f'].max()],
        [windows['Reads_b'].min(), windows['Reads_f'].max()],
        linewidth=2,
        color='black')
    return


def get_windows_scatterplot(df, fdr_cutoff):
    """Scatterplot color according to hue"""
    fdr_cutoff = float(fdr_cutoff)
    df['significant'][df['FDR'] <= fdr_cutoff] = 'yes'
    df['significant'][df['FDR'] > fdr_cutoff] = 'no'
    pal = dict(yes="red", no="gray")
    windows = pd.read_csv(
        df,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')
    windows = windows.replace(0, np.nan)
    windows.dropna(inplace=True)
    windows['Reads_f'] = np.log(windows['Reads_f'] / windows['Reads_f'].sum())
    windows['Reads_b'] = np.log(windows['Reads_b'] / windows['Reads_b'].sum())
    f, ax = plt.subplots(figsize=(15, 10))
    ax.scatter(
        windows['Reads_f'],
        windows['Reads_b'],
        hue=['significant'],
        s=0.7,
        palette=pal)
    # add a diagonal line
    plt.plot(
        [windows['Reads_b'].min(), windows['Reads_f'].max()],
        [windows['Reads_b'].min(), windows['Reads_f'].max()],
        linewidth=2,
        color='black')
    return


def get_zscore_distribution(df):
    windows = pd.read_csv(
        df,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')
    df_plots = windows.copy(deep=True)
    df_plots.index = np.arange(len(df_plots))
    m, std = norm.fit(df_plots['z_score'])
    x = np.random.normal(m, std, len(df_plots['z_score']))

    f, ax = plt.subplots(figsize=(8, 10))
    sns.kdeplot(
        df_plots['z_score'],
        bw=0.08,
        shade=True,
        color="red")

    sns.kdeplot(
        x,
        bw=0.08,
        shade=True,
        color="black")
    plt.ylabel('Frequency (log)')
    plt.xlabel('z_score')
    plt.yscale('log')
    plt.show()
    plt.clf()
    return


def get_zscore_rev_cum_distribution(df):
    windows = pd.read_csv(
        df,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')
    df_plots = windows.copy(deep=True)
    df_plots.index = np.arange(len(df_plots))
    f, ax = plt.subplots(figsize=(8, 10))
    ax.scatter(df_plots['z_score'], df_plots['FDR'], s=20, color='red')
    plt.clf()
    plt.close()
    return


def get_topology(peaks, annotation):
    peaks.sort_values(by=['z_score'], ascending=False, inplace=True)
    peaks = peaks.head(1000)
    peaks.index = np.arange(len(peaks))
    peaks['new_start'] = peaks['start'] + peaks['mi'] - peaks['sigma']
    peaks['new_end'] = peaks['start'] + peaks['mi'] + peaks['sigma']
    peaks['new_mi'] = peaks['start'] + peaks['mi']

    peak_annotation = pd.DataFrame()
    counter = 0
    names = []
    part = pd.DataFrame()
    for index, row in peaks.iterrows():
        if 'ENS' in row['chromosome']:
            transcript_id = row['chromosome'].split('.')[0]
            tr_annotation = annotation['seqname'][(
                annotation['transcript_id'] == transcript_id) & (
                annotation['feature'] == 'transcript')]
            chromosome = tr_annotation['seqname'].values[0]
            transcript_start = tr_annotation['start'].values[0]
            gene_id = tr_annotation['gene_id'].values[0]

            strand = annotation['strand'][(
                annotation['gene_id'] == gene_id) & (
                annotation['feature'] == 'gene')].values[0]
            start = int(row['new_start'] + transcript_start)
            end = int(row['new_end'] + transcript_start)
            mi = int(row['new_mi'] + transcript_start)

        else:
            strand = row['strand']
            chromosome = row['chromosome']
            start = row['new_start']
            end = row['new_end']
            mi = int(row['new_mi'])

        part = annotation[(
            annotation['seqname'] == chromosome) & (
            annotation['strand'] == strand) & (
            (annotation['start'] <= mi) & (mi <= annotation['end']))]

        if part.empty:
            sys.stderr.write('Z-score FDR scatterplot...\n')
            sys.stderr.flush()
            continue

        if 'five_prime_utr' in part['feature'].values:
            for index2, row2 in part.iterrows():
                if row2['feature'] == 'five_prime_utr':
                    part = part.loc[[index2]]
                    break
        elif 'three_prime_utr' in part['feature'].values:
            for index2, row2 in part.iterrows():
                if row2['feature'] == 'three_prime_utr':
                    part = part.loc[[index2]]
                    break
        elif 'CDS' in part['feature'].values:
            for index2, row2 in part.iterrows():
                if row2['feature'] == 'CDS':
                    part = part.loc[[index2]]
                    break
        elif 'transcript' in part['feature'].values:
            for index2, row2 in part.iterrows():
                if row2['feature'] == 'transcript':
                    part = part.loc[[index2]]
                    tr_info = annotation[(
                        annotation['feature'] == 'exon') & (
                        annotation['transcript_id'] == row2['transcript_id'])]
                    max_exon = tr_info['exon_number'].astype('int').max()
                    if tr_info.empty:
                        break
                    if max_exon == 1:
                        break
                    elif max_exon == 2:
                        intron_start = tr_info['end'][
                            tr_info['exon_number'] == 1].values[0] + 1
                        intron_distance = tr_info['start'][
                            tr_info['exon_number'] == 2].values[0] - \
                            tr_info['end'][
                                tr_info['exon_number'] == 1].values[0] - 2
                        rel_intron_distance = (
                            mi - intron_start) / intron_distance
                    else:
                        for exon in np.arange(2, max_exon + 1):
                            if ((mi > tr_info['end'][
                                tr_info['exon_number'] == exon - 1].values[0]
                            ) & (mi < tr_info['start'][
                                tr_info['exon_number'] == exon].values[0])
                            ):
                                intron_start = tr_info['end'][
                                    tr_info['exon_number'] == exon -
                                    1].values[0] + 1
                                intron_distance = tr_info['start'][
                                    tr_info['exon_number'] == exon
                                ].values[0] - tr_info['end'][
                                    tr_info['exon_number'] == exon - 1
                                ].values[0] - 2
                                rel_intron_distance = (
                                    mi - intron_start) / intron_distance
                                break
                    break
        if part.empty:
            continue

        for index1, row1 in part.iterrows():
            if row1['feature'] == 'gene':
                continue
            names.append(row1['gene_id'])
            peak_annotation.at[counter, 'peak'] = str(chromosome) + ':' + str(
                start) + '-' + str(end) + ':' + row['strand']

            if row1['transcript_biotype'] == 'protein_coding':
                peak_annotation.at[counter, 'type'] = 'protein_coding'

                if row1['feature'] == 'three_prime_utr':
                    peak_annotation.at[counter, 'topology'] = 'three_prime_utr'
                    peak_annotation.at[counter, 'distance'] = (
                        mi - row1['start']) / (row1['end'] - row1['start'])

                elif row1['feature'] == 'five_prime_utr':
                    peak_annotation.at[counter, 'topology'] = 'five_prime_utr'
                    peak_annotation.at[counter, 'distance'] = (
                        mi - row1['start']) / (row1['end'] - row1['start'])

                elif row1['feature'] == 'CDS':
                    peak_annotation.at[counter, 'topology'] = 'CDS'
                    peak_annotation.at[counter, 'distance'] = (
                        mi - row1['start']) / (row1['end'] - row1['start'])

                else:
                    peak_annotation.at[counter, 'topology'] = 'intron'
                    peak_annotation.at[counter, 'distance'] = \
                        rel_intron_distance

                peak_annotation.at[counter, 'zscore'] = row['z_score']
            else:
                peak_annotation.at[counter, 'type'] = row1[
                    'transcript_biotype']
                if row1['feature'] == 'three_prime_utr':
                    peak_annotation.at[counter, 'topology'] = \
                        'nc_three_prime_utr'
                    peak_annotation.at[counter, 'distance'] = (
                        mi - row1['start']) / (row1['end'] - row1['start'])

                elif row1['feature'] == 'five_prime_utr':
                    peak_annotation.at[counter, 'topology'] = \
                        'nc_five_prime_utr'
                    peak_annotation.at[counter, 'distance'] = (
                        mi - row1['start']) / (row1['end'] - row1['start'])

                elif row1['feature'] == 'CDS':
                    peak_annotation.at[counter, 'topology'] = 'nc_CDS'
                    peak_annotation.at[counter, 'distance'] = (
                        mi - row1['start']) / (row1['end'] - row1['start'])

                else:
                    peak_annotation.at[counter, 'topology'] = 'nc_intron'
                    peak_annotation.at[counter, 'distance'] = \
                        rel_intron_distance
            counter += 1
    return peak_annotation


def get_topology_plot(df):
    sns.set_style(style='white')
    sns.set_context("talk")
    f, ax = plt.subplots(figsize=(6, 8))
    sns.boxplot(
        x="distance",
        y="topology",
        data=df,
        palette="RdBu")
    sns.swarmplot(
        x="distance",
        y="topology",
        data=df,
        color='black',
        size=2,
        linewidth=0)
    plt.show()
    plt.clf()

    f, ax = plt.subplots(figsize=(8, 10))
    sns.swarmplot(
        x="distance",
        y="topology",
        palette=["r", "c", "y"],
        data=df)
    plt.show()
    plt.clf()
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

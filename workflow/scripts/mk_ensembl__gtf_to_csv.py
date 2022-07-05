# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import sys
import multiprocessing as mp
import numpy as np
import os
import pysam
import logging
from io import StringIO
from contextlib import closing
from csv import writer
logger = mp.log_to_stderr(logging.DEBUG)


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

def main():
    """ Main function """

    __doc__ = "Convert ensembl gtf file to csv format."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--gtf_file",
        dest="gtf_file",
        help="Gtf file",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--flag",
        dest="flag",
        help="finish flag",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--out_file",
        dest="out_file",
        help="Output file csv formatted ensembl gtf",
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

    reader = pd.read_table(
        options.gtf_file,
        chunksize=2000,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')

    try:
        pool = mp.Pool(processes=12)
        attrs = pool.map(
            get_attributes,
            reader)
    finally:
        pool.close()
        pool.join()

    attributes = []
    for attr in attrs:
        attributes.extend(attr)
    attributes = list(set(attributes))

    attributes = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame'] + attributes
    parameters = []
    for i in pd.read_table(
            options.gtf_file,
            chunksize=2000,
            sep='\t',
            index_col=None,
            comment='#',
            engine='python'):
        parameters.append([i, attributes])

    with mp.Pool() as pool:
        with open(options.out_file, 'w') as output:
            writer_ens = writer(output, delimiter='\t')
            writer_ens.writerow(attributes)
            for rows in pool.imap_unordered(ensembl_table, parameters):
                writer_ens.writerows(rows)
    out = open(options.flag, 'w')
    out.close()
    return

# _________________________________________________________________________
# -------------------------------------------------------------------------
# Functions
# -------------------------------------------------------------------------


# def df_chunking(df, chunksize):
#     """Splits df into chunks, drops data of original df inplace"""
#     while len(df):
#         yield df.iloc[:chunksize].copy()


# def get_parameters(df, chunksize, attribute):
#     while len(df):
#         yield [df.iloc[:chunksize].copy(), attribute]
#         df.drop(df.index[:chunksize], inplace=True)


def get_attributes(gtf):
    '''get a table with the info from ensembl'''
    gtf.columns = [
        'seqname',
        'source',
        'feature',
        'start',
        'end',
        'score',
        'strand',
        'frame',
        'attribute']
    gtf.index = np.arange(len(gtf))
    df1 = gtf['attribute'].str.split('; ', expand=True)
    gtf.drop(['attribute'], axis=1, inplace=True)
    attr = []
    for column in df1:
        column = df1[column].str.strip()
        column.dropna(axis=0, how='any', inplace=True)
        df2 = column.str.split(' ', 1, expand=True)
        attr.extend(list(set(df2[0].values)))
    return attr


def ensembl_table(params):
    gtf = params[0]
    attr = params[1]
    '''get a table with the info from ensembl'''
    gtf.columns = [
        'seqname',
        'source',
        'feature',
        'start',
        'end',
        'score',
        'strand',
        'frame',
        'attribute']

    gtf.index = np.arange(len(gtf))
    df1 = gtf['attribute'].str.split('; ', expand=True)
    gtf.drop(['attribute'], axis=1, inplace=True)
    for column in df1:
        column = df1[column].str.strip()
        column.dropna(axis=0, how='any', inplace=True)
        df2 = column.str.split(' ', 1, expand=True)
        for index, row in df2.iterrows():
            if not row[0]:
                continue
            row_name = str(row[0]).replace(';', '').strip()
            row_value = str(row[1]).replace('"', '').replace(';', '').strip()
            gtf.loc[index, row_name] = row_value
    gtf = gtf.reindex(columns=attr)
    return gtf.values


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

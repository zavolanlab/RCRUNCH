# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# import needed (external) modules
# ----------------------------------------------------------------------------------------------------------------------
import sys
import pandas as pd
import re
import io
import os
import numpy as np
from argparse import ArgumentParser, RawTextHelpFormatter
import shutil
from pathlib import Path


# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# Main function
# ----------------------------------------------------------------------------------------------------------------------


def main():
    """ Convert the two motifs detected by phylogibbs into motevo
    compatible motif format. These will be refined in the refinement
    mode of Motevo."""

    __doc__ = "Phylogibbs motif format to Motevo compatible format."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--pwms",
        dest="pwms",
        help="motifs provided by attract",
        nargs='?',
        const='',
        required=False)

    parser.add_argument(
        "--names",
        dest="names",
        help="names used by attract",
        nargs='?',
        const='',
        required=False)

    parser.add_argument(
        "--rbp_name",
        dest="rbp_name",
        help="rbp name",
        required=False,
        default="")

    parser.add_argument(
        "--de_novo",
        dest="de_novo",
        nargs='*',
        help="motifs predicted by phylogibbs",
        required=False,
        default='')

    parser.add_argument(
        "--organism",
        dest="organism",
        help="names used by attract",
        nargs='?',
        const='',
        required=False,
        default="0")

    parser.add_argument(
        "--outdir",
        dest="outdir",
        help="output directory for motifs",
        required=True,
        default="0")

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

    # -------------------------------------------------------
    if options.de_novo:
        paths = list(options.de_novo)
        for path in paths:
            my_file = Path(path)
            if my_file.is_file():
                new_name = os.path.join(
                    options.outdir, 'motif_' + path.split('/')[-1] + '.pwm')
                shutil.copy(path, new_name)
            else:
                sys.stderr.write('File {0} not found!\n'.format(path))
    if options.pwms:
        get_motifs(options.names, options.pwms, options.outdir,
                    options.organism, options.rbp_name)
    else:
        pass
    return


def get_motifs(names, pwm, outdir, organism, rbp_name):
    attract_info = pd.read_csv(
        names,
        header=0,
        sep='\t',
        index_col=None,
        comment='#',
        engine='python')
    attract_info['Matrix_id'] = attract_info['Matrix_id'].astype('str')
    if organism != "0":
        attract_info = attract_info[
            attract_info['Organism'] == organism]
    if rbp_name:
        attract_info = attract_info[attract_info['Gene_name'] == rbp_name]

    motifs = open(pwm, 'r')
    each_motif = pd.DataFrame()
    info = ''
    counter_motif = 0
    for line in motifs.readlines():
        if line.startswith('>'):
            if counter_motif >= 1:
                each_motif = each_motif.round(3)
                each_motif.index = each_motif.index.map("{:02}".format)
                # outdirectory = os.path.join(outdir, str(info))
                # if not os.path.exists(outdirectory):
                #     os.makedirs(outdirectory)
                outfile = os.path.join(outdir, 'motif_' + str(info) + '.pwm')
                template = """//\nNA """ + str(info) + """\n{}//"""
                with open(outfile, 'w') as fp:
                    fp.write(template.format(each_motif.to_csv(
                        sep='\t',
                        index=True, header=True)))
            count = 1
            each_motif = pd.DataFrame()
            name = line.strip('\n').strip('>').split('\t')
            each_motif = pd.DataFrame()

            if name[0] in attract_info['Matrix_id'].tolist():
                info = attract_info[['Gene_name', 'Matrix_id']][
                        attract_info['Matrix_id'] == name[0]].values[0]
            else:
                counter_motif = 0
                continue
            info = '_'.join(info)
            info = info.replace(' ', '')
            info = info.replace('(', '')
            info = info.replace(')', '')
            info = info.replace(']', '')
            info = info.replace('[', '')
            info = info.replace('-', '_')
            counter_motif += 1
        else:
            line = line.strip('\n')
            values = line.split('\t')
            each_motif.loc[count, 'A'] = float(values[0]) * 100
            each_motif.loc[count, 'C'] = float(values[1]) * 100
            each_motif.loc[count, 'G'] = float(values[2]) * 100
            each_motif.loc[count, 'T'] = float(values[3]) * 100
            count += 1
    if counter_motif >= 1:
        each_motif.index = each_motif.index.map("{:02}".format)
        each_motif.index.name = 'P0'
        outfile = os.path.join(outdir, 'motif_' + str(info)+ '.pwm')
        template = """//\nNA """ + str(info) + """\n{}//"""
        with open(outfile, 'w') as fp:
            fp.write(template.format(each_motif.to_csv(
                sep='\t',
                index=True, header=True)))
    motifs.close()
    return

def read_pwm(path):
    """convert the motif from trasfac to a matrix """
    pwm = re.compile(r"^(\/\/\nNA(.*)\n([^\/]+))", re.MULTILINE)
    with open(path, 'r') as myfile:
        print(path)
        filein = str(myfile.read())
    for match in pwm.finditer(filein):
        name = match.group(2).replace(' ', '')
        found = match.group(3)
        motif = pd.read_csv(
            StringIO(found), sep='\s+', comment='#', engine='python')
        try:
            motif.drop(['cons', 'inf'], axis=1, inplace=True)
        except KeyError:
            pass
        try:
            motif.drop(['P0'], axis=1, inplace=True)
        except KeyError:
            pass
        try:
            motif.drop(['PO'], axis=1, inplace=True)
        except KeyError:
            pass
        motif = motif.astype('float64')
        motif = motif + 0.01
        a = pd.Series(motif.sum(axis='columns'))
        motif = motif.divide(a, axis='rows')
        motif = motif.astype('float64')
        motif.index = np.arange(1, len(motif) + 1)
        motif.index.name = 'pos'
        motif.index = motif.index.astype('int')

        motif1 = pd.DataFrame()
        # convert the pwm to bits info for the logo
        for index, row in motif.iterrows():
            row_ft = 0
            for nt in ['A', 'C', 'G', 'T']:
                row_ft += -row.loc[nt] * np.log2(row.loc[nt])
            for nt in ['A', 'C', 'G', 'T']:
                motif1.loc[index, nt] = row.loc[nt] * (2 - row_ft)
        motif = motif1
        motif_logo = logomaker.Logo(
            motif,
            font_name='Stencil Std',
            color_scheme='classic',
            vpad=0,
            width=0.9)

        # # style using Logo methods
        motif_logo.style_xticks(anchor=0, spacing=1, rotation=0)
        # ww_logo.highlight_position(p=1, color='gold', alpha=.5)
        # ww_logo.highlight_position(p=26, color='gold', alpha=.5)
        # style using Axes methods
        motif_logo.ax.set_ylabel('information (bits)')
        motif_logo.ax.set_xlim([0.5, len(motif)+ 0.5])
        motif_logo.ax.set_ylim([0, 2])
        return motif, name


# ___________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)

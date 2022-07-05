# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# import needed (external) modules
# ----------------------------------------------------------------------------------------------------------------------
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import os
from io import StringIO
from csv import writer
import re
import io
import multiprocessing as mp
import logging
import logomaker
import matplotlib.pyplot as plt
import matplotlib
from subprocess import Popen, PIPE, CalledProcessError
logger = mp.log_to_stderr(logging.DEBUG)
matplotlib.rc('figure', max_open_warning=0)

# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# Main function
# ----------------------------------------------------------------------------------------------------------------------
def main():
    """Assess enrichment of motifs in binding sites, using Motevo """

    __doc__ = "Assess enrichment of motifs in binding sites, using Motevo."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)

    parser.add_argument(
        "--phylogenetic_tree",
        dest="phylogenetic_tree",
        help="phylogenetic tree of the organism at hand",
        required=False,
        metavar="FILE")
    
    parser.add_argument(
        "--test_pool",
        dest="test_pool",
        help="path of alignments",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--outpath",
        dest="outpath",
        help="output file containing the parameters \
            needed for running motevo in refinement mode",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--wms",
        dest="wms",
        help="motif paths",
        nargs='+',
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--genome_tag",
        dest="genome_tag",
        help="genome name (eg:hg38)",
        required=True)

    parser.add_argument(
        "--outfile",
        dest="outfile",
        help="output with final enrichement values",
        required=True)

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

    # ------------------------------------------------------------------------------------------------------------------
    # this part is due to Motevo which runs in the local folder
    # and does not take input path as an argument
    initial_path = os.getcwd()
    if not os.path.exists(options.outpath):
        os.makedirs(options.outpath)
    os.chdir(options.outpath)
    outpath = os.getcwd()
    # ------------------------------------------------------------------------------------------------------------------

    # training and test set
    # contains half of the significant peaks
    # contains also random sequences to serve as background
    test_pool = os.path.join(initial_path, options.test_pool)

    phylogenetic_type = check_phylogenetic(test_pool)
    if phylogenetic_type == 'phylogenetic':
        phylogenetic_tree = open(options.phylogenetic_tree, 'r')
        phylogenetic_tree = str(phylogenetic_tree.read())
    else:
        phylogenetic_tree = ''

    parameters = []
    for wm in options.wms:
        wm = os.path.join(initial_path, wm)
        parameters.append([wm, test_pool,
                           options.genome_tag, phylogenetic_tree, outpath])

    output = StringIO()
    csv_writer = writer(output)
    pool = mp.Pool(processes=2)
    for motif_enr in pool.imap_unordered(motevo_process, parameters):
        csv_writer.writerow(motif_enr)
    pool.close()
    pool.join()
    output.seek(0)
    total_enrichment = pd.read_csv(output, names=[
        'motif', 'mean_enrichment', 'std_enrichment',
        'll_ratio_sum', 'beta', 'prior', 'n_bg'])
    os.chdir(initial_path)
    total_enrichment.to_csv(options.outfile, index=True, header=True, sep='\t')
    return


# ________________________________________________________________________
# -------------------------------------------------------------------------
# Functions
# -------------------------------------------------------------------------
def motevo_process(params):
    '''prior and posterior calls to Motevo and beta optimisation'''
    wm = params[0]
    test_pool = params[1]
    genome_tag = params[2]
    phylogenetic_tree = params[3]
    outpath = params[4]

    wmname = wm.split('/')[-1]
    sys.stdout.write('start {0}\n'.format(wmname))
    sys.stdout.flush()

    # specify the wm path
    each_path = os.path.join(outpath, 'motevo', wmname)
    if not os.path.exists(each_path):
        os.makedirs(each_path)
    os.chdir(each_path)

    # ------------------------------------------------------------
    # Get motif logo
    each_logo = read_pwm(wm)
    if not each_logo:
        return pd.Series([])
    each_logo.fig.savefig(f'{wmname}.pdf', transparent=True)
    plt.clf()

    bgprior = 0.99

    # ----------------------------------------------------------------
    # calculate motevo posteriors on training set - used for beta fitting

    beta = 0.0001
    # ------------------------------------------------------------------------------------------------------------------
    # Calculate the posteriors for each of the motifs
    posterior_report = os.path.join(each_path, 'posterior_report')
    parameters_posterior = os.path.join(each_path, 'posterior_params')
    posterior_sites = os.path.join(each_path, 'posterior_sites')
    posteriors = os.path.join(each_path, 'posteriors')

    # make the motevo parameter file for the posteriors
    motevo_parameters_posterior(bgprior, posteriors, posterior_sites,
                                parameters_posterior, genome_tag,
                                phylogenetic_tree)
    posterior_motevo = ['motevo', test_pool, parameters_posterior, wm]

    # call Motevo
    with open(posterior_report, "wb") as outfile:
        # subprocess.call(posterior_motevo, stdout=outfile)
        try:
            posteriorcall = Popen(posterior_motevo,
                                  stdout=outfile, stderr=PIPE)
            posteriorcall.wait()
            stdout, stderr = posteriorcall.communicate()
        except CalledProcessError as err:
            print(f"Error ocurred: posterior {wmname} {err}")

    # ------------------------------------------------------------------------------------------------------------------
    # Fit the beta values
    window_nr_test = get_sequence_length(test_pool, wm)
    enr_avg, enr_std, ll_ratio_sum, n_bg, ll_ratio_p = calculate_enrichment(
        posterior_sites, beta, window_nr_test, wmname)
    if not enr_avg:
        enr_avg = 0
    if not enr_std:
        enr_std = 0
    if not ll_ratio_sum:
        ll_ratio_sum = 0
    if not n_bg:
        n_bg = 0
    if not bgprior:
        bgprior = 0
    enrichment = pd.Series([
        str(wmname), str(enr_avg), str(enr_std),
        str(ll_ratio_sum), str(beta), str(bgprior), str(n_bg)])
    return enrichment


def check_phylogenetic(sequences):
    '''check if alignments are phylogenetic or not'''
    double = 0
    single = 0
    with open(sequences, 'r') as myfile:
        counter = 0
        for line in myfile.readlines():
            counter += 1
            quote = 0
            for character in line:
                if character == '>':
                    quote += 1
            if quote == 1:
                single += 1
            elif quote == 2:
                double += 1
            if counter >= 30:
                break
    if ((single > 0) and (double > 0)):
        return 'phylogenetic'
    else:
        return 'simple'


def motevo_parameters_posterior(bgprior, posterior, sites, parampath,
                                genome_tag, phylogenetic_tree, ufewmprior=0,
                                ufewmfile=None, ufewmlen=8, ufeprint=0,
                                markovorderbg=1, bga=0.25, bgt=0.25, bgc=0.25,
                                bgg=0.25, restrictparses=0, printsiteals=1,
                                minposterior=0.01):
    if phylogenetic_tree:
        tree = str('TREE ' + phylogenetic_tree)
    else:
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
            'sitefile ' + sites + '\n' +
            'priorfile ' + posterior + '\n' +
            'printsiteals ' + str(printsiteals) + '\n' +
            'minposterior ' + str(minposterior) + '\n')
    return


def get_sequence_length(pool, motif_paths):
    '''get effective length of the sequence '''
    # (crunch_version) every sequence has plus and minus strand,
    # thus l_i = (len_i - len_wm)*2 - the values say on how many locations
    # a WM can sit on a sequence. window_nr contains fg and bg seqs
    # (rcrunch_version) we consider only the strand at hand since the binding
    # happens at the mRNA level, so l_i = (len_i - len_wm)
    wm_length = get_wmlength(motif_paths)
    length = {}
    with open(pool) as myfile:
        for line in myfile.readlines():
            if '>' in line:
                record_id = line.replace('>', '')
                record_id = record_id.replace('\n', '')

            elif len(line) > 1:
                seq_length = len(line.replace('\n', ''))
                if seq_length - wm_length > 0:
                    length[record_id] = seq_length - wm_length
            else:
                continue
    return length


def get_wmlength(motif_path):
    ''' get length of the pwm'''
    minlength = np.inf
    title = re.compile(r"^(\/\/\nNA.*\n([^\/]+))", re.MULTILINE)
    length = 0
    try:
        with open(motif_path, 'r') as myfile:
            filein = str(myfile.read())
        for match in title.finditer(filein):
            found = match.group(2)
            motif = pd.read_csv(
                io.StringIO(found), sep='\s+',
                comment='#', index_col=0, engine='python'
                )
            length = len(motif)
            if length < minlength:
                minlength = length
    except IOError:
        sys.stderr.write(
            f'Motif {motif_path} not found or incorrect format!\n')
    return minlength



def sum_of_posteriors(fname):
    fg_posteriors = {}
    bg_posterior = 0
    with open(fname, 'r') as myfile:
        for line in myfile.readlines():
            row = line.split(' ')
            if len(row) > 3:
                posterior = float(row[2])
                name = row[-1].strip()
                if re.search('_random_sequence', name):
                    # this means it is a bg sequence
                    bg_posterior += posterior
                else:
                    fg_posteriors.setdefault(name, 0.0)
                    fg_posteriors[name] += posterior
    return fg_posteriors, bg_posterior


def beta_derivative(fg_sites_norm, sites_norm, beta_i):
    derivative = -np.sum(sites_norm / (fg_sites_norm + beta_i))
    return derivative


def fit_beta(prior_sites, window_nr):
    """
    window_nr is a dictionary where the region ID is the key
    (the same ID as in the sitesfile) and the value is the length,
    """
    # binding_regions_sitecount, n_bg = sum_of_posteriors(sites)
    fg_posteriors, bg_posterior = sum_of_posteriors(prior_sites)
    if float(bg_posterior) == 0:
        bg_posterior = 1
    # since bg seqs have the same lengths as the fg peaks,
    # I can just take the average over all seqs
    bg_mean_length = np.mean(np.array(list(window_nr.values())))

    fg_sites_norm = []
    n_bg = 0
    for seq_id, length in window_nr.items():
        try:
            # normalise the number of detected sites by the length
            fg_sites_norm.append(fg_posteriors[seq_id] / length)
        except KeyError:
            # check if the sequence is randomised(bg)
            if re.search('_random_sequence', seq_id):
                n_bg += 1
            else:
                # fg sequence with no detected site
                fg_sites_norm.append(0.0)
    fg_sites_norm = np.array(fg_sites_norm)

    # the number of detected sites by motevo
    # divided by number of bg seqs (randomised)
    n_bg_avg = bg_posterior / n_bg

    # normalise by the average length
    bg_sites_norm = n_bg_avg / bg_mean_length
    sites_norm = fg_sites_norm - bg_sites_norm

    # find range for beta. Go up until negative:
    beta_max_cutoff = 10
    beta_min = 1.0e-12
    beta_max = 1.0e-12
    while beta_max <= beta_max_cutoff:
        if (
            beta_derivative(fg_sites_norm, sites_norm, beta_min) *
                beta_derivative(fg_sites_norm, sites_norm, beta_max) > 0):
            beta_max = beta_max * 1000
        else:
            break

    if beta_max == beta_min:
        sys.stderr.write(
            'Deriv is already negative at beta_min {0} \n'.format(beta_max))
        return beta_max

    elif beta_max >= beta_max_cutoff:
        sys.stderr.write(
            'Deriv does not become positive at beta {0} \n'.format(beta_max))
        return beta_max_cutoff

    der_max = beta_derivative(fg_sites_norm, sites_norm, beta_max)
    der_min = beta_derivative(fg_sites_norm, sites_norm, beta_min)
    cutoff = 1.0e-6
    beta_mid = np.mean([beta_min, beta_max])
    while abs(beta_max - beta_min) > cutoff:
        beta_mid = np.mean([beta_min, beta_max])
        der_mid = beta_derivative(fg_sites_norm, sites_norm, beta_mid)
        if der_mid * der_max < 0:
            beta_min = beta_mid
            der_min = der_mid
        elif der_mid * der_min < 0:
            beta_max = beta_mid
            der_max = der_mid
        elif der_mid == 0:
            return beta_mid
    return beta_mid


def calculate_enrichment(sitefile, beta, window_nr, wmname):
    """
    \\begin{equation}
        E{w} = \\exp[ \frac{1}{F} *
        \\sum_{p\\epsilon\\{Ptest\\}} \\log \\frac{n_{p,w} + \\beta l_p}
        { \\langle n_{bg,w} \\rangle + \\beta \\langle l \\rangle}]
    \\end{equation}
    """

    # siteFile = sites file of motevo on testPool
    # length = dict of lengths of fg and bg seqs (testPool)
    # beta = fitted on trainingPool
    # calculate how much more likely it is to immuno-precipitate
    # a true binding peak as opposed to a background sequence.

    fg_sites, tot_bg_sites = sum_of_posteriors(sitefile)

    # total number of background samples
    B = 0
    F = 0
    tot_bg_length = 0
    # nominator : n_{p,w} + \beta l_p for all peaks
    nominator = []
    for peak_id, peak_l in window_nr.items():
        try:
            nominator.append(np.log(fg_sites[peak_id] + (beta * peak_l)))
            F += 1
        except KeyError:
            if re.search('_random_sequence', peak_id):
                B += 1
                tot_bg_length += peak_l
            else:
                F += 1
                nominator.append(np.log(beta * peak_l))
    n_bg = tot_bg_sites / B
    bg_l = tot_bg_length / B
    # print(
    # 'name: ', wmname, '\n', '<nb>: ', tot_bg_sites, '\n', 'n_bg',
    # n_bg, '\n', 'F:', F, '\n', 'B', B, '\n', fg_sites, '\n' )
    nominator = np.array(nominator)
    denominator = np.log(n_bg + (beta * bg_l))
    ll_ratio_p = nominator - denominator

    mean_enrichment = np.exp(np.mean(ll_ratio_p))
    std_enrichment = np.exp(np.std(ll_ratio_p))
    ll_ratio_sum = np.sum(ll_ratio_p)
    sys.stdout.write('done!\n')
    sys.stdout.flush()
    return mean_enrichment, std_enrichment, ll_ratio_sum, n_bg, ll_ratio_p


def read_pwm(path):
    """convert the motif from trasfac to a matrix """
    pwm = re.compile(r"^(\/\/\nNA(.*)\n([^\/]+))", re.MULTILINE)
    with open(path, 'r') as myfile:
        filein = str(myfile.read())
    for match in pwm.finditer(filein):
        found = match.group(3)
        motif = pd.read_csv(
            io.StringIO(found), sep='\s+',
            comment='#', engine='python')
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
        if len(motif) == 0:
            return None
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
        motif_logo.ax.set_xlim([0.5, len(motif) + 0.5])
        motif_logo.ax.set_ylim([0, 2])
        return motif_logo


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

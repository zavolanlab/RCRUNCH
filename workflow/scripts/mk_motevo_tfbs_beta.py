
# -----------------------------------------------------------------------------
# Author : Severin Berger
# Company: Erik van Nimwegen, Biozentrum, Basel
# Part of the CRUNCH pipeline

# Modified and adapted for RCRUNCH by: Maria Katsantoni
# Company: Mihaela Zavolan, Biozentrum, Basel
# 11_07_2019
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# This script takes as input the output of Motevo run with TFBS mode and EM: 1
# These values are calculated by Motevo using the training set. Here, the final
# prior and the beta are calculated and plugged in so they are considered by
# Motevo when running in TFBS mode and EM: 0 for the test sequences.
# -----------------------------------------------------------------------------
# ______________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# import needed (external) modules
# ----------------------------------------------------------------------------------------------------------------------
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter
import os
import io
import re


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
        "--motif",
        dest="motif",
        help="motif in transfac format",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--training_pool",
        dest="training_pool",
        help="training pool set of sequences",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--test_pool",
        dest="test_pool",
        help="test pool set of sequences",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--prior_sites",
        dest="prior_sites",
        help="sites predicted in the training pool set",
        required=True)

    parser.add_argument(
        "--priorfile",
        dest="priorfile",
        help="priors predicted by the training pool set",
        required=True)

    parser.add_argument(
        "--posterior_sites",
        dest="posterior_sites",
        help="sites predicted for the test set",
        required=True)

    parser.add_argument(
        "--motif_name",
        dest="motif_name",
        help="name of the motif",
        required=True)

    parser.add_argument(
        "--outfile",
        dest="outfile",
        help="beta output file",
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

    # -------------------------------------------------------
    motif = options.motif
    if motif == ['']:
        sys.stderr.write("No motif provided!\n")
        sys.exit(0)

    window_nr_training = get_sequence_length(
        options.training_pool,
        options.motif)

    window_nr_test = get_sequence_length(
        options.test_pool,
        options.motif)

    prior = get_prior(options.priorfile)

    beta = fit_beta(options.prior_sites, window_nr_training)

    mean_enrichment, std_enrichment, ll_ratio_sum, n_bg = \
        calculate_enrichment_scores(
            options.posterior_sites,
            beta,
            window_nr_test)

    final_df = pd.DataFrame()
    final_df.loc[options.motif_name, 'mean_enrichment'] = str(mean_enrichment)
    final_df.loc[options.motif_name, 'std_enrichment'] = str(std_enrichment)
    final_df.loc[options.motif_name, 'll_ratio_sum'] = str(ll_ratio_sum)
    final_df.loc[options.motif_name, 'beta'] = str(beta)
    final_df.loc[options.motif_name, 'prior'] = str(prior)
    final_df.loc[options.motif_name, 'n_bg'] = str(n_bg)

    final_df.to_csv(
        options.outfile,
        index=False,
        header=True,
        sep='\t')


def get_sequence_length(pool, motif_paths):
    # (crunch_version) every sequence has plus and minus strand,
    # thus l_i = (len_i - len_wm)*2 thus, the values say on how many locations
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
    # print(length)
    return length


def get_wmlength(motif_path):
    minlength = np.inf
    title = re.compile(r"^(\/\/\nNA.*\n([^\/]+))", re.MULTILINE)
    length = 0
    try:
        with open(motif_path, 'r') as myfile:
            filein = str(myfile.read())
        for match in title.finditer(filein):
            found = match.group(2)
            motif = pd.read_csv(
                io.StringIO(found),
                sep='\s+',
                comment='#',
                index_col=0,
                engine='python')
            length = len(motif)
            if length < minlength:
                minlength = length
    except IOError:
        sys.stderr.write(
            'Motif {0} not found or incorrect format!\n'.format(
                motif_path))
    print('minlength', minlength)
    return minlength


def get_prior(priorfile):
    with open(priorfile, 'r') as myfile:
        for line in myfile.readlines():
            if line.startswith('background'):
                prior = str(line.split(' ')[1])
    print('prior', prior)
    return prior


def beta_derivative(fg_sites_norm, sites_norm, beta_i):
    derivative = -np.sum(sites_norm / (fg_sites_norm + beta_i))
    return derivative


def sum_of_posteriors(fname):
    fg_posteriors = {}
    bg_posterior = 0
    with open(fname, 'r') as myfile:
        for line in myfile.readlines():
            row = line.split(' ')
            print('row', row, row[2])
            posterior = float(row[2])
            if re.search('_random_sequence', row[-1]):
                # i.e. it is a bg seq
                bg_posterior += posterior
            else:
                fg_posteriors.setdefault(row[-1].strip(), 0.0)
                fg_posteriors[row[-1].strip()] += posterior

    return fg_posteriors, bg_posterior


def fit_beta(prior_sites, window_nr):
    """
    window_nr is a dictionary where the region ID is the key
    (the same ID as in the sitesfile) and the value is the length,
    """
    print('prior_sites', prior_sites)
    # binding_regions_sitecount, n_bg = sum_of_posteriors(sites)
    fg_posteriors, bg_posterior = sum_of_posteriors(prior_sites)
    print(fg_posteriors, bg_posterior)
    # since bg seqs have the same lengths as the fg peaks,
    # I can just take the average over all seqsPUM2_87
    bg_mean_length = np.mean(np.array(list(window_nr.values())))
    fg_sites_norm = []
    n_bg = 0
    for seq_id, length in window_nr.items():
        # print(seq_id)
        try:
            # normalise the number of detected sites by the length
            fg_sites_norm.append(fg_posteriors[seq_id] / length)
        except KeyError:
            # check if the sequence is a randomised(bg)
            if re.search('_random_sequence', seq_id):
                n_bg += 1
            else:
                # fg sequence with no detected site
                fg_sites_norm.append(0.0)

    fg_sites_norm = np.array(fg_sites_norm)
    # the number of detected sites by motevo
    # with EM:1 divided by number of bg seqs (randomised)
    n_bg_avg = bg_posterior / n_bg

    # normalise by the average length
    bg_sites_norm = n_bg_avg / bg_mean_length

    sites_norm = fg_sites_norm - bg_sites_norm

    # find range for beta. Go up until negative:
    beta_max_cutoff = 200

    beta_min = 1.0e-12
    beta_max = 1.0e-12
    while (beta_max <= beta_max_cutoff):
        if (beta_derivative(
                fg_sites_norm, sites_norm, beta_min) * beta_derivative(
                fg_sites_norm, sites_norm, beta_max)) > .0:
            beta_max = beta_max * 2
        elif (
            beta_derivative(fg_sites_norm, sites_norm, beta_min) *
                beta_derivative(fg_sites_norm, sites_norm, beta_max)) < .0:
            break
        else:
            break
    if beta_max >= beta_max_cutoff:
        sys.stderr.write(
            'For beta_max = {0} still no negative derivative.\
            Using beta_max.\n'.format(beta_max))
        return beta_max

    if beta_max == 1.0e-12:
        sys.stderr.write(
            'Deriv is already negative at beta_min {0} \n'.format(beta_max))
        return beta_max

    der_max = beta_derivative(fg_sites_norm, sites_norm, beta_max)
    der_min = beta_derivative(fg_sites_norm, sites_norm, beta_min)
    cutoff = 1.0e-12
    while abs(der_max - der_min) > cutoff:
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


def calculate_enrichment_scores(sitefile, beta, window_nr):
    """
    siteFile = sites file of motevo on testPool
    length = dict of lengths of fg and bg seqs (testPool)
    beta = fitted on trainingPool
    Calculate how much more likely it is to immuno-precipitate
    a true binding peak as opposed to a background sequence.

    \\begin{equation}
        E{w} = \\exp[ \frac{1}{F} *
        \\sum_{p\\epsilon\\{Ptest\\}} \\log \\frac{n_{p,w} + \\beta l_p}
        { \\langle n_{bg,w} \\rangle + \\beta \\langle l \\rangle}]
    \\end{equation}
    """

    fg_sites, tot_bg_sites = sum_of_posteriors(sitefile)
    total_n_p = 0

    # total number of background samples
    B = 0
    tot_bg_length = 0
    # nominator : n_{p,w} + \beta l_p for all peaks
    nominator = []
    for peak_id, peak_l in window_nr.items():
        try:
            n_p = fg_sites[peak_id]
            # print('np:', n_p)
            nominator.append(np.log(n_p + beta * peak_l))
            # print('nominator:', n_p + beta * peak_l)
            total_n_p += n_p
        except KeyError:
            if re.search('_random_sequence', peak_id):
                B += 1
                tot_bg_length += peak_l
            else:
                nominator.append(np.log(beta * peak_l))
                if peak_l == 0:
                    print('***nominator:', peak_id, peak_l)

    n_bg = tot_bg_sites / B
    bg_l = tot_bg_length / B

    nominator = np.array(nominator)
    denominator = np.log(n_bg + beta * bg_l)
    ll_ratio_p = nominator - denominator

    # print(
    #     '\nFG_sites:', fg_sites,
    #     '\ntotalbg sites:', tot_bg_sites,
    #     '\ntotal bg test seqs:', B,
    #     '\nnominator:', nominator,
    #     '\ndenominator:', denominator,
    #     '\nn_bg:', n_bg,
    #     '\nbg_l:', bg_l,
    #     '\nll_ratio_p:', ll_ratio_p)

    mean_enrichment = np.exp(np.mean(ll_ratio_p))
    std_enrichment = np.exp(np.std(ll_ratio_p))
    ll_ratio_sum = np.sum(ll_ratio_p)

    return mean_enrichment, std_enrichment, ll_ratio_sum, n_bg


# __________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!\n")
        sys.exit(0)

# -----------------------------------------------------------------------------
# Author : Mihaela Zavolan, Anneke
# Company: Mihaela Zavolan, Biozentrum, Basel
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# This script is takes as input a specific frequency matrix and a set of
# sequences and calculates the enrichment of this motif in those sequences.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import sys
from argparse import ArgumentParser, RawTextHelpFormatter


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------
def main():
    """ Main function """

    __doc__ = "Obtain enrichment of PSSMs on a set of sequences."

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter)
    # ----------------------------------------------------------
    parser.add_argument(
        "--sequences",
        dest="sequences",
        help="Sequences in fasta format",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--motif",
        dest="motif",
        help="PSSM motif",
        required=True,
        metavar="FILE")

    parser.add_argument(
        "--prefix",
        dest="prefix",
        help="Prefix for output files",
        required=False,
        default='')

    parser.add_argument(
        "--fg_prior",
        dest="fg_prior",
        help="Foreground prior",
        required=False,
        default=0.1)

    parser.add_argument(
        "--pseudocount",
        dest="pseudocount",
        help="Pseudocount added to the motifs",
        required=False,
        default=0.01)

    parser.add_argument(
        "--accuracy",
        dest="accuracy",
        help="Accuracy for prior convergence",
        required=False,
        default=0.001)

    parser.add_argument(
        "--maxiterations",
        dest="maxiterations",
        help="Maximum number of iterations for convergence",
        required=False,
        default=10)

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

    seqNames, seqs = input_seq(options.sequences)
    pseudo = options.pseudocount
    fg_wm, bg_wm = input_wm(options.motif, pseudo)

    fg_prior = options.fg_prior
    bg_prior = 1 - fg_prior
    # for prior optimization
    accuracy = options.accuracy
    maxIterations = options.maxiterations
    prefix = options.prefix
    # first: prior optimization
    fg_prior = all_prior_update(
        seqs, accuracy, fg_prior, maxIterations,
        fg_wm, bg_wm, prefix)
    bg_prior = 1 - fg_prior

    # scoring of the sites in each input sequence
    outposterior = open(
        prefix + '_PosteriorProbabilities_perSite_prior.tab' % fg_prior, 'w')
    outsumposterior = open(
        prefix +
        '_SumOfPosteriorProbabilities_perSequence_prior.tab' % fg_prior, 'w')
    for s in range(len(seqs)):
        seqname = seqNames[s]
        score = num_sites(seqs[s], fg_prior, bg_prior, fg_wm, bg_wm)
        posts = post_sites(seqs[s], fg_prior, bg_prior, fg_wm, bg_wm)

        outposterior.write(seqname + '\t' + '\t'.join(
            ['%s' % posts[m] for m in range(len(posts))]) + '\n')
        outsumposterior.write(seqname + '\t%1.5f\n' % score)

    outposterior.close()
    outsumposterior.close()
    return


def input_seq(path):
    '''get the name and seq from a fasta formated file'''
    with open(path, 'r') as fasta:
        seqnames = []
        seqs = []
        for line in fasta.readlines():
            line = line.rstrip()
            if line.startswith('>'):
                seqnames.append(line)
            elif line:
                seqs.append(line.upper())
    return seqnames, seqs


def input_wm(infile, pseudo):
    '''read input weight matrix'''
    fg_wm = pd.read_csv(
        infile,
        header=0,
        sep='\t',
        index_col=0,
        comment='#',
        engine='python')
    fg_wm = (fg_wm + pseudo)
    fg_wm = fg_wm.divide(fg_wm.sum(axis=1), axis=0)
    fg_wm = fg_wm.astype('float64')
    bg_wm = pd.DataFrame()
    for i in fg_wm.columns.values:
        bg_wm.loc[1, i] = 1 / len(fg_wm.columns)
    return fg_wm, bg_wm


def inverse_wm(wm):
    wm = wm.iloc[::-1]
    wm.index = np.arange(1, len(wm) + 1)
    return wm


def site_sum(s, i, lwm, fg_wm):
    '''Likelihood of the sub-sequence starting at i
    in sequence s under the WM model'''
    fg_sum = 1.0
    counter = 0
    for j in list(fg_wm.index.values):
        fg_sum = fg_sum * float(fg_wm.at[j, s[i + counter]])
        counter += 1
    return fg_sum


def bgsite_sum(s, i, l, bg_wm):
    '''Likelihood of the sub-sequence starting at i
    in sequence s under the background model'''
    bg_sum = 1.0
    for j in range(l):
        bg_sum = bg_sum * bg_wm.loc[1, s[i + j]]
    return bg_sum


def fw_xlog_sum(seq, fg_prior, bg_prior, fg_wm, bg_wm):
    '''Forward recurrence formulas for the ratios F_n/F_n-1 and R_m/R_m+1'''
    # parameters:
    # seq - sequence
    # fg_prior - prior probability of a site
    # bg_prior - prior probability of background
    # fg_wm - weight matrix
    # bg_wm - background model

    lwm = len(fg_wm)
    if len(seq) < lwm:
        # there is no place to put a WM
        # only 0th order background
        # initialize the forward sum vector
        fw_sum = [0 for i in range(len(seq))]
        for i in range(len(seq)):
            fw_sum[i] = bgsite_sum(seq, i, 1, bg_wm)
        x_sum = np.multiply(fw_sum, bg_prior)

    else:
        fw_sum = [0 for i in range(len(fg_wm))]
        x_sum = [0 for i in range(len(fg_wm))]
        # initialize the forward sum for values where only background
        # can be placed set values for all bases in the sequence
        # that are at an index smaller than the length of the WM
        for i in range(len(fg_wm)):
            fw_sum[i] = bgsite_sum(seq, i, 1, bg_wm)
        fw_sum = np.multiply(fw_sum, bg_prior)
        # cumulative product of the columns (multiply by the previous value)
        fw_sum = np.cumprod(fw_sum)

        # calculation for the first position in the sequence where the WM can
        # end (includes both foreground and background)
        fw_sum[lwm - 1] = fw_sum[lwm - 1] + fg_prior * site_sum(
            seq, 0, lwm, fg_wm)

        # to avoid underflows, we now work with ratios of forward sums
        # first calculate x_sum - ratios of forward sums
        # r(i) = Z(1,i)/Z(1,i-1)
        # update rule r(i) = p(bg) + p(fg)/prod_i-l+1^i-1(r(k))
        # cf. Bussemaker, Li, SIggia - ISMB 2000
        # positions at the beginning of the sequence, up until we can fit a WM
        x_sum[0] = fw_sum[0]
        for i in range(lwm - 1):
            x_sum[i + 1] = fw_sum[i + 1] / fw_sum[i]

        # positions further down the sequence where a WM can be fit
        for i in range(len(seq) - lwm):
            # partition sum at i + len(WM) - background base at i+len(WM)
            # or WM at i+1..i+len(WM) product of partition function ratios
            # has to be taken up to 1 nc before the end of the WM
            x_sum.append(
                bg_prior * bgsite_sum(seq, lwm + i, 1, bg_wm) +
                (fg_prior * site_sum(seq, i + 1, lwm, fg_wm)) /
                np.prod(x_sum[i + 1: lwm + i]))

    # now build the log-forward partition sums
    x_log = np.log(x_sum)
    x_logsum = np.cumsum(x_log)
    return x_logsum


def bw_ylog_sum(seq, fg_prior, bg_prior, fg_wm, bg_wm):
    seq = seq[::-1]
    fg_wm = inverse_wm(fg_wm)
    x_logsum = fw_xlog_sum(seq, fg_prior, bg_prior, fg_wm, bg_wm)
    x_logsum = x_logsum[::-1]
    return x_logsum


def post_sites(seq, fg_prior, bg_prior, fg_wm, bg_wm):
    '''Posterior for a site starting at each nucleotide'''
    # initialize vectors of probability and log-probability
    lwm = len(fg_wm)
    log_post = [0.0 for i in range(len(seq) - lwm + 1)]

    if(len(seq) < lwm):
        # sequence too short to put a site
        return []

    # sequence long enough to hold a site
    # calculate forward partition sum
    fw = fw_xlog_sum(seq, fg_prior, bg_prior, fg_wm, bg_wm)
    # Z(1, L)
    FL = fw[len(seq) - 1]
    # calculate the backward partition sum
    bw = bw_ylog_sum(seq, fg_prior, bg_prior, fg_wm, bg_wm)

    # calculate posteriors per site
    if (len(seq) == lwm):
        # only one nucleotide in the sequence where we can place a site
        log_post[0] = np.log(fg_prior) + np.log(
            site_sum(seq, 0, lwm, fg_wm)) - FL

    else:
        # there are more nucleotides in the sequence where a site can be placed
        # for the first nucleotide we do not have a stored forward sum
        log_post[0] = np.log(fg_prior) + np.log(
            site_sum(seq, 0, lwm, fg_wm)) + bw[lwm] - FL

        # for all the others we do, so we have to include it in the calculation
        for i in range(len(seq) - lwm - 1):
            log_post[i + 1] = fw[i] + np.log(fg_prior) + np.log(
                site_sum(seq, i + 1, lwm, fg_wm)) + bw[i + lwm + 1] - FL

        # for the last nucleotide we do not have a stored backward sum
        log_post[len(seq) - lwm] = fw[len(seq) - lwm - 1] + np.log(
            fg_prior) + np.log(site_sum(seq, len(seq) - lwm, lwm, fg_wm)) - FL

    # exponentiate all log-probabilities
    post = np.exp(log_post)

    return post


def post_bgsites(seq, fg_prior, bg_prior, fg_wm, bg_wm):
    '''Posterior for a background site to start at a given position'''

    # initiate vectors of probability and log-probability
    log_bgpost = [0.0 for i in range(len(seq))]
    # calculate forward partition sum
    fw = fw_xlog_sum(seq, fg_prior, bg_prior, fg_wm, bg_wm)
    # Z(1, L)
    FL = fw[len(seq) - 1]
    # calculate backward partition sum
    bw = bw_ylog_sum(seq, fg_prior, bg_prior, fg_wm, bg_wm)

    # base case: only one nucleotide in the sequence
    if(len(seq) == 1):
        log_bgpost[0] = np.log(bg_prior) + np.log(
            bgsite_sum(seq, 0, 1, bg_wm)) - FL
    else:
        # first position in the sequence, there is no stored forward sum
        log_bgpost[0] = np.log(bg_prior) + np.log(
            bgsite_sum(seq, 0, 1, bg_wm)) + bw[1] - FL

        # all other nucleotides, we have both forward and backward sums
        for i in range(len(seq) - 2):
            log_bgpost[i + 1] = fw[i] + np.log(bg_prior) + np.log(
                bgsite_sum(seq, i + 1, 1, bg_wm)) + bw[i + 2] - FL

        # last position does not have a stored backward sum
        log_bgpost[len(seq) - 1] = fw[len(seq) - 2] + np.log(
            bg_prior) + np.log(bgsite_sum(seq, len(seq) - 1, 1, bg_wm)) - FL
    # exponentiate to get posterior probabilities
    bgpost = np.exp(log_bgpost)
    return bgpost


def num_sites(seq, fg_prior, bg_prior, fg_wm, bg_wm):
    '''Calculate posteriors per fg site and expectation value'''
    posterior = post_sites(seq, fg_prior, bg_prior, fg_wm, bg_wm)
    numSites = np.sum(posterior)
    return numSites


def bgnum_sites(seq, fg_prior, bg_prior, fg_wm, bg_wm):
    '''Calculate posteriors per bg site and expectation value'''
    bgpost = post_bgsites(seq, fg_prior, bg_prior, fg_wm, bg_wm)
    bgnumSites = np.sum(bgpost)
    return bgnumSites


def all_prior_update(seqs, accuracy, initialSitePrior, maxIterations, fg_wm, bg_wm, prefix):
    # define an accuracy for convergence
    site_prior = initialSitePrior
    site_delta = 1

    # keep track of the convergence in a log file
    logfile = open(prefix + "UpdatedPriors.log", 'w')
    logfile.write('%1.15f\n' % site_prior)
    logfile.close()

    # keep track of how many iterations are doing
    # in case we don't converge (fast enough) we want to quit
    numIterations = 0

    # write iterations to log file
    logfile = open(prefix + "UpdatedPriors.log", 'a')

    while(abs(site_delta) > accuracy and numIterations < maxIterations):
        numSites = 0.0
        bgNumSites = 0.0
        for i in range(len(seqs)):
            numSites = numSites + num_sites(
                seqs[i], site_prior, 1 - site_prior, fg_wm, bg_wm)
            bgNumSites = bgNumSites + bgnum_sites(
                seqs[i], site_prior, 1 - site_prior, fg_wm, bg_wm)
        new_site_prior = numSites / (numSites + bgNumSites)
        logfile.write('%d\t%d\t%1.15f\n' % (
            numSites, bgNumSites, new_site_prior))
        site_delta = new_site_prior - site_prior
        site_prior = new_site_prior

    if(numIterations < maxIterations):
        return site_prior

    logfile.write(
        'Prior update did not converge after %d iterations\n') % numIterations
    logfile.close()

    return site_prior


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

# -----------------------------------------------------------------------------
# Author : Severin Berger
# Company: Erik van Nimwegen, Biozentrum, Basel
# Part of the CRUNCH pipeline

# Modified and adapted for RCRUNCH by: Maria Katsantoni
# Company: Mihaela Zavolan, Biozentrum, Basel
# 16_07_2019
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# This script takes as input the final enrichment predictions for all the
# motifs (known and de novo) and remove the redunandant ones
# (by keeping the one with the) highest enrichment.
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
from math import isnan


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
        "--enrichments",
        dest="enrichments",
        nargs='+',
        help="enrichment of motifs",
        required=True)

    parser.add_argument(
        "--motifs",
        dest="motifs",
        help="motif in transfac format",
        required=True,
        metavar="DIR")

    parser.add_argument(
        "--output",
        dest="output",
        help="output file",
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
    dfs = []
    enrichments = list(options.enrichments)
    for enrichment in enrichments:
        try:
            enr_df = pd.read_csv(
                enrichment,
                header=None,
                index_col=0,
                sep='\t')
            dfs.append(enr_df)
        except:
            continue

    total_df = pd.concat(dfs)

    total_df.sort_values(by=[1], ascending=False, inplace=True)
    print(total_df)
    wm_reduced = []
    counter = 0
    for index, row in total_df.iterrows():
        counter += 1
        if counter == 1:
            wm_reduced.append(index)
            continue
        wm2_name = index
        wm2_path = os.path.join(options.motifs, wm2_name)
        wm2 = get_wm(wm2_path)
        
        similarity = 0
        for wm1_name in wm_reduced:
            wm1_path = os.path.join(options.motifs, wm1_name)
            wm1 = get_wm(wm1_path)

            dissimilarity = 1 - (
                (2 * getSimilarityScore(wm1, wm2)) / (
                    getSimilarityScore(wm1, wm1) + getSimilarityScore(wm2, wm2)))

            if dissimilarity < 0.2:
                similarity += 1
                break
        if similarity == 0:
            wm_reduced.append(wm2_name)
        else:
            print(index, wm1_name, dissimilarity)
    total_df = total_df[total_df.index.isin(wm_reduced)]
    total_df.to_csv(
        options.output,
        sep='\t',
        header=False,
        index=True)
    return

def get_wm(path):
    with open(path, 'r') as myfile:
        title = re.compile(r"^(\/\/\nNA.*\n([^\/]+))", re.MULTILINE)
        filein = str(myfile.read())
        for match in title.finditer(filein):
            found = match.group(2)
            motif = pd.read_csv(
                io.StringIO(found),
                sep='\s+',
                comment='#',
                index_col=0,
                engine='python')
            if 'cons' in motif.columns.values:
                motif.drop(['cons', 'inf'], axis=1, inplace=True)
    return motif


def mysimilarityscore(wm1, wm2):
    inner_product = []
    lenwm1 = len(wm1)
    lenwm2 = len(wm2)

    wm1 = wm1.T
    wm2 = wm2.T

    for s in np.arange(- lenwm1 + 2, lenwm2 + 1):
        wm1_inv = wm1.copy(deep=True)
        wm2_inv = wm2.copy(deep=True)
        new = np.arange(s, s + lenwm1)
        wm1_inv.columns = new
        new_keep = list(set(new) & set(wm2.columns.values))
        new_keep = np.sort(new_keep)
        wm1_inv = wm1_inv[new_keep]
        wm2_inv = wm2_inv[new_keep]
        inner_product.append(np.sum(np.diag(np.dot(wm1_inv, wm2_inv.T))))
    return np.max(inner_product)

def getSimilarityScore(wm1, wm2):
    wm1 = np.array(wm1)
    wm2 = np.array(wm2)
    s1 = np.zeros((wm1.shape[0] + wm2.shape[0] - 1, ))
    s2 = np.zeros((wm1.shape[0] + wm2.shape[0] - 1, ))
    for n in np.arange(wm1.shape[1]):
        s1 += np.convolve(wm1[:, n], wm2[:: - 1, n])
    score = np.vstack((s1)) #, s2))
    idx = np.argmax(score, axis=0)
    max_score = np.max(score, axis=0)
    return max_score


# def rename_and_copy_WM(wm, outwm, wmname):

#     wmlines = open(wm).readlines()
#     #NA Logo
#     wmlines[1] = 'NA %s\n' %wmname

#     o = open(outwm, 'w')
#     for i in wmlines:
#         o.write(i)


# def filterQuarterWMs(wmdict):
#     """
#     Sometimes WMs just contain 1 everywhere (I produce them in RunMotevo).
#     This function filters those out.
#     """

#     badwms = []

#     for wm in wmdict:
#         name, mat, matrev = wm2mat(open(wm))
#         if mat.shape[0] * mat.shape[1] == len(where(mat == 0.25)[0]):
#             badwms.append(wm)

#     for wm in badwms:
#         del wmdict[wm]

#     return wmdict


# def wm2mat(wmFile,ps=0.5):
#     """ Read swiss regulon style WM file """
#     M = False
#     ok = False
#     lines = wmFile.read().splitlines()
#     for line in lines:
#         if not (line.startswith("//") or line==""):
#             if ok:
#                 elm = reshape(array(map(float,line.split()[1:5])),(1,4))
#                 if M is False:
#                     M = elm
#                 else:
#                     M = vstack((M,elm))
#             else:
#                 if line[0:2]=='NA':
#                     na,name = line.split()
#                 elif line[0:2]=='PO' or line[0:2]=='P0': # WM start after PO/0 line
#                     ok = True
#                 else:
#                     pass
#     # automatic detection if WM is already normalized
#     if all(M.sum(axis=1)-1.0 < 0.0001):
#         M += 0.0001 # to avoid the case that we have completely polarized column
#     else: # matrix is in counts; add pseudo count
#         M += 0.5

#     # convert counts to frequencies
#     M = M/M.sum(axis=1)[:,newaxis]
#     # also build the reverse complement matrix
#     Mrev = vstack((M[::-1,3],M[::-1,2],M[::-1,1],M[::-1,0])).T

#     return(name,M,Mrev)


# def getSimilarityScore(wm1, wm1rev, wm2, wm2rev):

#     s1 = zeros((wm1.shape[0]+wm2.shape[0]-1,))
#     s2 = zeros((wm1.shape[0]+wm2.shape[0]-1,))
#     for n in arange(wm1.shape[1]): # over A,T,G,C                                                                                      
#         s1 += convolve(wm1[:,n],wm2[::-1,n])
#         s2 += convolve(wm1[:,n],wm2rev[::-1,n])
#     score = vstack((s1,s2))
#     idx = argmax(score,axis=0)
#     idx2 = argmax(score[idx,arange(score.shape[1])])
#     max_score = score[idx,arange(score.shape[1])][idx2] # idx[idx2] = 0 origianl, idx[idx2] = 1, reverse complement matches

#     return max_score


# def reduceWMs(wmdict, dist_co):
#     """
#     This function takes a dictionary with WM filepaths as keys and a score (AUCs or likelihoods) as values.
#     The function goes down the sorted dictionary (by score), computes which WMs from the same dictionary 
#     are similar (distance smaller than dist_co) and removes them from the dictionary.
#     It does this as long as there is no WM left in the given wmdict.
#     """

#     final_wms = [] #list containing names of WMs to give out. WMs that have high likelihood and are not similar to any other in this list

#     while len(wmdict.keys()) != 0:

#         print 'Remaining candidate WMs:', len(wmdict.keys())
#         #reformat wmdict to a matrix
#         wmmat1 = []
#         for k in wmdict:
#             wmmat1.append([k, wmdict[k]])

#         wmmat = sorted(wmmat1, key = lambda k: k[1], reverse=True)

#         refWM = wmmat[0][0]
#         print '-----------------'
#         print refWM
#         final_wms.append(refWM)
#         ignore, rWM, rWMrev = wm2mat(open(refWM))
#         queries = array(wmmat).T[0] #wmdict.keys()

#         #compute score between reference WM and itself to be able to compute distance later
#         rscore = getSimilarityScore(rWM, rWMrev, rWM, rWMrev)

#         for wm in queries:

#             ignore, qWM, qWMrev = wm2mat(open(wm))

#             #compute score between query WM and itself to  be able to compute distance later
#             qscore = getSimilarityScore(qWM, qWMrev, qWM, qWMrev)

#             #compute score between reference WM and query WM
#             rqscore = getSimilarityScore(rWM, rWMrev, qWM, qWMrev)

#             # normalize score to a distance between 0 (=both WMs are the same) and 1                                                                  
#             dist = 1 - 2*rqscore / (rscore + qscore)

#             # compute divergence. Measures how much the longer is away from the shorter! I take this instead of symmetrical distance to detect sub matrices as sub matrices!
#             rWMlen = rWM.shape[0]
#             qWMlen = qWM.shape[0]
#             if rWMlen < qWMlen:
#                 diverg = 1 - (2*rqscore) / (2*rscore)
#             elif qWMlen < rWMlen:
#                 diverg = 1 - (2*rqscore) / (2*qscore)
#             else:
#                 diverg = dist

#             print wm, dist, diverg
#             #if diverg <= dist_co:
#             if dist <= dist_co:
#                 del wmdict[wm]
#     return final_wms


# def convertFloat(x):
#     res = float(x)
#     if isnan(res):
#         raise Exception
#     return res 


# def execute(cf):
#     """
#     This component reduces candidate WMs.
#     """

#     ##Ports and parameters
#     infile = cf.get_input("infile")
#     outdir = cf.get_output("WMdir")
#     log_file = cf.get_output("log_file")
#     dist_co = cf.get_parameter("distance_cutoff", "float")
#     minscore = cf.get_parameter("minscore", "float")
    
#     os.mkdir(outdir)    
#     wmdict = {} #filename: [AUC, Likelihood]

#     badwms = []

#     for i, line in enumerate(open(infile)):
#         if i == 0:
#             continue
#         t = line.strip().split()
#         try:
#             score = convertFloat(t[1])
#             if score > minscore:
#                 wmdict[t[0]] = score
#             else:
#                 print "Motif %s wasn't included to the initial list (too low score)" % t[0]
#         except:
#             print "Motif %s wasn't included to the initial list (score probably nan)" % t[0]
#             continue 

#     wmdict = filterQuarterWMs(wmdict)

#     final_wms = reduceWMs(wmdict, dist_co)
#     print final_wms

#     l = open(log_file, 'w')

#     l.write('Passed WMs:\n')
#     for i, wm in enumerate(final_wms):
#         os.system('cp \'%s\' %s' %(wm, outdir))
#         #rename_and_copy_WM(wm, '%s/WM_%i' %(outdir, i+1), 'WM_%i' %(i+1))
#         l.write('WM_%i\t%s\n' %(i+1, wm))

#     l.close()


#     return 0



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

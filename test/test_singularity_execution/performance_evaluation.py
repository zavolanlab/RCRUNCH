import pandas as pd
import numpy as np
import os
import sys
import glob
import shutil
import re
import io

def main():
    """ Main function """
    __doc__ = "Evaluate the run based on expected values"
    outpath = os.path.dirname(os.path.realpath(__file__))
    os.chdir(outpath)
    a = pd.read_csv(
        'results/GN/PUM2_K562_ENCSR661ICQ_2/PUM2_K562_ENCSR661ICQ_2_total_peaks.csv',
        sep='\t', header=0, index_col=None)
    expected = pd.read_csv(
        'means_crosslink.tsv',
        sep='\t', header=0, index_col=None)
    a['center'] = a['start'] + a['crosslink']

    for index, row in a.iterrows():
        distance = (
            row['center'] - expected['crosslink'][expected['strand'] == row['strand']]).abs().min()
        if distance <= 50:
            a.loc[index, 'distance'] = distance
        else:
            a.loc[index, 'distance'] = np.nan
    rmsd = np.sqrt(np.sum((a["distance"]**2)) / len(a[~a.distance.isna()]))
    outlier_percentage = len(a[a.distance.isna()])/len(a) * 100


    # expected  = pd.read_csv(
    #     'means_peak_center.tsv',
    #     sep='\t', header=0, index_col=None)
    # a['center'] = a['start'] +  a['mi']

    # for index, row in a.iterrows():
    #     distance = (row['center'] - expected['peak_center'][expected['strand'] == row['strand']]).abs().min()
    #     if distance <= 50:
    #         a.loc[index, 'distance'] = distance
    #     else:
    #         a.loc[index, 'distance'] = np.nan
    # sys.stdout.write(f'rmsd of peak_center as center: {np.sqrt(np.sum((a["distance"]**2)) / len(a[~a.distance.isna()]))}\n'
    # )
    # sys.stdout.write(f'percentage of outlier peaks: {len(a[a.distance.isna()])/len(a) * 100} %\n'
    # )

    # add the new test best pwm
    path = f'results/GN/PUM2_K562_ENCSR661ICQ_2/motif_analysis_final/peak_center/motif_enrichment.tsv'
    a = pd.read_csv(path, sep='\t', header=0, index_col=0)
    b = a.loc[a['mean_enrichment'][a['motif'].str.contains('trimmed')].idxmax(), 'motif']
    run = b.split('_')[-1].split('.')[0]
    name = path.split('/')[-4]
    path_new = path.replace(
        'motif_analysis_final',
        f'motif_analysis_crossvalidation/{run}').replace('motif_enrichment.tsv', f'wlib/{b}')
    shutil.copy(path_new, f'wlib/{name}')

    df = pd.DataFrame()
    for path1 in glob.glob('wlib/*'):
        name1 = path1.split('/')[-1]
        for path2 in glob.glob('wlib/*'):
            name2 = path2.split('/')[-1]
            sim = getsimilarity(path1, path2)
            df.loc[name1, name2] = sim
    mean_sim = np.mean([df.loc['PUM2_K562_ENCSR661ICQ_2', :].mean(), df.loc[:, 'PUM2_K562_ENCSR661ICQ_2'].mean()])
    
    sys.stdout.write(f'rmsd of crosslink centers: {rmsd}\n')
    sys.stdout.write(f'percentage of outlier peaks: {outlier_percentage} % \n')
    sys.stdout.write(f'motif similarity: {mean_sim}\n')
    if rmsd < 1.5:
        sys.stdout.write('Rmsd is low. Test passed. 1/3\n')
    else:
        sys.stdout.write('Rmsd seems to be too high 1/3\n')
    if outlier_percentage < 5:
        sys.stdout.write('Few outliers detected. Test passed. 2/3\n')
    else:
        sys.stdout.write('Rmsd seems to be too high 2/3\n')
    if mean_sim > 0.75:
        sys.stdout.write('Motif similarity is high.Test passed 3/3\n')
    else:
        sys.stdout.write('Similarity seems to be low 3/3\n')
    return



def getSimilarityScore(wm1, wm2):
    wm1 = np.array(wm1)
    wm2 = np.array(wm2)
    s1 = np.zeros((wm1.shape[0] + wm2.shape[0] - 1, ))
    s2 = np.zeros((wm1.shape[0] + wm2.shape[0] - 1, ))
    for n in np.arange(wm1.shape[1]):
        s1 += np.convolve(wm1[:, n], wm2[:: - 1, n])
    score = np.vstack((s1)) #, s2))
    idx = np.argmax(score, axis=0)
    max_score = np.max(score, axis=0)[0]
    return max_score


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


def getsimilarity(wm1_path, wm2_path):
    wm1 = get_wm(wm1_path)
    wm2 = get_wm(wm2_path)
    similarity =  (
        (2 * getSimilarityScore(wm1, wm2)) / (
            getSimilarityScore(wm1, wm1) + getSimilarityScore(wm2, wm2)))
    return similarity


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


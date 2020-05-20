import re
import sys
import os
import pkg_resources
from io import TextIOWrapper
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
import numpy as np
from argparse import Namespace

from .protein_functions import is_phage_func, is_unknown_func

def find_training_genome(training_flag):
    try:
        f = pkg_resources.resource_stream('PhiSpyModules', 'data/trainingGenome_list.txt')
    except IOError as e:
        sys.stderr.write(f"There was an error opening data/trainingGenome_list.txt: {e}\n")
        return ''

    for line in f:
        temp = re.split('\t', line.decode().strip())
        if int(temp[0]) == training_flag:
            f.close()
            return temp[1].strip()
    return ''


def call_randomforest(**kwargs):
    output_dir = kwargs['output_dir']
    training_file = kwargs['training_set']
    infile = os.path.join(output_dir, "testSet.txt")
    outfile = os.path.join(output_dir, "classify.tsv")

    if not pkg_resources.resource_exists('PhiSpyModules', training_file):
        sys.stderr.write("FATAL: Can not find data file {}\n".format(training_file))
        sys.exit(-1)
    strm = pkg_resources.resource_stream('PhiSpyModules', training_file)
    train_data = np.genfromtxt(TextIOWrapper(strm), delimiter="\t", skip_header=1, filling_values=1)
    test_data = np.genfromtxt(fname=infile, delimiter="\t", skip_header=1, filling_values=1)
    if 'phmms' not in kwargs:
        train_data = np.delete(train_data, 5, 1)
        test_data = np.delete(train_data, 5, 1)
    """
    Przemek's comment
    by default 10 until version 0.22 where default is 100
    number of estimators also implies the precision of probabilities, generally 1/n_estimators
    in R's randomForest it's 500 and the usage note regarding number of trees to grow says:
    "This should not be set to too small a number, to ensure that every input row gets predicted at least a few times."
    """
    clf = RandomForestClassifier(n_estimators=kwargs['randomforest_trees'], n_jobs=kwargs['threads'])
    clf.fit(train_data[:, :-1], train_data[:, -1].astype('int'))
    np.savetxt(outfile, clf.predict_proba(test_data)[:,1])


def my_sort(orf_list):
    n = len(orf_list)
    i = 1
    while i <= n:
        j = i + 1
        while j < n:
            flag = 0
            # direction for both
            if orf_list[i]['start'] < orf_list[i]['stop']:
                dir_i = 1
            else:
                dir_i = -1
            if orf_list[j]['start'] < orf_list[j]['stop']:
                dir_j = 1
            else:
                dir_j = -1

            # check whether swap need or not
            if dir_i == dir_j:
                if orf_list[i]['start'] > orf_list[j]['start']:
                    flag = 1
            else:
                if dir_i == 1:
                    if orf_list[i]['start'] > orf_list[j]['stop']:
                        flag = 1
                else:
                    if orf_list[i]['stop'] > orf_list[j]['start']:
                        flag = 1
            # swap
            if flag == 1:
                temp = orf_list[i]
                orf_list[i] = orf_list[j]
                orf_list[j] = temp
            j = j+1
        i = i+1
    return orf_list


def find_mean(all_len):
    s = 0.0
    for i in all_len:
        s = s + i
    return float(s)/len(all_len)


def calc_pp(func):
    func = func.replace('-', ' ')
    func = func.replace(',', ' ')
    x = 0
    if is_phage_func(func):
        x = 1
    elif is_unknown_func(func):
        x = 0.5

    # a few special cases
    if 'recombinase' in func or 'integrase' in func:
        x = 1.5
    if ('phage' in func) and ('shock' in func):
        x = 0
    if "dna binding domain" in func:
        x = 0
    return x


def make_initial_tbl(**kwargs):
    self = Namespace(**kwargs)
    x = []
    for entry in self.record:
        for feature in entry.get_features('CDS'):
            all = {}
            all['fig'] = feature.id
            all['function'] = feature.function
            all['contig'] = entry.id
            all['start'] = feature.start
            all['stop'] = feature.stop
            all['rank'] = 0.0
            all['status'] = 0
            all['pp'] = calc_pp(feature.function)
            x.append(all)
    try:
        infile = open(os.path.join(self.output_dir, 'classify.tsv'), 'r')
        outfile = open(os.path.join(self.output_dir, 'initial_tbl.tsv'), 'w')
    except:
        sys.exit('ERROR: Cannot open classify.tsv in make_initial_tbl')

    j = 0
    ranks = [[] for n in range(len(x))]
    for line in infile:
        val = float(line.strip())
        for k in range(j-int(self.window_size/2), j+int(self.window_size/2)):
            if k < 0 or k >= len(x) or j >= len(x) or x[k]['contig'] != x[j]['contig']:
                continue
            ranks[k].append(val)
        j += 1
    infile.close()
    if not self.keep:
        os.remove(os.path.join(self.output_dir, 'classify.tsv'))

    #calculate threshold
    y = []
    j = 0
    while j < len(x):
        x[j]['rank'] = sum(ranks[j]) / len(ranks[j])
        x[j]['extra'] =  ranks[j]
        y.append(x[j]['rank'])
        j = j+1

    y2 = np.array(y).reshape(-1, 1)
    km = KMeans(n_clusters = 2)
    km.fit(y2)
    centers = km.cluster_centers_
    threshold = max(centers[0][0], centers[1][0])
    """
    Note added by Rob:
    At this point we have the classifications for each ORF and we want to take a sliding window and decide where the phage should
    start. We have two calculations for a threshold for the rank: either the kmeans centers and finding things above the larger center
    or just a plain threshold.
    
    """
    j = 0
    outfile.write('fig_no\tfunction\tcontig\tstart\tstop\tposition\trank\tmy_status\tpp\tFinal_status\tstart of attL\tend of attL\tstart of attR\tend of attR\tsequence of attL\tsequence of attR\tReason for att site\n')
    while j < len(x):
        if x[j]['rank'] > threshold:
            x[j]['status'] = 1
        outfile.write(str(x[j]['fig']))
        outfile.write('\t')
        outfile.write(str(x[j]['function']))
        outfile.write('\t')
        outfile.write(str(x[j]['contig']))
        outfile.write('\t')
        outfile.write(str(x[j]['start']))
        outfile.write('\t')
        outfile.write(str(x[j]['stop']))
        outfile.write('\t')
        outfile.write(str(j))
        outfile.write('\t')
        outfile.write(str(x[j]['rank']))
        outfile.write('\t')
        outfile.write(str(x[j]['status']))
        outfile.write('\t')
        outfile.write(str(x[j]['pp']))
        outfile.write('\n')
        j = j+1
    outfile.close()

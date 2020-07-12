import re
import sys
import os
import pkg_resources
from io import TextIOWrapper
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
import numpy as np
from argparse import Namespace

from .protein_functions import is_phage_func, is_unknown_func, is_not_phage_func
from .genbank_accessory_functions import feature_id
from .log_and_message import log_and_message

def find_training_genome(training_flag):
    try:
        f = pkg_resources.resource_stream('PhiSpyModules', 'data/trainingGenome_list.txt')
    except IOError as e:
        log_and_message(f"There was an error opening data/trainingGenome_list.txt: {e}\n", c="RED", stderr=True,
                        loglevel="CRITICAL")
        return ''

    for line in f:
        temp = re.split('\t', line.decode().strip())
        if int(temp[0]) == training_flag:
            f.close()
            return temp[1].strip()
    return ''


def call_randomforest(**kwargs):
    training_file = kwargs['training_set']
    test_data = kwargs['test_data']

    if not pkg_resources.resource_exists('PhiSpyModules', training_file):
        log_and_message(f"FATAL: Can not find data file {training_file}\n", c="RED", stderr=True, loglevel="CRITICAL")
        sys.exit(11)
    strm = pkg_resources.resource_stream('PhiSpyModules', training_file)
    log_and_message(f"Using training set in {training_file}")
    train_data = np.genfromtxt(TextIOWrapper(strm), delimiter="\t", skip_header=1, filling_values=1)

    all_metrics = ['orf_length_med', 'shannon_slope', 'at_skew', 'gc_skew', 'max_direction', 'phmms']
    if kwargs['phmms']:
        kwargs['metrics'].append('phmms')

    if len(kwargs['metrics']) < len(all_metrics):
        log_and_message(f"Using the following metric(s): {', '.join(sorted(kwargs['metrics']))}.", c='GREEN', stderr=True, quiet=kwargs['quiet'])
        skip_metrics = [all_metrics.index(x) for x in set(all_metrics) - set(kwargs['metrics'])]
        train_data = np.delete(train_data, skip_metrics, 1)
        test_data = np.delete(test_data, skip_metrics, 1)
    else:
        log_and_message(f"Using all metrics: {', '.join(all_metrics)}.", c='GREEN', stderr=True, quiet=kwargs['quiet'])

    """
    Przemek's comment
    by default 10 until version 0.22 where default is 100
    number of estimators also implies the precision of probabilities, generally 1/n_estimators
    in R's randomForest it's 500 and the usage note regarding number of trees to grow says:
    "This should not be set to too small a number, to ensure that every input row gets predicted at least a few times."
    """
    log_and_message(f"Running the random forest classifier with {kwargs['randomforest_trees']} trees and {kwargs['threads']} threads", c='GREEN', stderr=True, quiet=kwargs['quiet'])
    if train_data.shape[1] > 1:
        clf = RandomForestClassifier(n_estimators=kwargs['randomforest_trees'], n_jobs=kwargs['threads'])
        clf.fit(train_data[:, :-1], train_data[:, -1].astype('int'))
        return clf.predict_proba(test_data)[:,1]
    else:
        log_and_message("None of the metrics were requested and so we did not run the random forest. Results may be variable.", c='PINK', stderr=True, quiet=kwargs['quiet'])
        return np.zeros(len(kwargs['test_data']))

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
    func = func.replace('/', ' ')
    func = func.lower()
    x = 0
    if is_phage_func(func):
        x = 1
    elif is_unknown_func(func):
        x = 0.5
    elif is_not_phage_func(func):
        x = 0

    # a few special cases
    if 'recombinase' in func or 'integrase' in func:
        x = 1.5

    return x


"""
Artemis colour codes

These come from the Setting Colours section of
https://sanger-pathogens.github.io/Artemis/Artemis/artemis-manual.pdf

0 white (RGB values: 255 255 255)
1 dark grey (RGB values: 100 100 100)
2 red (RGB values: 255 0 0)
3 green (RGB values: 0 255 0)
4 blue (RGB values: 0 0 255)
5 cyan (RGB values: 0 255 255)
6 magenta (RGB values: 255 0 255)
7 yellow (RGB values: 255 255 0)
8 pale green (RGB values: 152 251 152)
9 light sky blue (RGB values: 135 206 250)
10 orange (RGB values: 255 165 0)
11 brown (RGB values: 200 150 100)
12 pale pink (RGB values: 255 200 200)
13 light grey (RGB values: 170 170 170)
14 black (RGB values: 0 0 0)
15 mid red: (RGB values: 255 63 63)
16 light red (RGB values: 255 127 127)
17 pink (RGB values: 255 191 191)

"""


def make_initial_tbl(**kwargs):
    self = Namespace(**kwargs)
    data = []
    x = []
    for entry in self.record:
        for feature in entry.get_features('CDS'):
            pp_score = calc_pp(feature.function)
            if self.color:
                if pp_score > 1:
                    feature.qualifiers['colour'] = 2 # red
                elif 'mobile element protein' in feature.function.lower():
                    feature.qualifiers['colour'] = 6 # pink
                elif pp_score == 1:
                    feature.qualifiers['colour'] = 4 # blue
                elif is_unknown_func(feature.function):
                    feature.qualifiers['colour'] = 13  # light grey
                elif pp_score == 0.5:
                    feature.qualifiers['colour'] = 9 # light blue
            fid = feature_id(entry, feature)
            ft = {
                'fig': fid,
                'function': feature.function,
                'contig': entry.id,
                'start': feature.start,
                'stop': feature.stop,
                'rank': 0.0,
                'status': 0,
                'pp': pp_score,
            }


            x.append(ft)

    ranks = [[] for n in range(len(x))]
    for j, val in enumerate(self.rfdata):
        for k in range(j-int(self.window_size/2), j+int(self.window_size/2)):
            if k < 0 or k >= len(x) or j >= len(x) or x[k]['contig'] != x[j]['contig']:
                continue
            ranks[k].append(val)


    # calculate threshold
    y = []
    j = 0
    while j < len(x):
        x[j]['rank'] = sum(ranks[j]) / len(ranks[j])
        x[j]['extra'] = ranks[j]
        y.append(x[j]['rank'])
        j = j+1

    y2 = np.array(y).reshape(-1, 1)
    km = KMeans(n_clusters=2)
    km.fit(y2)
    centers = km.cluster_centers_
    threshold = max(centers[0][0], centers[1][0])
    log_and_message(f"Rank threshold calculated as {threshold}")
    """
    Note added by Rob:
    At this point we have the classifications for each ORF and we want to take a sliding window and decide where
    the phage should start. We have two calculations for a threshold for the rank: either the kmeans centers and
    finding things above the larger center or just a plain threshold.

    """

    """
    Eventually each row has:
        0. gene id           6. rank                12. start of attR
        1. function          7. my status           13. end of attR
        2. contig            8. pp                  14. sequence of attL
        3. start             9. final status        15. sequence of attR
        4. stop             10. start of attL       16. Reason for att site choice
        5. position         11. end of attL

    However, at this point we only have  0 .. 8
    """

    for i in range(len(x)):
        status = 1 if x[i]['rank'] > threshold else 0
        # for regions that are absolutely not prophages (esp. ribosomal proteins)
        # set the rank to 0. We should, perhaps, mark these as potentially mobile?
        if is_not_phage_func(x[i]['function']):
            x[i]['rank'] = 0
            status = 0
            x[i]['pp'] = 0
        thisrow = [
            x[i]['fig'],
            x[i]['function'],
            x[i]['contig'],
            x[i]['start'],
            x[i]['stop'],
            i,
            x[i]['rank'],
            status,
            x[i]['pp']
        ]
        data.append(thisrow)
    return data

"""
Compare the regions predicted to be prophages to the regions that are marked as prophages in our testing set

Probably the hardest part of this is the identifiers!
"""

import os
import sys
import argparse
import gzip
from Bio import SeqIO, BiopythonWarning
from PhiSpyModules import message, is_gzip_file

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def genbank_seqio(gbkf):
    if is_gzip_file(gbkf):
        handle = gzip.open(gbkf, 'rt')
    else:
        handle = open(gbkf, 'r')

    return SeqIO.parse(handle, "genbank")

def actual_phage_cds(gbkf, verbose=False):
    """
    Read the genbank file and return a list of features that are actually phage regions
    :param gbkf: the test genbank file with CDS marked with is_phage
    :param verbose: more output
    :return: a set of phage features
    """

    if verbose:
        message(f"Reading {gbkf}", "GREEN", "stderr")

    phage = {}
    nonphage = {}
    for seq in genbank_seqio(gbkf):
        for feat in seq.features:
            if feat.type == 'CDS':
                if 'product' not in feat.qualifiers:
                    feat.qualifiers['product'] = [f"Hypothetical protein (not annotated in {gbkf})"]
                if 'is_phage' in feat.qualifiers:
                    phage[str(feat.translate(seq, cds=False).seq).upper()] = feat.qualifiers['product'][0]
                else:
                    nonphage[str(feat.translate(seq, cds=False).seq).upper()] = feat.qualifiers['product'][0]

    return phage, nonphage

def predicted_genbank(predf, verbose=False):
    """
    Read the predictions from the genbank file and return
    a set of features
    :param predf: the predictions file
    :param verbose: more output
    :return: a set of predicted phage genes
    """

    if verbose:
        message(f"Reading {predf}", "GREEN", "stderr")

    predicted = {}
    for seq in genbank_seqio(predf):
        for feat in seq.features:
            if feat.type == 'CDS':
                if 'product' in feat.qualifiers:
                    predicted[str(feat.translate(seq, cds=False).seq).upper()] = feat.qualifiers['product'][0]
                else:
                    predicted[str(feat.translate(seq, cds=False).seq).upper()] = f"Hypothetical protein (not annotated in {predf})"

    if verbose:
        message(f"Found {len(predicted)} predicted prophage features", "BLUE", "stderr")

    return predicted

def predicted_regions(regf, gbkf, verbose):
    """
    Pull the phage genes from the regions
    :param regf: the regions file with contigs/start/stop
    :param gbkf: the genbank file used to make those predictions
    :param verbose: more output
    :return: a set of predicted phage genes
    """

    regions = {}
    if verbose:
        message(f"Reading {regf}", "GREEN", "stderr")
    with open(regf, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            assert(len(p) == 3), f"Expected a tple of [contig, start, stop] in {regf}"
            p[1] = int(p[1])
            p[2] = int(p[2])
            if p[0] not in regions:
                regions[p[0]] = []
            if p[2] < p[1]:
                regions[p[0]].append([p[2], p[1]])
            else:
                regions[p[0]].append([p[1], p[2]])

    if verbose:
        message(f"Reading {gbkf} again to get the phage regions", "GREEN", "stderr")

    predicted = {}
    for seq in genbank_seqio(gbkf):
        if seq.id in regions:
            for loc in regions[seq.id]:
                if verbose:
                    message(f"Getting from {loc[0]} to {loc[1]}", "PINK", "stderr")
                for feat in seq[loc[0]:loc[1]].features:
                    if feat.type == 'CDS':
                        if 'product' in feat.qualifiers:
                            predicted[str(feat.translate(seq[loc[0]:loc[1]], cds=False).seq).upper()] = feat.qualifiers['product'][0]
                        else:
                            predicted[str(feat.translate(seq[loc[0]:loc[1]], cds=False).seq).upper()] = f"Hypothetical protein (not annotated in {gbkf})"

    if verbose:
        message(f"Found {len(predicted)} predicted prophage features", "BLUE", "stderr")

    return predicted

def compare_real_predicted(phage: dict, nonphage: dict, predicted: dict, print_fp: bool, print_fn: bool, verbose: bool):
    """
    Compare the features that are real and predicted
    :param print_fn: print out the false negative matches
    :param print_fp: print out the false positive matches
    :param phage: actual phage features
    :param nonphage: actual non phage features
    :param predicted: predicted phage features
    :param verbose: more output
    :return:
    """

    if verbose:
        message(f"Comparing real and predicted", "GREEN", "stderr")

    # TP = phage intersection predicted
    # TN = nonphage intersection [not in predicted]
    # FP = nonphage intersection predicted
    # FN = phage intersection [not in predicted]

    # convert the keys to sets
    phage_set = set(phage.keys())
    nonphage_set = set(nonphage.keys())
    predicted_set = set(predicted.keys())

    # calculate not in predicted
    not_predicted = set()
    for s in phage_set.union(nonphage):
        if s not in predicted:
            not_predicted.add(s)

    print(f"Found:\nTest set:\n\tPhage: {len(phage)} Not phage: {len(nonphage)}")
    print(f"Predictions:\n\tPhage: {len(predicted)} Not phage: {len(not_predicted)}")

    tp = len(phage_set.intersection(predicted_set))
    tn = len(nonphage_set.intersection(not_predicted))
    fp = len(nonphage_set.intersection(predicted_set))
    fn = len(phage_set.intersection(not_predicted))

    print(f"TP: {tp}  FP: {fp}  TN: {tn}  FN: {fn}")
    try:
        accuracy = (tp+tn)/(tp + tn + fp + fn)
    except ZeroDivisionError:
        accuracy = "NaN"

    try:
        precision = tp/(tp+fp)
    except ZeroDivisionError:
        precision = "NaN"

    try:
        recall = tp/(tp+fn)
    except ZeroDivisionError:
        recall = "NaN"

    try:
        specificity = tn/(tn+fp)
    except ZeroDivisionError:
        specificity = "NaN"

    f1_score = "NaN"
    if accuracy != "NaN" and precision != "NaN" and recall != "NaN" and specificity != "NaN":
        try:
            f1_score = 2*(recall * precision) / (recall + precision)
        except ZeroDivisionError:
            f1_score = "NaN"

    if accuracy != "NaN":
        print(f"Accuracy:    {accuracy:.3f}\t(this is the ratio of the correctly labeled phage genes to the whole pool of genes")
    else:
        print("Accuracy: NaN")

    if precision != "NaN":
        print(f"Precision:   {precision:.3f}\t(This is the ratio of correctly labeled phage genes to all predictions)")
    else:
        print("Precision: NaN")

    if recall != "NaN":
        print(f"Recall:      {recall:.3f}\t(This is the fraction of actual phage genes we got right)")
    else:
        print("Recall: NaN")

    if specificity != "NaN":
        print(f"Specificity: {specificity:.3f}\t(This is the fraction of non phage genes we got right)")
    else:
        print("Specificity: NaN")
    
    if f1_score != "NaN":
        print(f"f1 score: {f1_score:.3f}\t(this is the harmonic mean of precision and recall, and is the best measure when, as in this case, there is a big difference between the number of phage and non-phage genes)")
    else:
        print("f1 score: NaN")

    if print_fp:
        for i in nonphage_set.intersection(predicted_set):
            print(f"FP\t{i}\t{nonphage[i]}")

    if print_fn:
        for i in phage_set.intersection(not_predicted):
            print(f"FN\t{i}\t{[phage[i]]}")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compare predictions to reality")
    parser.add_argument('-t', '--testfile', help='test file that has phage proteins marked with the is_phage qualifier', required=True)
    parser.add_argument('-p', '--predictfile', help='predictions genbank file that has each prophage as a sequence entry')
    parser.add_argument('-r', '--regionsfile', help='predictions regions file that has tuple of [contig, start, end]')
    parser.add_argument('--fp', help='print out the false positives', action='store_true')
    parser.add_argument('--fn', help='print out the false negatives', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    pred = None
    if args.predictfile:
        pred = predicted_genbank(args.predictfile, args.v)
    elif args.regionsfile:
        pred = predicted_regions(args.regionsfile, args.testfile, args.v)
    else:
        message("FATAL: Please provide either a predictions genbank or tsv file", "RED", "stderr")

    phage, nonphage = actual_phage_cds(args.testfile, args.v)

    compare_real_predicted(phage, nonphage, pred, args.fp, args.fn, args.v)

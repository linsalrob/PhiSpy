import os
import sys
import argparse
import re
from argparse import RawTextHelpFormatter
from argparse import ArgumentTypeError as err
from modules.pathtype import PathType


def print_list():
    f = None
    try:
        f = open(os.path.join(os.path.dirname(os.path.dirname(os.path.relpath(__file__))),'data/trainingGenome_list.txt'), 'r')
    except:
        sys.stderr.write('cannot find list')
        sys.exit(-1)
    for line in f:
        line = line.strip()
        temp = re.split('\t', line)
        if int(temp[3]) == 1:
            print("{}\t{}".format(temp[2], os.path.join(os.path.dirname(os.path.dirname(os.path.relpath(__file__))),'data/' + temp[1])))
    f.close()

def is_valid_file(x):
    if not x:
        x = os.path.join(os.path.dirname(os.path.dirname(os.path.relpath(__file__))),'data/genericAll.txt')
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

def get_args():
    usage = 'python3 PhiSpy.py [-opt1, [-opt2, ...]] infile'
    parser = argparse.ArgumentParser(
        description="phiSpy is a program for identifying prophages from among microbial genome sequences",
        epilog="(c) 2008-2018 Sajia Akhter, Katelyn McNair, Przemysław Decewicz, Rob Edwards, San Diego State University, San Diego, CA")
    parser.add_argument('infile', type=is_valid_file, help='Input file in genbank format', nargs='?')
    parser.add_argument('-m', '--make_training_data', type=str,
                             help='Create training data from a set of annotated genome files. Requires is_phage=1 qualifier in prophage\'s CDSs')
    parser.add_argument('-t', '--training_set', action='store', type=is_valid_file, default='',
                             help='Choose the most closely related set to your genome. [Default: 0]')
    parser.add_argument('-l', '--list', action='store_true', default=False,
                             help='List the available training sets and exit')
    #parser.add_argument('-c', '--choose', type=bool, default=False, const=True, nargs='?',
    #                         help='Choose a training set from a list (overrides -t)')
    parser.add_argument('-e', '--evaluate', type=bool, default=False, const=True, nargs='?',
                             help='Run in evaluation mode -- does not generate new data, but reruns the evaluation')
    parser.add_argument('-n', '--number', default=5, type=int,
                             help='Number of consecutive genes in a region of window size that must be prophage genes to be called. [Default: 5]')
    parser.add_argument('-w', '--window_size', default=30, type=int,
                             help='Window size of consecutive genes to look through to find phages. [Default: 30]')
    parser.add_argument('-g', '--nonprophage_genegaps', default=10, type=int,
                             help='The number of non phage genes betweeen prophages.  [Default: 10]')
    parser.add_argument('--expand_slope', action='store_true', default=False, 
                             help='Use the product of the slope of the Shannon scores in making test sets')
    parser.add_argument('--kmers_type', default='all', choices=['all', 'codon', 'simple'], type=str,
                             help='Type of kmers used for calculating Shannon scores. [Default: all]')
    #parser.add_argument('-i', '--input_dir', help='The input directory that holds the genome')
    parser.add_argument('-o', '--output_dir', help='The output directory to write the results')
    parser.add_argument('-qt', '--quiet', type=bool, default=False, const=True, nargs='?',
                             help='Run in quiet mode')
    parser.add_argument('-k', '--keep', type=bool, default=False, const=True, nargs='?',
                             help='Do not delete temp files')
    return parser.parse_args()

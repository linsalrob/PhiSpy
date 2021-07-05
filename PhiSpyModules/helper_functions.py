# -*- coding: utf-8 -*-
import os
import sys
import argparse
import re
import pkg_resources
import binascii
import logging

from .log_and_message import message
from .version import __version__

def print_list(list_type):
    f = None
    try:
        # with pip we use resource streams that may be files or from archives
        f = pkg_resources.resource_stream('PhiSpyModules', 'data/trainingGenome_list.txt')
    except:
        message('Cannot find the list of training sets. It should be in data/trainingGenome_list.txt', "RED", 'stderr')
        sys.exit(10)
    if list_type == 'short':
        message("Training Set\t# of genomes used\t1st genome", "GREEN", "stderr")
        for line in f:
            line = line.decode().strip()
            temp = re.split('\t', line)
            genomes = temp[2].split(';') if ';' in temp[2] else [temp[2]]
            message(f"data/{temp[1]}\t{temp[3]}\t{genomes[0]}", "PINK", "stderr")
    elif list_type == 'long':
        message("Training Set\t# of genomes used", "GREEN", "stderr")
        for line in f:
            line = line.decode().strip()
            temp = re.split('\t', line)
            genomes = temp[2].split(';') if ';' in temp[2] else [temp[2]]
            message(f"data/{temp[1]}\t{temp[3]}", "PINK", "stderr")
            for genome in genomes:
                message(f"- {genome}", "YELLOW", "stderr")
    f.close()


def is_valid_file(x):
    if not x or not os.path.exists(x):
        raise argparse.ArgumentTypeError("Checking for validity: {0} does not exist".format(x))
    return x


def is_gzip_file(f):
    """
    This is an elegant solution to test whether a file is gzipped by reading the first two characters.
    I also use a version of this in fastq_pair if you want a C version :)
    See https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed for inspiration
    :param f: the file to test
    :return: True if the file is gzip compressed else false
    """
    with open(f, 'rb') as i:
        return binascii.hexlify(i.read(2)) == b'1f8b'

def create_logger(self):
    """
    Add a logger to self
    :param self: the args parsed object
    :return: a logger for the args parsed object
    """
    if not self.log:
        self.log = os.path.join(self.output_dir, self.file_prefix + 'phispy.log')

    logger = logging.getLogger('PhiSpy')
    logger.setLevel(5)
    hdlr = logging.FileHandler(self.log)
    fmt = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    hdlr.setFormatter(fmt)
    logger.addHandler(hdlr)
    return logger

def get_args():
    parser = argparse.ArgumentParser(
        description="phiSpy is a program for identifying prophages from among microbial genome sequences",
        epilog="(c) 2008-2020 Sajia Akhter, Katelyn McNair, Przemysław Decewicz, Rob Edwards")
    parser.add_argument('infile', type=is_valid_file, help='Input file in genbank format', nargs='?')
    parser.add_argument('-o', '--output_dir', help='The output directory to write the results')
    parser.add_argument('-m', '--make_training_data', type=str,
                        help='Create training data from a set of annotated genome files. Requires is_phage=1 ' +
                             'qualifier in prophage\'s CDSs')
    parser.add_argument('-t', '--training_set', action='store', default='data/trainSet_genericAll.txt',
                        help='Choose the most closely related set to your genome. [Default %(default)s]')
    parser.add_argument('-l', '--list', choices=['short', 'long'], default=False,
                        help='List the available training sets in one of the formats [%(choices)s] and exit.')
    parser.add_argument('-p', '--file_prefix', default="",
                        help='An optional prefix to prepend to all of the output files')
    parser.add_argument('-e', '--evaluate', type=bool, default=False, const=True, nargs='?',
                        help='Run in evaluation mode -- does not generate new data, but reruns the evaluation')
    parser.add_argument('-n', '--number', default=5, type=int,
                        help='Number of consecutive genes in a region of window size that must be prophage genes' +
                             ' to be called. [Default: %(default)d]')
    parser.add_argument('-u', '--min_contig_size', default=5000, type=int,
                        help='Minimum contig size (in bp) to be included in the analysis. Smaller contigs will ' +
                             'be dropped. [Default: %(default)d]')
    parser.add_argument('-w', '--window_size', default=30, type=int,
                        help='Window size of consecutive genes to look through to find phages. [Default: %(default)d]')
    parser.add_argument('-g', '--nonprophage_genegaps', default=10, type=int,
                        help='The number of non phage genes betweeen prophages. [Default: %(default)d]')
    parser.add_argument('--phage_genes', default=1, type=int,
                        help='The minimum number of genes that must be identified as belonging to a phage for the ' +
                             'region to be included. The default is %(default)d or more genes.')
    parser.add_argument('--metrics', nargs='+', type=str, action='append',
                        help='The set of metrics to consider during classification. If not set, all metrics (orf_length_med, shannon_slope, at_skew, gc_skew, max_direction) will be considered.')
    parser.add_argument('-r', '--randomforest_trees', default=500, type=int,
                        help='Number of trees generated by Random Forest classifier. [Default: %(default)d]')
    parser.add_argument('--expand_slope', action='store_true', default=False,
                        help='Use the product of the slope of the Shannon scores in making test sets')
    parser.add_argument('--kmers_type', default='all', choices=['all', 'codon', 'simple'], type=str,
                        help='Type of kmers used for calculating Shannon scores. [Default: %(default)s]')
    parser.add_argument('--phmms', type=str,
                        help='Phage HMM profile database (like pVOGs) will be mapped against the genome of ' +
                            'interest and used as additional feature to identify prophages.\nNote that this ' +
                            'is experimental at the moment')
    incl_ann_group = parser.add_mutually_exclusive_group()
    incl_ann_group.add_argument('--include_annotations', action='store_true',
                        help='Use the annotations in the GenBank file in the phage predictions (default)')
    incl_ann_group.add_argument('--ignore_annotations', action='store_true',
                                help='Ignore the annotations in the GenBank file during predictions')
    parser.add_argument('--color', action='store_true',
                        help='If set, within the output GenBank file CDSs with phmms hits will be ' +
                             'colored (for viewing in Artemis).')
    parser.add_argument('--threads', type=int, default=2,
                        help='Number of threads to use while searching with phmms and the random forest. [default: %(default)d]')
    parser.add_argument('--output_choice', type=int, default=3,
                        help='Sum of codes for files to output. For more details see the README.md file at ' +
                            'https://github.com/linsalrob/PhiSpy#choosing-which-output-files-are-created. ' +
                            '[default: %(default)d]')
    parser.add_argument('--include_all_repeats', help='include all repeats in the genbank output', action='store_true')
    parser.add_argument('--keep_dropped_predictions', help='Write regions that might be prophages to the output files',
                        action='store_true')
    parser.add_argument('--extra_dna', type=int, default=2000,
                        help='additional DNA flanking the predicted prophage to test for repeats [Default %(default)d]')
    parser.add_argument('--min_repeat_len', type=int, default=10,
                        help='minimum repeat length to search for [Default: %(default)d]')
    parser.add_argument('--log', type=str,
                        help="Name of the log file to write details to [Default: phispy.log]")
    parser.add_argument('--quiet', action='store_true',
                        help='Run in quiet mode')
    parser.add_argument('-k', '--keep', type=bool, default=False, const=True, nargs='?',
                        help='Do not delete temp files')
    parser.add_argument('-v', '--version', action='version', version=__version__)
    args = parser.parse_args()


    ######################################
    #   list the training sets and exit  #
    ######################################
    if args.list:
        print_list(args.list)
        exit(0)

    if args.file_prefix != "":
        args.file_prefix += "_"
        args.file_prefix = args.file_prefix.replace(' ', '_')

    # check whether output directory was provided
    if not args.output_dir and not args.make_training_data:
        message("ERROR: Output directory (-o) is required. Use -h for more options\n", "RED", "stderr")
        sys.exit(3)
    elif args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
    else:
        args.output_dir=""


    # check whether any metrics were provided
    if not args.metrics:
        args.metrics = ['orf_length_med', 'shannon_slope', 'at_skew', 'gc_skew', 'max_direction']
    else:
        # flatten the list of lists caused by append. Note we should use extend but that is only python >=3.8
        m = [i for ml in args.metrics for i in ml]
        args.metrics = m

    # do we include the annotations?
    if args.ignore_annotations:
        args.include_annotations = False
    else:
        args.include_annotations = True

    args.logger = create_logger(args)
    return args

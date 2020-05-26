#!/usr/bin/env python
import os
import sys
import gzip
from functools import reduce

from Bio import SeqIO

INSTALLATION_DIR = os.path.dirname(os.path.realpath(__file__)) + '/'
sys.path.append(INSTALLATION_DIR)

import PhiSpyModules


def main(argv):  #organismPath, output_dir, trainingFlag, INSTALLATION_DIR, evaluateOnly, threshold_for_FN, phageWindowSize, quietMode, keep):

    ######################################
    #         parse the options          #
    ######################################
    args_parser = PhiSpyModules.get_args()

    ######################################
    #   list the training sets and exit  #
    ######################################
    if args_parser.list:
        PhiSpyModules.print_list()
        exit(0)

    if args_parser.version:
        print("PhiSpy version: {}".format(PhiSpyModules.version.__version__))
        sys.exit(0)

    # if we get here we need an input file
    if not args_parser.infile:
        PhiSpyModules.message("ERROR: Please provide an input file. Use -h for more options\n", "RED", 'stderr')
        sys.exit(-1)

    # check whether output directory was provided
    if not args_parser.output_dir and not args_parser.make_training_data:
        PhiSpyModules.message("ERROR: Output directory (-o) is required. Use -h for more options\n", "RED", 'stderr')
        sys.exit(-1)
    elif args_parser.output_dir:
        os.makedirs(args_parser.output_dir, exist_ok=True)

    ######################################
    #       add HMM search signal        #
    ######################################
    # if phmm search is required
    if args_parser.phmms:
        PhiSpyModules.message('Performing HMM search.\n', "GREEN", 'stderr')
        args_parser.infile = PhiSpyModules.search_phmms(**vars(args_parser))

    ######################################
    #        process input file          #
    ######################################
    # in future support other types
    # added a filter to remove short contigs. These break everything and we can't predict them to be 
    # phages anyway

    # RAE: Add support for gzipped files
    try:
        if PhiSpyModules.is_gzip_file(args_parser.infile):
            handle = gzip.open(args_parser.infile, 'rt')
        else:
            handle = open(args_parser.infile, 'r')
    except IOError as e:
        PhiSpyModules.message(f"There was an error reading {args_parser.infile}: {e}", "RED",'stderr')
        sys.exit(20)

    args_parser.record = PhiSpyModules.SeqioFilter(filter(lambda x: len(x.seq) > args_parser.min_contig_size, SeqIO.parse(handle, "genbank")))
    handle.close()

    # do we have any records left. Yes, this bug caught me out
    ncontigs = reduce(lambda sum, element: sum + 1, args_parser.record, 0)
    if ncontigs == 0:
        msg = f"Sorry, all of the contigs in {args_parser.infile} are less than {args_parser.min_contig_size}bp.\n"
        msg += "There is no data to process\n"
        PhiSpyModules.message(msg, "RED", "stderr")
        sys.exit(20)

    PhiSpyModules.message(f"Processing {ncontigs} contigs \n", "GREEN", 'stderr')

    ######################################
    #         make training set          #
    ######################################
    if args_parser.make_training_data:
        PhiSpyModules.message('Making Training Set...\n', "GREEN", 'stderr')
        my_make_train_flag = PhiSpyModules.make_set_train(**vars(args_parser))
        exit()

    ######################################
    #         make testing set           #
    ######################################
    PhiSpyModules.message('Making Testing Set...\n', "GREEN", 'stderr')
    args_parser.test_data = PhiSpyModules.measure_features(**vars(args_parser))

    ######################################
    #         do classification          #
    ######################################
    PhiSpyModules.message('Start Classification Algorithm...\n', "GREEN", 'stderr')
    args_parser.rfdata = PhiSpyModules.call_randomforest(**vars(args_parser))
    args_parser.initial_tbl = PhiSpyModules.make_initial_tbl(**vars(args_parser))

    ######################################
    #         i dont know what           #
    ######################################
    ###### added in this version 2.2 #####
    if (args_parser.training_set == 'data/trainSet_genericAll.txt'):
        PhiSpyModules.message('As the training flag is zero, down-weighting unknown functions\n', "RED", 'stderr')
        args_parser.initial_tbl = PhiSpyModules.downweighting_unknown_functions(args_parser)

    ######################################
    #         do evaluation              #
    ######################################
    PhiSpyModules.message('Evaluating...\n', "GREEN", 'stderr')
    PhiSpyModules.fixing_start_end(**vars(args_parser))
    PhiSpyModules.message('Done!!!\n', "GREEN", 'stderr')



if __name__== "__main__":
    main(sys.argv)

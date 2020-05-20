#!/usr/bin/env python
import os
import sys
import subprocess
import types
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
        sys.stderr.write("ERROR: Please provide an input file. Use -h for more options\n")
        sys.exit(-1)

    # check whether output directory was provided
    if not args_parser.output_dir:
        sys.stderr.write("ERROR: Output directory (-o) is required\n")
        sys.exit(-1)
    os.makedirs(args_parser.output_dir, exist_ok=True)

    ######################################
    #       add HMM search signal        #
    ######################################
    # if phmm search is required
    if args_parser.phmms:
        sys.stderr.write('Performing HMM search.\n')
        args_parser.infile = PhiSpyModules.search_phmms(**vars(args_parser))

    ######################################
    #        process input file          #
    ######################################
    # in future support other types
    # added a filter to remove short contigs. These break everything and we can't predict them to be 
    # phages anyway
    args_parser.record = PhiSpyModules.SeqioFilter(filter(lambda x: len(x.seq) > args_parser.min_contig_size, SeqIO.parse(args_parser.infile, "genbank")))
    # args_parser.record = input_file
    
    # do we have any records left. Yes, this bug caught me out
    ncontigs = reduce(lambda sum, element: sum + 1, args_parser.record, 0)
    if ncontigs == 0:
        sys.stderr.write(f"Sorry, all of the contigs in {args_parser.infile} are less than {args_parser.min_contig_size}bp.\n")
        sys.stderr.write("There is no data to process\n")
        sys.exit(20)

    sys.stderr.write(f"Processing {ncontigs} contigs \n")

    ######################################
    #         make training set          #
    ######################################
    if args_parser.make_training_data:
        sys.stderr.write('Making Train Set...\n')
        my_make_train_flag = PhiSpyModules.make_set_train(**vars(args_parser))
        exit()

    ######################################
    #         make testing set           #
    ######################################
    sys.stderr.write('Making Test Set...\n')
    my_make_test_flag = PhiSpyModules.make_test_set(**vars(args_parser))
    # check file im,plement later
    #if (my_make_test_flag == 0):
    #    print('The input organism is too small to predict prophages. Please consider large contig (having at least 40 genes) to use PhiSpy.')
    #    return

    ######################################
    #         do classification          #
    ######################################
    sys.stderr.write('Start Classification Algorithm...\n')
    PhiSpyModules.call_randomforest(**vars(args_parser))
    PhiSpyModules.make_initial_tbl(**vars(args_parser))

    ######################################
    #         i dont know what           #
    ######################################
    ###### added in this version 2.2 #####
    if (args_parser.training_set == 'data/trainSet_genericAll.txt'):
        sys.stderr.write('As the training flag is zero, considering unknown functions\n')
        PhiSpyModules.consider_unknown(args_parser.output_dir)

    ######################################
    #         do evaluation              #
    ######################################
    sys.stderr.write('Evaluating...')
    PhiSpyModules.fixing_start_end(**vars(args_parser))
    sys.stderr.write('Done!!!')



if __name__== "__main__":
    main(sys.argv)

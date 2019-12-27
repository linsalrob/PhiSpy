#!/usr/bin/env python
import os
import sys
import subprocess
import types

from Bio import SeqIO

INSTALLATION_DIR = os.path.dirname(os.path.realpath(__file__)) + '/'
sys.path.append(INSTALLATION_DIR)

from modules import seqio_filter
from modules import makeTrain
from modules import makeTest
from modules import classification
from modules import evaluation
from modules import unknownFunction
import modules.helper_functions as helpers


def main(argv):  #organismPath, output_dir, trainingFlag, INSTALLATION_DIR, evaluateOnly, threshold_for_FN, phageWindowSize, quietMode, keep):

    ######################################
    #         parse the options          #
    ######################################
    args_parser = helpers.get_args()

    ######################################
    #   list the training sets and exit  #
    ######################################
    if args_parser.list:
        helpers.print_list()
        exit(0)

    if args_parser.version:
        print(helpers.get_version())
        sys.exit(0)

    # if we get here we need an input file
    if not args_parser.infile:
        sys.stderr.write("ERROR: Please provide an input file. Use -h for more options\n")
        sys.exit(-1)

    # in future support other types
    input_file = seqio_filter.SeqioFilter(SeqIO.parse(args_parser.infile, "genbank"))
    args_parser.record = input_file
    if not args_parser.output_dir:
        sys.stderr.write("ERROR: Output directory (-o) is required\n")
        sys.exit(-1)
    os.makedirs(args_parser.output_dir, exist_ok=True)

    ######################################
    #         make training set          #
    ######################################
    if args_parser.make_training_data:
        print('Making Train Set...')
        my_make_train_flag = makeTrain.make_set_train(**vars(args_parser))
        exit()

    ######################################
    #         make testing set           #
    ######################################
    print('Making Test Set...')
    my_make_test_flag = makeTest.make_test_set(**vars(args_parser))
    # check file im,plement later
    #if (my_make_test_flag == 0):
    #    print('The input organism is too small to predict prophages. Please consider large contig (having at least 40 genes) to use PhiSpy.')
    #    return

    ######################################
    #         do classification          #
    ######################################
    print('Start Classification Algorithm...')
    classification.call_randomforest(**vars(args_parser))
    classification.make_initial_tbl(**vars(args_parser))

    ######################################
    #         i dont know what           #
    ######################################
    ###### added in this version 2.2 #####
    if (args_parser.training_set == 'data/trainSet_genericAll.txt'):
        print('As training flag is zero, considering unknown functions')
        unknownFunction.consider_unknown(args_parser.output_dir)

    ######################################
    #         do evaluation              #
    ######################################
    print('Evaluating...')
    evaluation.fixing_start_end(**vars(args_parser))
    print('Done!!!')

   

if __name__== "__main__":
    main(sys.argv)

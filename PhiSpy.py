#!/usr/bin/env python
import os
import sys
import subprocess

from Bio import SeqIO

INSTALLATION_DIR = os.path.dirname(os.path.realpath(__file__)) + '/'
sys.path.append(INSTALLATION_DIR)

from modules import makeTrain
from modules import makeTest
from modules import classification
from modules import evaluation
from modules import unknownFunction
import modules.helper_functions as helpers


def main(argv):  #organismPath, output_dir, trainingFlag, INSTALLATION_DIR, evaluateOnly, threshold_for_FN, phageWindowSize, quietMode, keep):
    ######################################
    #         check R install            #
    ######################################
    try:
        subprocess.call("type Rscript", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
    except OSError:
        sys.exit("The R programming language is not installed")

    ######################################
    #         parse the options          #
    ######################################
    args_parser = helpers.get_args()
    # in future support other types
    input_file = SeqIO.parse(args_parser.infile, "genbank")
    args_parser.records = input_file

    ######################################
    #         make training set          #
    ######################################
    if 0:
        print('Making Train Set... (need couple of minutes)')
        my_make_train_flag = makeTrain.call_make_train_set(**vars(args_parser))
        exit()

    ######################################
    #         make testing set           #
    ######################################
    print('Making Test Set... (need couple of minutes)')
    my_make_test_flag = makeTest.make_test_set(**vars(args_parser))
    exit()
    # check file im,plement later
    #if (my_make_test_flag == 0):
    #    print('The input organism is too small to predict prophages. Please consider large contig (having at least 40 genes) to use PhiSpy.')
    #    return
    ######################################
    #         do classification          #
    ######################################
    print('Start Classification Algorithm')
    classification.call_classification(**vars(args_parser))
    print('Done with classification Algorithm')

    ######################################
    #         i dont know what           #
    ######################################
    ###### added in this version 2.2 ##### 
    if (trainingFlag == 0):
        print('As training flag is zero, considering unknown functions')
        unknownFunction.consider_unknown(output_dir)

    ######################################
    #         do evaluation              #
    ######################################
    print('Start evaluation...')
    evaluation.call_start_end_fix(**vars(args_parser))
    print('Done!!!')

    ######################################
    #                                    #
    ######################################

if __name__== "__main__":
    main(sys.argv)

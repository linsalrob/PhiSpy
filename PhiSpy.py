#!/usr/bin/env python
import os
import sys
import re
import subprocess


try:
    import argparse
except ImportError:
    sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/argparse-1.2.1')
    import argparse


def call_phiSpy(organismPath, output_dir, trainingFlag, INSTALLATION_DIR, evaluateOnly, threshold_for_FN,
                phageWindowSize, quietMode, keep):
    sys.path.append(INSTALLATION_DIR + 'source/')
    import makeTest
    import classification
    import evaluation
    import unknownFunction

    sys.stderr.write("Running PhiSpy on " + organismPath + "\n")
    if (not evaluateOnly):
        if (quietMode == 0):
            print 'Making Test Set... (need couple of minutes)'

        my_make_test_flag = makeTest.call_make_test_set(organismPath, output_dir, INSTALLATION_DIR)
        if (my_make_test_flag == 0):
            print 'The input organism is too small to predict prophages. Please consider large contig (having at least 40 genes) to use PhiSpy.'
            return
        if (quietMode == 0):
            print 'Start Classification Algorithm'
        classification.call_classification(organismPath, output_dir, trainingFlag, phageWindowSize, INSTALLATION_DIR)

        if (quietMode == 0):
            print 'Done with classification Algorithm'

        ###### added in this version 2.2 ##### 
        if (trainingFlag == 0):
            if (quietMode == 0):
                print 'As training flag is zero, considering unknown functions'
            unknownFunction.consider_unknown(output_dir)
        ######################################

    if (quietMode == 0):
        print 'Start evaluation...'
    evaluation.call_start_end_fix(output_dir, organismPath, INSTALLATION_DIR, threshold_for_FN, phageWindowSize)
    if (quietMode == 0):
        print 'Done!!!'


def print_list(INSTALLATION_DIR):
    printstr = ''
    try:
        f = open(INSTALLATION_DIR + "data/trainingGenome_list.txt", "r")
    except:
        print'cannot find list'
    for line in f:
        line = line.strip()
        temp = re.split('\t', line)
        if int(temp[3]) == 1:
            printstr = printstr + temp[0] + ' ' + temp[2] + '\n'
    print printstr
    f.close()


def start_propgram(argv):
    try:
        subprocess.call("type Rscript", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
    except OSError:
        sys.exit("The R programming language is not installed")

    INSTALLATION_DIR = argv[0]
    if '/' in argv[0]:
        INSTALLATION_DIR = INSTALLATION_DIR[0:len(INSTALLATION_DIR) - INSTALLATION_DIR[::-1].find('/')]
    else:
        INSTALLATION_DIR = './'

    args_parser = argparse.ArgumentParser(
        description="phiSpy is a program for identifying prophages from among microbial genome sequences",
        epilog="(c) 2008-2018 Sajia Akhter, Katelyn McNair, Rob Edwards, San Diego State University, San Diego, CA")
    args_parser.add_argument('-t', '--training_set', default=0, type=int,
                             help='Choose a training set from the list of training sets.')
    args_parser.add_argument('-l', '--list', type=bool, default=False, const=True, nargs='?',
                             help='List the available training sets and exit')
    args_parser.add_argument('-c', '--choose', type=bool, default=False, const=True, nargs='?',
                             help='Choose a training set from a list (overrides -t)')
    args_parser.add_argument('-e', '--evaluate', type=bool, default=False, const=True, nargs='?',
                             help='Run in evaluation mode -- does not generate new data, but reruns the evaluation')
    args_parser.add_argument('-n', '--number', default=5, type=int,
                             help='Number of consecutive genes in a region of window size that must be prophage genes to be called')
    args_parser.add_argument('-w', '--window_size', default=30, type=int,
                             help='Window size of consecutive genes to look through to find phages')
    args_parser.add_argument('-i', '--input_dir', help='The input directory that holds the genome')
    args_parser.add_argument('-o', '--output_dir', help='The output directory to write the results')
    args_parser.add_argument('-qt', '--quiet', type=bool, default=False, const=True, nargs='?',
                             help='Run in quiet mode')
    args_parser.add_argument('-k', '--keep', type=bool, default=False, const=True, nargs='?',
                             help='Do not delete temp files')

    args_parser = args_parser.parse_args()

    if (args_parser.list):
        print_list(INSTALLATION_DIR)
        sys.exit(0)

    if not args_parser.input_dir and not args_parser.output_dir:
        print(sys.argv[0] + " [-h for help] [-l to list training sets] OPTIONS")
        print("Input and output directories are required")
        sys.exit(0)

    output_dir = args_parser.output_dir
    organismPath = args_parser.input_dir
    trainingFlag = args_parser.training_set

    output_dir = output_dir.strip()
    if output_dir[len(output_dir) - 1] != '/':
        output_dir = output_dir + '/'
    try:
        f = open(output_dir + 'testing.txt', 'w')
    except:
        try:
            os.makedirs(output_dir)
            f = open(output_dir + 'testing.txt', 'w')
        except:
            print "Cannot create the output directory or write file in the output directory", output_dir
            return
    f.close()
    os.system("rm " + output_dir + 'testing.txt')

    organismPath = organismPath.strip()
    if organismPath[len(organismPath) - 1] == '/':
        organismPath = organismPath[0:len(organismPath) - 1]

    try:
        f_dna = open(organismPath + '/contigs', 'r')
        f_dna.close()
    except:
        print "Cannot open", organismPath + '/contigs'
        return
    try:
        f = open(organismPath + '/Features/peg/tbl', 'r')
        f.close()
    except:
        print "Cannot open", organismPath + '/Features/peg/tbl'
        return
    try:
        f = open(organismPath + '/assigned_functions', 'r')
        f.close()
    except:
        print "Cannot open", organismPath + '/assigned_functions'
        # return
    try:
        f = open(organismPath + '/Features/rna/tbl', 'r')
        f.close()
    except:
        print "Cannot open", organismPath + '/Features/rna/tbl'
        # return

    if (args_parser.choose):
        while (1):
            print_list(INSTALLATION_DIR)
            temp = raw_input(
                "Please choose the number for a closely related organism we can use for training, or choose 0 if you don't know: ")
            try:
                trainingFlag = int(temp)
            except:
                continue
            if trainingFlag < 0 or trainingFlag > 30:
                continue
            break
        print ''

    call_phiSpy(organismPath, output_dir, trainingFlag, INSTALLATION_DIR, args_parser.evaluate, args_parser.number,
                args_parser.window_size, args_parser.quiet, args_parser.keep)


start_propgram(sys.argv)

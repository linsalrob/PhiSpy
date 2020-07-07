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
    PhiSpyModules.log_and_message(f"Starting PhiSpy.py with the following arguments\n{args_parser}")

    # if we get here we need an input file
    if not args_parser.infile:
        PhiSpyModules.log_and_message("ERROR: Please provide an input file. Use -h for more options", c="RED",
                                      stderr=True, stdout=False)
        sys.exit(2)

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
        PhiSpyModules.log_and_message(f"There was an error reading {args_parser.infile}: {e}", c="RED",
                                      stderr=True, stdout=False, quiet=args_parser.quiet)
        sys.exit(20)

    args_parser.record = PhiSpyModules.SeqioFilter(filter(lambda x: len(x.seq) > args_parser.min_contig_size, SeqIO.parse(handle, "genbank")))
    handle.close()

    # do we have any records left. Yes, this bug caught me out
    ncontigs = reduce(lambda sum, element: sum + 1, args_parser.record, 0)
    if ncontigs == 0:
        msg = f"Sorry, all of the contigs in {args_parser.infile} are less than {args_parser.min_contig_size}bp.\n"
        msg += "There is no data to process"
        PhiSpyModules.log_and_message(msg, c="RED", stderr=True, stdout=False, quiet=args_parser.quiet)
        sys.exit(30)

    PhiSpyModules.log_and_message(f"Processing {ncontigs} contigs", c="GREEN", stderr=True, stdout=False,
                                  quiet=args_parser.quiet)

    ######################################
    #       add HMM search signal        #
    ######################################
    # if phmm search is required
    if args_parser.phmms:
        args_parser.record = PhiSpyModules.search_phmms(**vars(args_parser))

    ######################################
    #         make training set          #
    ######################################
    if args_parser.make_training_data:
        if not args_parser.quiet:
            PhiSpyModules.log_and_message('Making Training Set...', c="GREEN", stderr=True, stdout=False,
                                          quiet=args_parser.quiet)
        my_make_train_flag = PhiSpyModules.make_set_train(**vars(args_parser))
        exit()

    ######################################
    #         make testing set           #
    ######################################
    PhiSpyModules.log_and_message('Making Testing Set...', c="GREEN", stderr=True, stdout=False,
                                  quiet=args_parser.quiet)
    args_parser.test_data = PhiSpyModules.measure_features(**vars(args_parser))
    if len(args_parser.test_data) < 100:
        m = f"Your genome only contains {len(args_parser.test_data)} ORFs. This is not enough to identify prophages.\n"
        m += "You might consider reannotating your genome with PROKKA or RAST"
        PhiSpyModules.log_and_message(m, c='RED', stderr=True, quiet=args_parser.quiet)
        sys.exit(41)

    ######################################
    #         do classification          #
    ######################################
    PhiSpyModules.log_and_message('Start Classification Algorithm...', c="GREEN", stderr=True, stdout=False,
                                  quiet=args_parser.quiet)
    args_parser.rfdata = PhiSpyModules.call_randomforest(**vars(args_parser))
    args_parser.initial_tbl = PhiSpyModules.make_initial_tbl(**vars(args_parser))

    ######################################
    #         i dont know what           #
    ######################################
    ###### added in this version 2.2 #####
    if (args_parser.training_set == 'data/trainSet_genericAll.txt'):
        PhiSpyModules.log_and_message('As the training flag is zero, down-weighting unknown functions', c="RED",
                                      stderr=True, stdout=False, quiet=args_parser.quiet)
        args_parser.initial_tbl = PhiSpyModules.downweighting_unknown_functions(args_parser)

    ######################################
    #         do evaluation              #
    ######################################
    PhiSpyModules.log_and_message('Evaluating...', c="GREEN", stderr=True, stdout=False, quiet=args_parser.quiet)
    args_parser.pp = PhiSpyModules.fixing_start_end(**vars(args_parser))

    ######################################
    #         write outputs              #
    ######################################
    PhiSpyModules.write_all_outputs(**vars(args_parser))

    # don't forget to close the log!
    PhiSpyModules.log_and_message(f'Done!!! Output is in {args_parser.output_dir} and the log is in {args_parser.log}',
                                  c="WHITE", stderr=True, stdout=False, quiet=args_parser.quiet)


if __name__== "__main__":
    main(sys.argv)

#!/usr/bin/python3
__author__ = 'Przemek Decewicz'

from argparse import ArgumentParser
from Bio import SeqIO
from glob import glob
from os import makedirs, path
from sys import argv


def read_genbank(infile):

    """
    Parses GenBank file's CDSs and groups them into host or phage groups based on the '/is_phage' qualifier.
    Return lists of sequences of both of these groups.
    """

    ref_orfs_list = []
    target_orf_list = []

    for record in SeqIO.parse(infile, 'genbank'):
        for f in record.features:
            if f.type == 'CDS':
                start = f.location.start
                end = f.location.end
                if f.location.strand == 1:
                    dna = str(record.seq[start : end])
                else:
                    dna = str(record.seq[start : end].reverse_complement())
                try:
                    status = f.qualifiers['is_phage'][0]
                    if status == '1':
                        target_orf_list.append(dna)
                    else:
                        ref_orfs_list.append(dna)
                except KeyError:
                    ref_orfs_list.append(dna)

    # print('Reference ORFs:', len(ref_orfs_list))
    # print('Targer ORFs:', len(target_orf_list))

    return ref_orfs_list, target_orf_list


def read_groups(infile):

    """
    Reads tab-delimited input file with the path to input file to use for training in the first column and the name of the group to put if afterwards while creating trainingGenome_list.txt.
    Returns a dict of groups and their input files.
    """

    groups = {}
    with open(infile) as inf:
        for line in inf:
            line = line.strip().split('\t')
            try:
                groups[line[1]].append(line[0])
            except KeyError:
                groups[line[1]] = [line[0]]

    return groups


def kmerize_orf(orf, k, t):

    """
    Creates kmers of size k and all, codon or simple type  from a single input sequence.
    Return a set of identified unique kmers.
    """

    kmers = []
    if t == 'simple':
        stop = len(orf) - (len(orf) % k)
        for i in range(0, stop, k):
            kmers.append(orf[i : i + k])
    elif t == 'all':
        for j in range(0, k):
            stop = len(orf) - ((len(orf) - j) % k)
            for i in range(j, stop, k):
                kmers.append(orf[i : i + k])
    elif t == 'codon':
        for j in range(0, k, 3):
            stop = len(orf) - ((len(orf) - j) % k)
            for i in range(j, stop, k):
                kmers.append(orf[i : i + k])

    return set(kmers)


def main():
    args = ArgumentParser(prog = 'make_training_sets.py', description = 'Automates making new or extending current PhiSpy\'s training sets.')

    args.add_argument('-i', '--infile',
                      type = str,
                      help = 'Path to input GenBank file.')

    args.add_argument('-d', '--indir',
                      type = str,
                      help = 'Path to input directory with multiple GenBank files.')

    args.add_argument('-g', '--groups',
                      type = str,
                      help = 'Path to file with path to input file and its group name in two tab-delimited columns. Otherwise each file will have its own training set.')

    args.add_argument('-o', '--outdir',
                      type = str,
                      help = 'Path to output directory. For each kmer creation approach subdirectory will be created.',
                      required = True)

    args.add_argument('-k', '--kmer_size',
                      type = int,
                      help = 'The size of required kmers. For codon approach use multiplicity of 3. [Default: 12]',
                      default = 12)

    args.add_argument('-t', '--type',
                      type = str,
                      help = 'Approach for creating kmers. Options are: simple (just slicing the sequence from the first position), all (all possible kmers), codon (all possible kmers made with step of 3 nts to get kmers corresponding translated aas). [Default: all]',
                      default = 'all')

    args.add_argument('--retrain',
                      type = str,
                      choices = ['yes', 'no'],
                      help = 'Retrain original training sets [yes] or extend them [no]. [Default: no]',
                      default = 'no')

    if len(argv[1:]) == 0:
        args.print_help()
        args.exit()

    try:
        args = args.parse_args()
    except:
        args.exit()

    if not args.infile and not args.indir: 
        print('You have to provide input data by either --infile or --indir.')
        exit(1)

    # Create output directory
    if not path.isdir(args.outdir): makedirs(args.outdir)

    # Create kmers for all input files
    infiles = glob(path.join(args.indir, '*')) if args.indir else [args.infile]
    infiles = infiles[:3]

    print('Making kmers from input file(s).')
    if not path.isdir(args.outdir): makedirs(args.outdir)

    kmers = {'PHAGE': {}, 'HOST': {}}
    for infile in infiles:
        print('  Processing: ', path.basename(infile))
        ref_orfs_list, target_orf_list = read_genbank(infile)

        ref_kmers = set()
        target_kmers = set()

        for i in ref_orfs_list:
            ref_kmers.update(kmerize_orf(i, args.kmer_size, args.type))
        for i in target_orf_list:
            target_kmers.update(kmerize_orf(i, args.kmer_size, args.type))

        # ref_minus_target = set()
        # target_minus_ref = set()
        # ref_minus_target = ref_kmers - target_kmers
        # target_minus_ref = target_kmers - ref_kmers
        # print('kmer_size\ttarget_kmers\treference_kmers\ttarget_specific\treference_specific\tintersection')
        # print(args.kmer_size, len(target_kmers), len(ref_kmers), len(ref_minus_target), len(target_minus_ref), len(ref_kmers.intersection(target_kmers)), sep = '\t')

        kmers['HOST'][infile] = ref_kmers
        kmers['PHAGE'][infile] = target_kmers

    # groups of resulting training sets
    if args.groups:
        groups = read_groups(args.groups)
    else:
        groups = {'group%05d' % (i + 1): [f] for i, f in enumerate(infiles)}


    # combine kmers into groups
    print('Making groups')
    host_kmers = set()
    phage_kmers = set()
    for g, i in groups.items():

        # HOST
        group_kmers = set()
        for infile in i:
            group_kmers.update(kmers['HOST'].pop(infile))

        # with open(path.join(args.outdir, '%s_HOST.txt' % g), 'w') as outf:
        #     outf.write('\n'.join(list(ref_kmers)))

        # kmers['HOST'][g] = group_kmers
        host_kmers.update(group_kmers)

        # PHAGE
        group_kmers = set()
        for infile in i:
            group_kmers.update(kmers['PHAGE'].pop(infile))

        # with open(path.join(args.outdir, '%s_PHAGE.txt' % g), 'w') as outf:
        #     outf.write('\n'.join(list(target_kmers)))

        # kmers['PHAGE'][g] = group_kmers
        phage_kmers.update(group_kmers)

    # write unique phage kmers
    print('Writing phage_kmers_' + args.type + '_wohost.txt') 
    with open(path.join(args.outdir, 'phage_kmers_' + args.type + '_wohost.txt'), 'w') as outf:
        outf.write('\n'.join(list(phage_kmers - host_kmers)))

    # make all training groups

    # make genericAll.txt

    # create trainingGenome_list.txt

if __name__ == '__main__':
    main()
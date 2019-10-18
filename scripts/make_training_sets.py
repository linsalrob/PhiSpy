#!/usr/bin/python3
__author__ = 'Przemek Decewicz'

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from Bio import SeqIO
from glob import glob
from os import makedirs, path
from subprocess import call
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


def read_groups(infile, indir):

    """
    Reads tab-delimited input file with the path to input file to use for training in the first column and the name of the group to put it afterwards while creating trainingGenome_list.txt. If --indir flag provided, then it also adds that path to input file name.
    Returns a dict of groups and their input files.
    """

    if not indir: indir = ''
    groups = {}
    with open(infile) as inf:
        for line in inf:
            line = line.strip().split('\t')
            try:
                groups[line[1]].append(path.join(indir, line[0]))
            except KeyError:
                groups[line[1]] = [path.join(indir, line[0])]

    return groups


def read_kmers(infile):

    """
    Simply reads input file with kmers and returns set of read kmers.
    """

    kmers = []
    with open(infile) as inf:
        kmers = [x.strip() for x in inf.readlines()]

    return set(kmers)


def kmerize_orf(orf, k, t):

    """
    Creates kmers of size k and all, codon or simple type from a single input sequence.
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

    return kmers


def main():
    args = ArgumentParser(prog = 'make_training_sets.py', 
                          description = 'Automates making new or extending current PhiSpy\'s training sets.',
                          epilog = 'Example usage:\npython3 scripts/make_training_sets.py -d tests -o data -g tests/groups.txt --retrain',
                          formatter_class = RawDescriptionHelpFormatter)

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
                      action = 'store_true',
                      help = 'If set, retrains original training sets, otherwise it extends what it finds in output directory.')

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


    # groups of resulting training sets
    if args.groups:
        groups = read_groups(args.groups, args.indir)
        infiles = set()
        for i in groups.values():
            infiles.update(set(i))
        infiles = sorted(list(infiles))
        print('Working on %i input files based on group file.' % len(infiles))
    else:
        if args.indir:
            infiles = glob(path.join(args.indir, r'*.gb'))
            infiles += glob(path.join(args.indir, r'*.gb[kf]'))
            infiles += glob(path.join(args.indir, r'*.gbff'))
            infiles = sorted(infiles)
        else:
            infiles = [args.infile]
        groups = {'group%05d' % (i + 1): [f] for i, f in enumerate(infiles)}


    # Create kmers for all input files
    print('Making kmers from input file(s).')
    if not path.isdir(args.outdir): makedirs(args.outdir)

    kmers = {'PHAGE': {}, 'HOST': {}}
    for i, infile in enumerate(infiles):
        print('  Processing %i/%i: %s' % (i + 1, len(infiles), path.basename(infile)))
        ref_orfs_list, target_orf_list = read_genbank(infile)

        ref_kmers = []
        target_kmers = []

        for i in ref_orfs_list:
            ref_kmers.extend(kmerize_orf(i, args.kmer_size, args.type))
        for i in target_orf_list:
            target_kmers.extend(kmerize_orf(i, args.kmer_size, args.type))

        kmers['HOST'][infile] = set(ref_kmers)
        kmers['PHAGE'][infile] = set(target_kmers)


    # combine kmers into groups
    print('Combining kmers into HOST and PHAGE groups.')
    host_kmers = set()
    phage_kmers = set()
    for infile in infiles:
        host_kmers.update(kmers['HOST'].pop(infile))
        phage_kmers.update(kmers['PHAGE'].pop(infile))


    # write unique phage kmers
    print('Writing phage_kmers_' + args.type + '_wohost.txt.')
    kmers_file = path.join(args.outdir, 'phage_kmers_' + args.type + '_wohost.txt')

    if path.isfile(kmers_file):
        if not args.retrain:
            print('  Reading %s.' % kmers_file)
            prev_kmers = read_kmers(kmers_file)
            phage_kmers.update(prev_kmers)
    else:
        print('  %s is missing - just making a new one.' % kmers_file)

    with open(kmers_file, 'w') as outf:
        outf.write('\n'.join(list(phage_kmers - host_kmers)))


    # make all training groups
    print('Making trainSets for each input file.')
    trainsets_outdir = path.join(args.outdir, 'trainSets')
    if not path.isdir(trainsets_outdir): makedirs(trainsets_outdir)
    phispy = path.join(path.dirname(path.dirname(path.realpath(__file__))), 'PhiSpy.py')
    for infile in infiles:
        print('  Processing %s' % infile)
        cmd = ['python3', phispy, infile, '-o', trainsets_outdir, '-m', path.basename(infile) + '.trainSet']
        call(cmd)


    # create trainingGenome_list.txt
    print('Writing trainingGenome_list.txt.')
    tg_file = path.join(args.outdir, 'trainingGenome_list.txt')
    if args.retrain or not path.isfile(tg_file):
        with open(tg_file, 'w') as outf:
            outf.write('0\ttestSet_genericAll.txt\tGeneric Test Set\t%i\n' % len(infiles))
            gcnt = 1
            for g, i in groups.items():
                outf.write('%i\ttrainSet_%s.txt\t%s\t%i\n' % (gcnt, g, ';'.join([path.basename(x) for x in i]), len(i)))
                gcnt += 1
    else:
        with open(tg_file) as inf:
            train_sets = {}
            inf.readline()
            for line in inf:
                line = line.strip().split('\t')
                line[0] = int(line[0])
                line[-1] = int(line[-1])
                train_sets[line[1].rsplit('.', 1)[0]] = line

        for t, l in train_sets.items():
            t = t.split('_')[1]
            if t in groups and l[3] == len(groups[t]):
                groups.pop(t)
            elif t in groups:
                infiles = [path.basename(i) for i in groups.pop(t)]
                train_sets = [l[0], l[1], ';'.join(infiles), len(infiles)]

        gcnt = len(train_sets) + 1
        if len(groups) > 0:
            for g, infiles in groups.items():
                infiles = [path.basename(i) for i in infiles]
                train_sets[g] = [gcnt, 'trainSet_%s.txt' % g, ';'.join(infiles), len(infiles)]
                gcnt += 1

        gsize = 0
        for t, l in train_sets.items():
            if l[3] == 1:
                gsize += 1

        with open(tg_file, 'w') as outf:
            outf.write('0\ttestSet_genericAll.txt\tGeneric Test Set\t%i\n' % (gsize))
            for t in sorted(train_sets.values()):
                outf.write('\t'.join([str(x) for x in t]) + '\n')


    # make or extend genericAll.txt
    print('Making training sets.')
    for g, i in groups.items():
        with open(path.join(args.outdir, 'trainSet_%s.txt' % g), 'w') as outf:
            for infile in i:
                trainset = path.join(trainsets_outdir, path.basename(infile) + '.trainSet')
                with open(trainset) as inf:
                    outf.write(inf.read())

    if args.retrain:
        with open(path.join(args.outdir, 'trainSet_genericAll.txt'), 'w') as outf:
            for infile in infiles:
                trainset = path.join(trainsets_outdir, path.basename(infile) + '.trainSet')
                with open(trainset) as inf:
                    outf.write(inf.read())
    else:
        # read all testSets in directory and combine them into generic test set
        print('*WARNING* - for updating generic train set only trainSets from single reference files are considered!')
        with open(path.join(args.outdir, 'trainingGenome_list.txt')) as inf:
            with open(path.join(args.outdir, 'trainSet_genericAll.txt'), 'w') as outf:
                for line in inf:
                    line = line.split()
                    if int(line[-1]) == 1: 
                        trainset = path.join(args.outdir, path.basename(line[1]))
                        with open(trainset) as infts:
                            outf.write(infts.read())



    print('Done!')

if __name__ == '__main__':
    main()
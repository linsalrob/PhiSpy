#!/usr/bin/python3
__author__ = 'Przemek Decewicz'

from argparse import ArgumentParser
from Bio import SeqIO
from glob import glob
from os import makedirs, path
from sys import argv


def read_genbank(infile):

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
                    dna = str(record.seq[end : start])
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


def kmerize_orf(orf, k, t):

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
    args = ArgumentParser()

    args.add_argument('-i', '--infile',
                      type = str,
                      help = 'Path to input GenBank file.')

    args.add_argument('-d', '--indir',
                      type = str,
                      help = 'Path to input directory with multiple GenBank files.')

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
                      help = 'Approach for creating kmers. Options are: simple (just slicing the sequence from the first position), all (all possible kmers), codon (all possible kmers made with step of 3 nts to get kmers corresponding translated aas). Multiple comma-separated types can be provided. [Default: simple,codon,all]',
                      default = 'simple,codon,all')

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

    # Read and process input files
    infiles = glob(path.join(args.indir, '*')) if args.indir else [args.infile]

    types = args.type.split(',') if ',' in args.type else [args.type]
    for t in types:
        print('Approach:', t)
        outdir = path.join(args.outdir, t)
        if not path.isdir(outdir): makedirs(outdir)
        for infile in infiles:
            print('  Processing: ', path.basename(infile))
            ref_orfs_list, target_orf_list = read_genbank(infile)

            ref_kmers = set()
            target_kmers = set()
            ref_minus_target = set()
            target_minus_ref = set()

            for i in ref_orfs_list:
                ref_kmers.update(kmerize_orf(i, args.kmer_size, t))
            for i in target_orf_list:
                target_kmers.update(kmerize_orf(i, args.kmer_size, t))

            ref_minus_target = ref_kmers - target_kmers
            target_minus_ref = target_kmers - ref_kmers
            # print('kmer_size\ttarget_kmers\treference_kmers\ttarget_specific\treference_specific\tintersection')
            # print(args.kmer_size, len(target_kmers), len(ref_kmers), len(ref_minus_target), len(target_minus_ref), len(ref_kmers.intersection(target_kmers)), sep = '\t')

            outfile_base = path.join(outdir, path.basename(infile).rsplit('.', 1)[0])
            with open(outfile_base + '_' + str(args.kmer_size) + '_HOST.txt', 'w') as outf:
                outf.write('\n'.join(list(ref_minus_target)))
            with open(outfile_base + '_' + str(args.kmer_size) + '_PHAGE.txt', 'w') as outf:
                outf.write('\n'.join(list(target_minus_ref)))


if __name__ == '__main__':
    main()
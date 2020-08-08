"""
Convert a phispy coordinates file and a phispy genbank file to fasta file(s) of bacterial and/or phage sequences
"""

import os
import sys
import argparse
from Bio import SeqIO
import gzip
import PhiSpyModules

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-c', help='PhiSpy coordinates file', required=True)
    parser.add_argument('-g', help='GenBank sequence file', required=True)
    parser.add_argument('-b', help='bacterial fraction output file')
    parser.add_argument('-p', help='phage fraction output file')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not args.p and not args.b:
        PhiSpyModules.message("FATAL: Please provide either -p or -b (or both) so we can make at least one output. Use -h for more help",
                "RED", "stderr")
        sys.exit(0)

    phageout = None
    bactout = None
    if args.p:
        phageout = open(args.p, 'w')
    if args.b:
        bactout = open(args.b, 'w')

    # read the coordinates
    phages = {}
    with open(args.c, 'r') as f:
        for l in f:
            p = l.strip().split("\t")
            if p[1] not in phages:
                phages[p[1]] = []
            phages[p[1]].append([p[2], p[3], p[0]])


    try:
        if PhiSpyModules.is_gzip_file(args.g):
            handle = gzip.open(args.g, 'rt')
        else:
            handle = open(args.g, 'r')
    except IOError as e:
        PhiSpyModules.message(f"There was an error reading {args.g}: {e}", "RED", "stderr")
        sys.exit(20)

    s = SeqIO.parse(handle, "genbank")


    if not phages:
        PhiSpyModules.message("No phages were found in args.c", "BLUE", "stderr")
        if args.b:
            for seq in s:
                bactout.write(f">{seq.id}\n{seq.seq}\n")
    else:
        for seq in s:
            if seq.id in phages:
                start = 0
                for p in phages[seq.id]:
                    bactend = p[0]
                    if args.b:
                        bactout.write(f">{seq.id}_{start}_{bactend}\n{seq.seq[start:bactend]}")
                    if args.p:
                        phageout.write(f">{seq.id}_{p[0]}_{p[1]} [prophage {p[2]}]\n{seq.seq[p[0]:p[1]]}")
                    start = p[1]+1
                bactout.write(f">{seq.id}_{start}_{len(seq.seq)}\n{seq.seq[start:]}")
            else:
                bactout.write(f">{seq.id}\n{seq.seq}\n")


    handle.close()
    if args.p:
        phageout.close()
    if args.b:
        bactout.close()
"""
Use biopython to read a gzip compressed genbank file, and then for 
each CDS in each record, use the is_phage_func method 
to count which proteins are phage proteins
"""

import os
import sys
import gzip
import re
import argparse
from Bio import SeqIO
from roblib import bcolors
__author__ = 'Rob Edwards'

def is_phage_func(func):
    func = func.lower()
    func = func.replace('-', ' ')
    func = func.replace(',', ' ')
    func = func.replace('/', ' ')
    a = re.split(' ', func)
    if (
            'integrase'     in a or
            'phage'         in a or
            'lysin'         in a or
            'endolysin'     in a or
            'holin'         in a or
            'capsid'        in a or
            'tail'          in a or
            'bacteriophage' in a or
            'prophage'      in a or
            'portal'        in a or
            'terminase'     in a or
            'tapemeasure'   in a or
            'baseplate'     in a or
            'virion'        in a or
            'antirepressor' in a or
            'excisionase'   in a or
            re.search(r"\b%s\b" % "tape measure", func) or
            re.search(r"\b%s\b" % "Cro-like repressor", func) or
            re.search(r"\b%s\b" % "CI-like repressor", func) or
            re.search(r"\b%s\b" % "rIIA lysis", func) or
            re.search(r"\b%s\b" % "rI lysis", func) or
            re.search(r"\b%s\b" % "rIIB lysis", func) or
            re.search(r"\b%s\b" % "base plate", func) or
            ("head" in a and "decoration" in a) or
            ("helix" in a and "turn" in a) or
            "HNH endonuclease" in func or
            "single stranded dna binding protein" in func
    ):
        return True
    return False





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', help='genbank file with prophage sequences', required=True)
    parser.add_argument('-d', help='include genbank definition line in outptu', action='store_true')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()


    opener = open
    if args.f.endswith('.gz'):
        opener = gzip.open

    handle = opener(args.f, 'rt')
    for record in SeqIO.parse(handle, "genbank"):
        count = 0
        phage_count = 0
        if args.v:
            print(f"Processing record: {record.id}", file=sys.stderr)
        for feature in record.features:
            if feature.type == 'CDS':
                count += 1
                if 'product' in feature.qualifiers:
                    func = feature.qualifiers['product'][0]
                    if is_phage_func(func):
                        phage_count += 1
                    if args.v:
                        print(f"function: |{func}| phage: {is_phage_func(func)}", file=sys.stderr)
                elif 'function' in feature.qualifiers:
                    func = feature.qualifiers['function'][0]
                    if is_phage_func(func):
                        phage_count += 1
                        if args.v:
                            print(f"Phage function found: {func}")
        if args.d:
            print(f"{bcolors.OKBLUE}{record.description}{bcolors.ENDC}", file=sys.stderr)
            print(f"{args.f}\t{record.description}\t{record.id}\t{phage_count}\t{count}")
        else:
            print(f"{args.f}\t{record.id}\t{phage_count}\t{count}")
    handle.close()

"""
Merge the information in prophage_information.tsv and test_data.tsv

We keep these separate on purpose, but sometimes you want to merge them all
"""

import os
import sys
import argparse
from PhiSpyModules.log_and_message import message

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-i', help='prophage information tsv file (default="prophage_information.tsv")', default='prophage_information.tsv')
    parser.add_argument('-t', help='test data file (default="test_data.tsv")', default='test_data.tsv')
    parser.add_argument('-o', help='output file (default="merged_prophage_information.tsv")', default='merged_prophage_information.tsv')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    if not os.path.exists(args.i):
        message(f"FATAL: {args.i} does not exist", "RED", "stderr")
        sys.exit(-1)
    if not os.path.exists(args.t):
        message(f"FATAL: {args.t} does not exist", "RED", "stderr")
        sys.exit(-1)
    with open(args.i, 'r') as i, open(args.t, 'r') as t, open(args.o, 'w') as out:
        for x,y in zip(i,t):
            p = x.split("\t")
            q = y.split("\t")
            assert(p[0] == q[0])
            p.insert(6, y.strip())
            out.write("\t".join(p))

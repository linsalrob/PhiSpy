"""
Create a prophage.tbl file from a phispy directory that does not contain one.
"""

import os
import sys
import argparse
from evaluation import  make_prophage_tbl

def make_new_prophage_tbl(phispydir):
    """
    Make a new prophage table
    :param phispydir: the directory to read the input and create the output
    :return: nothing
    """

    make_prophage_tbl(os.path.join(phispydir, 'prophage_tbl.txt'), os.path.join(phispydir, 'prophage.tbl'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="create a prophage.tbl file for a PhiSpy directory")
    parser.add_argument('-d', help='phispy directory')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    make_new_prophage_tbl(args.d)

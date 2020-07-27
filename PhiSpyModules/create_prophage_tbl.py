#!/usr/bin/env python
"""
Create a prophage.tbl file from a phispy directory that does not contain one.
"""

import os
import sys
import argparse

INSTALLATION_DIR = os.path.dirname(os.path.dirname(__file__))
sys.path.append(INSTALLATION_DIR)

from .writers import  prophage_measurements_to_tbl

def main(phispydir):
    """
    Make a new prophage table
    :param phispydir: the directory to read the input and create the output
    :return: nothing
    """
    prophage_measurements_to_tbl(os.path.join(phispydir, 'prophage_tbl.txt'), os.path.join(phispydir, 'prophage.tbl'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="create a prophage.tbl file for a PhiSpy directory")
    parser.add_argument('-d', help='phispy directory')
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()
    main(args.d)

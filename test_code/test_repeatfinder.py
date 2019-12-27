"""
Test the implementation of the repeatfinder extension
"""

import os
import sys
import argparse


import PhiSpyRepeatFinder
import pprint

s = "TTTTTTTTTTTTagcaTTTTTTTTTTTT"
print(f"s: {s}")
r = PhiSpyRepeatFinder.repeatFinder(s, 0)
pp = pprint.PrettyPrinter(indent=4)
pp.pprint(r)

with open("test.fna", 'r') as f:
    seqid = f.readline().strip()
    fna = f.readline().strip()
    r = PhiSpyRepeatFinder.repeatFinder(fna, 0)
    pp.pprint(r)


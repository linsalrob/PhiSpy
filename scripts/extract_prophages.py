#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import io
import sys
import re
import argparse
import textwrap
from argparse import RawTextHelpFormatter

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def main(args):
	with open(args.infile) as fp:
		locus = None
		dna = ''
		in_pp = False
		in_features = False
		for line in fp:
			if line.startswith('LOCUS'):
				locus = line.split()[1]
				dna = ''
				in_pp = False
				in_features = False
			elif line.startswith('FEATURES'):
				in_features = True
				for pp in args.coords.get(locus, []):
					args.coords[locus][pp].file.write(line)
					args.coords[locus][pp].file.write('     source          1..')
					args.coords[locus][pp].file.write(str(1+args.coords[locus][pp].right-args.coords[locus][pp].left))
					args.coords[locus][pp].file.write('\n')
			elif line.startswith('ORIGIN'):
				in_features = False
				dna = '\n'
			elif line.startswith('//'):
				dna = dna.replace('\n', '')
				for pp in args.coords.get(locus, []):
						args.coords[locus][pp].file.write('ORIGIN')
						i = 0
						for block in textwrap.wrap(dna[ args.coords[locus][pp].left-1 : args.coords[locus][pp].right ], 10):
							if(i%60 == 0):
								args.coords[locus][pp].file.write('\n')
								args.coords[locus][pp].file.write(str(i+1).rjust(9))
								args.coords[locus][pp].file.write(' ')
								args.coords[locus][pp].file.write(block.lower())
							else:
								args.coords[locus][pp].file.write(' ')
								args.coords[locus][pp].file.write(block.lower())
							i += 10
						args.coords[locus][pp].file.write('\n')
						args.coords[locus][pp].file.write('//\n')
			elif in_features:
				if not line.startswith('      '):
					in_pp = False
					left = min(map(int, re.findall(r"\d+", line)))
					right = max(map(int, re.findall(r"\d+", line)))
					for pp in args.coords.get(locus, []):
						if args.coords[locus][pp].left <= left and left <= args.coords[locus][pp].right and right <= args.coords[locus][pp].right:
							offset = args.coords[locus][pp].left - 1
							for match in re.findall(r"\d+", line):
								line = line.replace(match, str(int(match)-offset), 1)
							args.coords[locus][pp].file.write(line)
							in_pp = pp
				elif in_pp:
					args.coords[locus][in_pp].file.write(line)
			elif dna:
				line = line[10:].replace(' ','')
				dna += line.upper()
			if not in_features and not dna:
				for pp in args.coords.get(locus, []):
					if line.startswith('LOCUS'):
						match = re.findall(r"\d+ bp", line)[0]
						replacement = str( 1 + args.coords[locus][pp].right - args.coords[locus][pp].left ) + " bp"
						line = line.replace(match, replacement.rjust(len(match)))
						args.coords[locus][pp].file.write(line)
					else:
						args.coords[locus][pp].file.write(line)



if __name__ == '__main__':
	usage = 'extract_prophages.py [-opt1, [-opt2, ...]] -c coordinates infile'
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file in genbank format')
	parser.add_argument('-c', '--coordinates', type=is_valid_file, required=True, help='prophage coordinates')
	#parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
	parser.add_argument('--ids', action="store", help=argparse.SUPPRESS)
	args = parser.parse_args()

	args.coords = dict()
	with open(args.coordinates) as fp:
		for line in fp:
			pp, locus, left, right, *_ = line.split('\t')
			obj = lambda: None
			obj.left = int(left)
			obj.right = int(right)
			obj.file = open(os.path.splitext(args.infile)[0] + "_" + pp + '.gb', 'w')
			args.coords.setdefault(locus,{})[pp] = obj
	main(args)




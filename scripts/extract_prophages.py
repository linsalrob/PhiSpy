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
				for pp in args.coords:
					if args.coords[pp][0] == locus:
						args.out[pp].write(line)
						args.out[pp].write('     source          1..')
						args.out[pp].write(str(1+args.coords[pp][2]-args.coords[pp][1]))
						args.out[pp].write('\n')
			elif line.startswith('ORIGIN'):
				in_features = False
				dna = '\n'
			elif line.startswith('//'):
				dna = dna.replace('\n', '')
				for pp in args.coords:
					if args.coords[pp][0] == locus:
						args.out[pp].write('ORIGIN\n')
						i = 0
						for block in textwrap.wrap(dna[ args.coords[pp][1]-1 : args.coords[pp][2] ], 10):
							if(i%60 == 0):
								args.out[pp].write('\n')
								args.out[pp].write(str(i+1).rjust(9))
								args.out[pp].write(' ')
								args.out[pp].write(block.lower())
							else:
								args.out[pp].write(' ')
								args.out[pp].write(block.lower())
							i += 10
						args.out[pp].write('\n')
						args.out[pp].write('//\n')
			elif in_features:
				if not line.startswith('      '):
					in_pp = False
					left = min(map(int, re.findall(r"\d+", line)))
					for pp in args.coords:
						if args.coords[pp][0] == locus and args.coords[pp][1] <= left and left <= args.coords[pp][2]:
							offset = args.coords[pp][1]-1
							for match in re.findall(r"\d+", line):
								line = line.replace(match, str(int(match)-offset))
							args.out[pp].write(line)
							in_pp = pp
				elif in_pp:
					args.out[in_pp].write(line)
			elif dna:
				line = line[10:].replace(' ','')
				dna += line.upper()
			if not in_features and not dna:
				for pp in args.out:
					if args.coords[pp][0] != locus:
						pass
					elif line.startswith('LOCUS'):
						match = re.findall(r"\d+ bp", line)[0]
						replacement = str(1+args.coords[pp][2]-args.coords[pp][1]) + " bp"
						line = line.replace(match, replacement.rjust(len(match)))
						args.out[pp].write(line)
					else:
						args.out[pp].write(line)




if __name__ == '__main__':
	usage = 'extract_prophages.py [-opt1, [-opt2, ...]] -c coordinates infile'
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file in genbank format')
	parser.add_argument('-c', '--coordinates', type=is_valid_file, required=True, help='prophage coordinates')
	#parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
	parser.add_argument('--ids', action="store", help=argparse.SUPPRESS)
	args = parser.parse_args()

	args.coords = dict()
	args.out = dict()
	with open(args.coordinates) as fp:
		for line in fp:
			cols = line.split('\t')
			args.coords[cols[0]] = tuple([cols[1], int(cols[2]), int(cols[3])])
			args.out[cols[0]] = open(os.path.splitext(args.infile)[0] + "_" + cols[0] + '.gb', 'w')
	main(args)




#!/usr/bin/python3
__author__ = 'Przemek Decewicz'

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from glob import glob
from os import makedirs, path
from sys import argv
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np
from io import TextIOWrapper
import math


def plot_stats(infile, outfile):

    """Reads test/training set and plots all identified stats.

    Stats are slightly transformed to retained a visible scale.
    Two types of plots are provided:
    - transformed stats

    Parameters
    ----------
    infile: str
        Path to input file.
    outfile: str
        Path to resulting PNG file with plots.
    """

    # read input file
    with open(infile) as inf:
        colnames = inf.readline().strip().split('\t')
    data = np.genfromtxt(infile, delimiter="\t", filling_values=1, dtype=np.float64, skip_header=1)

    for i, name in enumerate(colnames):
        if name == 'orf_length_med':
            data[:, i] = data[:, i] / 50
        elif name == 'shannon_slope':
            data[:, i] = data[:, i] * 200
        elif name == 'at_skew':
            data[:, i] = data[:, i] * 2
        elif name == 'gc_skew':
            data[:, i] = data[:, i] * 2
        elif name == 'max_direction':
            data[:, i] = data[:, i] / 3
        elif name == 'phmms':
            data[:, i] = data[:, i] * 2
        elif name == 'status':
            data[:, i] = data[:, i] * 20

    # make a plot
    fig, ax = plt.subplots(figsize=(18, 4.5), dpi = 150)

    plt.plot(data, '-', linewidth=.8, alpha = 0.9)
    plt.legend(colnames, loc='lower center', bbox_to_anchor=(0.5,-0.17), ncol = len(colnames))
    plt.margins(x=0.01)
    plt.subplots_adjust(left=0.03, right=0.99, top=0.9, bottom=0.15)
    plt.title(path.basename(infile))

    plt.savefig(outfile)
    plt.close()



def main():
    args = ArgumentParser(prog = 'plot_trainSets_stats.py',
                          description = 'Plots PhiSpy\'s training/test sets statistics.',
                          epilog = 'Example usage:\npython3 scripts/plot_trainSets_stats.py -d PhiSpyModules/data -o PhiSpyModules/data/trainSets_stats ',
                          formatter_class = RawDescriptionHelpFormatter)

    args.add_argument('-i', '--infile',
                      type = str,
                      help = 'Path to input GenBank file.')

    args.add_argument('-d', '--indir',
                      type = str,
                      help = 'Path to input directory with multiple GenBank files.')

    args.add_argument('-s', '--suffix',
                      type = str,
                      help = 'Suffix that will be added to input file name.')

    args.add_argument('-o', '--outdir',
                      type = str,
                      help = 'Path to output directory.',
                      required = True)

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
    elif args.indir:
        infiles = glob(path.join(args.indir, '*.txt'))
    else:
        infiles = [args.infile]

    # Create output directory
    if not path.isdir(args.outdir): makedirs(args.outdir)

    # Process all input files
    for infile in infiles:
        plot_file_name = f'{path.basename(infile).rsplit(".", 1)[0]}.{args.suffix}.png'
        plot_file = path.join(args.outdir, plot_file_name)
        plot_stats(infile, plot_file)
        print(f'Done with plot: {plot_file}')

if __name__ == '__main__':
    main()

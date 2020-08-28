#!/usr/bin/python3
__author__ = 'Przemek Decewicz'

import gzip
import sys
import pkg_resources
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from Bio import SeqIO
from compare_predictions_to_phages import genbank_seqio
from os import makedirs, path
from PhiSpyModules import log_and_message, is_gzip_file

from argparse import ArgumentParser
from os import makedirs, path
from sys import argv

def read_prophage_coordinates(infile):
    """
    Reads PhiSpy's prophage coordinates file.
    :param infile: path to the file
    :return prophages: dictionary of replicons ids and prophage coordinates
    """

    prophages = {}

    with open(infile) as inf:
        for line in inf:
            pp, replicon_id, start, end, r1_start, r1_end, r2_start, r2_end, r1_seq, r2_seq, note = line.strip().split('\t')

            # adjust coords to Python's indexing
            coords = (int(start) - 1, int(end) - 1)

            try:
                prophages[replicon_id][coords] = pp
            except KeyError:
                prophages[replicon_id] = {coords: pp}

    return prophages

def read_prophage_table(infile):
    """
    Reads tab-delimited table with columns for:
    1 - path to GenBank file,
    2 - replicon id,
    3 - prophage start coordinate,
    4 - prophage end coordinate,
    5 (optional) - prophage name (if not provided pp1, pp2, etc.
    will be assigned for each file)
    :param infile: path to the file
    :return prophages: dictionary of GenBank file(s) paths, replicons and prophages coordinates
    """

    prophages = {}
    with open(infile) as inf:
        for line in inf:
            if line.startswith('#'): continue

            line = line.strip().split('\t')
            file_path, replicon_id, start, end = line[:4]

            # adjust coords to Python's indexing
            coords = (int(start) - 1, int(end) - 1)

            if start > end:
                start, end = end, start

            if file_path in prophages:
                pp_cnt = len(prophages[file_path]) + 1
            else:
                pp_cnt = 1

            pp = line[4] if len(line) == 5 else  f"pp{pp_cnt}"

            try:
                prophages[file_path][replicon_id][coords] = pp
            except KeyError:
                try:
                    prophages[file_path][replicon_id] = {coords: pp}
                except KeyError:
                    prophages[file_path] = {replicon_id: {coords: pp}}

    return prophages

def feature_within_prophage(feature, pps):
    """
    Checks whether a feature is located within any of prophage regions.
    :param feature: SeqFeature object
    :param pps: dictionary of prophage coordinates tuples
    :return pp: prophage name
    """

    f_start = int(feature.location.start)
    f_end = int(feature.location.end)

    for pp_coords in pps:
        if f_start >= pp_coords[0] and f_end <= pp_coords[1]:
            return pps[pp_coords]

    return None


def mark_prophage_features(prophages, outdir, gzip_files, ungzip_files):
    """
    Marks features within prophage regions with "is_phage='1'" qualifier.
    :param prophages: dictionary of GenBank file(s) paths, replicons and prophages coordinates
    :param outdir: output directory to write updated GenBank file(s)
    :param gzip_files: gzip all processed input files
    :param ungzip_files: ungzip all processed input files
    :return:
    """

    for gbkf, pps in prophages.items():
        log_and_message(f"Updating {path.basename(gbkf)}.", c="YELLOW", stderr = True)
        records = list(genbank_seqio(gbkf))
        cnt_all = 0
        cnt_cds = 0
        for record in records:
            if record.id in pps:
                for feature in record.features:
                    pp = feature_within_prophage(feature, pps[record.id])
                    if pp:
                        feature.qualifiers['is_phage'] = ['1']
                        feature.qualifiers['color'] = ['6']
                        try:
                            feature.qualifiers['note'].append(pp)
                        except KeyError:
                            feature.qualifiers['note'] = [pp]

                        cnt_all += 1
                        if feature.type == 'CDS': cnt_cds += 1

        log_and_message(f"-- Marked {cnt_all} features, including {cnt_cds} CDSs, with 'is_phage' qualifier.", stderr=True)

        outfile = path.join(outdir, path.basename(gbkf))
        if gzip_files:
            if not outfile.endswith('.gz'):
                outfile += '.gz'
            handle = gzip.open(outfile, 'wt')
        elif ungzip_files:
            if outfile.endswith('.gz'):
                outfile = outfile[:-3]
            handle = open(outfile, 'w')
        elif is_gzip_file(gbkf):
            handle = gzip.open(outfile, 'wt')
        else:
            handle = open(outfile, 'w')

        log_and_message(f"Writing {outfile}", c="GREEN", stderr=True)
        SeqIO.write(records, handle, 'genbank')


def main():
    args = ArgumentParser(prog = "mark_prophage_features.py",
                          description = "Updates GenBank files based on prophage_coordinates.tsv file or other tab-delimited table and marks features within indicated coordinates with 'is_phage=1' qualifier.",
                          epilog = "Example usages:\npython3 scripts/mark_prophage_features.py --genbank infile.gbk --outdir updated_genbanks --ppcoords PhiSpy_output/prophage_coordinates.tsv\npython3 scripts/mark_prophage_features.py --outdir updated_genbanks --table my_confirmed_predictions.tsv",
                          formatter_class = RawDescriptionHelpFormatter)

    args.add_argument('-g', '--genbank',
                      type = str,
                      help = 'Path to input GenBank file.')

    args.add_argument('-o', '--outdir',
                      type = str,
                      help = 'Path to output directory to write updated GenBank(s).',
                      required = True)

    args.add_argument('-c', '--ppcoords',
                      type = str,
                      help = 'Path to prophage_coordinates.tsv file.')

    args.add_argument('-t', '--table',
                      type = str,
                      help = 'Path to tab-delimited file with confirmed prophage regions to mark. The file has to have the following columns: 1 - path to GenBank file, 2 - replicon id, 3 - prophage start coordinate, 4 - prophage end coordinate, 5 (optional) - prophage name (if not provided pp1, pp2, etc. will be assigned for each file)')

    args.add_argument('--gzip_files',
                      action = 'store_true',
                      default = False,
                      help = 'Gzip all output files. \'.gz\' extension will be added if missing. [Default: %(default)s]')

    args.add_argument('--ungzip_files',
                      action = 'store_true',
                      default = False,
                      help = 'Ungzip all output files. \'.gz\' extension will be removed if present. [Default: %(default)s]')

    if len(argv[1:]) == 0:
        args.print_help()
        args.exit()

    try:
        args = args.parse_args()
    except:
        args.exit()

    if not path.isdir(args.outdir): makedirs(args.outdir)

    if args.genbank and args.outdir and args.ppcoords:
        log_and_message(f"Marking prophage features for a single GenBank file based on PhiSpy's prophage coordinates file.", c="GREEN", stderr=True)

        prophages = {args.genbank: read_prophage_coordinates(args.ppcoords)}

    elif args.outdir and args.table:
        log_and_message(f"Marking prophage features for GenBank file(s) based on prophage table.", c="GREEN", stderr=True)

        prophages = read_prophage_table(args.table)

    else:
        log_and_message(f"Incorrect input data.", c="RED", stderr=True)
        log_and_message(f"Use --genbank, --outdir and --ppcords or --table and --outdir.", c="PINK", stderr=True)

    mark_prophage_features(prophages, args.outdir, args.gzip_files, args.ungzip_files)

    log_and_message(f"Done. You can now use output file(s) for training PhiSpy.", c="GREEN", stderr=True)

if __name__ == '__main__':
    main()

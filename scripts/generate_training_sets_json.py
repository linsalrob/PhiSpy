#!/usr/bin/env python3
"""
Generate a JSON file describing the available PhiSpy training sets.

This script reads the training genome list used by `phispy --list short`
and produces a structured JSON file suitable for use by PhiSpyWeb.
"""

import argparse
import json
import os
import re
import sys

try:
    from importlib.resources import files, as_file
except ImportError:
    from importlib_resources import files, as_file


def label_from_genome_file(genome_file):
    """
    Derive a human-readable label from a genome filename.

    Removes the `.gb.gz` extension and replaces underscores with spaces.

    :param genome_file: filename such as ``Escherichia_coli_O157-H7_EDL933.gb.gz``
    :return: label such as ``Escherichia coli O157-H7 EDL933``
    """
    label = re.sub(r'\.gb\.gz$', '', genome_file)
    label = label.replace('_', ' ')
    return label


def read_training_sets():
    """
    Read training set metadata from the package data file.

    :return: list of dicts with keys ``value``, ``count``, ``genomeFile``, ``label``
    """
    try:
        data_file = files('PhiSpyModules').joinpath('data/trainingGenome_list.txt')
        file_ctx = as_file(data_file)
    except (FileNotFoundError, TypeError, AttributeError) as exc:
        print(f"ERROR: Cannot find data/trainingGenome_list.txt: {exc}", file=sys.stderr)
        sys.exit(10)

    entries = []
    with file_ctx as path, open(path, 'rb') as f:
        for raw_line in f:
            line = raw_line.decode().strip()
            if not line:
                continue
            parts = re.split(r'\t', line)
            # columns: index, filename, genome(s), count
            if len(parts) < 4:
                continue
            train_file = parts[1]
            genomes_field = parts[2]
            try:
                count = int(parts[3])
            except ValueError:
                continue

            # Skip the generic test set (index 0)
            if train_file.startswith('testSet_'):
                continue

            first_genome = genomes_field.split(';')[0]
            entries.append({
                'value': f'data/{train_file}',
                'count': count,
                'genomeFile': first_genome,
                'label': label_from_genome_file(first_genome),
            })

    entries.sort(key=lambda e: e['label'])
    return entries


def get_phispy_version():
    """Return the installed PhiSpy version string."""
    try:
        from PhiSpyModules.version import __version__
        return __version__
    except (ImportError, AttributeError):
        return 'unknown'


def generate(output_path):
    """
    Generate the training-sets JSON file and write it to *output_path*.

    :param output_path: destination file path (created if necessary)
    """
    training_sets = read_training_sets()
    payload = {
        'phispyVersion': get_phispy_version(),
        'schemaVersion': 1,
        'trainingSets': training_sets,
    }

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=2)
        fh.write('\n')

    print(f"Wrote {len(training_sets)} training sets to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate a JSON file describing the available PhiSpy training sets.'
    )
    parser.add_argument(
        '--output', '-o',
        default='release-assets/training-sets.json',
        help='Output file path [default: %(default)s]',
    )
    args = parser.parse_args()
    generate(args.output)


if __name__ == '__main__':
    main()

"""
Tests for scripts/generate_training_sets_json.py
"""

import json
import os
import sys
import tempfile
import unittest

# Ensure the scripts directory is importable
_SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'scripts')
sys.path.insert(0, os.path.abspath(_SCRIPTS_DIR))

from generate_training_sets_json import label_from_genome_file, read_training_sets, generate


class LabelTransformTest(unittest.TestCase):
    """Tests for label_from_genome_file()."""

    def test_removes_gb_gz_extension(self):
        self.assertEqual(
            label_from_genome_file('Escherichia_coli_O157-H7_EDL933.gb.gz'),
            'Escherichia coli O157-H7 EDL933',
        )

    def test_replaces_underscores_with_spaces(self):
        self.assertEqual(
            label_from_genome_file('Bacillus_subtilis_str_168.gb.gz'),
            'Bacillus subtilis str 168',
        )

    def test_hyphen_preserved(self):
        label = label_from_genome_file('Staphylococcus_aureus_subsp._aureus_Mu50.gb.gz')
        self.assertNotIn('.gb.gz', label)
        self.assertNotIn('_', label)
        self.assertIn('subsp.', label)

    def test_no_trailing_extension_left(self):
        label = label_from_genome_file('Some_Organism_Name.gb.gz')
        self.assertFalse(label.endswith('.gb.gz'))
        self.assertFalse(label.endswith('.gz'))

    def test_file_without_gb_gz(self):
        """A filename that doesn't end in .gb.gz should have underscores replaced."""
        label = label_from_genome_file('My_Genome.fasta')
        self.assertEqual(label, 'My Genome.fasta')


class ReadTrainingSetsTest(unittest.TestCase):
    """Tests for read_training_sets()."""

    def setUp(self):
        self.entries = read_training_sets()

    def test_returns_list(self):
        self.assertIsInstance(self.entries, list)

    def test_at_least_one_entry(self):
        self.assertGreater(len(self.entries), 0)

    def test_each_entry_has_required_keys(self):
        required = {'value', 'count', 'genomeFile', 'label'}
        for entry in self.entries:
            with self.subTest(entry=entry['value']):
                self.assertEqual(set(entry.keys()), required)

    def test_value_starts_with_data_prefix(self):
        for entry in self.entries:
            with self.subTest(value=entry['value']):
                self.assertTrue(entry['value'].startswith('data/'))

    def test_count_is_integer(self):
        for entry in self.entries:
            with self.subTest(value=entry['value']):
                self.assertIsInstance(entry['count'], int)

    def test_genome_file_is_string(self):
        for entry in self.entries:
            with self.subTest(value=entry['value']):
                self.assertIsInstance(entry['genomeFile'], str)
                self.assertTrue(entry['genomeFile'])

    def test_label_has_no_gb_gz(self):
        for entry in self.entries:
            with self.subTest(value=entry['value']):
                self.assertNotIn('.gb.gz', entry['label'])

    def test_label_has_no_underscores(self):
        for entry in self.entries:
            with self.subTest(value=entry['value']):
                self.assertNotIn('_', entry['label'])

    def test_sorted_by_label(self):
        labels = [e['label'] for e in self.entries]
        self.assertEqual(labels, sorted(labels))

    def test_test_set_excluded(self):
        for entry in self.entries:
            with self.subTest(value=entry['value']):
                self.assertNotIn('testSet_', entry['value'])

    def test_ecoli_entry_present(self):
        ecoli_values = [e['value'] for e in self.entries if 'Ecoli' in e['value']]
        self.assertTrue(ecoli_values, "Expected at least one Ecoli training set")

    def test_ecoli_count(self):
        ecoli = next((e for e in self.entries if e['value'] == 'data/trainSet_Ecoli.txt'), None)
        self.assertIsNotNone(ecoli, "trainSet_Ecoli.txt not found")
        self.assertEqual(ecoli['count'], 4)
        self.assertEqual(ecoli['genomeFile'], 'Escherichia_coli_O157-H7_EDL933.gb.gz')
        self.assertEqual(ecoli['label'], 'Escherichia coli O157-H7 EDL933')


class GenerateJsonFileTest(unittest.TestCase):
    """Tests for generate() – JSON file structure and content."""

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.output_path = os.path.join(self.tmp_dir, 'training-sets.json')
        generate(self.output_path)
        with open(self.output_path, encoding='utf-8') as fh:
            self.data = json.load(fh)

    def test_file_is_valid_json(self):
        self.assertIsNotNone(self.data)

    def test_top_level_keys(self):
        self.assertIn('phispyVersion', self.data)
        self.assertIn('schemaVersion', self.data)
        self.assertIn('trainingSets', self.data)

    def test_schema_version_is_one(self):
        self.assertEqual(self.data['schemaVersion'], 1)

    def test_phispy_version_is_string(self):
        self.assertIsInstance(self.data['phispyVersion'], str)
        self.assertTrue(self.data['phispyVersion'])

    def test_training_sets_is_list(self):
        self.assertIsInstance(self.data['trainingSets'], list)

    def test_training_sets_not_empty(self):
        self.assertGreater(len(self.data['trainingSets']), 0)

    def test_each_training_set_structure(self):
        for ts in self.data['trainingSets']:
            with self.subTest(value=ts.get('value')):
                self.assertIn('value', ts)
                self.assertIn('count', ts)
                self.assertIn('genomeFile', ts)
                self.assertIn('label', ts)
                self.assertIsInstance(ts['count'], int)

    def test_training_sets_sorted_by_label(self):
        labels = [ts['label'] for ts in self.data['trainingSets']]
        self.assertEqual(labels, sorted(labels))


if __name__ == '__main__':
    unittest.main()

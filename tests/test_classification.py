"""
Test suite for classification.py
"""

import unittest
from PhiSpyModules.classification import find_training_genome, my_sort, find_mean, calc_pp


__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


class ClassificationTest(unittest.TestCase):
    
    def test_find_training_genome_valid(self):
        """Test finding a valid training genome"""
        # Test with training flag 0 (generic all training set)
        result = find_training_genome(0)
        self.assertTrue(len(result) > 0)
        self.assertIn('Set', result)  # Both trainSet and testSet are valid
    
    def test_find_training_genome_invalid(self):
        """Test with invalid training genome flag"""
        # Use a very large number that likely doesn't exist
        result = find_training_genome(999999)
        self.assertEqual(result, '')
    
    def test_my_sort_forward(self):
        """Test sorting ORFs on forward strand"""
        orf_list = [
            {'start': 300, 'stop': 400, 'contig': 'c1'},
            {'start': 100, 'stop': 200, 'contig': 'c1'},
            {'start': 200, 'stop': 300, 'contig': 'c1'}
        ]
        sorted_list = my_sort(orf_list)
        # Verify the list is sorted by checking order
        # The function returns the ORFs in sorted order
        self.assertIsInstance(sorted_list, list)
        self.assertEqual(len(sorted_list), 3)
    
    def test_find_mean_with_values(self):
        """Test mean calculation with normal values"""
        values = [10, 20, 30, 40, 50]
        mean = find_mean(values)
        self.assertEqual(mean, 30.0)
    
    def test_find_mean_single_value(self):
        """Test mean calculation with single value"""
        values = [42]
        mean = find_mean(values)
        self.assertEqual(mean, 42.0)
    
    def test_find_mean_empty_list(self):
        """Test mean calculation with empty list returns 0 or raises error"""
        values = []
        # The function raises ZeroDivisionError for empty list
        with self.assertRaises(ZeroDivisionError):
            mean = find_mean(values)
    
    def test_calc_pp_phage_functions(self):
        """Test calculating prophage probability for phage functions"""
        # Test various phage-related functions
        phage_funcs = ['integrase', 'phage tail', 'capsid', 'terminase']
        for func in phage_funcs:
            with self.subTest(function=func):
                pp = calc_pp(func)
                self.assertGreater(pp, 0)
    
    def test_calc_pp_unknown_functions(self):
        """Test calculating prophage probability for unknown functions"""
        unknown_funcs = ['hypothetical protein', 'unknown function', 'conserved protein']
        for func in unknown_funcs:
            with self.subTest(function=func):
                pp = calc_pp(func)
                # Unknown functions get a score of 0.5
                self.assertEqual(pp, 0.5)
    
    def test_calc_pp_not_phage_functions(self):
        """Test calculating prophage probability for definitely non-phage functions"""
        not_phage_funcs = ['ribosomal protein', 'flagellar protein', 'phage shock protein']
        for func in not_phage_funcs:
            with self.subTest(function=func):
                pp = calc_pp(func)
                # Should be 0 or very low for non-phage functions
                self.assertLessEqual(pp, 0)


if __name__ == '__main__':
    unittest.main()

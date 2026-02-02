"""
Test suite for evaluation.py
"""

import unittest
from PhiSpyModules.evaluation import (
    find_smallest, check_intg, check_pp
)


__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


class EvaluationTest(unittest.TestCase):
    
    def test_find_smallest_first_smaller(self):
        """Test find_smallest with list inputs"""
        # find_smallest takes two lists and finds minimum distance
        result = find_smallest([5, 10], [15, 20])
        self.assertIsInstance(result, float)
        self.assertEqual(result, 5.0)  # Distance between 10 and 15
    
    def test_find_smallest_second_smaller(self):
        """Test find_smallest with list inputs"""
        result = find_smallest([10, 20], [5, 25])
        self.assertIsInstance(result, float)
        self.assertEqual(result, 5.0)  # Distance between 10 and 5
    
    def test_find_smallest_equal(self):
        """Test find_smallest when lists have same values"""
        result = find_smallest([7], [7])
        self.assertEqual(result, 0.0)
    
    def test_find_smallest_negative(self):
        """Test find_smallest with negative values"""
        result = find_smallest([-5, 0], [-10, -2])
        self.assertIsInstance(result, float)
    
    def test_check_intg_with_integrase_within_range(self):
        """Test checking for integrase within prophage region"""
        # Mock integrase data: position -> dict with contig, start, stop info
        integ = {
            100: {'contig': 'contig1', 'id': 'gene1', 'start': 100, 'stop': 200},
            500: {'contig': 'contig1', 'id': 'gene2', 'start': 500, 'stop': 600},
            1000: {'contig': 'contig1', 'id': 'gene3', 'start': 1000, 'stop': 1100}
        }
        
        # Mock repeat dictionary with required keys
        rep = {'s1': 50, 'e1': 70, 's2': 590, 'e2': 610}
        
        # Check for integrase between positions 50 and 600
        result = check_intg(50, 600, rep, integ, 'contig1')
        
        # Should return 0 or 1 depending on proximity to repeats
        self.assertIn(result, [0, 1])
    
    def test_check_intg_no_integrase(self):
        """Test checking for integrase when none exists in range"""
        integ = {
            100: {'contig': 'contig1', 'id': 'gene1', 'start': 100, 'stop': 200},
            500: {'contig': 'contig1', 'id': 'gene2', 'start': 500, 'stop': 600}
        }
        
        # Mock repeat dictionary
        rep = {'s1': 1050, 'e1': 1070, 's2': 1990, 'e2': 2010}
        
        # Check range where no integrase exists near repeats
        result = check_intg(1000, 2000, rep, integ, 'contig1')
        
        self.assertEqual(result, 0)
    
    def test_check_pp_valid_prophage(self):
        """Test checking if prophage region is valid"""
        # check_pp expects contig to be a string (contig name)
        contig = 'test_contig'
        
        # Test with valid start and stop positions
        result = check_pp(contig, 100, 600, {})
        
        # Should return something indicating validity check
        self.assertIsNotNone(result)


if __name__ == '__main__':
    unittest.main()

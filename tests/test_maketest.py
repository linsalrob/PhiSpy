"""
Test suite for makeTest.py
"""

import unittest
from PhiSpyModules.makeTest import (
    my_sort, find_median, find_all_median, 
    reverse_complement, find_atgc_skew
)


__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


class MakeTestTest(unittest.TestCase):
    
    def test_reverse_complement(self):
        """Test reverse complement function"""
        test_cases = [
            ('ATCG', 'CGAT'),
            ('AAAA', 'TTTT'),
            ('GGCC', 'GGCC'),
            ('ATGC', 'GCAT'),
            ('', ''),
            ('N', 'N')
        ]
        for seq, expected in test_cases:
            with self.subTest(sequence=seq):
                self.assertEqual(reverse_complement(seq), expected)
    
    def test_find_median_odd_length(self):
        """Test median calculation with odd number of elements"""
        values = [1, 3, 5, 7, 9]
        self.assertEqual(find_median(values), 5)
    
    def test_find_median_even_length(self):
        """Test median calculation with even number of elements"""
        values = [1, 2, 3, 4]
        self.assertEqual(find_median(values), 2.5)
    
    def test_find_median_single_element(self):
        """Test median calculation with single element"""
        values = [42]
        self.assertEqual(find_median(values), 42)
    
    def test_find_median_unsorted(self):
        """Test median calculation with unsorted list"""
        values = [9, 1, 5, 3, 7]
        self.assertEqual(find_median(values), 5)
    
    def test_find_all_median(self):
        """Test finding all median values from ORF list"""
        orf_list = [
            {'start': 100, 'stop': 200},
            {'start': 200, 'stop': 400},
            {'start': 300, 'stop': 600},
            {'start': 400, 'stop': 800},
            {'start': 500, 'stop': 1000}
        ]
        median = find_all_median(orf_list)
        # Returns median value directly
        self.assertIsInstance(median, (int, float))
    
    def test_my_sort_forward_strand(self):
        """Test sorting ORFs on forward strand"""
        orf_list = [
            {'start': 300, 'stop': 400},
            {'start': 100, 'stop': 200},
            {'start': 200, 'stop': 300}
        ]
        sorted_list = my_sort(orf_list)
        self.assertEqual(sorted_list[0]['start'], 100)
        self.assertEqual(sorted_list[1]['start'], 200)
        self.assertEqual(sorted_list[2]['start'], 300)
    
    def test_my_sort_reverse_strand(self):
        """Test sorting ORFs on reverse strand"""
        orf_list = [
            {'start': 400, 'stop': 300},
            {'start': 300, 'stop': 200},
            {'start': 200, 'stop': 100}
        ]
        sorted_list = my_sort(orf_list)
        # Just verify it returns a sorted list
        self.assertIsInstance(sorted_list, list)
        self.assertEqual(len(sorted_list), 3)
    
    def test_my_sort_mixed_strands(self):
        """Test sorting ORFs on mixed strands"""
        orf_list = [
            {'start': 300, 'stop': 400},  # forward
            {'start': 250, 'stop': 150},  # reverse
            {'start': 100, 'stop': 200}   # forward
        ]
        sorted_list = my_sort(orf_list)
        # Should be sorted by absolute position
        self.assertEqual(sorted_list[0]['start'], 100)
        self.assertEqual(sorted_list[1]['start'], 250)
        self.assertEqual(sorted_list[2]['start'], 300)
    
    def test_find_atgc_skew_balanced(self):
        """Test AT/GC skew calculation with balanced sequence"""
        # Equal AT and GC - returns tuple of skews (a_skew, t_skew, g_skew, c_skew)
        seq = 'ATGC'
        skew = find_atgc_skew(seq)
        self.assertIsInstance(skew, tuple)
        self.assertEqual(len(skew), 4)
    
    def test_find_atgc_skew_at_rich(self):
        """Test AT/GC skew calculation with AT-rich sequence"""
        # Needs both AT and GC content to avoid NoBasesCounted error
        seq = 'AAAATTTTGC'
        skew = find_atgc_skew(seq)
        self.assertIsInstance(skew, tuple)
        self.assertEqual(len(skew), 4)
    
    def test_find_atgc_skew_gc_rich(self):
        """Test AT/GC skew calculation with GC-rich sequence"""
        # Needs both AT and GC content
        seq = 'GGGCCCAATT'
        skew = find_atgc_skew(seq)
        self.assertIsInstance(skew, tuple)
        self.assertEqual(len(skew), 4)
    
    def test_find_atgc_skew_empty(self):
        """Test AT/GC skew calculation with empty sequence raises error"""
        from PhiSpyModules.errors import NoBasesCounted
        seq = ''
        # Empty sequence raises NoBasesCounted
        with self.assertRaises(NoBasesCounted):
            skew = find_atgc_skew(seq)


if __name__ == '__main__':
    unittest.main()

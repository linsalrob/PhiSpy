"""
Test suite for protein_functions.py
"""

import unittest
from PhiSpyModules.protein_functions import (
    is_phage_func, is_not_phage_func, is_unknown_func, 
    is_mobile_element, downweighting_unknown_functions
)


__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


class ProteinFunctionsTest(unittest.TestCase):
    
    def test_is_phage_func_positive_cases(self):
        """Test that phage-related functions are correctly identified"""
        phage_functions = [
            'integrase',
            'phage tail protein',
            'lysin domain',
            'endolysin protein',
            'holin family',
            'capsid protein',
            'tail fiber',
            'bacteriophage protein',
            'prophage integrase',
            'portal protein',
            'terminase large subunit',
            'tapemeasure',
            'baseplate protein',
            'virion morphogenesis',
            'antirepressor protein',
            'excisionase protein',
            'tape measure protein',
            'head decoration protein',
            'helix turn helix domain'
        ]
        for func in phage_functions:
            with self.subTest(function=func):
                self.assertTrue(is_phage_func(func), f"Failed to identify '{func}' as phage function")
    
    def test_is_phage_func_negative_cases(self):
        """Test that non-phage functions are correctly rejected"""
        non_phage_functions = [
            'ribosomal protein',
            'DNA polymerase',
            'ATP synthase',
            'transcription factor',
            'hypothetical protein'
        ]
        for func in non_phage_functions:
            with self.subTest(function=func):
                self.assertFalse(is_phage_func(func), f"Incorrectly identified '{func}' as phage function")
    
    def test_is_not_phage_func_positive_cases(self):
        """Test functions that are definitely NOT phage-related"""
        not_phage_functions = [
            'phage shock protein',
            'trna synthase complex',
            'conjugal transfer protein',
            'conjugative plasmid',
            'flagella biosynthesis',
            'flagellar hook protein',
            'flagellin domain',
            'flagellum assembly',
            'ribosomal protein L1',
            'translation elongation factor EF-Tu',
            'secy protein',
            'summary phrase annotation',
            'dna binding domain protein',
            'abortive infection bacteriophage resistance protein'
        ]
        for func in not_phage_functions:
            with self.subTest(function=func):
                self.assertTrue(is_not_phage_func(func), f"Failed to identify '{func}' as NOT phage function")
    
    def test_is_not_phage_func_negative_cases(self):
        """Test that regular functions are not marked as definitely non-phage"""
        regular_functions = [
            'DNA helicase',
            'protease',
            'kinase'
        ]
        for func in regular_functions:
            with self.subTest(function=func):
                self.assertFalse(is_not_phage_func(func), f"Incorrectly identified '{func}' as definitely NOT phage")
    
    def test_is_unknown_func_positive_cases(self):
        """Test that unknown/hypothetical functions are correctly identified"""
        unknown_functions = [
            '',
            'hypothetical protein',
            'mobile element protein',
            'conserved protein of unknown function',
            'gene product XYZ',
            'interpro domain',
            'uncharacterized protein',
            'pseudogene fragment',
            'similar to protein ABC',
            'similarity to known proteins',
            'glimmer predicted gene',
            'unknown function',
            'complete genome sequence',
            'ensangp00000012345',
            'unnamed protein product',
            'expressed protein',
            'ortholog of gene ABC',
            'structural feature only',
            'cds_12345',
            'predicted by psort analysis',
            'AGR_C_1234',
            'EG:12345',
            'RIKEN cDNA clone',
            'lmo1234 protein',
            'lmo1234protein',
            'U123',
            'ABC|123',
            'orf123',
            'orf1',
            'predicted gene',
            'bh1234',
            'yabc'
        ]
        for func in unknown_functions:
            with self.subTest(function=func):
                self.assertTrue(is_unknown_func(func), f"Failed to identify '{func}' as unknown function")
    
    def test_is_unknown_func_negative_cases(self):
        """Test that known functions are not marked as unknown"""
        known_functions = [
            'DNA polymerase III',
            'ATP synthase subunit alpha',
            'ribosomal protein S1',
            'integrase',
            'DNA helicase RecQ'
        ]
        for func in known_functions:
            with self.subTest(function=func):
                self.assertFalse(is_unknown_func(func), f"Incorrectly identified '{func}' as unknown function")
    
    def test_is_mobile_element_positive_cases(self):
        """Test that mobile elements are correctly identified"""
        mobile_elements = [
            'mobile element protein',
            'transposon Tn5',
            'transposase IS1',
            'insertion sequence IS1',
            'insertion element IS1'
        ]
        for func in mobile_elements:
            with self.subTest(function=func):
                self.assertTrue(is_mobile_element(func), f"Failed to identify '{func}' as mobile element")
    
    def test_is_mobile_element_negative_cases(self):
        """Test that non-mobile elements are not identified as mobile elements"""
        non_mobile_elements = [
            'DNA polymerase',
            'integrase',
            'recombinase'
        ]
        for func in non_mobile_elements:
            with self.subTest(function=func):
                self.assertFalse(is_mobile_element(func), f"Incorrectly identified '{func}' as mobile element")
    
    def test_downweighting_unknown_functions(self):
        """Test downweighting of unknown functions in initial table"""
        # Create mock object with initial_tbl
        class MockArgs:
            def __init__(self):
                # initial_tbl format: [contig, function, start, stop, strand, 
                #                      length, median, shannon, phmm_score]
                # Note: The actual structure has 9 elements (index 8 is the score)
                self.initial_tbl = [
                    ['contig1', 'DNA polymerase', 100, 200, 1, 100, 50, 1.5, 0.8],
                    ['contig1', 'hypothetical protein', 300, 400, 1, 100, 50, 1.5, 0.7],
                    ['contig1', 'integrase', 500, 600, 1, 100, 50, 1.5, 0.6]
                ]
        
        args = MockArgs()
        result = downweighting_unknown_functions(args)
        
        # First entry (known function) should keep original score
        self.assertEqual(result[0][8], 0.8)
        # Second entry (unknown function) should have score downweighted to 0.5
        self.assertEqual(result[1][8], 0.5)
        # Third entry (known function) should keep original score
        self.assertEqual(result[2][8], 0.6)


if __name__ == '__main__':
    unittest.main()

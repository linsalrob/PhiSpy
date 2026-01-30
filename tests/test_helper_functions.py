"""
Test suite for helper_functions.py
"""

import unittest
import os
import tempfile
import gzip
from PhiSpyModules.helper_functions import is_gzip_file


__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


class HelperFunctionsTest(unittest.TestCase):
    
    def setUp(self):
        """Create temporary files for testing"""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create a regular text file
        self.text_file = os.path.join(self.temp_dir, 'test.txt')
        with open(self.text_file, 'w') as f:
            f.write('This is a test file\n')
        
        # Create a gzipped file
        self.gzip_file = os.path.join(self.temp_dir, 'test.gz')
        with gzip.open(self.gzip_file, 'wt') as f:
            f.write('This is a gzipped test file\n')
    
    def tearDown(self):
        """Clean up temporary files"""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_is_gzip_file_with_gzip(self):
        """Test that gzipped files are correctly identified"""
        self.assertTrue(is_gzip_file(self.gzip_file))
    
    def test_is_gzip_file_with_text(self):
        """Test that regular text files are not identified as gzipped"""
        self.assertFalse(is_gzip_file(self.text_file))


if __name__ == '__main__':
    unittest.main()

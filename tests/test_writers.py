"""
Test suite for writers.py
"""

import unittest
import tempfile
import os
import shutil
from argparse import Namespace
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


class WritersTest(unittest.TestCase):
    
    def setUp(self):
        """Create temporary directory and mock data for testing"""
        self.temp_dir = tempfile.mkdtemp()
        
        # Create a simple mock SeqRecord for testing
        self.record = SeqRecord(
            Seq("ATGCATGCATGC"),
            id="test_contig",
            name="test",
            description="test sequence"
        )
        self.record.features.append(
            SeqFeature(FeatureLocation(0, 12, strand=1), type="CDS")
        )
        
        # Create mock prophage data
        self.mock_pp = {
            1: {
                'contig': 'test_contig',
                'start': 100,
                'stop': 500,
                'att': [90, 110, 490, 510]
            },
            2: {
                'contig': 'test_contig',
                'start': 1000,
                'stop': 1500
            }
        }
        
        # Create namespace object with required attributes
        self.args = Namespace(
            output_dir=self.temp_dir,
            file_prefix='test_',
            quiet=True,
            pp=self.mock_pp,
            record=[self.record]
        )
    
    def tearDown(self):
        """Clean up temporary files"""
        shutil.rmtree(self.temp_dir)
    
    def test_prophage_measurements_to_tbl(self):
        """Test converting prophage measurements to table format"""
        from PhiSpyModules.writers import prophage_measurements_to_tbl
        
        # The function expects output in a specific format
        # Just test that it runs without errors
        input_file = os.path.join(self.temp_dir, 'input.txt')
        output_file = os.path.join(self.temp_dir, 'output.tbl')
        
        # Create mock input file with complete columns (10 columns required)
        with open(input_file, 'w') as f:
            f.write("PP\tcontig\tstart\tstop\tatt_l_start\tatt_l_stop\tatt_r_start\tatt_r_stop\tlen\tnum_genes\n")
            f.write("1\tcontig1\t100\t500\t90\t110\t490\t510\t400\t5\n")
            f.write("2\tcontig1\t1000\t1500\t990\t1010\t1490\t1510\t500\t6\n")
        
        # Test the function
        prophage_measurements_to_tbl(input_file, output_file)
        
        # Verify output file exists
        self.assertTrue(os.path.exists(output_file))
        
        # Verify output file is not empty
        with open(output_file, 'r') as f:
            content = f.read()
            self.assertTrue(len(content) > 0)


if __name__ == '__main__':
    unittest.main()

"""
Test suite for seqio_filter.py

This is (currently) not complete, and is a work in progres.
"""

import unittest
import gzip

from Bio import SeqIO, SeqFeature
from Bio.SeqFeature import FeatureLocation
from PhiSpyModules import SeqioFilter


__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'





class SeqIO_FilterTest(unittest.TestCase):
    def test_merge_split(self):
        """
        Test the merge and split function

        Paracoccus_yeei_TT13.gb has the following compound locations:
        PYTT13_06780: join{[1366920:1367182](-), [1365992:1366921](-)}
                      [1365992:1367182](-)

        PYTT13_11395: join{[2283890:2284152](+), [2284151:2285080](+)}
                      [2283890:2285080](+)

        PYTT13_11465: join{[2301567:2301817](+), [2301816:2302655](+)}
                      [2301567:2302655](+)

        PYTT13_12460: join{[2495319:2495581](+), [2495580:2496509](+)}
                      [2495319:2496509](+)

        PYTT13_16505: join{[3331106:3331356](+), [3331355:3332194](+)}
                      [3331106:3332194](+)

        :return:
        """
        correct_locations = {
            'PYTT13_06780': (1365992, 1367182),
            'PYTT13_11395': (2283890, 2285080),
            'PYTT13_11465': (2301567, 2302655),
            'PYTT13_12460': (2495319, 2496509),
            'PYTT13_16505': (3331106, 3332194)
        }

        testgbk = "test_genbank_files/Paracoccus_yeei_TT13.gb.gz"
        handle = gzip.open(testgbk, 'rt')
        record = SeqioFilter(SeqIO.parse(handle, "genbank"))
        handle.close()
        for s in record:
            for f in s.get_features("CDS"):
                if 'locus_tag' in f.qualifiers and f.qualifiers['locus_tag'][0] in correct_locations:
                    lt = f.qualifiers['locus_tag'][0]
                    self.assertIsInstance(f.location,FeatureLocation)
                    self.assertEqual(f.location.start, correct_locations[lt][0])
                    self.assertEqual(f.location.end, correct_locations[lt][1])


if __name__ == '__main__':
    unittest.main()
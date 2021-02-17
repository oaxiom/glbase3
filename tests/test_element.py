"""
expression.py tester code.

Part of glbase

Tests for expression class

"""

import unittest

# get glbase
import sys, os, math
sys.path.append(os.path.realpath("../../"))

import glbase3 as gl

class Test_Motif_Element(unittest.TestCase):
    def setUp(self):
        self.sox2 = gl.motif("Sox2", "ywttswn")

    def test_motif_fasta_scan(self):
        fasta = [
            'ACcactcacccattgtaaAAAcCCCAaaaa', 
            'ACACCCATTGTc',
            'cagtgtCCtggcC',
            'CAGTCtCaTTgtC',
            'NNNNNNNNNNNNN',
            'GCGCGCGCGCGCGC']
        fasta = [{'seq': s} for s in fasta]
        fgl = gl.genelist()
        fgl.load_list(fasta)
        
        r = self.sox2.scan_sequences(fgl)
        #print(r.all())
        #r.saveTSV('')
        f = r.get(self.sox2.name, 'Found') # Found seqs only;
        #print(f)
        self.assertEqual(len(r), 6)
        self.assertEqual(len(f), 3)
        #print(f)
        self.assertDictEqual(f[0], {'seq': 'ACcactcacccattgtaaAAAcCCCAaaaa', 'Sox2': 'Found', 'Sox2_seqs': 'cattgta'})
        

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Motif_Element)
    unittest.TextTestRunner(verbosity=2).run(suite)

"""

Tester Suite:

**Purpose**

This one checks

"""

import unittest

# get glbase
import sys, os
sys.path.append(os.path.realpath("../../"))

import glbase

glbase.config.SILENT = True
glbase.config.set_log_level(None)

class Test_Genome(unittest.TestCase):
    def setUp(self):
        self.gsql = glbase.genome_sql(new=True, filename='/tmp/test_genome_sql.sql') # This is platform specific and breaks on Windows
        self.gsql.add_feature(glbase.location(chr='chr1', left=110, right=120), 
            glbase.location(chr='chr1', left=110, right=120), 
            10, [1,2,3,4], [5,6,7,8], 
            'Nanog', '+')

    def test_genome_sql(self):   
        # Test collision capabilities
        self.assertEqual(len(self.gsql.getFeatures(loc='chr1:100-130')), 1)
        self.assertEqual(len(self.gsql.getFeatures(loc='chr1:100-110')), 1)   # Should still return the entry
        self.assertEqual(len(self.gsql.getFeatures(loc='chr1:200-330')), 0)
        
        # Test key return formatting:
        self.assertEqual(self.gsql.getFeatures(loc='chr1:100-130')[0]['loc']['chr'], '1')
        self.assertTrue(isinstance(self.gsql.getFeatures(loc='chr1:100-130')[0]['exonStarts'], list))

    def test_get_sequences(self):
        genome_mm10 = glbase.genome()
        genome_mm10.bindSequence(os.path.join(os.path.expanduser("~"), "mm10/seq"))
        seq = genome_mm10.getSequence("chr1:10000000-10000200")
        self.assertEqual(seq, 'TTTTCAATGCAGGAAATGCAATTGTTCTGTAGGTACAAGTGGGTCAGATTTGTGGTGTAATTCAGGTTAGTGACTTGACTAATGCGATTATCATATAAATATAAAACACTCAGGTTTCTGCAAAGAGAGAGGTCATCCTGAAAAGTAAACAAAACAGGCCCTATTTAATTACCTCACAAGCTTACAAGTTGGATTTTAAGA')
        seq = genome_mm10.getSequence("chrX:10000000-10000200")
        self.assertEqual(seq, 'AGTAATGGTAGTCATATGGTCCTTTGACACTTCAAGATTTATATTTAATTGGAAAAGAGAAAGCCATAGAAAGAATAATGGTAAGCCCTATTTATAGGAATAAAGTTGGTAGAAATACCAAGTCCAAAAATCCTTTGAAACTGAAAATCTGTAGCTGTGTTGGGTTTTGTTTTTTCCATGAAATACACCACACTGATCTGA')       

    def test_save_fasta(self):
        genome_mm10 = glbase.genome()
        genome_mm10.bindSequence(os.path.join(os.path.expanduser("~"), "mm10/seq"))
        newl = [{"name": "A", "loc": glbase.location(loc="chr1:10000000-10000200")},
            {"name": "X", "loc": glbase.location(loc="chrX:10000000-10000200")}]
        newgl = glbase.genelist()
        newgl.load_list(newl)
        fasta = genome_mm10.getSequences(newgl)
        fasta.saveFASTA(filename="/tmp/test_fasta.fa", name=["loc", "name"])
        
        oh = open("/tmp/test_fasta.fa")
        self.assertEqual(oh.readline().strip(), '>chr1:10000000-10000200_A')
        self.assertEqual(oh.readline().strip(), 'TTTTCAATGCAGGAAATGCAATTGTTCTGTAGGTACAAGTGGGTCAGATTTGTGGTGTAATTCAGGTTAGTGACTTGACTAATGCGATTATCATATAAATATAAAACACTCAGGTTTCTGCAAAGAGAGAGGTCATCCTGAAAAGTAAACAAAACAGGCCCTATTTAATTACCTCACAAGCTTACAAGTTGGATTTTAAGA')
        self.assertEqual(oh.readline().strip(), '>chrX:10000000-10000200_X')
        self.assertEqual(oh.readline().strip(), 'AGTAATGGTAGTCATATGGTCCTTTGACACTTCAAGATTTATATTTAATTGGAAAAGAGAAAGCCATAGAAAGAATAATGGTAAGCCCTATTTATAGGAATAAAGTTGGTAGAAATACCAAGTCCAAAAATCCTTTGAAACTGAAAATCTGTAGCTGTGTTGGGTTTTGTTTTTTCCATGAAATACACCACACTGATCTGA')
        oh.close()
        
        fasta.saveFASTA(filename="/tmp/test_fasta.fa")       
        oh = open("/tmp/test_fasta.fa")
        self.assertEqual(oh.readline().strip(), '>chr1:10000000-10000200')
        self.assertEqual(oh.readline().strip(), 'TTTTCAATGCAGGAAATGCAATTGTTCTGTAGGTACAAGTGGGTCAGATTTGTGGTGTAATTCAGGTTAGTGACTTGACTAATGCGATTATCATATAAATATAAAACACTCAGGTTTCTGCAAAGAGAGAGGTCATCCTGAAAAGTAAACAAAACAGGCCCTATTTAATTACCTCACAAGCTTACAAGTTGGATTTTAAGA')
        self.assertEqual(oh.readline().strip(), '>chrX:10000000-10000200')
        self.assertEqual(oh.readline().strip(), 'AGTAATGGTAGTCATATGGTCCTTTGACACTTCAAGATTTATATTTAATTGGAAAAGAGAAAGCCATAGAAAGAATAATGGTAAGCCCTATTTATAGGAATAAAGTTGGTAGAAATACCAAGTCCAAAAATCCTTTGAAACTGAAAATCTGTAGCTGTGTTGGGTTTTGTTTTTTCCATGAAATACACCACACTGATCTGA')
        oh.close()

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Genome)
    unittest.TextTestRunner(verbosity=2).run(suite)

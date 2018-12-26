"""

Tester Suite:

**Purpose**

This one checks

"""

import unittest

# get glbase
import sys, os
sys.path.append(os.path.realpath("../../"))

import glbase3

glbase3.config.SILENT = True
glbase3.config.set_log_level(None)

class Test_Genome(unittest.TestCase):
    def setUp(self):
        self.gsql = glbase3.genome_sql(new=True, filename='/tmp/test_genome_sql.sql') # This is platform specific and breaks on Windows
        self.gsql.add_feature(glbase3.location(chr='chr1', left=110, right=120), 
            glbase3.location(chr='chr1', left=110, right=120), 
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
        genome_mm10 = glbase3.genome()
        genome_mm10.bindSequence("test_data/seq")
        seq = genome_mm10.getSequence("chr1:100-150")
        self.assertEqual(seq, 'ATCAGACAGGTAGATCATCTCGCTCCGAGCTTGCCACCAGCAAACCATTGC')
        seq = genome_mm10.getSequence("chrA:100-150")
        self.assertEqual(seq, 'GTAAAAACCCGATGGAATACTCATCCAGTAAGTCCGAACCACTTCAACATC')       

    def test_save_fasta(self):
        genome_mm10 = glbase3.genome()
        genome_mm10.bindSequence("test_data/seq")
        newl = [{"name": "A", "loc": glbase3.location(loc="chr1:100-150")},
            {"name": "X", "loc": glbase3.location(loc="chrA:100-150")}]
        newgl = glbase3.genelist()
        newgl.load_list(newl)
        fasta = genome_mm10.getSequences(newgl)
        fasta.saveFASTA(filename="/tmp/test_fasta.fa", name=["loc", "name"])
        
        oh = open("/tmp/test_fasta.fa")
        self.assertEqual(oh.readline().strip(), '>chr1:100-150_A')
        self.assertEqual(oh.readline().strip(), 'ATCAGACAGGTAGATCATCTCGCTCCGAGCTTGCCACCAGCAAACCATTGC')
        self.assertEqual(oh.readline().strip(), '>chrA:100-150_X')
        self.assertEqual(oh.readline().strip(), 'GTAAAAACCCGATGGAATACTCATCCAGTAAGTCCGAACCACTTCAACATC')
        oh.close()
        
        fasta.saveFASTA(filename="/tmp/test_fasta.fa")       
        oh = open("/tmp/test_fasta.fa")
        self.assertEqual(oh.readline().strip(), '>chr1:100-150')
        self.assertEqual(oh.readline().strip(), 'ATCAGACAGGTAGATCATCTCGCTCCGAGCTTGCCACCAGCAAACCATTGC')
        self.assertEqual(oh.readline().strip(), '>chrA:100-150')
        self.assertEqual(oh.readline().strip(), 'GTAAAAACCCGATGGAATACTCATCCAGTAAGTCCGAACCACTTCAACATC')
        oh.close()

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Genome)
    unittest.TextTestRunner(verbosity=2).run(suite)

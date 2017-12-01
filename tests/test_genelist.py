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

class Test_GeneList(unittest.TestCase):
    def setUp(self):
        self.a = glbase3.genelist(filename="test_data/array_data.csv", format=glbase3.format.sniffer)

    def test_get_by_slice(self):
        self.assertEqual(len(self.a[0:2]), 2)
        # notice that this arrangement is not a microarray and Mash and GFP keys are separate.
        # single slices now return a plain dictionary
        
        self.assertDictEqual(self.a[-1], {'name': 'Pdia4', 'GFP': 1.18, 'Mash': 0.6, 'array_systematic_name': 'scl29051.11.1_27-S', 'refseq': 'NM_009787', 'entrez': 12304})
        self.assertDictEqual(self.a[2], {'name': 'Srpr', 'GFP': 1, 'Mash': 0.77, 'array_systematic_name': 'scl0067398.1_126-S', 'refseq': 'NM_026130', 'entrez': 67398})   
        
    def test_get_by_name(self):
        self.assertEqual(len(self.a["name"]), len(self.a))
        self.assertEqual(self.a["name"][0], "Lypla1")
        self.assertEqual(self.a["name"][1], "Smc1a")

    def test_edit_item_persistence(self):
        get = self.a[2] # this returns a vanilla dictionary
        get["name"] = "Bleh"
        self.assertNotEqual(self.a[2]["name"], "Srpr")
        self.assertEqual(get["name"], "Bleh")
    
    def test_rename_key(self):
        new = self.a.renameKey("name", "meh")
                
        self.assertTrue("meh" in new.linearData[0])
        self.assertFalse("name" in new.linearData[0])
        
        self.assertTrue("name" in self.a.linearData[0]) # test original list intact  
        self.assertFalse("meh" in self.a.linearData[0]) 

    def test_rename_key_keep_old_key(self):
        new = self.a.renameKey("name", "meh", keep_old_key=True)
                
        self.assertTrue("meh" in new.linearData[0])
        self.assertTrue("name" in new.linearData[0])
        
        self.assertTrue("name" in self.a.linearData[0]) # test original list intact  
        self.assertFalse("meh" in self.a.linearData[0]) 
        
    def test_addfakekey(self):       
        b = self.a.addFakeKey("meh", "++")
        self.assertTrue("meh" in b.linearData[0])
        self.assertEqual(b.linearData[1]["meh"], "++")
        
    def test_splitbykey(self):
        newl = self.a.getRowsByKey(key="name", values="Ptp")
        self.assertTrue(len(newl) == 3) # A little basic?

    def test_renameKey(self):
        newl = self.a.renameKey("name", "other-name")
        self.assertTrue("name" in self.a[0])
        self.assertTrue("other-name" not in self.a[0])
        self.assertTrue("name" not in newl[0])
        self.assertTrue("other-name" in newl[0])

    def test_load_gzips(self):
        self.a = glbase3.genelist(filename="test_data/array_data.csv.gz", format=glbase3.format.sniffer, gzip=True)
        self.assertEqual(str(self.a[-1]), "{'name': 'Pdia4', 'GFP': 1.18, 'Mash': 0.6, 'array_systematic_name': 'scl29051.11.1_27-S', 'refseq': 'NM_009787', 'entrez': 12304}")
        self.assertEqual(str(self.a[2]),  "{'name': 'Srpr', 'GFP': 1, 'Mash': 0.77, 'array_systematic_name': 'scl0067398.1_126-S', 'refseq': 'NM_026130', 'entrez': 67398}")

    def test_load_FASTA_gzips(self):
        self.a = glbase3.genelist(filename="test_data/Fasta_file.fa.gz", format=glbase3.format.fasta, gzip=True)
        self.assertEqual(self.a[0]['seq'], 'AAATctggatacagtggcctttatttctagttccagtgactgggagactgaaacaagagagtcacttgagtacaggagtgcaaggctagcttgagcaatatagtaagactatctcaaaaTGTGAATTtagatcaacagaattgacatcaagaaaaatactgatatcactcaaagcaatctacagattcaacacaatctccatcaacatgacaatgacttccatcaGCATGACAATGACTCCATCAACATGCCAATGGGCCCCATCAACATAACAATGACCCCTATCATCATGACAATGATCCCCATCAACATGACAATGACCTCCATCAACATGACAATTACTCCTGTCAACATGCCAATtgttggggttcagaagtcaccctgcaaaccacaagaacact')

    def test_removeDuplicatesByLoc(self):
        a = [{'loc': glbase3.location(chr=1, left=100, right=200)},
            {'loc':  glbase3.location(chr=1, left=100, right=200)},
            {'loc':  glbase3.location(chr=1, left=100, right=200)},
            {'loc':  glbase3.location(chr=1, left=100, right=200)},
            {'loc':  glbase3.location(chr=1, left=100, right=200)},
            {'loc':  glbase3.location(chr=1, left=100, right=200)},
            {'loc':  glbase3.location(chr=1, left=100, right=200)},
            {'loc':  glbase3.location(chr=1, left=100, right=200)},
            {'loc':  glbase3.location(chr=1, left=130, right=230)},
            {'loc':  glbase3.location(chr=1, left=130, right=230)},
            {'loc':  glbase3.location(chr=1, left=9800, right=9990)}, # across bucket
            {'loc':  glbase3.location(chr=1, left=10001, right=10200)},]
        gl = glbase3.genelist()
        gl.load_list(a)
        dups = gl.removeDuplicatesByLoc('loc', 10)
        self.assertEqual(len(dups), 4)
        dups = gl.removeDuplicatesByLoc('loc', 200)
        self.assertEqual(len(dups), 2)
        
        #config.bucket_size = 10000 

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_GeneList)
    unittest.TextTestRunner(verbosity=2).run(suite)

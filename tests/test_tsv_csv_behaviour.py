"""

Tester Suite:

**Purpose**

Check that the force_tsv args and format specifier is working correctly

"""

import unittest

# get glbase
import sys, os
sys.path.append(os.path.realpath("../../"))

import glbase as gl
from glbase.utils import qcollide

gl.config.SILENT = True
gl.config.set_log_level(None)

class Test_TSVCSV(unittest.TestCase):  
    
    def test_force_tsvarg(self):
        form = dict(tss_loc=1, skiplines=0) # This loads tss_loc as strings
        form_delayed = dict(tss_loc=1, skiplines=0) # delayedlists must have skiplines
        a = gl.genelist(filename="../example/shared_raw_data/mm9_refGene.tsv", force_tsv=True, format=form)
        c = gl.delayedlist(filename="../example/shared_raw_data/mm9_refGene.tsv", format=form, force_tsv=True)
        d = gl.genome(filename="../example/shared_raw_data/mm9_refGene.tsv", format=form_delayed, force_tsv=True)
        e = gl.expression(filename="../example/shared_raw_data/mm9_refGene.tsv", format=form, force_tsv=True, expn="column[5:]") # fake array data # must go last as it modifies format

        # Make sure glbase is not just bodging it all in in one key:
        self.assertEqual("chr1:134212701-134212701", a[0]["tss_loc"])
        self.assertEqual("chr1:134212701-134212701", c[0]["tss_loc"])
        self.assertEqual("chr1:134212701-134212701", d[0]["tss_loc"]) # dls should work as __getitem__() will return the zeroth entry.
        self.assertEqual("chr1:134212701-134212701", e[0]["tss_loc"])
    
    def test_force_tsv_format(self):
        form = dict(tss_loc=1, force_tsv=True, chr=1)
        form_delayed = dict(tss_loc=1, force_tsv=True, skiplines=0) # delayedlists must have skiplines
        a = gl.genelist(filename="../example/shared_raw_data/mm9_refGene.tsv", format=form)
        c = gl.delayedlist(filename="../example/shared_raw_data/mm9_refGene.tsv", format=form_delayed)
        d = gl.genome(filename="../example/shared_raw_data/mm9_refGene.tsv", format=form)
        e = gl.expression(filename="../example/shared_raw_data/mm9_refGene.tsv", format=form, expn="column[5:]") # must go last as it modifies format
        
        # Make sure glbase is not just bodging it all in in one key:
        self.assertEqual("chr1:134212701-134212701", a[0]["tss_loc"])
        self.assertEqual("chr1:134212701-134212701", c[0]["tss_loc"])
        self.assertEqual("chr1:134212701-134212701", d[0]["tss_loc"])
        self.assertEqual("chr1:134212701-134212701", e[0]["tss_loc"])
    
    def test_tsv_sniffer_force(self):
        # These are all tsv files
        # sniffer correctly loads locations.
        a = gl.genelist(filename="../example/shared_raw_data/mm9_refGene.tsv", force_tsv=True, format=gl.format.sniffer)
        d = gl.genome(filename="../example/shared_raw_data/mm9_refGene.tsv", force_tsv=True, format=gl.format.sniffer)
        # Microarrays and delayedlists can't be sniffed
        
        # Make sure glbase is not just bodging it all in in one key:
        self.assertEqual("chr1:134212701-134212701", a[0]["tss_loc"])
        self.assertEqual("chr1:134212701-134212701", d[0]["tss_loc"])
    
    """
    def test_tsv_sniffer_auto_detect(self):
        # These are all tsv files
        a = gl.genelist(filename="../example/raw_data/mm9_refGene.tsv")
        b = gl.peaklist(filename="../example/raw_data/mm9_refGene.tsv")
        d = gl.genome(filename="../example/raw_data/mm9_refGene.tsv")
        # Microarrays and delayedlists can't be sniffed
        
        # Make sure glbase is not just bodging it all in in one key:
        self.assertEqual("chr1:134212701-134230065", a[0]["tss_loc"])
        self.assertEqual("chr1:134212701-134230065", b[0]["tss_loc"])
        self.assertEqual("chr1:134212701-134230065", d[0]["tss_loc"])
    """
        
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_TSVCSV)
    unittest.TextTestRunner(verbosity=2).run(suite)

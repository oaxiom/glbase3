"""
track.py tester code.

Part of glbase

Tests that track is performing accurately.

TODO:
-----
. some of the tests are not very extensive yet.

"""

import unittest

# get glbase
import sys, os
import glbase3
import numpy

class Test_Flat_Function(unittest.TestCase):
    def setUp(self):
        self.t = glbase3.flat_track(filename="/tmp/test.flat", new=True, name="Test Track", bin_format="f")
        for i in range(0, 100):
            self.t.add_score(chromosome=1, left=i, right=i+1, score=i)
        for i in range(100,110):
            self.t.add_score(chromosome=2, left=i, right=i+1, score=0)
        self.t.add_score(chromosome=2, left=111, right=112, score=2)
        self.t.add_score(chromosome=2, left=99, right=100, score=2)
        self.t.finalise()

    def test_get(self):
        a = self.t.get(glbase3.location(loc="chr1:10-20"))
        self.assertListEqual(list(a), [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0])

    def test_reload(self):
        t = glbase3.flat_track(filename="/tmp/test.flat", bin_format="f")
        a = self.t.get(glbase3.location(loc="chr1:0-100"))
        expected_result = [i for i in range(0, 100)] +[0]
        self.assertTrue(False not in [int(x) == int(y) for x, y in zip(a, expected_result)]) # all seqs.
    
    def test_meta_data(self):
        t = glbase3.flat_track(filename="/tmp/test.flat", new=False, name="Test Track", bin_format="f")
        self.assertEqual(t["name"], "Test Track")
        self.assertEqual(t["glbase_version"], glbase3.config.version)
        self.assertEqual(t["bin_format"], "f")

    def test_name_override(self):
        t =  glbase3.flat_track(filename="/tmp/test.flat", new=False, name="Changed Name", bin_format="f")
        self.assertEqual(t["name"], "Changed Name")

    def test_pileup(self):
        t = glbase3.flat_track(filename="/tmp/test.flat", bin_format="f")
        
        g = glbase3.genelist(filename="test_data/track_test.bed", format=glbase3.format.bed).pointify().expand('loc', 15)
        L, _ = t.pileup(genelists=g, filename="test_output.png", bandwidth=15, respect_strand=False)
        expected_result = numpy.array(list(range(0, 30)))
        self.assertListEqual(list(L[0]), list(expected_result))

    def test_pileup_respect_strand(self):
        t = glbase3.flat_track(filename="/tmp/test.flat", bin_format="f")
        
        g = glbase3.genelist(filename="test_data/track_test.bed", format=glbase3.format.bed)
        L, meh = t.pileup(genelists=g, filename="test_output.png", bandwidth=15, respect_strand=True)
        
        expected_result = numpy.zeros(31)
        expected_result += 14.5
        # IF respect strand is true then you will get 4 arrays --> --> and <-- <-- the average of the 4
        # will be 14.5 for all points. 
        # Note this is identical to the above test_draw_pielup()
        # except respect_strand=False
        print(L)
        print(expected_result)
        self.assertTrue(False not in [x == y for x, y in zip(L, expected_result)])

    def test_mask(self):
        t = glbase3.flat_track(filename="/tmp/test.flat", bin_format="f")
        a = t.get(glbase3.location(loc="chr2:99-111"))
        unmasked = [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,2]
        self.assertTrue(False not in [int(x) == int(y) for x, y in zip(a, unmasked)]) # all seqs.
        a = t.get(glbase3.location(loc="chr2:99-111"), mask_zero=True)
        expected = "[2.0 -- -- -- -- -- -- -- -- -- -- --]" # not sure how to test this apart from a string.
        self.assertEqual(str(a), expected)
        
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Flat_Function)
    unittest.TextTestRunner(verbosity=2).run(suite)

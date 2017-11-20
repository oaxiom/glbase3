"""
track.py tester code.

Part of glbase

Tests that track is performing accurately.

TODO:
-----
. some of the tests are not very extensive yet.
. There is a bug, whereby the order of the retrieved reads is random. Making the test erroneously fail.
. The /tmp/ stuff is platform specific

"""

import unittest

# get glbase
import sys, os
#sys.path.append(os.path.realpath("../"))

import glbase3 as gl
import numpy

gl.config.SILENT = True
gl.config.set_log_level(None)

# The format of the returned reads can sometimes change...
# track:
li = 0
ri = 1
st = 2
# track3:
#li = 2
#ri = 3
#st = 4

def generate_track(filename, norm_factor):
    t = gl.track(filename=filename, new=True, name="Test Track", norm_factor=norm_factor)
    t.add_location(gl.location(loc="chr1:10-20"))
    t.add_location(gl.location(loc="chr1:10-20")) # test duplicates
    t.add_location(gl.location(loc="chr1:21-22")) # test 1's
    t.add_location(gl.location(loc="chr1:23-23")) # test single outside
    t.add_location(gl.location(loc="chr1:1-100")) # test massive span
    t.add_location(gl.location(loc="chr1:9-19")) # inside test
    t.add_location(gl.location(loc="chr1:15-15")) # 1bp inside
    t.add_location(gl.location(loc="chr1:19-25")) # over right border
    t.add_location(gl.location(loc="chr1:5-11")) # over left border
    t.add_location(gl.location(loc="chrX:1-1000")) # letter chromsome
    t.add_location(gl.location(loc="chr2:2-2000")) # other numeric chr
    t.add_location(gl.location(loc="chr3:1-2"), strand="-") # test strand
    t.add_location(gl.location(loc="chr10:1-2"), strand="-") # test strand
    t.add_location(gl.location(loc="chr100:1-2"), strand="-") # test strand
    t.finalise()
    return(t)

def generate_track_frags_only(filename, norm_factor):
    t = gl.track(filename=filename, new=True, name="Test Track", norm_factor=norm_factor)
    t.add_location(gl.location(loc="chr1:10-20"))
    t.add_location(gl.location(loc="chr1:10-21")) # test duplicates
    t.add_location(gl.location(loc="chr1:10-22")) # test 1's
    t.add_location(gl.location(loc="chr1:9-23")) # test 1's
    t.finalise()
    return(t)

class Test_Track_Function(unittest.TestCase):
    def setUp(self):
        self.t = generate_track("/tmp/test.trk", 1.0)
        self.frags = generate_track_frags_only("/tmp/test_frags.trk", 1.0)

    def test_get_array(self):
        a = self.t.get(gl.location(loc="chr1:10-20"))
        self.assertEqual(str(a), "[ 5.  5.  4.  4.  4.  5.  4.  4.  4.  5.  4.]")

    def test_get_reads(self):
        # Test a range of reads and edges.
        expected_reads = [(10, 20, '+'), (10, 20, '+'), (21, 22, '+'), (23, 23, '+'), (1, 100, '+'), (9, 19, '+'), 
            (15, 15, '+'), (19, 25, '+'), (5, 11, '+')]
        reads = self.t.get_reads(gl.location(loc="chr1:1-100"))
        res = [(r[li], r[ri], r[st]) in expected_reads for r in reads]
        self.assertTrue(False not in res)
        self.assertTrue(len(expected_reads) == len(reads))
        
        # Test an inner get for a long read
        expected_reads = [(1, 1000, '+')]
        reads = self.t.get_reads(gl.location(loc="chrX:10-20"))
        res = [(r[li], r[ri], r[st]) in expected_reads for r in reads] # Can't test chrom, track1 does not return chrom
        self.assertTrue(False not in res)
        self.assertTrue(len(expected_reads) == len(reads))
        
        expected_reads = [(1, 2, '-')]
        reads = self.t.get_reads(gl.location(loc="chr3:1-5"))
        res = [(r[li], r[ri], r[st]) in expected_reads for r in reads] # Can't test chrom, track1 does not return chrom
        self.assertTrue(False not in res)
        self.assertTrue(len(expected_reads) == len(reads))

        reads = self.t.get_reads(gl.location(loc="chr3:50-100")) # empty region test
        self.assertTrue(len(reads) == 0)
        
        expected_reads = [(10, 20, '+'), (10, 20, '+'), (1, 100, '+'), (9, 19, '+'), 
            (5, 11, '+'), (19,25, '+'), (15, 15, '+')]
        reads = self.t.get_reads(gl.location(loc="chr1:10-20")) #  overspanning read test
        res = [(r[li], r[ri], r[st]) in expected_reads for r in reads] # Can't test chrom, track1 does not return chrom
        self.assertTrue(False not in res)
        self.assertTrue(len(expected_reads) == len(reads))
        
        expected_reads = [(15, 15, '+'), (10, 20, '+'), (10, 20, '+'), (1, 100, '+'), (9, 19, '+')]
        reads = self.t.get_reads(gl.location(loc="chr1:15-15")) #  single point location read test
        res = [(r[li], r[ri], r[st]) in expected_reads for r in reads] # Can't test chrom, track1 does not return chrom
        self.assertTrue(False not in res)
        self.assertTrue(len(expected_reads) == len(reads))

    def test_read_extend(self):
        a = self.t.get(gl.location(loc="chr1:10-20"), read_extend=1)
        self.assertEqual(str(a), "[ 5.  5.  5.  4.  4.  5.  5.  4.  4.  5.  5.]") # These are correct. Always returns floats now
        a = self.t.get(gl.location(loc="chr1:10-20"), read_extend=2)
        self.assertEqual(str(a), "[ 5.  5.  5.  5.  4.  5.  5.  5.  4.  5.  5.]") # .

    def test_read_extend_frags_only(self):
        a = self.frags.get(gl.location(loc="chr1:5-25"), read_extend=1)
        self.assertTrue(False not in [i == e for i, e in zip(a, [0,0,0,0,1,4,4,4,4,4,4,4,4,4,4,4,4,3,2,1,0])])

    def test_reload(self):
        t = gl.track(filename="/tmp/test.trk") # uses self.t
        # list of left and rights to compare
        expected_reads = [(10, 20, '+'), (10, 20, '+'), (21, 22, '+'), (23, 23, '+'), (1, 100, '+'), (9, 19, '+'), 
            (15, 15, '+'), (19, 25, '+'), (5, 11, '+')]
        reads = self.t.get_reads(gl.location(loc="chr1:1-100"))
        res = [(r[li], r[ri], r[st]) in expected_reads for r in reads]
        self.assertTrue(False not in res)
        self.assertTrue(len(expected_reads) == len(reads))
        
    def test_meta_data(self):
        t = gl.track(filename="/tmp/test.trk")
        self.assertEqual(t["name"], "Test Track")
        self.assertEqual(t["glbase_version"], gl.config.version)
        #self.assertEqual(t["bin_format"], None)

    def test_name_override(self):
        t = gl.track(filename="/tmp/test.trk", name="Changed Name")
        self.assertEqual(t["name"], "Changed Name")

    def test_pileup(self):
        t = gl.track(filename="/tmp/test_pileup.trk", new=True, name="Test Track")
        for i in [10, 10, 10, 10, 10, 10]:
            t.add_location(gl.location(chr="chr1", left=i, right=i+5))
        t.finalise()
        
        g = gl.genelist(filename="track_test.bed", format=gl.format.bed)
        L = t.pileup(genelist=g, filename="test_images/test_output.png", heatmap_filename="test_images/test_heatmap.png", window_size=15, 
            bin_size=1, respect_strand=True, normalise=False, read_extend=1, raw_tag_filename="test_images/test_tags.tsv")
        
        expected_result = numpy.array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  3.,  3.,  3.,  6.,  6.,  6.,  6.,  3., 3.,  3.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,])
        # units are now reads per item in genelist
        
        #print L['pileup']
        #print expected_result
        
        self.assertTrue(False not in [x == y for x, y in zip(L["pileup"], expected_result)])

    def test_saveBedGraph(self):
        self.t.saveBedGraph('/tmp/bedgraph.bg')
        oh = open('/tmp/bedgraph.bg', 'rU')
        l1 = oh.readline().strip().split('\t')
        oh.close()
        self.assertEqual(l1, ['chr1', '0', '100', '2.68'])
        
    '''def test_seqtotrk(self):
        gl.seqToTrk("track_test2.bed", "test_trk2.trk", name="Test Track", format=gl.format.bed)
        t = gl.track(filename="test_trk2.trk")
        b = t.get("chr16:23984998-23985025")
        
        expected_result = numpy.array([0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0])
        self.assertTrue(False not in [x == y for x, y in zip(b, expected_result)])
    '''
    def test_num_reads(self):
        num_reads = self.t.get_total_num_reads()
        self.assertEqual(num_reads, 14)

    def test_norm_factor(self):
        n1 = generate_track(filename="/tmp/test_1.0.trk", norm_factor=1.0) # Test silent argument
        actual = self.t.get(gl.location(loc="chr1:10-20"), read_extend=1)
        observed = n1.get(gl.location(loc="chr1:10-20"), read_extend=1)
        self.assertTrue(False not in [x == y for x, y in zip(actual, observed)])
        
        n2 = generate_track(filename="/tmp/test_2.0.trk", norm_factor=2.0)
        observed = n2.get(gl.location(loc="chr1:10-20"), read_extend=1)
        self.assertTrue(False not in [x == y for x, y in zip(actual*2.0, observed)])
        
        n3 = generate_track(filename="/tmp/test_0.5.trk", norm_factor=0.5)
        observed = n3.get(gl.location(loc="chr1:10-20"), read_extend=1)
        self.assertTrue(False not in [x == y for x, y in zip(actual*0.5, observed)])       

    def test_catch_incomplete_trk_assertion(self):
        t = gl.track(filename="/tmp/test_incomplete.trk", name="test", new=True)
        self.assertRaises(AssertionError, gl.track, "/tmp/test_incomplete.trk") # should throw an assertion error only
        
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Track_Function)
    unittest.TextTestRunner(verbosity=2).run(suite)

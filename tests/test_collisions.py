"""

Tester Suite:

**Purpose**

This one checks collisions and overlaps using the three toy testN.csvs.

"""

import unittest

# get glbase
import sys, os
sys.path.append(os.path.realpath("../../"))

import glbase3
from glbase3.utils import qcollide

#glbase3.config.SILENT = True
#glbase3.config.set_log_level(None)

class Test_Collisions_Overlaps(unittest.TestCase):
    def setUp(self):
        self.a = glbase3.genelist(filename="test_data/testA.csv", format=glbase3.format.sniffer)
        self.b = glbase3.genelist(filename="test_data/testB.csv", format=glbase3.format.sniffer)
        self.c = glbase3.genelist(filename="test_data/testC.csv", format=glbase3.format.sniffer)
        self.d = glbase3.genelist(filename="test_data/ccat_list.region", format=glbase3.format.ccat_output)
        self.e = glbase3.genelist(filename="test_data/macs_list.xls", format=glbase3.format.macs_output)       
        
        fake1 = [{"name": "gene1"}, 
            {"name": "gene2"}, 
            {"name": "gene3"}, 
            {"name": "gene4"}, 
            {"name": "gene5"}]
        self.f1 = glbase3.genelist()
        self.f1.load_list(fake1)
        
        fake2 = [{"name": "gene4", "alt_key": "meh"}, 
            {"name": "gene5"}, 
            {"name": "gene6"}, 
            {"name": "gene7"}, 
            {"name": "gene8"},
            {"name": "gene9"}] 
        self.f2 = glbase3.genelist()
        self.f2.load_list(fake2)
        
        fake3 = [{"name": "gene4", "alt_key": False}, 
            {"name": "gene4", "alt_key": True}, 
            {"name": "gene5", "alt_key": False}, 
            {"name": "gene6"}, 
            {"name": "gene7"}]
        self.f3 = glbase3.genelist()
        self.f3.load_list(fake3)

    def test_load(self):
        self.assertEqual(len(self.a), 10) # len() reports are not *zero* based
        self.assertEqual(len(self.b), 7)
        self.assertEqual(len(self.c), 6) # tests format=sniffer
        self.assertEqual(len(self.d), 114) # tests skiplines = -1
        self.assertEqual(len(self.e), 16) # tests a much more complex format.
    
    def test_collisions(self):
        # only identical points:
        self.assertEqual(len(self.a.collide(genelist=self.b, delta=0)), 5) # first three worked out by hand.
        self.assertEqual(len(self.a.collide(genelist=self.c, delta=0)), 5)
        self.assertFalse(self.b.collide(genelist=self.c, delta=0)) # no overlaps are now None
        self.assertFalse(self.d.collide(genelist=self.e, delta=0)) # no overlaps are now None

    def test_collisions_logic(self):
        self.assertEqual(len(self.a.collide(genelist=self.b, delta=0, logic="notinright")),5) # i.e. only if items in self.a are not in self.b
        self.assertEqual(len(self.a.collide(genelist=self.b, delta=0, logic="notinleft")), 2)
        self.assertEqual(len(self.b.collide(genelist=self.a, delta=0, logic="notinright")), 2)
        self.assertEqual(len(self.b.collide(genelist=self.a, delta=0, logic="notinleft")), 5)

    def test_collisions_delta(self):
        self.assertEqual(len(self.a.collide(genelist=self.b, delta=15)), 7) # first three worked out by hand.
        self.assertEqual(len(self.a.collide(genelist=self.c, delta=15)), 15)
        self.assertEqual(len(self.b.collide(genelist=self.c, delta=15)), 2)
        self.assertFalse(self.d.collide(genelist=self.e, delta=500)) # no overlaps are now None

    def test_overlap(self):
        self.assertEqual(len(self.a.overlap(genelist=self.b, delta=0)), 7) # first three worked out by hand.
        self.assertEqual(len(self.a.overlap(genelist=self.c, delta=0)), 6)
        self.assertEqual(len(self.b.overlap(genelist=self.c, delta=0)), 1)

    def test_overlap_delta(self):
        self.assertEqual(len(self.a.overlap(genelist=self.b, delta=15)), 9) # zero based
        self.assertEqual(len(self.a.overlap(genelist=self.c, delta=15)), 19) 
        self.assertEqual(len(self.b.overlap(genelist=self.c, delta=15)), 3) 
        self.assertFalse(self.d.overlap(genelist=self.e, delta=500))  # no overlaps are now None

		# funny mac bug/feature here where linux gets \n, mac gets \r\n...
		# this may fail on Windows too...
		# Currently it works on Mac...
        #self.assertEqual(str(self.a.overlap(genelist=self.b, delta=0)).replace("\n", "\r\n"), "0: loc: chr1:100-200, name: Loc5\r\n1: loc: chr1:150-200, name: Loc6\r\n2: loc: chr1:1000-2000, name: Loc7\r\n... truncated, showing 3/7\r\n6: loc: chr1:10000-20000, name: Loc9")

		#0: loc: chr1:100-200, name: Loc5\n1: loc: chr1:150-200, name: Loc6\n2: loc: chr1:1000-2000, name: Loc7\n... truncated, showing 3/7\n6: loc: chr1:10000-20000, name: Loc9'
		#0: loc: chr1:100-200, name: Loc5\r\n1: loc: chr1:150-200, name: Loc6\r\n2: loc: chr1:1000-2000, name: Loc7\r\n... truncated, showing 3/7\r\n6: loc: chr1:10000-20000, name: Loc9

    def test_raw_qcollide(self):
        # test the raw underlying utils.qcollide() for robustness
        self.assertEqual(qcollide(10, 20, 11, 12), True) # internal
        self.assertEqual(qcollide(10, 20, 15, 25), True) # right overhang
        self.assertEqual(qcollide(10, 20, 5, 15), True) # left overhang
        self.assertEqual(qcollide(10, 20, 5, 25), True) # outdies
        self.assertEqual(qcollide(10, 20, 10, 10), True) # edge left
        self.assertEqual(qcollide(10, 20, 20, 20), True) # edge right
        self.assertEqual(qcollide(10, 20, 21, 21), False) # edge right
        self.assertEqual(qcollide(10, 20, 9, 9), False) # edge right
        self.assertEqual(qcollide(10, 20, 10, 20), True) # perfect match
        
    def test_map(self):
        o = self.f1.map(genelist=self.f2, key="name")
        self.assertEqual(len(o), 2)
        
    def test_map_not(self):
        o = self.f1.map(genelist=self.f2, key="name", logic="notright")
        self.assertEqual(len(o), 4)
        o = self.f2.map(genelist=self.f1, key="name", logic="notright")
        self.assertEqual(len(o), 3)
        
    def test_map_lazygreedy(self):
        # default is greedy:
        o = self.f2.map(genelist=self.f3, key="name", greedy=False) # greedy=True is default behaviour
        self.assertEqual(len(o), 5)
        o = self.f2.map(genelist=self.f3, key="name", greedy=True)
        self.assertEqual(len(o), 5)

    def test_map_keymerge(self):
        # Keys must inherit from the right
        o = self.f2.map(genelist=self.f3, key="name")
        self.assertFalse(o[0]["alt_key"])

        o = self.f2.map(genelist=self.f3, key="name", greedy=True)
        self.assertTrue(o[1], 5)

    def test_collisions_add_tags(self):
        res = self.a.collide(genelist=self.b, delta=0, add_tags="score") # first three worked out by hand.
        
        self.assertEqual(res[0]["score"], 7.0)
        self.assertEqual(res[1]["score"], 9.0)
        self.assertEqual(res[2]["score"], 11.0)
    
    def test_buckets(self):
        glbase3.config.bucket_size = 100 # change to a smaller value for testing purposes.

        g = glbase3.genelist()
        
        data = [{"loc": glbase3.location(loc="chr1:1000-1200")},
            {"loc": glbase3.location(loc="chr1:1200-1300")},
            {"loc": glbase3.location(loc="chr1:1200-1201")},
            {"loc": glbase3.location(loc="chr1:1300-1400")},
            {"loc": glbase3.location(loc="chr1:1400-1500")},
            {"loc": glbase3.location(loc="chr1:1500-1600")},
            {"loc": glbase3.location(loc="chr1:1600-1600")}, # point locs on edges of buckets
            {"loc": glbase3.location(loc="chr1:1423-1423")}, # point locs in middle of buckets
            {"loc": glbase3.location(loc="chr1:0-1500")}] # span much larger than bucket 
        
        g.load_list(data)
       
        left_buck = int((1299-1)/glbase3.config.bucket_size)*glbase3.config.bucket_size
        right_buck = int((1788)/glbase3.config.bucket_size)*glbase3.config.bucket_size
        buckets_reqd = list(range(left_buck, right_buck+glbase3.config.bucket_size, glbase3.config.bucket_size)) # make sure to get the right spanning and left spanning sites
        
        loc_ids = set()
        if buckets_reqd:
            for buck in buckets_reqd:
                if buck in g.buckets["1"]:
                    loc_ids.update(g.buckets["1"][buck]) # unique ids
                
        self.assertSetEqual(loc_ids, set([0, 1, 2, 3, 4, 5, 6, 7, 8]))
        self.assertEqual(len(g.buckets), 1)
        self.assertEqual(len(g.buckets["1"]), 17)
        
        glbase3.config.bucket_size = 10000 # change it back

    def test_remove_dupes_by_loc(self):
        data = [{"loc": glbase3.location(loc="chr1:1000-1200")},
            {"loc": glbase3.location(loc="chr1:1000-1200")},
            {"loc": glbase3.location(loc="chr1:1100-1200")},
            {"loc": glbase3.location(loc="chr1:1300-1400")},
            {"loc": glbase3.location(loc="chr1:1300-1400")},
            {"loc": glbase3.location(loc="chr1:1300-1400")},
            {"loc": glbase3.location(loc="chr1:1600-1600")},
            {"loc": glbase3.location(loc="chr1:1423-1423")},
            {"loc": glbase3.location(loc="chr2:1000-1200")}]   
        
        g = glbase3.genelist()            
        g.load_list(data)
        
        newl = g.removeDuplicatesByLoc(delta=100, mode='pointify_expand')        
        
        self.assertEqual(len(newl), 4)
        
    
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Collisions_Overlaps)
    unittest.TextTestRunner(verbosity=2).run(suite)

"""

Tester Suite:

**Purpose**

This one checks collisions and overlaps using the three toy testN.csvs.

"""

import unittest

# get glbase
import sys, os
sys.path.append(os.path.realpath("../../"))

import glbase
from glbase.utils import qcollide

glbase.config.SILENT = True
glbase.config.set_log_level(None)

class Test_Annotate(unittest.TestCase):
    def setUp(self):
        self.a = glbase.genelist(filename="testA.csv", format=glbase.format.sniffer)
        self.g = glbase.genome(filename="test-genome.csv", format=glbase.format.sniffer)

    def test_annotate(self):
        ann = self.g.annotate(genelist=self.a, key_to_match="loc", closest_only=False,
            distance=10)
        self.assertEqual(len(ann), 6)
    
        ann = self.g.annotate(genelist=self.a, key_to_match="loc", closest_only=False,
            distance=5)
        self.assertEqual(len(ann), 2)
    
    def test_annotate_closest_only(self):
        ann = self.g.annotate(genelist=self.a, key_to_match="loc", closest_only=True, distance=10)
        self.assertEqual(len(ann), 4)
    
    def test_annotate_keymerge(self):
        ann = self.g.annotate(genelist=self.a, key_to_match="loc", closest_only=False, distance=10)
        self.assertEqual(ann[0]["score"], 1) # score from self.a, despite overlapping names
        self.assertEqual(ann[0]["extra_key"], "B")

        ann = self.g.annotate(genelist=self.a, key_to_match="loc", closest_only=True, distance=10)
        self.assertEqual(ann[0]["score"], 1) # score from self.a, despite overlapping names
        self.assertEqual(ann[0]["extra_key"], "C")

    def test_annotate_image_draw(self):
        ann = self.g.annotate(genelist=self.a, key_to_match="loc", distance=10, window=1,
            closest_only=False,
            image_filename="test_images/annotate-moving-average.png")
        ann = self.g.annotate(genelist=self.a, key_to_match="loc", distance=10, window=1, 
            closest_only=True,
            image_filename="test_images/annotate-moving-average.png")
        
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Annotate)
    unittest.TextTestRunner(verbosity=2).run(suite)

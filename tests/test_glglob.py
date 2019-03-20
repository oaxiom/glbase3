"""

Tester Suite:

**Purpose**

This one checks glglob (replaces glglob_test.py)

"""

import unittest, numpy

# get glbase
import sys, os
sys.path.append(os.path.realpath("../../"))

import glbase3

glbase3.config.SILENT = True
glbase3.config.set_log_level(None)

class Test_glglob(unittest.TestCase):
    def setUp(self):
        # get some data;
        self.data1 = glbase3.genelist(filename="test_data/testA.csv", format={'loc': 0, 'name':1, 'score': 2, 'skiplines': 0})
        self.data2 = glbase3.genelist(filename="test_data/testB.csv", format={'loc': 0, 'name':1})
        self.data3 = glbase3.genelist(filename="test_data/testC.csv", format={'loc': 0, 'name':1})
        #self.data4 = glbase3.genelist(filename="test_data/ccat_list.region", format=glbase3.format_ccat_output)
        print(self.data1)
        self.g = glbase3.glglob(self.data1, self.data2, self.data3, type="peaklist")

    def test_chip_seq_cluster_heatmap_error(self):
        no_loc_gl = glbase3.genelist()
        no_loc_gl.load_list([{'name': 'missing'}, {'name': 'a'}, {'name': 'loc'}, {'name': 'key'}])

        self.assertRaises(ValueError, self.g.chip_seq_cluster_heatmap, [self.data1, self.data2, self.data3], []) # Fails at a differnet stage, but passes the assertion
        self.assertRaises(glbase3.errors.AssertionError, self.g.chip_seq_cluster_heatmap, [self.data1, self.data2, no_loc_gl], [])
        self.assertRaises(glbase3.errors.AssertionError, self.g.chip_seq_cluster_heatmap, [self.data1, no_loc_gl, no_loc_gl], [])



if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_glglob)
    unittest.TextTestRunner(verbosity=2).run(suite)


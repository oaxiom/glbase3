"""

Tester Suite:

**Purpose**

This one checks glglob (replaces glglob_test.py)

"""

import unittest, numpy

# get glbase
import sys, os
sys.path.append(os.path.realpath("../../"))

import glbase

glbase.config.SILENT = True
glbase.config.set_log_level(None)

class Test_glglob(unittest.TestCase):
    def setUp(self):
        # get some data;
        self.data1 = glbase.peaklist(filename="testA.csv")
        self.data2 = glbase.peaklist(filename="testB.csv")
        self.data3 = glbase.peaklist(filename="testC.csv")
        self.data4 = glbase.peaklist(filename="../example/ccat_list.region", format=glbase.format_ccat_output)
        self.g = glbase.glglob(self.data1, self.data2, self.data3, self.data4, type="peaklist")

	"""
	# This test is obselete now as I use Scipy to calculate?
    def test_compare_euclidean_collide(self):
        r = self.g.compare(key="loc", filename="matrix_compare.png", method="collide", distance="euclidean")

        result_matrix = numpy.array(
            [[   0.       ,    21.9544984,    17.02938637,  230.24769271],
             [  21.9544984 ,    0.        ,   17.49285568,  228.45130772],
             [  17.02938637,   17.49285568,    0.        ,  229.63884689],
             [ 230.24769271,  228.45130772,  229.63884689,    0.        ]])# calculated from R

        array_equal = True
        for ri, row in enumerate(result_matrix):
            for index, item in enumerate(row):
                # have to use own compare here as the == compares bits and for some reason I am last three sig nums off.
                if abs(item - r[ri][index]) > 1e-6: # 6 sig figs precision
                    array_equal = False
                    break
        self.assert_(array_equal)
	"""
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_glglob)
    unittest.TextTestRunner(verbosity=2).run(suite)


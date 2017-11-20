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
        self.data1 = glbase3.peaklist(filename="test_data/testA.csv")
        self.data2 = glbase3.peaklist(filename="test_data/testB.csv")
        self.data3 = glbase3.peaklist(filename="test_data/testC.csv")
        self.data4 = glbase3.peaklist(filename="test_data/ccat_list.region", format=glbase3.format_ccat_output)
        self.g = glbase3.glglob(self.data1, self.data2, self.data3, self.data4, type="peaklist")

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


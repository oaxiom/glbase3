"""
track.py tester code.

Part of glbase

Tests that track is performing accurately.

TODO:
-----
. some of the tests are not very extensive yet.

"""

import unittest, numpy

# get glbase
import sys, os
import glbase3.utils as utils

class Test_Utils(unittest.TestCase):
    def test_scale_data(self):
        # equal, simple case:
        data = [1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0]
        res = utils.scale_data(data, range=(0,2))
        exp = [0.16666667, 0.33333333]
        self.assertTrue(False not in [int(x*10000) == int(y*10000) for x, y in zip(res, exp)])

    def test_movingAverage(self):
        data = list(range(21)) + list(range(20,0,-1))
        #print data
        ma = utils.movingAverage(data, window=4)
        expected = [6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 74, 77, 78, 77, 74, 70, 66, 62, 58, 54, 50, 46, 42, 38, 34, 30, 26, 22, 18, 14]
        self.assertTrue(False not in [x == y for x, y in zip(list(ma[1]), list(expected))])

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Utils)
    unittest.TextTestRunner(verbosity=2).run(suite)

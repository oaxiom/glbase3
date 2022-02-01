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

class Test_MassSpec(unittest.TestCase):
    def setUp(self):
        pass

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_GeneList)
    unittest.TextTestRunner(verbosity=2).run(suite)

"""

Tester Suite:

this one tests the draw class, makes sure all of the draw calls
are working correctly.

"""

import unittest

# get glbase
import sys, os
sys.path.append(os.path.realpath("../"))
from draw import draw

class Test_Draw(unittest.TestCase):
    def setUp(self):
        """
        start up and assign draw.
        """
        self.draw = draw(self)

    def test_vennDiagrams(self):
        """
        draw a series of venn diagrams.

        Not sure this is entirely satisfactory
        for a testcase as I can't test the output very
        easily...
        """
        labels = {"left": "Left label", "right": "Right label", "title": "This is the title"}
        res = self.draw._vennDiagram2(10,2,20, filename="test_venn.png", labels=labels)
        self.failIf(not res)
        # You have to test these by eye.
        self.draw._vennDiagram2(1000, 1000, 500, filename="messing_eq.png", proportional=True, labels=labels)
        self.draw._vennDiagram2(500, 1000, 1000, filename="messing_big_over.png", proportional=True, labels=labels)
        self.draw._vennDiagram2(500, 1000, 20, filename="messing_small_over.png", proportional=True, labels=labels)
        self.draw._vennDiagram2(1000, 10, 20, filename="messing_ludicrous_left.png", proportional=True, labels=labels)
        self.draw._vennDiagram2(100, 100, 1, filename="messing_small_centre.png", proportional=True, labels=labels)
        self.draw._vennDiagram2(10, 10, 100, filename="messing_big_centre.png", proportional=True, labels=labels)

        # all pass except ludicrous_left.

#-----------------------------------------------------------------------
# run

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Draw)
    unittest.TextTestRunner(verbosity=3).run(suite)


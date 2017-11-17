"""

Tester Suite:

**Purpose**

This one checks the delayedlist functions.
Mainly I'm testing to make sure it behaves similarly (within reason)
to genelist

"""

import unittest

# get glbase
import sys, os
sys.path.append(os.path.realpath("../../"))

import glbase

class Test_BarChart(unittest.TestCase):
    def setUp(self):
        data = [{"name": "Nanog", "items": 2},
            {"name": "Nanog", "items": 3},
            {"name": "Nanog", "items": 4},
            {"name": "Pou5f1", "items": 5},
            {"name": "Pou5f1", "items": 6},
            {"name": "Sox2", "items": 7}]
            
        self.a = glbase.genelist()
        self.a.load_list(data)
        
    def test_frequency_bar_chart(self):
        self.a.frequency_bar_chart(filename="test_images/frequency_bar_chart_test_percent.png", key="name")
        self.a.frequency_bar_chart(filename="test_images/image_frequency_bar_chart_test.png", key="name", percents=False)
        
    def test_bar_chart(self):
        self.a.bar_chart(filename="test_images/image_bar_chart_test_percents.png", labels="name", data="items", percents=True)
        self.a.bar_chart(filename="test_images/image_bar_chart_test.png", labels="name", data="items")

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_BarChart)
    unittest.TextTestRunner(verbosity=2).run(suite)

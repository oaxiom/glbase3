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

import glbase3

format = {"array_systematic_name": 0, "entrez": 1, "refseq": 2, "name": 3,
        "GFP": 4, "Mash": 5, "skiplines": 0}

class Test_Delayedlist(unittest.TestCase):
    def setUp(self):
        self.a = glbase3.delayedlist(filename="test_data/array_data.csv", format=format) # although I don't actually need this at all.
        spoof_gl = [{"name": "Lypla1"}, {"name": "Pdia4"}]
        self.b = glbase3.genelist()
        self.b.load_list(spoof_gl)

    def test_len(self):
        self.b = glbase3.genelist(filename="test_data/array_data.csv", format=format)
        self.a = glbase3.delayedlist(filename="test_data/array_data.csv", format=format)
        self.assertEqual(len(self.a), len(self.b))

    def test_gzipped_delayedlist(self):
        self.b = glbase3.genelist(filename="test_data/array_data.csv", format=format)
        self.a = glbase3.delayedlist(filename="test_data/array_data.csv.gz", format=format, gzip=True)
        self.assertEqual(len(self.a), len(self.b))

        for item in self.a:
            self.assertEqual(item["name"], "Lypla1")
            self.assertEqual(item["array_systematic_name"], 'scl000965.1_10-S')
            break

    def test_iteration(self):
        for item in self.a:
            self.assertEqual(item["name"], "Lypla1")
            self.assertEqual(item["array_systematic_name"], 'scl000965.1_10-S')
            break

        #self.assertEqual()
        self.a.reset()

        for item in self.a.linearData: # test internal access
            self.assertEqual(item["name"], "Lypla1")
            self.assertEqual(item["array_systematic_name"], 'scl000965.1_10-S')
            break

    def test_map_assertion(self):
        self.a.reset()
        res = self.b.map(genelist=self.a, key='name') # only support map() this way around.
        newl = []
        for i in res:
            newl.append(i)
        self.assertEqual(newl[0]["name"], "Lypla1")
        # delayedlist must be the right hand side
        # Will throw a qkeyfind error.

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Delayedlist)
    unittest.TextTestRunner(verbosity=2).run(suite)

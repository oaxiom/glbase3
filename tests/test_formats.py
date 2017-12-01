"""

Tester Suite:

this one tests the location class, makes sure it's working
as advertised.

"""

import sys, os, unittest

sys.path.append(os.path.realpath("../"))

import glbase3 as gl

class Test_Formats(unittest.TestCase):
    def setUp(self):
        # I spook a nice neat fc easier for testing.
        self.f = gl.fc("f",
            format={"name": 0, "description": 1, "strand": 2, "force_tsv": True},
            description="A simple format")
        
    def test_format_containers(self):
        self.assertDictEqual(dict(self.f), {'force_tsv': True, 'name': 0, 'strand': 2, 'description': 1})
               
        self.assertEqual(self.f.description, "A simple format")
        
        self.assertEqual(self.f["name"], 0)
        self.assertEqual(self.f["description"], 1)
        self.assertEqual(self.f["strand"], 2)
        
        # test updating
        self.f.update({"loc": "location()"})
        self.assertEqual(self.f["loc"],  "location()")
        
    def test_fc_iterator(self):
        expected_keys = set(["name", "description", "strand", "force_tsv"])
        expected_values = set([0, 1, 2, True])
    
        observed_keys = []
        observed_values = []
        for k in self.f:
            observed_keys.append(k)
            observed_values.append(self.f[k])
        
        self.assertEqual(set(observed_keys), expected_keys)
      
    def test_catalogue(self):
        self.assertEqual(len(gl.format.catalogue), 31) # will need to be updated each time I register a format

    def test_sam_tophat_xs(self):
        newgl = gl.genelist("test_data/test.sam", format=gl.format.sam_tophat_xs)        
        
        self.assertEqual(newgl[0]["loc"], "chr1:3035081-3035081")
        self.assertEqual(newgl[0]["seq"], "AAACATTCCTGGGAACATCTTGACCATAAGATAAAGGGGACTGTGAAGACATAGCAGGGCTATATTATCTAAGTCAACACCATCTGGCCG")
        self.assertEqual(newgl[0]["strand"], "+")
        self.assertEqual(newgl[1]["strand"], "-")
        
        # test it also works for a delayedlist, which is where I'd usually use it:
        newgl = gl.delayedlist("test_data/test.sam", format=gl.format.sam_tophat_xs)    
        
        for index, item in enumerate(newgl):    
            #print item
            if index == 0:
                self.assertEqual(item["loc"], "chr1:3035081-3035081")
                self.assertEqual(item["seq"], "AAACATTCCTGGGAACATCTTGACCATAAGATAAAGGGGACTGTGAAGACATAGCAGGGCTATATTATCTAAGTCAACACCATCTGGCCG")
                self.assertEqual(item["strand"], "+")
            elif index == 1:
                self.assertEqual(item["strand"], "-")
        
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Formats)
    unittest.TextTestRunner(verbosity=3).run(suite)

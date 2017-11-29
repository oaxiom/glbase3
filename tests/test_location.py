"""

Tester Suite:

this one tests the location class, makes sure it's working
as advertised.

"""

import sys, os, unittest

sys.path.append(os.path.realpath("../"))

from location import location

class Test_Location(unittest.TestCase):
    def test_loading(self):
        a = location(loc="Chr1:10-20")
        self.assertEqual(str(a), "chr1:10-20")
        
        a = location(chr="1", left=1000, right=2000)
        self.assertEqual(str(a), "chr1:1000-2000")
        
        a = location(chr=1, left=1000, right=2000)
        self.assertEqual(str(a), "chr1:1000-2000")
        
        a = location(chr="X", left=1000, right=2000)
        self.assertEqual(str(a), "chrX:1000-2000")
        
        a = location(loc="Chr2_RANDOM:100-200") # should still load...
        self.assertEqual(str(a), "chr2_RANDOM:100-200")

    def test_pass_throughs(self):
        # pass a var through and check it is properly copied/modifiable.
        a = location(loc="chr1:100-200")
        b = location(loc=a) # this will copy.
        c = location(loc=a)
        a = a.pointify()
        b = b.expand(5)
        self.assertEqual(str(a), "chr1:150-150")
        self.assertEqual(str(b), "chr1:95-205")
        self.assertEqual(str(c), "chr1:100-200") # should still be original.

        b = location(chr=c["chr"], left=c["left"], right=c["right"])
        b = b.expand(5)
        self.assertEqual(str(c), "chr1:100-200")
        self.assertEqual(str(b), "chr1:95-205")

        a = b.expand(10) # b should be left untouched.
        self.assertEqual(str(b),"chr1:95-205")
        self.assertEqual(str(a),"chr1:85-215")

    def test_location_modifications(self):
        a = location(loc="chr1:1000-2000")

        a = a.expand(100) # answer = chr1:900-2100
        self.assertEqual(str(a), "chr1:900-2100")

        a = a.expandLeft(10)# answer = chr1:890-2100
        self.assertEqual(str(a), "chr1:890-2100")

        a = a.expandRight(10) # answer = chr1:890-2110
        self.assertEqual(str(a), "chr1:890-2110")

        a = a.shrinkLeft(10) # answer = chr1:900-2110
        self.assertEqual(str(a), "chr1:900-2110")

        a = a.shrinkRight(10) # answer = chr1:900-2100
        self.assertEqual(str(a), "chr1:900-2100")

        a = a.shrink(100) # should be back where it started answer = chr1:1000-2000
        self.assertEqual(str(a), "chr1:1000-2000")

        a = a.pointify() # get the middle # answer = chr1:1500-1500
        self.assertEqual(str(a), "chr1:1500-1500")

        a = a.expand(100) # and make a 200bp window answer = chr1:1400-1600
        self.assertEqual(str(a), "chr1:1400-1600")

    def test_location_access(self):
        a = location(loc="chr1:1000-2000")

        t = a["chr"] # emulate print
        self.assertEqual(t, "1")

        t = a["left"]
        self.assertEqual(t, 1000)

        t = a["right"]
        self.assertEqual(t, 2000)

        t = a["string"]
        self.assertEqual(t, "chr1:1000-2000")

        t = str(sorted(a["dict"])) # emulate a print
        self.assertEqual(t, str(sorted({'chr': '1', 'right': 2000, 'left': 1000}))) # this is may be wrong as dicts can be unordered

        t = len(a) # should be the length of the location string. = span between the two elements
        self.assertEqual(t, 1000)

        t = a.split(":") # the value is actually ignored.
        self.assertEqual(t, ('1', 1000, 2000))

        t = a.split() # extra syntax, not available for a genuine string
        self.assertEqual(t, ('1', 1000, 2000))

        t = a.split("-") # also ignored
        self.assertEqual(t, ('1', 1000, 2000))

        t = repr(a) # the debug output
        self.assertEqual(t, "<location chr1:1000-2000>")

    def test_location_qcollide(self):
        a = location(chr=1, left=10, right=20)
        # test the raw underlying loc.qcollide() for robustness
        self.assertEqual(a.qcollide(location(loc="chr1:11-12")), True) # internal
        self.assertEqual(a.qcollide(location(loc="chr1:15-25")), True) # right overhang
        self.assertEqual(a.qcollide(location(loc="chr1:5-15")), True) # left overhang
        self.assertEqual(a.qcollide(location(loc="chr1:5-25")), True) # outdies
        self.assertEqual(a.qcollide(location(loc="chr1:10-10")), True) # edge left
        self.assertEqual(a.qcollide(location(loc="chr1:20-20")), True) # edge right
        self.assertEqual(a.qcollide(location(loc="chr1:21-21")), False) # edge right
        self.assertEqual(a.qcollide(location(loc="chr1:9-9")), False) # edge right
        
        # qcollide() is not the expected method. These should fail as qcollide
        # does not check chromosome.
        self.assertEqual(a.qcollide(location(loc="chr2:12-16")), True) # This will pass for qcollide()
        self.assertEqual(a.qcollide(location(loc="chr2:8-9")), False) # chromosome and loc fail

    def test_location_collide(self):
        a = location(chr=1, left=10, right=20)
        # test the raw underlying loc.qcollide() for robustness
        self.assertEqual(a.collide(location(loc="chr1:11-12")), True) # internal
        self.assertEqual(a.collide(location(loc="chr1:15-25")), True) # right overhang
        self.assertEqual(a.collide(location(loc="chr1:5-15")), True) # left overhang
        self.assertEqual(a.collide(location(loc="chr1:5-25")), True) # outdies
        self.assertEqual(a.collide(location(loc="chr1:10-10")), True) # edghe left
        self.assertEqual(a.collide(location(loc="chr1:20-20")), True) # edghe right
        self.assertEqual(a.collide(location(loc="chr1:21-21")), False) # edghe right
        self.assertEqual(a.collide(location(loc="chr1:9-9")), False) # edghe right
        self.assertEqual(a.collide(location(loc="chr2:12-16")), False) # chromosome fail
        self.assertEqual(a.collide(location(loc="chr2:8-9")), False) # chromosome and loc fail

    def test_location_equality_tests(self):
        a = location("chr1:1000-2000")
        b = location("chr1:1001-2000")
        c = location("chr1:1000-2001")
        d = location("chr2:1000-2000")
        e = location("chr1:1000-2000")

        self.assertEqual(a == b, False) # internal
        self.assertEqual(a == c, False) # internal
        self.assertEqual(a == d, False) # internal
        self.assertEqual(a == e, True) # internal
        self.assertEqual(a == b, False) # internal
        self.assertEqual(a == "chr1:1000-2000", True)
        self.assertEqual(a == "chr1:1000-2001", False)

    def left_right_pointify(self):
        a = location("chr1:1000-2000")
        l = a.pointLeft()
        r = a.pointRight()
        
        self.assertEqual(l, "chr1:1000-1000")
        self.assertEqual(r, "chr1:2000-2000")

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Location)
    unittest.TextTestRunner(verbosity=3).run(suite)

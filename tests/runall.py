"""

run all of the tests

"""
import unittest

# get glbase
import sys, os
sys.path.append(os.path.realpath("../../"))
import glbase3.config
print(sys.version)
glbase3.version()
glbase3.config.set_log_level(None) # Silence output

# please maintain these in alphabetical order.
from test_annotate import Test_Annotate
from test_collisions import Test_Collisions_Overlaps
from test_delayedlist import Test_Delayedlist
from test_draw import Test_Draw
from test_ecrbase import Test_EcrBase_Interface
from test_expression import Test_Expression
from test_fastq import Test_Fastq
from test_flats import Test_Flat_Function
from test_formats import Test_Formats
from test_genelist import Test_GeneList
from test_glglob import Test_glglob
from test_location import Test_Location
from test_track import Test_Track_Function
from test_tsv_csv_behaviour import Test_TSVCSV
from test_utils import Test_Utils
from test_genome import Test_Genome

def get_suite():
    suite = [
        unittest.TestLoader().loadTestsFromTestCase(Test_Collisions_Overlaps),
        unittest.TestLoader().loadTestsFromTestCase(Test_Location),
        unittest.TestLoader().loadTestsFromTestCase(Test_Draw),
        unittest.TestLoader().loadTestsFromTestCase(Test_glglob),
        unittest.TestLoader().loadTestsFromTestCase(Test_GeneList),
        unittest.TestLoader().loadTestsFromTestCase(Test_Track_Function),
        unittest.TestLoader().loadTestsFromTestCase(Test_Flat_Function),
        unittest.TestLoader().loadTestsFromTestCase(Test_Delayedlist),
        unittest.TestLoader().loadTestsFromTestCase(Test_Expression),
        unittest.TestLoader().loadTestsFromTestCase(Test_EcrBase_Interface),
        unittest.TestLoader().loadTestsFromTestCase(Test_Utils),
        unittest.TestLoader().loadTestsFromTestCase(Test_Annotate),
        unittest.TestLoader().loadTestsFromTestCase(Test_TSVCSV),
        unittest.TestLoader().loadTestsFromTestCase(Test_Formats),
        unittest.TestLoader().loadTestsFromTestCase(Test_Fastq),
        unittest.TestLoader().loadTestsFromTestCase(Test_Genome)
        ]

    return unittest.TestSuite(suite)

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=2).run(get_suite())

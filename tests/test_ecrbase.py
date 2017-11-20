"""

Tester Suite:

**Purpose**

Check the ecrbase readers

"""

import unittest

# get glbase
import sys, os
sys.path.append(os.path.realpath("../../"))

import glbase3
from glbase3.utils import qcollide

glbase3.config.SILENT = True
glbase3.config.set_log_level(None)

class Test_EcrBase_Interface(unittest.TestCase):
    def setup(self):
        self.e = ecrbase(path=".")

    def test_tfbs_iter(self):
        i = glbase3.tfbs_iter("test_data/tfbs_ecrs.mm9hg18.v102_top100_lines.txt")

        two_locs = {"chr1:3017653-3017842": "['CDXA_02', 'SRY_01', 'SRY_02', 'OCT1_06', 'NFY_C', 'EGR1_01', 'EGR2_01', 'XFD3_01', 'AML1_01', 'MIF1_01', 'TCF11MAFG_01', 'ATATA_B', 'PAX3_B', 'DTYPEPA_B', 'OCT1_B', 'GATA2_02', 'GATA2_03', 'GATA3_02', 'MEIS1AHOXA9_01', 'MEIS1BHOXA9_02', 'IPF1_Q4']",
            "chr1:3017653-3017842": "['E2F_02', 'COMP1_01', 'TAL1BETAE47_01', 'DELTAEF1_01', 'FOXD3_01', 'SRY_01', 'POLY_C', 'EGR1_01', 'NGFIC_01', 'RFX1_01', 'APOLYA_B', 'LDSPOLYA_B', 'PAX4_02', 'PAX4_04', 'MEIS1AHOXA9_01', 'NKX61_01', 'E2F_Q3', 'E2F_Q4', 'E2F_Q6', 'E2F1_Q3', 'E2F1_Q4', 'E2F1_Q6', 'FAC1_01', 'FOXO4_01', 'FOXO1_01', 'PAX2_02', 'CRX_Q4', 'HSF_Q6', 'SP3_Q3', 'E2F1DP1_01', 'E2F1DP2_01', 'E2F4DP1_01', 'E2F4DP2_01', 'E2F1DP1RB_01', 'HFH4_01', 'LEF1_Q2', 'E2F_Q3_01', 'E2F_Q4_01', 'E2F_Q6_01', 'E2F1_Q4_01', 'E2F1_Q6_01', 'DR3_Q4', 'TBX5_01', 'HSF1_Q6', 'CMAF_01', 'SZF11_01']"}

        for tfbs in i:
            if tfbs["motifs"]:
                if str(tfbs["loc"]) in two_locs:
                    self.assertEqual(str(tfbs["motifs"]), str(two_locs[str(tfbs["loc"])]))

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_EcrBase_Interface)
    unittest.TextTestRunner(verbosity=2).run(suite)

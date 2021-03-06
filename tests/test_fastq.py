"""
FASTQ SE and PE splitter loader tester

TODO:
-----
. not in test suite

"""

import unittest, numpy

# get glbase
import sys, os
import glbase3

class Test_Fastq(unittest.TestCase):
    def test_splitPE(self):
        """
        @HWI-ST507:76:A81MKNABXX:7:1:5257:1711:1
        CTCCTAGAGGGAAATATGGAGCAATTACATATTGTTCTCTAGGA
        +
        *625229@@@@@@@@@@@@7;@@;@@@@@;@@@@=@7@;@7@##
        @HWI-ST507:76:A81MKNABXX:7:1:5257:1711:2
        CAGGAGGGTCTGTGGTAGAAGGCTGTTACATACATAATAAA
        +
        HHHHHHHHCHBEEEE9EEBEEEDCECFBFBCB>?ACC>C##
        """
    
        fq = glbase3.fastq("test_data/fastq_typical_data.fq", "phred33")        
        fq.splitPE("/tmp/out1.fq", "/tmp/out2.fq")

        # Uh... It works, but is not a real test.
        oh1 = open("/tmp/out1.fq", "r")
        id1 = oh1.readline().strip()
        oh2 = open("/tmp/out2.fq", "r")
        id2 = oh2.readline().strip()
                
        self.assertEqual(id1[:-1], id2[:-1])

    def test_splitPE2(self):
        fq = glbase3.fastq("test_data/fastq_typical_data2.fq", "phred64")    
        """
        This interleaved file was also seen in the wild:
        
        @HWI-ST507:91:D062KACXX:1:1101:1242:2186:1:1:0:TAGCTT:1
        TTCAAAGGGGGACAGTCCTTGTGGAAGTTTTGAAAGTGGGTTTCATTCACTATGTGACCTTGGCTGTCTTGGAACCCACTATATAGACCAGGCTCGCCTT
        +
        bbbeeeeegfggghhhhhhhhfghgghefghhhhhhfghhagfghhhhhhhhhhdghhhhhgggeeceedecbdd`accccccdcccccccccccaccca
        @HWI-ST507:91:D062KACXX:1:1101:1242:2186:2:1:0:TAGCTT:2
        CAGAAACTAAGAGTGACAGATGCCATTCAGATACATCGTAAGGGTGTGTACATACCAGTAATCCCAACACTGGGAAGACAGAGCAGAAGGATGTCCTAAG
        +
        bbbeeeeeggfggfghhhhhhhhhhhhhhhhhhhhhhhfhhhhh^ccgfhcfffhhghfhhhhhhhhhhhhhgggfdeeeddddddcccccccbcccccb
        """
        
        fq.splitPE("/tmp/out1.fq", "/tmp/out2.fq")

        # Uh... It works, but is not a real test.
        oh1 = open("/tmp/out1.fq", "r")
        id1 = oh1.readline().strip()
        oh2 = open("/tmp/out2.fq", "r")
        id2 = oh2.readline().strip()
                
        self.assertEqual(id1.split(":")[0:6], id2.split(":")[0:6])

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Fastq)
    unittest.TextTestRunner(verbosity=2).run(suite)

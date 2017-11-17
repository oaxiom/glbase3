"""

Tester:

This tester sends a load of mangled csvs to the various loadCSV()s to
check on their error checking.

It has things like csvs in one format and another in another format,
genuinely mangled csvs, binary files, stuff like that.

It's to check that the error reporting makes sense and can correctly
handle bad stuff.

"""

import sys, os
sys.path.append(os.path.realpath("../..")) # only extra path required, adds the parent makes glbase importable from here.
from glbase import *

print(genelist.load.__gui__)
print(genelist.__gui__avail__)

chipseq_result = peaklist(filename="../example/chip_seq_data.csv")
print(chipseq_result)
print(chipseq_result[0]["loc"])
print(repr(chipseq_result[0]["loc"]))
print(chipseq_result[0]["loc"]["left"])

# lest the loading
data1 = peaklist(filename="../example/chip_seq_data.csv")
print(data1)

data2 = peaklist(filename="../example/sissrs_list.txt", format=format_bed)
print(data2)

data3 = peaklist(filename="../example/ccat_list.region")
print(data3)

data4 = peaklist(filename="../example/macs_list.xls")
print(data4)

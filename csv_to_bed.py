
import sys, os, csv, random, copy
sys.path.append(os.path.expanduser("~"))
from glbase import *
user_root = os.path.expanduser("~")

filename=sys.argv[1]
outname= "".join(filename.split(".")[:-1])

a = genelist(filename=filename)
a.saveBED("%s.bed" % outname)
print "%s -> %s.bed" % (filename, outname)

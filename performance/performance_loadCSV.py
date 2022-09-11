import cProfile as profile
import pstats

p = profile.Profile()
p.enable()

import os.path
from glbase3 import *

form = {'force_tsv': True,
    'name': 0, 'enst': 1, 'loc': 3}

try:
    # Need a big enough file to do this:
    f = genelist(os.path.expanduser('~/mm10/mm10_ensembl_v95_enst.tsv'), format=form)

finally:
    p.disable()
    pstats.Stats(p).sort_stats('tottime').print_stats(60)

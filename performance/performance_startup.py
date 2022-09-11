import cProfile as profile
import pstats

p = profile.Profile()
p.enable()

try:
    from glbase3 import *

finally:
    p.disable()
    pstats.Stats(p).sort_stats('tottime').print_stats(60)

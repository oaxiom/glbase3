"""

SOMS are now completely deterministic.

"""

import numpy, sys, os, math
from glbase3 import *
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
config.draw_mode = "svg"

sys.path.append("../sam_annotations/")
import sam_map

expn = glload("1.filter_data/test_data_set.glb")
#expn = glload("1.filter_data/real_data_set.glb") # These data are already log2 and sorted by germ layer
expn.abssub(math.log(10**1.6, 2)) # Without this the SOM seems to get confused by all the non-zero values

expn.som.config(nodenames="ensg", initmethod='fullpca',
    #threshold_value=sam_map.low_expressed_threshold_log2, digitize=12, # I don't use this system as I want to show the MDS on the PCA
    seed=123456, init_whiten=True,
    image_debug='debug/debug', components=13) # >2 will invoke MDS on the PCA
expn.som.train()

# 7 works very nicely, but 13 gives a sweet MDS. Try 13

expn.save('trained_som.glb')

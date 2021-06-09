
"""

SOMS are now completely deterministic.

"""


import numpy, sys, os
from glbase3 import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
config.draw_mode = "pdf"

sys.path.append("../sam_annotations/")
import sam_map
user_path = os.path.expanduser("~")
ensg = glload(os.path.join(user_path, "mm10", "mm10_ensembl_v79_ensg.glb"))

expn = glload('trained_som.glb')

# Build the lineage meaned SOMS
cond_names = expn.getConditionNames()

som_score = 0.2

# Extract all cell-type specific genes:
config.draw_mode = "png"
res = expn.som.threshold_SOM_nodes(som_score, filename="images/cell_lineages_thresholds.png", som_names=expn.getConditionNames(), text_size=8)
config.draw_mode = "pdf"
expn.som.threshold_SOM_nodes(som_score, filename="images/cell_lineages_thresholds.pdf", som_names=expn.getConditionNames())

for som_name in res:
    print()
    print(som_name)
    print(res[som_name])

for som_name in res:
    res[som_name] = ensg.map(genelist=res[som_name], key='ensg')
    res[som_name].saveTSV("cell_lineages/%s.tsv" % som_name.replace("/", "-"))



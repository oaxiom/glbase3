
"""

SOMS are now completely deterministic.

"""


import numpy, pickle, sys, os
from glbase3 import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
config.draw_mode = "pdf"

sys.path.append("../sam_annotations/")
import sam_map

som_score = 0.2

expn = glload('trained_som.glb')

# Build the lineage meaned SOMS
cond_names = expn.getConditionNames()
#expn_numpy_tab = expn.getExpressionTable()
merged_soms = {"max": None}
in_lay = {} # number
for lay in sam_map.gene_layer:
    if lay == "Unassigned" or lay == "Mixed":
        continue

    if lay not in merged_soms:
        merged_soms[lay] = numpy.zeros(expn.som.codebook.shape[0])
        in_lay[lay] = 0

    for ct in sam_map.gene_layer[lay]:
        if ct not in cond_names:
            continue
        cond_index = cond_names.index(ct)
        merged_soms[lay] += expn.som.codebook[:,cond_index]
        in_lay[lay] += 1
    merged_soms[lay] /= in_lay[lay] # normalised for number in layer

    if merged_soms['max'] is None:
        merged_soms["max"] = (merged_soms[lay] - merged_soms[lay].min()) / merged_soms[lay].max()
    else:
        normalised_minmax = (merged_soms[lay] - merged_soms[lay].min()) / merged_soms[lay].max()
        merged_soms["max"] = numpy.maximum(merged_soms["max"], normalised_minmax)

lays = list(merged_soms.keys())
lays.remove("max")
merged_soms["sum"] = numpy.array(sum([merged_soms[k] for k in lays[1:]], merged_soms[lays[0]]))
merged_soms["max"] = numpy.array(merged_soms["max"], dtype=numpy.float) / len(lay)

# Measure the Euclidean_dist between all merged_soms and soms:

merge_som_keys = [k for k in merged_soms.keys() if k != "max" and k != "sum"]
each_lineage = {}
for k in merge_som_keys:
    each_lineage[k] = []
    this_lineage = merged_soms[k]
    for cond_index, som in enumerate(cond_names):
        this_som = expn.som.codebook[:,cond_index]
        euclidean = numpy.sum(numpy.abs(this_lineage - this_som))
        each_lineage[k].append(euclidean)
    # normalise the distances
    each_lineage[k] = numpy.array(each_lineage[k])
    each_lineage[k] = (each_lineage[k] - each_lineage[k].min()) / (each_lineage[k].max() - each_lineage[k].min())

expn_like = []
for i, e in enumerate(cond_names):
    conds = [each_lineage[k][i] for k in merge_som_keys]
    for k in sam_map.gene_layer:
        if e in sam_map.gene_layer[k]:
            layer = k
            break
    lowest = None
    lowest_score = 999999999
    for si, k in enumerate(merge_som_keys):
        if conds[si] < lowest_score:
            lowest = k
            lowest_score = conds[si]
    expn_like.append({"cell_type": e, "conditions": conds, "annotated_layer": layer, "best_layer": lowest})
som_dists = expression(loadable_list=expn_like, cond_names=merge_som_keys)
som_dists.saveTSV("som_dists.tsv", key_order=["cell_type", "annotated_layer", 'best_layer'])
som_dists.save("som_dists.glb")
#som_dists.normalize_rows()
#som_dists.saveTSV("som_dists_row_norm.tsv")

# Extract lineage-type specific SOMs:
config.draw_mode = "png"
somres = expn.som.threshold_SOM_nodes(som_score, filename="images/germ_lineages_thresholds.png",
    alternate_soms=merged_soms, img_number_of_cols=4)

config.draw_mode = "svg"
masks = expn.som.threshold_SOM_nodes(som_score, filename="images/germ_lineages_thresholds.svg",
    alternate_soms=merged_soms, return_masks=True)

# Save the masks for script 4
oh = open('masks.pickle', 'wb')
pickle.dump(masks, oh, -1)
oh.close()
print("Saved pickle!")

# topological lineage-surfaces
X, Y = numpy.mgrid[0:expn.som.mapsize[0], 0:expn.som.mapsize[1]]
#Z = numpy.zeros((msz0,msz0))

for axisNum, k in enumerate(merged_soms):
    mp = merged_soms[k].reshape(expn.som.mapsize[0], expn.som.mapsize[1])#[::-1]
    mp -= mp.max()
    mp = numpy.abs(mp)
    # get the inverse

    fig = plt.figure()
    ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=75, azim=5)
    surf = ax.plot_surface(X, Y, mp,
        cmap=cm.PuOr_r, shade=True,
        cstride=1, rstride=1, linewidth=0, antialiased=False) # antialiased stops transparent edges
    #surf = ax.contourf(X, Y, mp, cmap=cm.PuOr_r)
    ax.set_xlim(0, expn.som.mapsize[0]-1)
    ax.set_ylim(0, expn.som.mapsize[1]-1)
    ax.set_zlim(mp.min(), mp.max())
    ax.set_xticklabels('')
    ax.set_yticklabels('')
    ax.set_zticklabels('')
    fig.savefig("topo/topo_%s.png" % k)
    fig.savefig("topo/topo_%s.pdf" % k)

for som_name in somres:
    somres[som_name].saveTSV("germ_lineages/%s.tsv" % som_name)


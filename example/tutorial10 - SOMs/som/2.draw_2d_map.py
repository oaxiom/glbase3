
"""

SOMS are now completely deterministic.

"""


import numpy, sys
from glbase3 import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

sys.path.append("../sam_annotations/")
import sam_map

expn = glload('trained_som.glb')

config.draw_mode = "png"
expn.som.tree(filename="images/som_tree.png", label_size=4, size=[3,10])
expn.som.plot_gene_density(filename="images/plot_gene_density.png", topological=False)
expn.som.plot_gene_density(filename="images/plot_gene_density_top.png", topological=True)
expn.som.view_map(which_dim='all', pack=True, text_size=7, filename='images/som_test.png')

config.draw_mode = "pdf"
expn.som.tree(filename="images/som_tree.png", label_size=4, size=[3,10])
expn.som.plot_gene_density(filename="images/plot_gene_density.pdf", topological=False)
expn.som.plot_gene_density(filename="images/plot_gene_density_top.pdf", topological=True)
expn.som.view_map(which_dim='all', pack=True, text_size=7, filename='images/som_test.pdf')

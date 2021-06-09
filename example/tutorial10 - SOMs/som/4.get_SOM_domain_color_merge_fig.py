
"""

SOMS are now completely deterministic.

"""

import sys, six, os
import numpy, pickle
from glbase3 import *
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from matplotlib import colors
config.draw_mode = "png"
sys.path.append("../sam_annotations/")
import sam_map

col_dict = {}
for name, hex in six.iteritems(colors.cnames):
    col_dict[name] = colors.colorConverter.to_rgb(hex)

oh = open('masks.pickle', "rb")
masks = pickle.load(oh)
oh.close()

# make an imshow array, with colors for each domain, and grey for shared between >1 domains
shape = masks['Embryonic'].shape[0]
all_mask = numpy.zeros((shape, shape, 3), dtype=float)
used_mask = numpy.zeros((shape, shape), dtype=int)

for k in masks:
    if k in ('max', 'sum'):
        continue
    c = col_dict[sam_map.colour_guide[k]] # int(col_dict[sam_map.colour_guide[k]].replace('#', ''), 16) # get the hex

    for x in range(shape):
        for y in range(shape):
            if masks[k][x,y] == 1.0:
                if used_mask[x, y] == 0:
                    all_mask[x,y,0] = c[0] # Bite me, I am in a rush
                    all_mask[x,y,1] = c[1]
                    all_mask[x,y,2] = c[2]
                    used_mask[x,y] = 1
                else: # Already in another domain
                    all_mask[x,y,0] = 0.7 # grey
                    all_mask[x,y,1] = 0.7 # grey
                    all_mask[x,y,2] = 0.7 # grey

# Draw the color-merged plot:
fig = plot.figure(figsize=[6,6])
ax = fig.add_subplot(111) # Top is nromal SOM, bottom is the threshold
ax.imshow(all_mask, interpolation='nearest',
    extent=[0, shape, 0, shape],
    aspect="auto", origin='lower')
ax.set_xticklabels('')
ax.set_yticklabels('')
ax.set_xlim(0, shape-1)
ax.set_ylim(0, shape-1)
ax.tick_params(top="off", bottom="off", left="off", right="off")
fig.savefig('col_merge.png')
fig.savefig('col_merge.svg')


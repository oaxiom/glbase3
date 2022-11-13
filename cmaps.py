
"""

A small collection of matplotlib colourmaps, primarily segemented maps for discreet heatmaps

From the matplotlib docs:

row i:   x  y0  y1
               /
              /
row i+1: x  y0  y1

Each row in the table for a given color is a sequence of x, y0, y1 tuples.
In each sequence, x must increase monotonically from 0 to 1. For any input
value z falling between x[i] and x[i+1], the output value of a given color
will be linearly interpolated between y1[i] and y0[i+1]:

"""

import matplotlib.cm as cm
import numpy
from matplotlib.colors import LinearSegmentedColormap

def discretize(cmap, N):
    """
    From:
    http://www.scipy.org/Cookbook/Matplotlib/ColormapTransformations

    Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if isinstance(cmap, str):
        cmap = cm.get_cmap(cmap)

    colors_i = numpy.concatenate((numpy.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = numpy.linspace(0, 1., N+1)
    cdict = {
        key: [
            (indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki])
            for i in range(N + 1)
        ]
        for ki, key in enumerate(('red', 'green', 'blue'))
    }

    return LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)

# Do the reverse ones here:
# Um...

if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    a = np.linspace(0, 1, 256).reshape(1,-1)
    a = np.vstack((a,a))

    maps = [discretize("RdBu", 6), discretize("PuBu", 6), discretize("Accent", 12)]
    nmaps = len(maps)

    fig = plt.figure(figsize=(5,10))
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.2, right=0.99)
    for i,m in enumerate(maps):
        ax = plt.subplot(nmaps, 1, i+1)
        plt.axis("off")
        plt.imshow(a, aspect='auto', cmap=plt.get_cmap(m), origin='lower')
        pos = list(ax.get_position().bounds)
        fig.text(pos[0] - 0.01, pos[1], m.name, fontsize=10, horizontalalignment='right')

    fig.savefig("colourmaps.png")

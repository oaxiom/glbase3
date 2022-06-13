"""
draw class for glbase

this is a static class containing various generic methods for drawing etc.

**TODO**

There's a change in direction here. Instead of draw containing lots of generic draw functions
instead its more like a set of wrappers around common ways to do matplotlib stuff.
Drawing inside genelists is fine, as long as it follows this paradigm:

fig = self.draw.getfigure(**kargs)

ax = ...

etc ...

self.draw.do_common_args(fig, **kargs)
filename = fig.savefigure(fig, filename)

It would probably be an improvement if the class part was removed.

Instead a series of methods, exposed by draw.method() at the module level would be better.

This makes them more like helpers for matplotlib than a full fledged object.

This could be easily refactored by changing lines like::

class genelist:
    ... init
        self.draw = draw()

        to

        self.draw = draw

For now, until I refactor the code to remove lines like that.
Also I want to rename this file gldraw to remove name clashes.

Then it can go::

    gldraw.heatmap()
    gldraw.scatter()

"""

import sys, os, copy, random, numpy, math
from collections.abc import Iterable

from numpy import array, arange, mean, max, min, std, float32
from scipy.cluster.hierarchy import distance, linkage, dendrogram
from scipy.spatial.distance import pdist # not in scipy.cluster.hierarchy.distance as you might expect :(
from scipy import polyfit, polyval
from scipy.stats import linregress
import scipy.stats
import numpy as np
import matplotlib
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from cycler import cycler # colours
from matplotlib.colors import ColorConverter, rgb2hex, ListedColormap
import matplotlib.colors as matplotlib_colors
import matplotlib.mlab as mlab
from matplotlib.patches import Ellipse, Circle
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
from .adjustText import adjust_text

from . import config, cmaps, utils
from .errors import AssertionError

# This helps AI recognise the text as text:
matplotlib.rcParams['pdf.fonttype']=42

# this is a work around in the implementation of
# scipy.cluster.hierarchy. It does some heavy
# recursion and even with relatively small samples quickly eats the
# available stack.
# I may need to implement this myself later.
# This can deal with ~23,000 x 14 at least.
# No idea on the upper limit.
sys.setrecursionlimit(5000) # 5x larger recursion.

# define static class here.
class draw:
    def __init__(self, bad_arg=None, **kargs):
        """please deprecate me"""
        pass

    def bracket_data(self,
        data,
        min:int,
        max:int):
        """
        brackets the data between min and max (ie. bounds the data with no scaling)

        This should be a helper?
        """
        ran = max - min
        newd = copy.deepcopy(data)
        for x, row in enumerate(data):
            for y, value in enumerate(row):
                if value < min:
                    newd[x][y] = min
                elif value > max:
                    newd[x][y] = max
        return(newd)

    def heatmap(self,
        filename:str = None,
        cluster_mode:str = "euclidean",
        row_cluster:bool = True,
        col_cluster:bool = True,
        vmin = 0,
        vmax = None,
        colour_map=cm.RdBu_r,
        col_norm:bool = False,
        row_norm:bool = False,
        heat_wid = 0.25,
        heat_hei = 0.85,
        highlights = None,
        digitize:bool = False,
        border:bool = False,
        draw_numbers:bool = False,
        draw_numbers_threshold = -9e14,
        draw_numbers_fmt = '{:.1f}',
        draw_numbers_font_size = 6,
        grid:bool = False,
        row_color_threshold:bool = None,
        col_names:bool = None,
        row_colbar:bool = None,
        col_colbar:bool = None,
        optimal_ordering:bool = True,
        dpi:int = 300,
        _draw_supplied_cell_labels = False,
        **kargs):
        """
        my own version of heatmap.

        This will draw a dendrogram ... etc...

        See the inplace variants as to how to use.

        row_names is very important as it describes the order of the data.
        cluster_mode = pdist method. = ["euclidean"] ??????!

        **Arguments**
            data (Required)
                the data to use. Should be a 2D array for the heatmap.

            filename (Required)
                The filename to save the heatmap to.

            col_norm (Optional, default=False)
                normalise each column of data between 0 .. max => 0.0 .. 1.0

            row_norm (Optional, default=False)
                similar to the defauly output of heatmap.2 in R, rows are normalised 0 .. 1

            row_tree (Optional, default=False)
                provide your own tree for drawing. Should be a Scipy tree. row_labels and the data
                will be rearranged based on the tree, so don't rearrnge the data yourself.
                i.e. the data should be unclustered. Use tree() to get a suitable tree for loading here

            col_tree (Optional, default=False)
                provide your own tree for ordering the data by. See row_tree for details.
                This one is applied to the columns.

            row_font_size or yticklabel_fontsize (Optional, default=guess suitable size)
                the size of the row labels (in points). If set this will also override the hiding of
                labels if there are too many elements.

            col_font_size or xticklabel_fontsize (Optional, default=8)
                the size of the column labels (in points)

            heat_wid (Optional, default=0.25)
                The width of the heatmap panel. The image goes from 0..1 and the left most
                side of the heatmap begins at 0.3 (making the heatmap span from 0.3 -> 0.55).
                You can expand or shrink this value depending wether you want it a bit larger
                or smaller.

            heat_hei (Optional, default=0.85)
                The height of the heatmap. Heatmap runs from 0.1 to heat_hei, with a maximum of 0.9 (i.e. a total of 1.0)
                value is a fraction of the entire figure size.

            colbar_label (Optional, default=None)
                the label to place beneath the colour scale bar

            highlights (Optional, default=None)
                sometimes the row_labels will be suppressed as there is too many labels on the plot.
                But you still want to highlight a few specific genes/rows on the plot.
                Send a list to highlights that matches entries in the row_names.

            digitize (Optional, default=False)
                change the colourmap (either supplied in cmap or the default) into a 'discretized' version
                that has large blocks of colours, defined by the number you send to discretize.

                Note that disctretize only colorises the comlourmap and the data is still clustered on the underlying numeric data.

                You probably want to use expression.digitize() for that.

            imshow (Optional, default=False)
                optional ability to use images for the heatmap. Currently experimental it is
                not always supported in the vector output files.

            draw_numbers (Optional, default=False)
                draw the values of the heatmaps in each cell see also draw_numbers_threshold

            draw_numbers_threshold (Optional, default=-9e14)
                draw the values in the cell if > draw_numbers_threshold

            draw_numbers_fmt (Optional, default= '{:.1f}')
                string formatting for the displayed values

            draw_numbers_font_size (Optional, default=6)
                the font size for the numbers in each cell

            _draw_supplied_cell_labels (Optional, default=False)
                semi-undocumented function to draw text in each cell.

                Please provide a 2D list, with the same dimensions as the heatmap, and this text
                will be drawn in each cell. Useful for tings like drawing a heatmap of expression
                and then overlaying p-values on top of all significant cells.

            col_colbar (Optional, default=None)
                add a colourbar for the samples names. This is designed for when you have too many
                conditions, and just want to show the different samples as colours

                Should be a list of colours in the same order as the condition names

            row_colbar (Optional, default=None)
                add a colourbar for the samples names. This is designed for when you have too many
                conditions, and just want to show the different samples as colours

                Should be a list of colours in the same order as the row names.

                Note that unclustered data goes from the bottom to the top!

            optimal_ordering (Optional, default=True)
                See https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html

        **Returns**
            The actual filename used to save the image.
        """
        assert filename, "heatmap() - no specified filename"

        # preprocess data
        if isinstance(kargs["data"], dict):
            # The data key should be a serialised Dict, I need to make an array.
            data = array([kargs["data"][key] for key in col_names]).T
            # If the lists are not square then this makes a numpy array of lists.
            # Then it will fail below with a strange error.
            # Let's check to make sure its square:
            ls = [len(kargs["data"][key]) for key in col_names]

            if not all(x == ls[0] for x in ls):
                raise Exception("Heatmap data not Square")
        else:
            # the default is a numpy like array object which can be passed right through.
            data = array(kargs["data"], dtype=float32)

        if col_colbar:
            assert len(col_colbar) == data.shape[1], "col_colbar not the same length as data.shape[1]"
        if row_colbar:
            assert len(row_colbar) == data.shape[0], "row_colbar not the same length as data.shape[0]"

        if col_norm:
            for col in range(data.shape[1]):
                data[:,col] /= float(data[:,col].max())

        if row_norm:
            for row in range(data.shape[0]):
                mi = min(data[row,:])
                ma = max(data[row,:])
                data[row,:] = (data[row,:]-mi) / (ma-mi)

        if "square" in kargs and kargs["square"]:
            # make the heatmap square, for e.g. comparison plots
            left_side_tree =          [0.15,    0.15,   0.10,   0.75]
            top_side_tree =           [0.25,    0.90,  0.55,   0.08]
            heatmap_location =        [0.25,    0.15,   0.55,   0.75]
            loc_col_colbar = [] # not supported?
        else:
            # positions of the items in the plot:
            # heat_hei needs to be adjusted. as 0.1 is the bottom edge. User wants the bottom
            # edge to move up rather than the top edge to move down.
            if row_cluster:
                mmheat_hei = 0.93 - heat_hei # this is also the maximal value (heamap edge is against the bottom)

                left_side_tree =         [0.01,  mmheat_hei,   0.186,      heat_hei]
                top_side_tree =          [0.198,   0.932,        heat_wid,   0.044]
                heatmap_location =       [0.198,   mmheat_hei,   heat_wid,   heat_hei]

                if col_colbar: # Slice a little out of the tree
                    top_side_tree =          [0.198,   0.946,  heat_wid,   0.040]
                    loc_col_colbar =         [0.198,   mmheat_hei+heat_hei+0.002,   heat_wid,  0.012]

                if row_colbar: # Slice a little out of the tree
                    left_side_tree =         [0.01,  mmheat_hei,   0.186-0.018, heat_hei]
                    loc_row_colbar =         [0.198-0.016,   mmheat_hei,   0.014,  heat_hei]

            else:
                # If no row cluster take advantage of the extra width available, but shift down to accomodate the scalebar
                mmheat_hei = 0.89 - heat_hei # this is also the maximal value (heamap edge is against the bottom)
                #top_side_tree =          [0.03,   0.852,  heat_wid,   0.044]
                top_side_tree =          [0.03,   0.891,  heat_wid,   0.020]
                heatmap_location =       [0.03,   mmheat_hei,   heat_wid,  heat_hei]
                loc_row_colbar =         [0.03-0.016,   mmheat_hei,   0.014,  heat_hei] # No need to cut the tree, just squeeze i into the left edge

                if col_colbar:
                    top_side_tree =          [0.03,   0.906,  heat_wid,   0.025] # squeeze up the colbar
                    loc_col_colbar =         [0.03,   0.892,   heat_wid,  0.012] #

        scalebar_location = [0.01,  0.98,   0.14,   0.015]

        # set size of the row text depending upon the number of items:
        row_font_size = 0
        if "row_font_size" in kargs:
            row_font_size = kargs["row_font_size"]
        elif "yticklabel_fontsize" in kargs:
            row_font_size = kargs["yticklabel_fontsize"]
        else:
            if "row_names" in kargs and kargs["row_names"]:
                if len(kargs["row_names"]) <= 100:
                    row_font_size = 6
                elif len(kargs["row_names"]) <= 150:
                    row_font_size = 4
                elif len(kargs["row_names"]) <= 200:
                    row_font_size = 3
                elif len(kargs["row_names"]) <= 300:
                    row_font_size = 2
                else:
                    if not highlights: # if highlights, don't kill kargs["row_names"] and don't print warning.
                        config.log.warning("heatmap has too many row labels to be visible. Suppressing row_labels")
                        kargs["row_names"] = None
                        row_font_size = 1
            else:
                row_font_size = 0

        if highlights:
            highlights = list(set(highlights)) # Make it unique, because sometimes duplicates get through and it erroneously reports failure.
            found = [False for i in highlights]
            if row_font_size == 0:
                row_font_size = 5 # IF the above sets to zero, reset to a reasonable value.
            if "row_font_size" in kargs: # but override ir row_font_size is being used.
                row_font_size = kargs["row_font_size"] # override size if highlights == True
                # I suppose this means you could do row_font_size = 0 if you wanted.

            # blank out anything not in row_names:
            new_row_names = []
            for item in kargs["row_names"]:
                if item in highlights:
                    new_row_names.append(item)
                    index = highlights.index(item)

                    found[index] = True
                else:
                    new_row_names.append("")
            kargs["row_names"] = new_row_names

            for i, e in enumerate(found): # check all highlights found:
                if not e:
                    config.log.warning("highlight: '%s' not found" % highlights[i])

        col_font_size = 6
        if "col_font_size" in kargs:
            col_font_size = kargs["col_font_size"]
        elif "xticklabel_fontsize" in kargs:
            col_font_size = kargs["xticklabel_fontsize"]

        if "bracket" in kargs: # done here so clustering is performed on bracketed data
            data = self.bracket_data(data, kargs["bracket"][0], kargs["bracket"][1])
            vmin = kargs["bracket"][0]
            vmax = kargs["bracket"][1]

        if not vmax:
            """
            I must guess the vmax value. I will do this by working out the
            mean then determining a symmetric colour distribution
            """
            me = mean(data)
            ma = abs(me - max(data))
            mi = abs(min(data) + me)
            if ma > mi:
                vmin = me - ma
                vmax = me + ma
            else:
                vmin = me - mi
                vmax = me + mi

        if "row_cluster" in kargs:
            row_cluster = kargs["row_cluster"]
        if "col_cluster" in kargs:
            col_cluster = kargs["col_cluster"]
        if not "colbar_label" in kargs:
            kargs["colbar_label"] = ""
        if "cmap" in kargs:
            colour_map = kargs["cmap"]
        if digitize:
            colour_map = cmaps.discretize(colour_map, digitize)

        # a few grace and sanity checks here;
        if len(data) <= 1: row_cluster = False # clustering with a single point?
        if len(data[0]) <= 1: col_cluster = False # ditto.

        if not "aspect" in kargs:
            kargs["aspect"] = "long"
        fig = self.getfigure(**kargs)

        row_order = False
        if row_cluster:
            # ---------------- Left side plot (tree) -------------------
            ax1 = fig.add_subplot(141)

            # from scipy;
            # generate the dendrogram
            if "row_tree" in kargs:
                assert "dendrogram" in kargs["row_tree"], "row_tree appears to be improperly formed ('dendrogram' is missing)"
                Z = kargs["row_tree"]["linkage"]
            else:
                Y = pdist(data, metric=cluster_mode)
                Z = linkage(Y, method='complete', metric=cluster_mode, optimal_ordering=optimal_ordering)

            if row_color_threshold:
                row_color_threshold = row_color_threshold*((Y.max()-Y.min())+Y.min()) # Convert to local threshold.
                a = dendrogram(Z, orientation='left', color_threshold=row_color_threshold, ax=ax1)
                ax1.axvline(row_color_threshold, color="grey", ls=":")
            else:
                a = dendrogram(Z, orientation='left', ax=ax1)

            ax1.set_position(left_side_tree)
            ax1.set_frame_on(False)
            ax1.set_xticklabels("")
            ax1.set_yticklabels("")
            ax1.set_ylabel("")
            # clear the ticks.
            ax1.tick_params(top=False, bottom=False, left=False, right=False)

            # Use the tree to reorder the data.
            row_order = [int(v) for v in a['ivl']]
            # resort the data by order;
            if "row_names" in kargs and kargs["row_names"]: # make it possible to cluster without names
                newd = []
                new_row_names = []
                for index in row_order:
                    newd.append(data[index])
                    new_row_names.append(kargs["row_names"][index])
                data = array(newd)
                kargs["row_names"] = new_row_names
            else: # no row_names, I still want to cluster
                newd = []
                for index in row_order:
                    newd.append(data[index])
                data = array(newd)

            if row_colbar:
                row_colbar = [row_colbar[index] for index in row_order]

        col_order = False
        if col_cluster:
            # ---------------- top side plot (tree) --------------------
            transposed_data = data.T

            ax2 = fig.add_subplot(142)
            ax2.set_frame_on(False)
            ax2.set_position(top_side_tree)
            if "col_tree" in kargs and kargs["col_tree"]:
                assert "dendrogram" in kargs["col_tree"], "col_tree appears to be improperly formed ('dendrogram' is missing)"
                #if kargs["col_names"] and kargs["col_names"]:
                #    assert len(kargs["col_tree"]["Z"]) == len(kargs["col_names"]), "tree is not the same size as the column labels"
                Z = kargs["col_tree"]["linkage"]
            else:
                Y = pdist(transposed_data, metric=cluster_mode)
                Z = linkage(Y, method='complete', metric=cluster_mode, optimal_ordering=optimal_ordering)
            a = dendrogram(Z, orientation='top', ax=ax2)

            ax2.tick_params(top=False, bottom=False, left=False, right=False)

            ax2.set_xticklabels("")
            ax2.set_yticklabels("")

            col_order = [int(v) for v in a["ivl"]]
            # resort the data by order;
            if col_names: # make it possible to cluster without names
                newd = []
                new_col_names = []
                for index in col_order:
                    newd.append(transposed_data[index])
                    new_col_names.append(col_names[index])
                data = array(newd).T # transpose back orientation
                col_names = new_col_names

            if col_colbar:
                col_colbar = [col_colbar[index] for index in col_order]

        # ---------------- Second plot (heatmap) -----------------------
        ax3 = fig.add_subplot(143)
        if 'imshow' in kargs and kargs['imshow']:
            ax3.set_position(heatmap_location) # must be done early for imshow
            hm = ax3.imshow(data, cmap=colour_map, vmin=vmin, vmax=vmax, aspect="auto",
                origin='lower', extent=[0, data.shape[1], 0, data.shape[0]],
                interpolation=config.get_interpolation_mode(filename)) # Yes, it really is nearest. Otherwise it will go to something like bilinear

        else:
            edgecolors = 'none'
            if grid:
                edgecolors = 'black'
            hm = ax3.pcolormesh(data, cmap=colour_map, vmin=vmin, vmax=vmax, antialiased=False, edgecolors=edgecolors, lw=0.4)

        if col_colbar:
            # Must be reordered by the col_cluster if present, done above;
            newd = {}
            colors_ = dict(list(matplotlib_colors.cnames.items()))
            for c in colors_:
                newd[c] = matplotlib_colors.hex2color(colors_[c])

            new_colbar = []
            for c in col_colbar:
                if '#' in c:
                    new_colbar.append([utils.hex_to_rgb(c)]) # needs to be tupled?
                else: # must be a named color:
                    new_colbar.append([newd[c]])

            col_colbar = numpy.array(new_colbar)#.transpose(1,0,2)

            ax4 = fig.add_axes(loc_col_colbar)
            if 'imshow' in kargs and kargs['imshow']:
                col_colbar = numpy.array(new_colbar).transpose(1,0,2)
                ax4.imshow(col_colbar, aspect="auto",
                    origin='lower', extent=[0, len(col_colbar),  0, 1],
                    interpolation=config.get_interpolation_mode(filename))

            else:
                col_colbar = numpy.array(new_colbar)
                # unpack the oddly contained data:
                col_colbar = [tuple(i[0]) for i in col_colbar]
                cols = list(set(col_colbar))
                lcmap = ListedColormap(cols)
                col_colbar_as_col_indeces = [cols.index(i) for i in col_colbar]

                ax4.pcolormesh(numpy.array([col_colbar_as_col_indeces,]), cmap=lcmap,
                    vmin=min(col_colbar_as_col_indeces), vmax=max(col_colbar_as_col_indeces),
                    antialiased=False, edgecolors=edgecolors, lw=0.4)

            ax4.set_frame_on(False)
            ax4.tick_params(top=False, bottom=False, left=False, right=False)
            ax4.set_xticklabels("")
            ax4.set_yticklabels("")

        if row_colbar:
            # Must be reordered by the row_cluster if present, done above;
            newd = {}
            colors_ = dict(list(matplotlib_colors.cnames.items()))
            for c in colors_:
                newd[c] = matplotlib_colors.hex2color(colors_[c])

            new_colbar = []
            for c in row_colbar:
                if '#' in c:
                    new_colbar.append([utils.hex_to_rgb(c)]) # needs to be tupled?
                else: # must be a named color:
                    new_colbar.append([newd[c]])

            row_colbar = numpy.array(new_colbar)

            ax4 = fig.add_axes(loc_row_colbar)
            if 'imshow' in kargs and kargs['imshow']:
                ax4.imshow(row_colbar, aspect="auto",
                    origin='lower', extent=[0, len(row_colbar),  0, 1],
                    interpolation=config.get_interpolation_mode(filename))
            else:
                # unpack the oddly contained data:
                row_colbar = [tuple(i[0]) for i in row_colbar]
                cols = list(set(row_colbar))
                lcmap = ListedColormap(cols)
                row_colbar_as_col_indeces = [cols.index(i) for i in row_colbar]
                ax4.pcolormesh(numpy.array([row_colbar_as_col_indeces,]).T, cmap=lcmap,
                    vmin=min(row_colbar_as_col_indeces), vmax=max(row_colbar_as_col_indeces),
                    antialiased=False, edgecolors=edgecolors, lw=0.4)

            ax4.set_frame_on(False)
            ax4.tick_params(top=False, bottom=False, left=False, right=False)
            ax4.set_xticklabels("")
            ax4.set_yticklabels("")

        if draw_numbers:
            for x in range(data.shape[0]):
                for y in range(data.shape[1]):
                    if data[x, y] >= draw_numbers_threshold:
                        if '%' in draw_numbers_fmt:
                            ax3.text(y+0.5, x+0.5, draw_numbers_fmt, size=draw_numbers_font_size,
                                ha='center', va='center')
                        else:
                            ax3.text(y+0.5, x+0.5, draw_numbers_fmt.format(data[x, y]), size=draw_numbers_font_size,
                                ha='center', va='center')

        if _draw_supplied_cell_labels:
            assert len(_draw_supplied_cell_labels) == data.shape[0], '_draw_supplied_cell_labels X does not equal shape[1]'
            assert len(_draw_supplied_cell_labels[0]) == data.shape[1], '_draw_supplied_cell_labels Y does not equal shape[0]'

            # This hack will go wrong if col_cluster or row_cluster is True
            # So fix the ordering:
            x_order = row_order if row_order else range(data.shape[0])
            y_order = col_order if col_order else range(data.shape[1])

            print(x_order, y_order)

            for xp, x in zip(range(data.shape[0]), x_order):
                for yp, y in zip(range(data.shape[1]), y_order):
                    val = _draw_supplied_cell_labels[x][y]
                    if draw_numbers_threshold and val < draw_numbers_threshold:
                        ax3.text(yp+0.5, xp+0.5, draw_numbers_fmt.format(val), size=draw_numbers_font_size,
                            ha='center', va='center')


        ax3.set_frame_on(border)
        ax3.set_position(heatmap_location)
        if col_names:
            ax3.set_xticks(arange(len(col_names))+0.5)
            ax3.set_xticklabels(col_names, rotation="vertical")
            ax3.set_xlim([0, len(col_names)])
            if "square" in kargs and kargs["square"]:
                ax3.set_xticklabels(col_names, rotation="vertical")
        else:
            ax3.set_xlim([0,data.shape[1]])

        if "row_names" in kargs and kargs["row_names"]:
            ax3.set_yticks(arange(len(kargs["row_names"]))+0.5)
            ax3.set_ylim([0, len(kargs["row_names"])])
            ax3.set_yticklabels(kargs["row_names"])
        else:
            ax3.set_ylim([0,data.shape[0]])
            ax3.set_yticklabels("")

        ax3.yaxis.tick_right()
        ax3.tick_params(top=False, bottom=False, left=False, right=False)
        [t.set_fontsize(row_font_size) for t in ax3.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(col_font_size) for t in ax3.get_xticklabels()]

        # Make it possible to blank with x/yticklabels
        if "xticklabels" in kargs:
            ax3.set_xticklabels(kargs["xticklabels"])
        if "yticklabels" in kargs:
            ax3.set_yticklabels(kargs["yticklabels"])

        ax0 = fig.add_subplot(144)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)

        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0)
        cb.set_label(kargs["colbar_label"], fontsize=6)
        cb.ax.tick_params(labelsize=4)

        return{
            "real_filename": self.savefigure(fig, filename, dpi=dpi),
            "reordered_cols": col_names,
            "reordered_rows": kargs["row_names"],
            "reordered_data": data
            }

    def heatmap2(self, filename=None, cluster_mode="euclidean", row_cluster=True, col_cluster=True,
        vmin=0, vmax=None, colour_map=cm.plasma, col_norm=False, row_norm=False, heat_wid=0.25,
        imshow=False,
        **kargs):
        """
        **Purpose**
            This version of heatmap is a simplified heatmap. It does not accept colnames, row_names
            and it outputs the heatmap better centred and expanded to fill the available space
            it does not draw a tree and (unlike normal heatmap) draws a black border
            around. Also, the scale-bar is optional, and by default is switched off.

            It is ideal for drawing sequence tag pileup heatmaps. For that is what it was originally made for
            (see track.heatmap())

        **Arguments**
            data (Required)
                the data to use. Should be a 2D array for the heatmap.

            filename (Required)
                The filename to save the heatmap to.

            col_norm (Optional, default=False)
                normalise each column of data between 0 .. max => 0.0 .. 1.0

            row_norm (Optional, default=False)
                similar to the defauly output of heatmap.2 in R, rows are normalised 0 .. 1

            colbar_label (Optional, default="expression")
                the label to place beneath the colour scale bar

            colour_map (Optional, default=afmhot)
                a matplotlib cmap for colour

            imshow (Optional, default=False)
                optional ability to use images for the heatmap. Currently experimental it is
                not always supported in the vector output files.

        **Returns**
            The actual filename used to save the image.
        """
        assert filename, "heatmap() - no specified filename"

        data = array(kargs["data"], dtype=float32) # heatmap2 can only accept a numpy array

        if col_norm:
            for col in range(data.shape[1]):
                data[:,col] /= float(data[:,col].max())

        if row_norm:
            for row in range(data.shape[0]):
                mi = min(data[row,:])
                ma = max(data[row,:])
                data[row,:] = (data[row,:]-mi) / (ma-mi)

        # positions of the items in the plot:
        heatmap_location =  [0.05,   0.01,   0.90,   0.90]
        scalebar_location = [0.05,  0.97,   0.90,   0.02]

        if "bracket" in kargs: # done here so clustering is performed on bracketed data
            data = self.bracket_data(data, kargs["bracket"][0], kargs["bracket"][1])
            vmin = kargs["bracket"][0]
            vmax = kargs["bracket"][1]
        else:
            vmin = data.min()
            vmax = data.max()

        if not "colbar_label" in kargs:
            kargs["colbar_label"] = "density"

        if "cmap" in kargs: colour_map = kargs["cmap"]

        # a few grace and sanity checks here;
        if len(data) <= 1: row_cluster = False # clustering with a single point?
        if len(data[0]) <= 1: col_cluster = False # ditto.

        if "size" not in kargs:
            kargs["size"] = (3,6)
        fig = self.getfigure(**kargs)

        # ---------------- (heatmap) -----------------------
        ax3 = fig.add_subplot(111)

        if imshow:
            ax3.set_position(heatmap_location) # must be done early for imshow
            hm = ax3.imshow(data, cmap=colour_map, vmin=vmin, vmax=vmax, aspect="auto",
                origin='lower', extent=[0, data.shape[1], 0, data.shape[0]],
                interpolation=config.get_interpolation_mode(filename))
        else:
            hm = ax3.pcolormesh(data, cmap=colour_map, vmin=vmin, vmax=vmax, antialiased=False)

        #ax3.set_frame_on(True)
        ax3.set_position(heatmap_location)

        ax3.set_xlim([0,data.shape[1]])

        ax3.set_ylim([0,data.shape[0]])
        ax3.set_yticklabels("")

        ax3.yaxis.tick_right()
        ax3.tick_params(top=False, bottom=False, left=False, right=False)
        [t.set_fontsize(1) for t in ax3.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(1) for t in ax3.get_xticklabels()]

        ax0 = fig.add_subplot(144)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)

        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0, cmap=colour_map)
        cb.set_label(kargs["colbar_label"])
        [label.set_fontsize(5) for label in ax0.get_xticklabels()]

        return(self.savefigure(fig, filename))

    def _heatmap_and_plot(self, peak_data=None, match_key=None,
        arraydata=None, peakdata=None, bin=None, draw_frames=False,
        imshow=False, **kargs):
        """
        Required:

        filename

            the filename to save the png to.

        array_data

            the data, usually serialisedArrayDataDict

        peak_data

            locations of the TF binding sites.

        match_key

            the key to match between the array and the peaklist.

        draw_frames (Optional, default=False)
            draw a frame around each of the elements in the figure.

        imshow (Optional, default=False)


        Optional:

        window

            the size of the moving average window (defaults to 10% of the list)
            If you specify a number it is the number of elements in the list
            and not a percentage.

        use_tag_score (defaults to False)

            use the key set by "use_tag_score" to determine the plot intesity for
            the peakdata. By default if the array and the peaklist
            match then it simply uses the
        """
        # ----------------------defaults:
        vmin = 0.0
        vmax = 1.0

        # ----------------------modify defaults:
        if "window" in kargs: moving_window = kargs["window"]
        if "match_key" in kargs: match_key = kargs["match_key"]
        if "bracket" in kargs:
            vmin = kargs["bracket"][0]
            vmax = kargs["bracket"][1]
        if "cmap" in kargs:
            cmap = kargs["cmap"]
        else:
            cmap = cm.RdBu_r

        # Positions of the items in the figure:
        left_heatmap = [0.10,  0.05,  0.20,  0.85]
        scale_bar =    [0.10,  0.97,  0.30,  0.02]
        binding_map =  [0.32,  0.05,  0.08,  0.85]
        freq_plot =    [0.42,   0.05,  0.4,   0.85]

        # Now do the plots:
        #plot.subplot(111)
        #plot.cla()

        fig = self.getfigure(**kargs)

        # heatmap ------------------------------------------------------

        ax0 = fig.add_subplot(141) # colour bar goes in here.
        ax0.set_frame_on(draw_frames)
        ax0.set_position(scale_bar)
        ax0.tick_params(left=False, right=False)

        ax1 = fig.add_subplot(142)
        plot_data = arraydata.T
        if imshow:
            ax1.set_position(left_heatmap)
            #hm = ax1.imshow(plot_data, cmap=cmap, vmin=vmin, vmax=vmax,
            #    interpolation=config.get_interpolation_mode(kargs["filename"]))

            hm = ax1.imshow(
                plot_data,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                aspect="auto",
                origin='lower',
                extent=[0, plot_data.shape[1], 0, plot_data.shape[0]],
                interpolation=config.get_interpolation_mode(kargs["filename"])
                )

        else:
            hm = ax1.pcolormesh(plot_data, cmap=cmap, vmin=vmin, vmax=vmax, antialiased=False)

        ax1.set_frame_on(draw_frames)
        ax1.set_position(left_heatmap)
        if "col_names" in kargs and kargs["col_names"]:
            ax1.set_xticks(arange(len(kargs["col_names"]))+0.5)
            ax1.set_xticklabels(kargs["col_names"])
            ax1.set_xlim([0, len(kargs["col_names"])])
        else:
            ax1.set_xlim([0,plot_data.shape[1]])
            ax1.set_xticklabels("")

        if "row_names" in kargs and kargs["row_names"] and len(kargs["row_names"]) <= 200: # you can't meaningfully see >200 labels. So suppress them:
            ax1.set_yticks(arange(len(kargs["row_names"]))+0.5)
            ax1.set_ylim([0, len(kargs["row_names"])])
            ax1.set_yticklabels(kargs["row_names"])
        else:
            ax1.set_ylim([0,plot_data.shape[0]])
            ax1.set_yticklabels("")

        ax1.yaxis.tick_left()
        ax1.tick_params(top=False, bottom=False, left=False, right=False)
        [t.set_fontsize(6) for t in ax1.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(6) for t in ax1.get_xticklabels()]

        fig.colorbar(hm, cax=ax0, orientation="horizontal")
        for label in ax0.get_xticklabels():
            label.set_fontsize(6)

        # binding map --------------------------------------------------

        ax2 = fig.add_subplot(143)

        a = array(bin) # reshape the bin array
        a.shape = 1,len(bin)

        if imshow:
            ax2.set_position(binding_map)
            #hm = ax1.imshow(plot_data, cmap=cmap, vmin=vmin, vmax=vmax,
            #    interpolation=config.get_interpolation_mode(kargs["filename"]))
            hm = ax2.imshow(
                a.T,
                cmap=cm.binary,
                aspect="auto",
                origin='lower',
                extent=[0, a.T.shape[1], 0, a.T.shape[0]],
                interpolation=config.get_interpolation_mode(kargs["filename"])
                )

        else:
            hm = ax2.pcolormesh(a.T, cmap=cm.binary, antialiased=True)

        ax2.set_frame_on(draw_frames)
        ax2.set_position(binding_map)
        ax2.set_yticks(arange(len(kargs["row_names"]))+0.5)
        ax2.set_yticklabels("")
        ax2.set_xticklabels("")
        ax2.set_xlim([0,1])
        ax2.set_ylim([0,len(kargs["row_names"])])
        ax2.yaxis.tick_left()
        ax2.tick_params(top=False, bottom=False, left=False, right=False)

        # linegraph -----------------------------------------------------

        ax3 = fig.add_subplot(144)
        ax3.plot(peakdata, arange(len(peakdata))) # doesn't use the movingAverage generated x, scale it across the entire graph.
        ax3.set_frame_on(draw_frames)
        ax3.set_position(freq_plot)
        ax3.set_yticklabels("")
        ax3.set_ylim([0, len(peakdata)])
        ax3.set_xlim([min(peakdata), (max(peakdata))+(max(peakdata)/10.0)])
        ax3.tick_params(left=False, right=False)
        [item.set_markeredgewidth(0.2) for item in ax3.xaxis.get_ticklines()]
        [t.set_fontsize(6) for t in ax3.get_xticklabels()]

        m = utils.mean(peakdata)
        s = utils.std(peakdata)

        ax3.axvline(x=m, color='black', linestyle=":", linewidth=1)
        ax3.axvline(x=(m+s), color='r', linestyle=":", linewidth=0.5)
        ax3.axvline(x=(m-s), color='r', linestyle=":", linewidth=0.5)

        return self.savefigure(fig, kargs["filename"], dpi=600)

    def multi_heatmap(self,
        list_of_data=None,
        filename=None,
        groups=None,
        titles=None,
        vmin=0, vmax=None,
        colour_map=cm.YlOrRd,
        col_norm=False,
        row_norm=False,
        heat_wid=0.25,
        frames=True,
        imshow=False,
        size=None,
        dpi:int = 80,
        **kargs):
        """
        **Purpose**
            Draw a multi-heatmap figure, i.e. containing multiple heatmaps. And also supports a
            last column indicating the different groups the blocks belong to.

        **Arguments**
            list_of_data (Required)
                Should be a list of 2D arrays ALL OF THE SAME SIZE! and all must be numpy arrays.

            groups (Optional)
                A list indicating which group each row belongs to.

            filename (Required)
                The filename to save the heatmap to.

            colbar_label (Optional, default="expression")
                the label to place beneath the colour scale bar

            size (Optional, default=None)
                override the guessed figure size with your own dimensions.

        **Returns**
            The actual filename used to save the image.
        """
        assert filename, "heatmap() - no specified filename"

        # work out a suitable size for the figure.
        num_heatmaps = len(list_of_data)

        if size:
            fig = self.getfigure(size=size, figsize=size)
        else: # guess:
            fig = self.getfigure(size=(3*num_heatmaps, 10))

        pad = 1.0 / (num_heatmaps+1)
        # item positions:
        heatmap_locations = [(0.01+(i*(pad+0.005)), 0.01, pad, 0.90) for i in range(num_heatmaps)]
        scalebar_location = [0.05,  0.97,   0.90,   0.02]

        if not "colbar_label" in kargs:
            kargs["colbar_label"] = "density"

        if "cmap" in kargs:
            colour_map = kargs["cmap"]

        # ---------------- (heatmap) -----------------------
        for index, data in enumerate(list_of_data):
            ax = fig.add_subplot(1, len(heatmap_locations)+1, index+1)

            if "bracket" in kargs: # done here so clustering is performed on bracketed data
                list_of_data[index] = self.bracket_data(data, kargs["bracket"][0], kargs["bracket"][1])
                vmin = kargs["bracket"][0]
                vmax = kargs["bracket"][1]
            else:
                vmin = list_of_data[index].min()
                vmax = list_of_data[index].max()

            if imshow:
                hm = ax.imshow(list_of_data[index], cmap=colour_map,
                    vmin=vmin, vmax=vmax, aspect="auto",
                    origin='lower',
                    extent=[0, list_of_data[index].shape[1], 0, list_of_data[index].shape[0]],
                    interpolation=config.get_interpolation_mode(filename))
            else:
                hm = ax.pcolormesh(list_of_data[index], cmap=colour_map, vmin=vmin, vmax=vmax, antialiased=False)

            ax.set_frame_on(frames)
            ax.set_position(heatmap_locations[index])
            if titles:
                ax.set_title(titles[index], fontdict={'fontsize': 6})

            ax.set_xlim([0,list_of_data[index].shape[1]])
            ax.set_ylim([0,list_of_data[index].shape[0]])

            #ax.yaxis.tick_right()
            ax.tick_params(top=False, bottom=False, left=False, right=False)
            [t.set_visible(False) for t in ax.get_yticklabels()] # generally has to go last.
            [t.set_visible(False) for t in ax.get_xticklabels()]

        ax0 = fig.add_subplot(1, len(heatmap_locations)+2, len(heatmap_locations)+1)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)
        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0, cmap=colour_map)

        # group membership
        if groups:
            ax = fig.add_subplot(1, len(heatmap_locations)+2, len(heatmap_locations)+2)
            right_most = list(heatmap_locations[-1])
            right_most[2] = 0.10
            right_most[0] += pad+0.02
            ax.set_position(right_most)

            dd = np.vstack((np.array(groups), np.array(groups))).T
            if imshow:
                hm = ax.imshow(dd, cmap=cm.Paired, vmin=min(groups), vmax=max(groups), aspect="auto",
                    origin='lower', extent=[0, dd.shape[1], 0, dd.shape[0]],
                    interpolation=config.get_interpolation_mode(filename))
            else:
                ax.pcolormesh(dd, vmin=min(groups), vmax=max(groups), antialiased=False, cmap=cm.Paired)

            ax.set_frame_on(frames)

            ax.set_xlim([0,dd.shape[1]])
            ax.set_ylim([0,dd.shape[0]])

            ax.tick_params(top=False, bottom=False, left=False, right=False)
            [t.set_visible(False) for t in ax.get_yticklabels()] # generally has to go last.
            [t.set_visible(False) for t in ax.get_xticklabels()]

            # add group labels:
            # get the upper and lower y axis coords for each group:
            last_item = None
            res = {}
            hist = {}
            for index, item in enumerate(groups):
                if item not in res:
                    res[item] = [index, None]
                if last_item and last_item != item:
                    res[last_item][1] = index-1
                last_item = item
                if item not in hist:
                    hist[item] = 0
                hist[item] += 1

            res[item][1] = len(groups) # fill in last item

            for r in res:
                ax.text(1, (res[r][0] + res[r][1]) / 2, "%s (%s)" % (r, hist[r]), ha="center", va="center")

        cb.set_label(kargs["colbar_label"])
        [label.set_fontsize(5) for label in ax0.get_xticklabels()]

        return self.savefigure(fig, filename, dpi=dpi)

    def grid_heatmap(self,
        data_dict_grid:dict = None,
        filename:str = None,
        row_labels=None,
        col_labels=None,
        titles=None,
        vmin:int = 0,
        vmax=None,
        colour_map=cm.YlOrRd,
        col_norm=False,
        row_norm=False,
        heat_wid=0.25,
        frames=False,
        imshow=False,
        size=None,
        dpi:int = 80,
        **kargs):
        """
        **Purpose**
            Draw a grid-based-heatmap figure, i.e. containing multiple heatmaps, with row and column labels

        **Arguments**
            data_dict_grid (Required)
                Should be adict in the form (for a 2x3 grid):
                    {
                        0: {0: numpy.array, 1: numpy.array},
                        1: {0: numpy.array, 1: numpy.array},
                        2: {0: numpy.array, 1: numpy.array}
                    }

            row_labels (Required)
                Row labels for the heatmaps;

            col_labels (Required)
                Col labels for the heatmaps;

            filename (Required)
                The filename to save the heatmap to.

            size (Optional, default=None)
                override the guessed figure size with your own dimensions.

        **Returns**
            The actual filename used to save the image.
        """
        assert filename, "heatmap() - no specified filename"

        num_cols = len(data_dict_grid)
        num_rows = len(data_dict_grid[0])

        hei_rats = []
        for pindex, _ in enumerate(data_dict_grid[0]): # must be all the same size...
            hei_rats.append(data_dict_grid[0][pindex].shape[0])

        # work out a suitable size for the figure.
        if size:
            fig = self.getfigure(size=size)
        else: # guess:
            fig = self.getfigure(size=(0.1+(2*num_cols), 1+num_rows)) # simpler jsut to use the number of rows

        gs = gridspec.GridSpec(num_rows, num_cols, height_ratios=hei_rats,
            top=0.95, bottom=0.08, left=0.05, right=0.98,
            hspace=0.03, wspace=0.05)

        gs2 = gridspec.GridSpec(1, num_cols,
            top=0.07, bottom=0.05, left=0.05, right=0.98,
            hspace=0.05, wspace=0.05)

        if not "colbar_label" in kargs:
            kargs["colbar_label"] = "density"

        if "cmap" in kargs:
            colour_map = kargs["cmap"]

        # ---------------- (heatmap) -----------------------
        for col in range(num_cols):
            for row in range(num_rows):
                data = data_dict_grid[col][row] # Yeah... That weird way around;
                ax = fig.add_subplot(gs[row, col])

                if "brackets" in kargs and kargs['brackets']: # done here so clustering is performed on bracketed data
                    bracket = kargs['brackets'][col]
                    data = self.bracket_data(data, bracket[0], bracket[1])
                    vmin = bracket[0]
                    vmax = bracket[1]
                elif "bracket" in kargs and kargs['bracket']:
                    bracket = kargs['bracket']
                    data = self.bracket_data(data, bracket[0], bracket[1])
                    vmin = bracket[0]
                    vmax = bracket[1]
                else:
                    vmin = data.min()
                    vmax = data.max()

                if imshow:
                    hm = ax.imshow(data,
                        cmap=colour_map,
                        vmin=vmin, vmax=vmax,
                        aspect="auto",
                        origin='lower',
                        extent=[0, data.shape[1], 0, data.shape[0]],
                        interpolation=config.get_interpolation_mode(filename))
                else:
                    hm = ax.pcolormesh(data, cmap=colour_map, vmin=vmin, vmax=vmax, antialiased=False)

                if col_labels and row == 0:
                    ax.set_title(col_labels[col], fontsize=6)

                if row_labels and col == 0:
                    ax.set_ylabel(row_labels[row], fontsize=6)

                ax.set_frame_on(frames)

                ax.set_xlim([0,data.shape[1]])
                ax.set_ylim([0,data.shape[0]])

                ax.tick_params(top=False, bottom=False, left=False, right=False)
                [t.set_visible(False) for t in ax.get_yticklabels()] # generally has to go last.
                [t.set_visible(False) for t in ax.get_xticklabels()]


            ax = fig.add_subplot(gs2[0, col])

            #ax0 = fig.add_subplot()
            #ax0.set_position(scalebar_location)
            #ax0.set_frame_on(False)
            cb = fig.colorbar(hm, cax=ax, orientation="horizontal") # cmap is from hm now;
            cb.set_label(kargs["colbar_label"], fontsize=6)
            [label.set_fontsize(6) for label in ax.get_xticklabels()]

        return self.savefigure(fig, filename, dpi=dpi)

    def boxplot(self,
        data=None,
        filename=None,
        labels=None,
        showfliers=True,
        whis=1.5,
        showmeans=False,
        meanline=False,
        tight_layout=False,
        grid=True,
        facecolors=None,
        **kargs):
        """
        wrapper around matplotlib's boxplot
        """
        assert filename, "no filename specified"
        assert labels, "boxplots must have labels"

        fig = self.getfigure(**kargs)

        if showmeans:
            meanline = True

        ax = fig.add_subplot(111)
        #ax.axhline(0, ls=":", color="grey") # add a grey line at zero for better orientation
        if grid:
            ax.grid(axis="y", ls=":", color="grey", lw=0.5, zorder=1000000)

        r = ax.boxplot(data, showfliers=showfliers, whis=whis, widths=0.5,
            patch_artist=True,
            showmeans=showmeans, meanline=meanline)

        plot.setp(r['medians'], color='green') # set nicer colours
        plot.setp(r['whiskers'], color='grey', lw=0.5)
        plot.setp(r['boxes'], color='black', lw=0.5)
        if facecolors:
            for patch, color in zip(r['boxes'], facecolors):
                patch.set_facecolor(color)
        else:
            for patch in r['boxes']:
                patch.set_facecolor('lightgrey')

        plot.setp(r['caps'], color='grey', lw=0.5)
        plot.setp(r['fliers'], color="grey", lw=0.5)

        #print(r.keys())

        ax.set_xticklabels(labels)

        if tight_layout:
            fig.tight_layout()

        fig.autofmt_xdate() # autorotate labels

        ax.tick_params(axis='x', labelsize=6)
        ax.tick_params(axis='y', labelsize=6)

        self.do_common_args(ax, **kargs)

        return self.savefigure(fig, filename)

    def cute_boxplot(self,
        filename:str,
        data,
        qs=None,
        title:str = None,
        xlims=None,
        sizer:float = 0.022,
        vert_height:int = 4,
        colors = None,
        bot_pad = 0.1,
        vlines = [0], # override the do_comon_args
        hlines = None, # override the do_comon_args
        showmeans = False,
        showfliers = False,
        **kargs):

        assert filename, 'You must specify a filename'
        if isinstance(colors, list): assert len(colors) == len(data.keys()), 'colors is not the same length as the data'
        if isinstance(qs, list): assert len(colors) == len(data.keys()), 'qs is not the same length as the data'

        # Because it draws bottom to top, almost always the user will prefer if the
        # data is reversed.

        keys = list(data.keys())
        keys.reverse()
        dats = list(data.values())
        dats.reverse()
        if qs: qs.reverse()
        if colors: colors.reverse()

        mmheat_hei = 0.1+(sizer*len(data))

        if 'figsize' not in kargs:
            kargs['figsize'] = [2.8,vert_height]

        fig = self.getfigure(**kargs)

        fig.subplots_adjust(left=0.4, right=0.8, top=mmheat_hei, bottom=bot_pad)
        ax = fig.add_subplot(111)

        ax.tick_params(right=True)

        m = 0
        if vlines:
            [ax.axvline(vli, ls=":", lw=0.5, color="grey") for vli in vlines]
        if hlines:
            [ax.axhline(vli, ls=":", lw=0.5, color="grey") for vli in hlines]

        r = ax.boxplot(dats,
            showfliers=showfliers,
            whis=True,
            patch_artist=True,
            widths=0.5,
            vert=False,
            showmeans=showmeans)

        plot.setp(r['medians'], color='black', lw=2) # set nicer colours
        plot.setp(r['boxes'], color='black', lw=0.5)
        plot.setp(r['caps'], color="grey", lw=0.5)
        plot.setp(r['whiskers'], color="grey", lw=0.5)

        ax.set_yticks(numpy.arange(len(data.keys()))+1)
        ax.set_yticklabels(keys)

        xlim = ax.get_xlim()[1]
        if xlims:
            ax.set_xlim(xlims)
            xlim = xlims[1]

        if qs:
            for i, q, p in zip(range(0, len(keys)+1), qs, r['boxes']):
                ax.text(xlim+(xlim/8), i+1, '{:.1e}'.format(q), ha='left', va='center', fontsize=6,)

        if isinstance(colors, list):
            [p.set_facecolor(k) for k, p in zip(colors, r['boxes'])]

        if title:
            ax.set_title(title, fontsize=6)

        [t.set_fontsize(6) for t in ax.get_yticklabels()]
        [t.set_fontsize(6) for t in ax.get_xticklabels()]

        self.do_common_args(ax, **kargs)
        return self.savefigure(fig, filename)

    def _scatter(self, x=None, y=None, filename=None, **kargs):
        """
        super thin wrapper aroung matplotlib's scatter
        """
        assert len(x), "x data missing"
        assert len(y), "y data missing"
        assert filename, "filename missing"

        fig = self.getfigure(**kargs)
        axis = fig.add_subplot(111)
        axis.scatter(x,y)

        if "logx" in kargs and kargs["logx"]: axis.set_xscale("log")
        if "logy" in kargs and kargs["logy"]: axis.set_yscale("log")

        if "title" in kargs: axis.set_title(kargs["title"])
        if "xlabel" in kargs: axis.set_xlabel(kargs["xlabel"])
        if "ylabel" in kargs: axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs: axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs: axis.set_yticks(kargs["yaxis"])

        return(self.savefigure(fig, filename))

    def _qhist(self, filename=None, data=None, bins=60, **kargs):
        """
        Very thin wrapper around matplotlibs' hist
        """
        assert filename, "Internal Error: _qhist no filename"

        log = False
        if "log" in kargs:
            log = kargs["log"]

        fig = self.getfigure(**kargs)
        axis = fig.add_subplot(111)
        axis.hist(data, bins=bins, facecolor='orange', ec="none", alpha=1.0, log=log)

        if "title" in kargs:
            axis.set_title(kargs["title"])
        if "xlabel" in kargs:
            axis.set_xlabel(kargs["xlabel"])
        if "ylabel" in kargs:
            axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs:
            axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs:
            axis.set_yticks(kargs["yaxis"])

        return self.savefigure(fig, filename)

    def _plot(self, filename=None, data=None, **kargs):
        """
        Internal very very thin wrapper around matplotlib's plot
        """

        fig = self.getfigure(**kargs)
        axis = fig.add_subplot(111)

        if "x" in kargs:
            axis.plot(kargs["x"], data)
        else:
            axis.plot(data)

        if "title" in kargs: axis.set_title(kargs["title"])
        if "xlabel" in kargs: axis.set_xlabel(kargs["xlabel"])
        if "ylabel" in kargs: axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs: axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs: axis.set_yticks(kargs["yaxis"])
        if "xlims" in kargs: axis.set_xlim(kargs["xlims"])
        if "ylims" in kargs: axis.set_ylim(kargs["ylims"])

        return self.savefigure(fig, filename)

    def _genome_segment(self, figure_axis, loc, feature_list):
        """
        draw a representation of the genome using the axis provided by figure_axis
        loc = the genomic location
        feature_list = a list of refseq features to draw on the graph.
        """
        ax = figure_axis
        # set the x axis to match.
        ax.set_xlim([0, len(loc)]) # 1 = 1bp scale.
        ax.set_ylim([0, 10]) # arbitary.

        ax.set_xticks([0, len(loc)])
        ax.set_yticks([0, 10])

        ax.tick_params(top=False, bottom=False, left=False, right=False)
        left_base = loc["left"]
        for item in feature_list:
            if item["type"] == "gene":
                left_most_base = item["loc"]["left"] - loc["left"]
                right_most_base = item["loc"]["right"] - loc["left"]

                if item["strand"] in positive_strand_labels:
                    ax.text(left_most_base, 7, item["name"], size=10, ha="left", va="center")
                    ax.arrow(left_most_base, 6, right_most_base-left_most_base, 0,
                        alpha=0.5, fc=(0,0,0), width=0.02)
                elif item["strand"] in negative_strand_labels:
                    ax.text(right_most_base, 3, item["name"], size=10, ha="right", va="center")
                    ax.arrow(left_most_base, 4, right_most_base-left_most_base, 0,
                        alpha=0.5, fc=(0,0,0), width=0.02)
        ax.axhline(y=5, color='gray', linestyle=":", linewidth=0.5, alpha=0.5)
        # tidy up stuff:
        ax.set_frame_on(False)
        ax.set_yticklabels("")
        ax.set_xticklabels("")
        # append the chromosome coords to the figure.
        ax.text(100, 8, str(loc), size=10, ha="left", va="center")

    def _labeled_figure(self, data=None, axissize=None, filename=None,
        horizontal_line=True, **kargs):
        """
        **Purpose**
            Draw a figure with a set of labels.

        **Arguments**
            filename
                filename to save to. May get modified depending upon the current
                draw mode.

            data
                A set of values if this form: {"pos": (x,y), "label": label}

            horizontal_line (Optional, default=True)
                draw  horizontal line at y axis 0.

            axissize (Required)
                the axis dimensions (x and y maximal limits).

            Common Arguments also supported:
                figsize - tuple specifying the figure aspect
                dpi - the dpi (only supported for ps outputs)
        """
        position_plot = [0.02, 0.05, 0.96, 0.9]

        fig = self.getfigure(**kargs)

        ax1 = fig.add_subplot(131)
        ax1.set_position(position_plot)
        ax1.set_xlim([0, axissize[0]])
        ax1.set_ylim([0, axissize[1]])
        if "ylim" and kargs["ylim"]:
            ax1.set_ylim(kargs["ylim"])
        if horizontal_line:
            ax1.axhline(y=4.5, color='gray', linestyle=":", linewidth=1.5)

        for l in data:
            # draw an arrow.
            ax1.text(l["pos"][0], l["pos"][1], l["label"], size=8, ha="left", va="bottom",
                rotation="vertical")

        if "genomic_features" in kargs:
            # I've got some genomic features, I want to draw them.
            position_genomic = [0.02, 0.05, 0.96, 0.1]
            ax3 = fig.add_subplot(133)
            ax3.set_position(position_genomic)
            self._genome_segment(ax3, kargs["loc"], kargs["genomic_features"])

        return(self.savefigure(fig, filename))

    def _plot_and_histogram(self, filename=None, data=None, figsize=(5,5), **kargs):
        """
        Draw a graph plot and a histogram on the right hand side.
        """
        position_plot = [0.02, 0.05, 0.83, 0.9]
        position_histogram = [0.87, 0.05, 0.12, 0.9]


        fig = self.getfigure(**kargs)
        ax1 = fig.add_subplot(131)

        ax1.set_position(position_plot)
        if "x" in kargs: # fake some x data.
            ax1.plot(kargs["x"], data)
        else:
            ax1.plot(data)

        # histogram plot.
        ax2 = fig.add_subplot(132)
        ax2.set_position(position_histogram)
        n, bins, patches = ax2.hist(data, bins=20, orientation='horizontal', histtype="stepfilled", color=(0,0,0))

        m = mean(data) # hehe, numpy pawns my homemade routine for nontrivial samples.
        s = std(data)

        y = mlab.normpdf( bins, m, s)
        l = ax2.plot(bins, y, 'r--', linewidth=1)
        # add mean, std lines to the first plot.
        # add mean line
        ax1.axvline(x=m, color='black', linestyle="-", linewidth=0.5)
        # add std lines
        for z in [2,3,4]:
            zz = s*z # num stds away
            ax1.axhline(y=(m+zz), color='gray', linestyle=":", linewidth=1.5)
            ax1.axhline(y=(m-zz), color='gray', linestyle=":", linewidth=1)

        self.do_common_args(ax1, **kargs)
        if "xlims" in kargs: ax1.set_xlim(kargs["xlims"])
        if "ylims" in kargs: ax1.set_ylim(kargs["ylims"])

        #if "xlims" in kargs: ax2.set_xlim(kargs["xlims"])
        if "ylims" in kargs: ax2.set_ylim(kargs["ylims"])

        if "genomic_features" in kargs:
            # I've got some genomic features, I want to draw them.
            # I have to generate my own axis for draw._genome_segment
            position_genomic = [0.02, 0.05, 0.83, 0.1]
            ax3 = fig.add_subplot(133)
            ax3.set_position(position_genomic)
            self._genome_segment(ax3, kargs["loc"], kargs["genomic_features"])

        return(self.savefigure(fig, filename))

    # Salted for deprecation
    def _qplotxy(self, list_of_tuples_data, filename=None, labels=None, **kargs):
        self.qplotxy(list_of_tuples_data, filename=filename, labels=labels, **kargs)

    def qplotxy(self, list_of_tuples_data, filename=None, labels=None, **kargs):
        """
        thin wrapper around plot.
        expects list_of_tuples_data to be a list of tuples containing X and Y data:

        [ ([0..x], [0..y]), ([],[]) .. ([],[]) ] ???
        # valid kargs:
        labels = labels for each line (a list or iterable)
        dpi = dpi. Uses config.DEFAULT_DPI if not specified
        """
        assert filename, "Internal Error: _qplot no filename"
        assert labels, "Internal Error: _qplot, labels must be a list"

        # set up figure
        fig = self.getfigure(**kargs)
        axis = fig.add_subplot(111)

        if "title" in kargs: axis.set_title(kargs["title"])
        if "xlabel" in kargs: axis.set_xlabel(kargs["xlabel"])
        if "ylabel" in kargs: axis.set_ylabel(kargs["ylabel"])
        if "xaxis" in kargs: axis.set_xticks(kargs["xaxis"])
        if "yaxis" in kargs: axis.set_yticks(kargs["yaxis"])
        if "xlim" in kargs: axis.set_xlim(kargs["xlim"])
        if "ylim" in kargs: axis.set_ylim(kargs["ylim"])

        for index, item in enumerate(list_of_tuples_data):
            axis.plot(item[0], item[1], label=labels[index])

        #self.draw.py.yaxis([0,max(item[1])])
        if labels:
            leg = axis.legend()#bbox_to_anchor=(0., 1.02, 1., .102), loc=3)
            leg.get_frame().set_alpha(0.5)
            [t.set_fontsize(6) for t in leg.get_texts()]

        return(self.savefigure(fig, filename))

    def _vennDiagram2(self, left, right, overlap, filename=None,
        proportional=False, plot_sig=False, labels=None,
        **kargs):
        """
        draw a venn Diagram.

        **Arguments**

            (Required Arguments)

            left (relative)
                left - overlap. _vennDiagram2() does not correct for
                overlap.

            right (relative)
                right - overlap. _vennDiagram2() does not correct for
                overlap.

            overlap
                number of sites in common between left and right

            filename
                save name

            labels
                A dictionary in this form:
                {"left": "", "right": "", "title": ""}

            (Optional Arguments)

            proportional (True|False, default False)
                Use proportional circles suggesting the sizes of the overlaps

            plot_sig (True|False, default False)
                plot a binomial p-value test result for the overlap.

            simulations (100 by default)
                number of simulations to run

            world_size (defaults to sum of left, right, overlap)
                If your world_size is bigger you need to send it as the
                simulations will overestimate the overlaps.

            (Supported generic kargs)

            dpi
                dpi of the output file.

        **Result**
            saves a venn Diagram and returns the actual path to the saved file.

        **todo**
            support multiple overlaps.
            support for weighted circles.

        """
        # no valid_args specified here. You should know what you are doing...
        assert filename, "No filename specified!"
        assert labels, "No labels!"

        fig = self.getfigure(**kargs)
        axis = fig.add_subplot(111)
        axis.set_position([0.02, 0.02, 0.96, 0.96])

        artists = []
        if not proportional:
            artists.append(Circle((7, 10), 5.5, alpha=0.5, facecolor="#FFA200")) # (loc), size
            artists.append(Circle((13, 10), 5.5, alpha=0.5, facecolor="#104BA9"))

            artists.append(Circle((7, 10), 5.5, alpha=1.0, fill=False, lw=2)) # I draw the circle twice to give bold outer lines.
            artists.append(Circle((13, 10), 5.5, alpha=1.0, fill=False, lw=2))

            axis.text(5, 10, str(left), size=30, ha="center", va="center")
            axis.text(10, 10, str(overlap), size=25, ha="center", va="center")
            axis.text(15, 10, str(right), size=30, ha="center", va="center")
        else:
            max_size = 5.5
            max_attraction = max_size * 2
            min_size = 1.0
            total_score = float(left + right + overlap)
            left_w =  ((left + overlap) / total_score) * max_size
            left_w2 =  ((left) / total_score) * max_attraction
            right_w = ((right + overlap) / total_score) * max_size
            right_w2 = ((right)/ total_score) * max_attraction

            centre_w = (overlap / total_score)

            delta_dist = centre_w

            # sanity checking:
            #if delta_dist > max_attraction: delta_dist = max_attraction
            #if left_w < min_size: left_w = min_size
            #if right_w < min_size: right_w = min_size

            #print delta_dist

            artists.append(Circle((4.5+left_w2, 10), left_w, alpha=0.5, facecolor="#FFA200")) # (loc), size
            artists.append(Circle((15.5-right_w2, 10), right_w, alpha=0.5, facecolor="#104BA9"))

            #artists.append(Circle((6, 10), left_weight, alpha=1.0, fill=False, lw=2)) # I draw the circle twice to give bold outer lines.
            #artists.append(Circle((14, 10), right_weight, alpha=1.0, fill=False, lw=2))

        for a in artists: # add all artists...
            axis.add_artist(a)

        axis.set_xticks([0,20])
        axis.set_yticks([0,20])
        # clear frame and axis markers:
        axis.set_frame_on(False)
        axis.set_yticklabels("")
        axis.set_xticklabels("")
        axis.tick_params(top=False, bottom=False, left=False, right=False)

        # add the labels:
        axis.text(5, 17, labels["left"], size=15, ha="center", va="center")
        axis.text(15, 17, labels["right"], size=15, ha="center", va="center")
        axis.text(10, 19, labels["title"], size=28, ha="center", va="center")

        return self.savefigure(fig, filename)

    def venn2(self, A, B, AB, labelA, labelB, filename, **kargs):
        """
        same as _vennDiagram2 except it performs a triple overlap.
        """
        if not "aspect" in kargs:
            kargs["aspect"] = "square"
        if not "size" in kargs:
            kargs["size"] = "small"

        fig = self.getfigure(**kargs)
        ax = fig.add_subplot(111)
        ax.set_position([0.02, 0.02, 0.96, 0.96])
        ax.set_xlim([0,30])
        ax.set_ylim([0,30])

        artists = []
        artists.append(Circle((10, 15), 8, alpha=1, lw=1, edgecolor='grey', facecolor="none"))
        artists.append(Circle((20, 15), 8, alpha=1, lw=1, edgecolor='grey', facecolor="none"))

        for a in artists: # add all artists...
            ax.add_artist(a)

        #ax.set_xticks([0,20])
        #ax.set_yticks([0,20])
        # clear frame and axis markers:
        ax.set_frame_on(False)
        ax.set_yticklabels("")
        ax.set_xticklabels("")

        ax.tick_params(top=False, bottom=False, left=False, right=False)

        # add the labels:
        ax.text(7.5, 15, A-AB, size=6, ha="center", va="center")
        ax.text(22.5, 15, B-AB, size=6, ha="center", va="center")

        ax.text(15, 15, AB, size=6, ha="center", va="center")

        ax.text(7.5,  25, labelA, size=6, ha="center", va="center")
        ax.text(22.5,  25, labelB, size=6, ha="center", va="center")

        self.do_common_args(ax, **kargs)
        return self.savefigure(fig, filename)

    def venn3(self, A, B, C, AB, AC, BC, ABC, labelA, labelB, labelC, filename, **kargs):
        """
        same as _vennDiagram2 except it performs a triple overlap.
        """
        if not "aspect" in kargs:
            kargs["aspect"] = "square"
        if not "size" in kargs:
            kargs["size"] = "small"

        fig = self.getfigure(**kargs)
        ax = fig.add_subplot(111)
        ax.set_position([0.02, 0.02, 0.96, 0.96])
        ax.set_xlim([0,30])
        ax.set_ylim([0,30])

        artists = []
        artists.append(Circle((10, 20), 8, alpha=1, lw=1, edgecolor='grey', facecolor="none"))
        artists.append(Circle((20, 20), 8, alpha=1, lw=1, edgecolor='grey', facecolor="none"))
        artists.append(Circle((15, 10), 8, alpha=1, lw=1, edgecolor='grey', facecolor="none"))

        for a in artists: # add all artists...
            ax.add_artist(a)

        #ax.set_xticks([0,20])
        #ax.set_yticks([0,20])
        # clear frame and axis markers:
        ax.set_frame_on(False)
        ax.set_yticklabels("")
        ax.set_xticklabels("")

        ax.tick_params(top=False, bottom=False, left=False, right=False)

        # add the labels:
        ax.text(7.5, 20, A-AB-AC+ABC, size=15, ha="center", va="center")
        ax.text(22.5, 20, B-AB-BC+ABC, size=15, ha="center", va="center")
        ax.text(15, 7.5, C-AC-BC+ABC, size=15, ha="center", va="center")

        ax.text(15, 21, AB-ABC, size=15, ha="center", va="center")
        ax.text(11, 14, AC-ABC, size=15, ha="center", va="center")
        ax.text(19, 14, BC-ABC, size=15, ha="center", va="center")

        ax.text(15, 16.5, ABC, size=15, ha="center", va="center")

        ax.text(7.5,  29, labelA, size=16, ha="center", va="center")
        ax.text(22.5,  29, labelB, size=16, ha="center", va="center")
        ax.text(15, 1, labelC, size=16, ha="center", va="center")

        return self.savefigure(fig, filename)

    def venn4(self, lists, labels, filename=None, **kargs):
        """
        Draw a 4-way venn Diagram.

        lists should be in the order:

        A, B, C, D, AB, AC, AD, BC, BD, CD, ABC, ABD, ACD, BCD, ABCD

        labels should be in the order:

        A, B, C, D
        """
        # no valid_args specified here. You should know what you are doing...
        assert filename, "No filename specified!"

        # Verbosity for clarity
        A = lists[0]
        B = lists[1]
        C = lists[2]
        D = lists[3]
        AB = lists[4]
        AC = lists[5]
        AD = lists[6]
        BC = lists[7]
        BD = lists[8]
        CD = lists[9]
        ABC = lists[10]
        ABD = lists[11]
        ACD = lists[12]
        BCD = lists[13]
        ABCD = lists[14]
        labelA = labels[0]
        labelB = labels[1]
        labelC = labels[2]
        labelD = labels[3]

        fig = self.getfigure(**kargs)
        ax = fig.add_subplot(111)
        ax.set_position([0.02, 0.02, 0.96, 0.96])
        ax.set_xlim([0,10])
        ax.set_ylim([0,10])

        artists = []
        artists.append(Ellipse((7,7), width=3, height=7, angle=45, edgecolor='black', facecolor="none"))
        artists.append(Ellipse((6.5,6), width=3, height=7, angle=45, edgecolor="black", facecolor="none"))
        artists.append(Ellipse((6.5,4), width=3, height=7, angle=135, edgecolor="black", facecolor="none"))
        artists.append(Ellipse((7,3), width=3, height=7, angle=135, edgecolor="black", facecolor="none"))

        for a in artists: # add all artists...
            ax.add_artist(a)

        ax.set_xticks([0,10])
        ax.set_yticks([0,10])
        # clear frame and axis markers:
        ax.set_frame_on(False)
        ax.set_yticklabels("")
        ax.set_xticklabels("")
        ax.tick_params(top=False, bottom=False, left=False, right=False)

        # New version using sets to do everything
        # add the labels:
        ax.text(7, 8.4, A, size=12, ha="center", va="center")
        ax.text(5, 6.3, B, size=12, ha="center", va="center")
        ax.text(5, 3.7, C, size=12, ha="center", va="center")
        ax.text(7, 1.6, D, size=12, ha="center", va="center")

        ax.text(6, 7,     AB, size=12, ha="center", va="center")
        ax.text(8.7, 6.2, AC, size=12, ha="center", va="center")
        ax.text(9.4, 5,   AD, size=12, ha="center", va="center")
        ax.text(6.1, 5,   BC, size=12, ha="center", va="center")
        ax.text(8.7, 3.8, BD, size=12, ha="center", va="center")
        ax.text(6, 3,     CD, size=12, ha="center", va="center")

        ax.text(7.2, 5.9, ABC, size=12, ha="center", va="center")
        ax.text(9.0, 4.6, ABD, size=12, ha="center", va="center")
        ax.text(9.0, 5.4, ACD, size=12, ha="center", va="center")
        ax.text(7.2, 4.1, BCD, size=12, ha="center", va="center")

        ax.text(8, 5, ABCD, size=12, ha="center", va="center")

        # labels
        ax.text(4.2,  9, labelA, size=13, ha="right", va="center")
        ax.text(3.7,  8, labelB, size=13, ha="right", va="center")
        ax.text(3.7,  2, labelC, size=13, ha="right", va="center")
        ax.text(4.2,  1, labelD, size=13, ha="right", va="center")

        return self.savefigure(fig, filename)

    def venn5(self, lists, labels, filename=None, **kargs):
        """
        Draw a 5-way venn Diagram.

        lists should be in the order:

        A, B, C, D, AB, AC, AD, BC, BD, CD, ABC, ABD, ACD, BCD, ABCD

        labels should be in the order:

        A, B, C, D
        """
        # no valid_args specified here. You should know what you are doing...
        assert filename, "No filename specified!"

        # Verbosity for clarity
        A = lists[0]
        B = lists[1]
        C = lists[2]
        D = lists[3]
        E = lists[4]
        AB = lists[5]
        AC = lists[6]
        AD = lists[7]
        AE = lists[8]
        BC = lists[9]
        BD = lists[10]
        BE = lists[11]
        CD = lists[12]
        CE = lists[13]
        DE = lists[14]
        ABC = lists[15]
        ABD = lists[16]
        ABE = lists[17]
        ACD = lists[18]
        ACE = lists[19]
        ADE = lists[20]
        BCD = lists[21]
        BCE = lists[22]
        BDE = lists[23]
        CDE = lists[24]
        ABCD = lists[25]
        ABCE = lists[26]
        ACDE = lists[27]
        BCDE = lists[28]
        ABCDE = lists[29]

        fig = self.getfigure(**kargs)
        ax = fig.add_subplot(111)
        #ax.set_position([0.02, 0.02, 0.96, 0.96])
        ax.set_xlim([0,10])
        ax.set_ylim([0,10])

        artists = []
        artists.append(Ellipse((5,       5-0.6), width=3, height=8, angle=0, facecolor="none"))
        artists.append(Ellipse((5+0.75,  5-0.3), width=3, height=8, angle=72, facecolor="none"))
        artists.append(Ellipse((5+0.45,  5+0.9), width=3, height=8, angle=144, facecolor="none"))
        artists.append(Ellipse((5-0.45,  5+0.9), width=3, height=8, angle=216, facecolor="none"))
        artists.append(Ellipse((5-0.75,  5-0.3), width=3, height=8, angle=288, facecolor="none"))

        for a in artists: # add all artists...
            ax.add_artist(a)

        #ax.set_xticks([0,10])
        #ax.set_yticks([0,10])
        # clear frame and axis markers:
        #ax.set_frame_on(False)
        #ax.set_yticklabels("")
        #ax.set_xticklabels("")
        #[i.set_markeredgewidth(0.0) for i in ax.yaxis.get_ticklines()]
        #[i.set_markeredgewidth(0.0) for i in ax.xaxis.get_ticklines()]

        # add the labels:
        ax.text(5, 5, 'ABCDE', size=12, ha="center", va="center")

        # labels
        ax.text(4.2, 9, labels[0], size=16, ha="right", va="center")
        ax.text(3.7, 8, labels[1], size=16, ha="right", va="center")
        ax.text(3.7, 2, labels[2], size=16, ha="right", va="center")
        ax.text(4.2, 1, labels[3], size=16, ha="right", va="center")

        return self.savefigure(fig, filename)

    def getfigure(self, size=None, aspect=None, **kargs):
        """
        **Purpose**
            setup a valid figure instance based on size.

        **Arguments**
            size or figsize (Optional, default="medium")
                if size is a tuple then that tuple is the specified size in inches (don't ask)
                You can also specify "small", "medium", "large" and "huge". Corrsponding to approximate pixel
                sizes of (with the "normal" aspect)
                small  :
                medium :
                large  :
                huge   :

                If size is specified then it takes preference over aspect and aspect will be ignored.

            aspect (Optional, default="normal")
                the aspect of the image.
                currently only "normal", "long" and "square" are respected).

        **Returns**
            A valid matplotlib figure object
        """
        if "figsize" in kargs and kargs["figsize"]:
            size = kargs["figsize"]

        if not size:
            size = config.draw_size
        elif len(size) == 2: # A tuple or list?
            size_in_in = (size[0], size[1])

        if not aspect:
            aspect = config.draw_aspect

        if len(size) == 2: # overrides aspect/size
            size_in_in = (size[0], size[1])
            return(plot.figure(figsize=size_in_in))

        data = {"normal": {"small": (5,4), "medium": (8,6), "large": (12,9), "huge": (16,12)},
                "square": {"small": (4,4), "medium": (7,7), "large": (9,9), "huge": (12,12)},
                "long": {"small": (4,5), "medium": (6,8), "large": (9,12), "huge": (12,16)},
                "wide": {"small": (7,4), "medium": (12,6), "large": (18,9), "huge": (24,12)}
                }
        dpi = {"small": 75, "medium": 150, "large": 300, "huge": 600} # This dpi doesn't actually work here...
        # See savefigure() for the actual specification
        return plot.figure(figsize=data[aspect][size])

    def savefigure(self, fig, filename, size=config.draw_size, bbox_inches=None, dpi=None):
        """
        **Purpose**
            Save the figure
            to filename, modifying the filename based on the current drawing mode
            (if required)
        **Arguments**
            fig
                the figure handle

            filename
                the filename to save the file to

        **Returns**
            the actual filename used to save the image
        """
        temp_draw_mode = config.draw_mode
        if isinstance(config.draw_mode, str):
            temp_draw_mode = [config.draw_mode] # for simple compat

        for mode in temp_draw_mode:
            assert mode in config.valid_draw_modes, "'%s' is not a supported drawing mode" % temp_draw_mode

            if mode == 'svg':
                matplotlib.rcParams["image.interpolation"] = 'nearest'
            # So that saving supports relative paths.
            path, head = os.path.split(filename)
            if "." in filename: # trust Ralf to send a filename without a . in it Now you get your own special exception!
                save_name = "%s.%s" % (".".join(head.split(".")[:-1]), mode) # this will delete .. in filename, e.g. file.meh.png
            else:
                save_name = "%s.%s" % (head, mode)

            fig.savefig(os.path.join(path, save_name), bbox_inches=bbox_inches, dpi=dpi)
            plot.close(fig) # Saves a huge amount of memory.
        return save_name

    def do_common_args(self, ax, **kargs):
        """
        **Purpose**
            deal with common arguments to matplotlib (may not always work, depending upon the figure type.

        **Arguments**
            ax
                an matplotlib axes object

            These are based loosly on the matplotlib versions
                xlabel - x-axis label
                ylabel - y-axis label
                title  - title
                xlims - x axis limits
                ylims - y-axis limits
                zlims - z-axis limits (For 3D plots only)
                xticklabels - list (or not) of labels for the x axis
                logx - set the x scale to a log scale argument should equal the base
                logy - set the y scale to a log scale
                legend_size - size of the legend, small, normal, medium
                xticklabel_fontsize - x tick labels fontsizes
                xticklabels - labels to draw on the x axis
                yticklabel_fontsize - y tick labels fontsizes
                yticklabels - labels to draw on the y axis
                vlines - A list of X points to draw a vertical line at
                hlines - A list of Y points to draw a horizontal line at
                ticks_top - True/False, display the axis ticks on the top
                ticks_bottom - True/False, display the axis ticks on the bottom
                ticks_left - True/False, display the axis ticks on the left
                ticks_right - True/False, display the axis ticks on the right
                ticks - True/False, display any ticks at all
                xticks - List of tick positions you want to draw
                yticks - List of tick positions you want to draw

        **Returns**
            None
        """
        legend = None
        try:
            legend = ax.get_legend()
        except AttributeError:
            pass
        if legend: # None in no legend on this plot
            legend.get_frame().set_alpha(0.5) # make the legend transparent

        if "xlabel" in kargs:
            ax.set_xlabel(kargs["xlabel"])
        if "ylabel" in kargs:
            ax.set_ylabel(kargs["ylabel"])
        if "title" in kargs:
            if "title_fontsize" in kargs:
                ax.set_title(kargs["title"], fontdict={'fontsize': kargs['title_fontsize']})
            else:
                ax.set_title(kargs["title"], fontdict={'fontsize': 6})
        if "xlims" in kargs:
            ax.set_xlim(kargs["xlims"])
        if "ylims" in kargs:
            ax.set_ylim(kargs["ylims"])
        if "zlims" in kargs: # For 3D plots
            ax.set_zlim([kargs["zlim"][0], kargs["zlim"][1]])
        if "logx" in kargs:
            ax.set_xscale("log", basex=kargs["logx"])
        if "logy" in kargs:
            ax.set_yscale("log", basey=kargs["logy"])
        if "log" in kargs and kargs["log"]:
            ax.set_xscale("log", basex=kargs["log"])
            ax.set_yscale("log", basey=kargs["log"])
        if "legend_size" in kargs:
            [t.set_fontsize(kargs["legend_size"]) for t in legend.get_texts()]
        if "xticklabel_fontsize" in kargs:
            ax.tick_params(axis='x', labelsize=kargs["xticklabel_fontsize"])
        if "yticklabel_fontsize" in kargs:
            ax.tick_params(axis='y', labelsize=kargs["yticklabel_fontsize"])
        if "xticklabels" in kargs:
            ax.set_xticklabels(kargs["xticklabels"])
        if "yticklabels" in kargs:
            ax.set_yticklabels(kargs["yticklabels"])
        if "vlines" in kargs and kargs["vlines"]:
            for l in kargs["vlines"]:
                ax.axvline(l, ls=":", color="grey")
        if "hlines" in kargs and kargs["hlines"]:
            for l in kargs["hlines"]:
                ax.axhline(l, ls=":", color="grey")
        if "alines" in kargs and kargs['alines']:
            for quple in kargs['alines']:
                # quples are interleaved, list of x's and list of y's
                ax.plot((quple[0], quple[2]), (quple[1], quple[3]), ls=':', color='grey', lw=0.5)
        if "grid" in kargs and kargs["grid"]:
            ax.grid()
        if "ticks" in kargs and not kargs["ticks"]:
            ax.tick_params(top="off", bottom="off", left="off", right="off")
        if "ticks_top" in kargs and not kargs["ticks_top"]:
            ax.tick_params(top="off")
        if "ticks_bottom" in kargs and not kargs["ticks_bottom"]:
            ax.tick_params(bottom="off")
        if "ticks_left" in kargs and not kargs["ticks_left"]:
            ax.tick_params(left="off")
        if "ticks_right" in kargs and not kargs["ticks_right"]:
            ax.tick_params(right="off")
        if 'xticks' in kargs and kargs["xticks"]:
            ax.set_xticks(kargs['xticks'])
        if 'yticks' in kargs and kargs["yticks"]:
            ax.set_yticks(kargs['yticks'])

    def _stacked_plots(self, filename, loc, features, graphs, merged=False, **kargs):
        """
        Draw a stacked set of line plots

        graph_data is 'stackable', the more graphs added the graphs
        will seperate into different graphs.

        See pwm.scan_seq_with_features() for details of usage.

        graphs should be a dictionary with of which the key will be used as
        a name
        """

        dpi = config.DEFAULT_DPI
        if "dpi" in kargs:
            dpi = kargs["dpi"]

        fig = self.getfigure(**kargs)

        spacer = 0.9 / len(graphs)

        for index, key in enumerate(graphs):
            ax2 = fig.add_subplot(1, len(graphs), index)
            ax2.set_position([0.02, spacer*index, 0.96, spacer])
            ax2.plot(graphs[key])

            #ax2.set_frame_on(False)
            ax2.annotate(key, xy=(0.02, spacer*index))
            ax2.set_xticklabels("")
            ax2.set_yticklabels("")
            ax2.set_ylabel("")
            ax2.tick_params(top=False, bottom=False, left=False, right=False)
            ax2.set_ylim([0, 1])
            ax2.set_xlim([0, len(graphs[key])])

        return self.savefigure(fig, filename)

    def _simple_heatmap(self, filename=None, colour_map = cm.Reds, vmin=0, vmax=None, symmetric=False, **kargs):
        """
        A simplified version of heatmap, with no clustering, and a simpler representation
        Also, you can change the size and aspect of the display.

        **Arguments**
            data
                an array of arrays or equivalent.

            colour_map
                Default is YlOrRd

            bracket
                specify a tuple for the min max values for the heatmap.

            symmetric (Optional, default=False)
                If set to true, find the mean, the max and the min n the data
                and then set the range of colours to span min .. mean .. max, so that
                the gap between mean and max and min and max are identical.
                If set to False (default behaviour) then simply range the colours from min(data) to
                max(data).

            fig_size (Optional, default=(6,6))
                change the figure size aspect.
        """
        # This should be a wrapper around draw.heatmap() to take advantage of heatmaps
        # better code.

        assert filename, "_heatmap() missing filename"

        data = kargs["data"]

        # positions of the items in the plot:
        heatmap_location =  [0.12,   0.01,   0.75,   0.98]
        scalebar_location = [0.01,  0.96,   0.10,   0.03]

        if "bracket" in kargs:
            data = self.bracket_data(data, kargs["bracket"][0], kargs["bracket"][1])
            vmin = kargs["bracket"][0]
            vmax = kargs["bracket"][1]

        if not vmax:
            """
            I must guess the vmax value. I will do this by working out the
            mean then determining a symmetric colour distribution
            """
            if symmetric:
                me = mean(data)
                ma = abs(me - max(data))
                mi = abs(min(data) + me)
                if ma > mi:
                    vmin = me - ma
                    vmax = me + ma
                else:
                    vmin = me - mi
                    vmax = me + mi
            else:
                vmax = max(data)
                vmin = min(data)

        if not "aspect" in kargs:
            kargs["aspect"] = "long"
        fig = self.getfigure(**kargs)

        # ---------------- Second plot (heatmap) -----------------------
        ax3 = fig.add_subplot(121)
        hm = ax3.pcolormesh(data, cmap=colour_map, vmin=vmin, vmax=vmax, antialiased=False)

        ax3.set_frame_on(True)
        ax3.set_position(heatmap_location)
        ax3.set_xlim([0,data.shape[1]])
        ax3.set_ylim([0,data.shape[0]])
        ax3.set_yticklabels("")
        ax3.set_xticklabels("")
        ax3.yaxis.tick_right()
        ax3.tick_params(top=False, bottom=False, left=False, right=False)
        #[t.set_fontsize(6) for t in ax3.get_yticklabels()] # generally has to go last.
        #[t.set_fontsize(6) for t in ax3.get_xticklabels()]

        ax0 = fig.add_subplot(122)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)

        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0, cmap=colour_map)
        cb.set_label("")
        [label.set_fontsize(5) for label in ax0.get_xticklabels()]

        return self.savefigure(fig, filename)

    def nice_scatter(self, x=None, y=None, filename=None, do_best_fit_line=False,
        print_correlation=False, spot_size=4, plot_diag_slope=False, label_fontsize=14,
        **kargs):
        """
        **Purpose**
            Draw a nice simple scatter plot

        **Arguments**
            x, y (Required)
                x and y values

            filename (Required)
                the filename to save as.

            spots (Optional, must be a 2-length tuple containing (x, y) data)
                These spots will be empahsised with whatever spots_cols is or an
                "Orange" colour by default

            spot_labels (Optional, labels to write on the spots)
                A list of labels to write over the spots.

            label_fontsize (Optional, default=14)
            	labels fontsize

            plot_diag_slope (Optional, default=False)
                Plot a diagonal line across the scatter plot.

            do_best_fit_line (Optional, default=False)
                Draw a line of best fit and the

            print_correlation (Optional, default=None)
                You have to spectify the type of correlation to print on the graph.
                valid are:
                    r = R (Correlation coefficient)
                    r2 = R^2.

            spot_size (Optional, default=5)
                The size of each dot.

            Supported keyword arguments:
                xlabel, ylabel, title, logx, logy

        **Returns**
            the real filename, which may get modified depending upon the current drawing mode
            (usually results in a png)
        """
        fig = self.getfigure(**kargs)
        ax = fig.add_subplot(111)

        ax.scatter(x, y, s=spot_size, c="grey", alpha=0.2, edgecolors="none")

        if "spots" in kargs and kargs["spots"]:
            if "spots_cols" in kargs and kargs["spots_cols"]:
                # Will recognise a string or sequence autmagivally.
                ax.scatter(kargs["spots"][0], kargs["spots"][1], s=spot_size*2, c=kargs["spots_cols"], alpha=0.8, edgecolor="none")
            else:
                ax.scatter(kargs["spots"][0], kargs["spots"][1], s=spot_size*2, c="orange", alpha=0.8, edgecolor="none")

            if "spot_labels" in kargs and kargs["spot_labels"]:
                # for the matplot.lib < 100: I want to label everything.
                for i, n in enumerate(kargs["spot_labels"]):
                    ax.annotate(n, (kargs["spots"][0][i], kargs["spots"][1][i]), size=label_fontsize, color="black", ha="center", va="center")

        if print_correlation or do_best_fit_line:
            # linear regression
            (ar, br) = polyfit(x, y, 1)
            xr = polyval([ar,br], x)
            slope, intercept, r_value, p_value, std_err = linregress(x,y)

            # I think this line is actually wrong?
            mx = [min(x), max(x)]
            my = [slope * min(x) + intercept, slope * max(x) + intercept]

            # Only draw if specified:
            if do_best_fit_line:
                ax.plot(mx, my, "r-.", lw=0.5)

            if print_correlation:
                if print_correlation == "r":
                    ax.set_title("R=%.4f" % r_value)
                elif print_correlation == "r2":
                    ax.set_title("R2=%.4f" % (r_value*r_value))
                elif print_correlation == "pearson":
                    ax.set_title("Pearson=%.4f" % scipy.stats.pearsonr(x,y)[0])
        if plot_diag_slope:
            ax.plot([min(x+y), max(x+y)], [min(x+y), max(x+y)], ":", color="grey")

        if "logx" in kargs and kargs["logx"]:
            ax.set_xscale("log", basex=kargs["logx"])
        if "logy" in kargs and kargs["logy"]:
            ax.set_yscale("log", basey=kargs["logy"])

        self.do_common_args(ax, **kargs)

        return self.savefigure(fig, filename)

    def bar_chart(self, filename=None, genelist=None, data=None, cols=None, **kargs):
        """
        **Purpose**
            draw a bar chart with error bars and interpret and package the data coming from a genelist-like object

        **Args**
            filename

            genelist
                a genelist-like object

            data
                the key to look for in the genelist for the data

            labels
                the key to look for in the genelist for labels

            title (Optional)
                the title

            err (Optional)
                the key to look for in the genelist for error bar values.
                This one assumes symmetric values +- around the data

            err_up (Optional)
                the key to look for errorbars going up

            err_dn (Optional)
                the key to look for error bars going down

            errs_are_absolute (Optional, default=False)
                error bars are not +- from the data, but are values that specify where the error
                bars extend to. This needs to be set to True commonly for confidence intervals
                and left as False for standard errors.

            cols (Optional, default=Use a default set from matplotlib)
                the colours to use for the bar charts, there should be one for each bar.

            Other kargs respected by bar_chart:
                aspect
                size
                xlabel - x-axis label
                ylabel - y-axis label
                title  - title
                xlims - x axis limits
                ylims - y-axis limits
                logx - set the x scale to a log scale argument should equal the base
                logy - set the y scale to a log scale

        **Returns**
            The actual_filename used to save the image
        """

        da = []
        err = []
        err_up = []
        err_dn = []
        for i in genelist:
            da.append(i[data])
            if "err" in kargs:
                err.append(i[kargs["errs"]])
            if "err_up" in kargs:
                err_up.append(i[kargs["err_up"]])
            if "err_dn" in kargs:
                err_dn.append(i[kargs["err_dn"]])

        if "errs_are_absolute" in kargs and kargs["errs_are_absolute"]:
            for i, n in enumerate(da): # normalise the values so matplotlib can underastand them
                if err_up:
                    err_up[i] = [a - b for a, b in zip(err_up[i], n)]
                if err_dn:
                    err_dn[i] = [b - a for a, b in zip(err_dn[i], n)]

        if "cond_names" in kargs:
            labs = kargs["cond_names"]
        else: # fake one
            labs = ["" for t in da]

        if not cols:
            # I need to generate a series of colours.
            cmap = cm.get_cmap(cm.Paired, len(labs))
            cols = []
            step = 256 // len(labs)
            for t in range(1, 256, step):
                cols.append(cmap(t))
            #print cols

        fig = self.getfigure(**kargs)
        ax = fig.add_subplot(111)
        ax.set_position([0.3, 0.1, 0.68, 0.8]) # plenty of space for labels

        # Convert to Numpy arrays for type laziness
        da = array(da).T
        err = array(err).T
        err_up = array(err_up).T
        err_dn = array(err_dn).T
        wid = (1.0 / len(da))-0.05
        x = arange(len(da[0]))

        general_args = {"ec": "black", "ecolor": "black"}

        for i, r in enumerate(da):
            if err:
                ax.barh(x+(wid*i), r, wid, xerr=err, label=labs[i], fc=cols[i], **general_args)
            elif "err_up" in kargs and "err_dn" in kargs:
                ax.barh(x+(wid*i), r, wid, xerr=(err_dn[i], err_up[i]), label=labs[i], fc=cols[i], **general_args)
            elif "err_up" in kargs:
                ax.barh(x+(wid*i), r, wid, xerr=err_up[i], label=labs[i], fc=cols[i], **general_args)
            else:
                ax.barh(x+(wid*i), r, wid, label=labs[i], fc=cols[i], **general_args)
            # I'm sure you don't mean just err_dn
        ax.set_ylim([0, x[-1]+(wid*2)])

        if "cond_names" in kargs:
            leg = ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
            leg.get_frame().set_alpha(0.5)

        if "labels" in kargs and kargs["labels"]:
            ax.set_yticklabels(genelist[kargs["labels"]], rotation="horizontal")
            ax.set_yticks(x)#+0.5)

        self.do_common_args(ax, **kargs)

        ax.tick_params(top=False, bottom=False, left=False, right=False)
        #[t.set_fontsize(6) for t in ax3.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(6) for t in ax.get_yticklabels()]

        return(self.savefigure(fig, filename))

    def pie(self, data, labels, filename, aspect="square", title=None, colours=None,
        draw_percents=False, cmap=None, font_size=12, **kargs):
        """
        Draw a PIE!
        """

        fig = self.getfigure(aspect=aspect, **kargs)
        ax = fig.add_subplot(111)

        extra_args = {}
        if draw_percents:
            extra_args["autopct"] = "%1.1f%%"

        if colours:
            texts = ax.pie(data, labels=labels, shadow=False, colors=colours, labeldistance=1.03, **extra_args)

        elif cmap:
            ld = len(data)
            cmap = cm.get_cmap(cmap, ld)
            ran = np.linspace(0.0, 1.0, num=ld)#/float(ld)
            colours = cmap(ran)[:,:-1]
            colours = [rgb2hex(col) for col in colours] # rgb2hex from matplotlib not needed, but easier to print

            texts = ax.pie(data, labels=labels, shadow=False, colors=colours, labeldistance=1.1, **extra_args)
        else:
            texts = ax.pie(data, labels=labels, shadow=False, labeldistance=1.03, **extra_args)

        plot.setp(texts[1], fontsize=font_size)

        if title:
            ax.set_title(title)

        return(self.savefigure(fig, filename))

    def beanplot(self, data, filename, violin=True, order=None, mean=False, median=True, alpha=0.2,
        beans=True, IQR=False, covariance=0.2, **kargs):
        """
        http://statsmodels.sourceforge.net/devel/_modules/statsmodels/graphics/boxplots.html#beanplot
        **Purpose**
            A beanplot/beeswarm/violin plot!

        **Arguments**
            data
                a dict of lists:
                {"condition A": [0, 1, ... n],
                "condition B": [0, 1, ... n]}
                These do not have to be the same length, if you are calling this function directly.

            filename (Required)
                filename to save the file to. modified based on current draw_mode

            violin (Optional, default=True)
                Draw a KDE violin around the points

            beans (Optional, default=True)
                draw the beans!

            mean (Optional, default=False)
                draw red lines for means on each plot

            median (Optional, default=True)
                draw a line at the median for each plot

            IQR (Optional, default=False)
                draw lines indicating the interquartile ranges at 75%

            order (Optional, default=None)
                list of the keys in the order you want them.
                Also doubles as labels.

            alpha (Optional, default=0.2)
                The alpha value to use for transparency of the spots (0...1)

            colors (Optional, default="blue")
                A color or list of colors to use to color the spots.

            Other standard supported args:
                aspect
                size
                xlabel - x-axis label
                ylabel - y-axis label
                title  - title
                xlims - x axis limits
                ylims - y-axis limits

            covariance (Optional, default=2)


        **Returns**
            None
        """
        fig = self.getfigure(**kargs)
        ax = fig.add_subplot(111)

        if not order:
            order = list(data.keys())

        num_cats = len(data)
        #xs = np.arange(num_cats)
        cmin = 0
        cmax = 0
        bins = 1000
        for x, d in enumerate(order):
            x_data = np.linspace(min(data[d]), max(data[d]), bins+2)
            # get the violin: required, even if not drawn.
            # Check that there is some variation. If no variation then utils.kde will break
            if numpy.around(numpy.std(data[d]), 3) > 0:
                bracket = (min(data[d]), max(data[d]))
                if 'ylims' in kargs and kargs['ylims']:
                    bracket = kargs['ylims']
                y_violin = utils.kde(data[d], range=bracket, bins=bins)
                y_violin = ((y_violin / max(y_violin)) * 0.4) # normalise
                y_violin = numpy.insert(y_violin, 0, 0)
                y_violin = numpy.append(y_violin, 0)

                if violin: # envelope
                    ax.plot(x+y_violin, x_data, color="grey", lw=0.75)
                    ax.plot(x-y_violin, x_data, color="grey", lw=0.75)

                jitter_envelope = np.interp(data[d], x_data, y_violin)
                jitter_coord = jitter_envelope * (np.random.uniform(low=0, high=1, size=len(data[d])) + np.random.uniform(low=-1, high=0, size=len(data[d])))
                if beans:
                    ax.scatter(jitter_coord+x, data[d], alpha=alpha, edgecolor="none", s=10)#, color="red")
            else:
                config.log.warning('bean/violinplot: No variation in sample "%s", not drawing' % d)

            if median:
                ax.plot([x-0.3, x+0.3], [numpy.median(data[d])]*2, c="green")
            if mean:
                ax.plot([x-0.3, x+0.3], [numpy.mean(data[d])]*2, c="red")
            if IQR:
                q75, q25 = numpy.percentile(data[d], [75 ,25])
                ax.plot([x-0.15, x+0.15], [q75]*2, c="green", lw=.5)
                ax.plot([x, x], [q75,q75-(q25*0.1)], c="green", lw=.5) # Yes q25
                ax.plot([x-0.15, x+0.15], [q25]*2, c="green", lw=.5)
                ax.plot([x, x], [q25,q25+(q25*0.1)], c="green", lw=.5)

            if max(data[d]) > cmax:
                cmax = max(data[d])
            if min(data[d]) < cmin:
                cmin = min(data[d])

        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order)

        ax.set_ylim([cmin, cmax])
        ax.set_xlim([-0.6, len(order)-0.4])
        ax.set_xticks([i for i in range(len(order))]) # xticks must be 1 separated to get all labels for line up

        fig.autofmt_xdate()

        ax.tick_params(axis='x', labelsize=6)
        ax.tick_params(axis='y', labelsize=6)

        self.do_common_args(ax, **kargs)
        return self.savefigure(fig, filename)

    def violinplot(self,
        data,
        filename:str,
        violin=True,
        order=None,
        mean=False,
        median=True,
        colors=None,
        **kargs):
        '''
        Uses the matplotlib implementation;

        '''
        fig = self.getfigure(**kargs)
        ax = fig.add_subplot(111)

        if not order:
            order = list(data.keys())

        pos = numpy.arange(len(order))

        r = ax.violinplot(
            [data[k] for k in order],
            pos,
            points=50,
            widths=0.5,
            showmeans=mean,
            showmedians=median,
            )

        for vidx, viol in enumerate(r['bodies']):
            if colors and isinstance(colors, str):
                viol.set_facecolor(colors)
            elif colors and isinstance(colors, Iterable):
                viol.set_facecolor(colors[vidx])

        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order)

        #ax.set_ylim([cmin, cmax])
        ax.set_xlim([-0.6, len(order)-0.4])
        ax.set_xticks([i for i in range(len(order))]) # xticks must be 1 separated to get all labels for line up

        fig.autofmt_xdate()

        ax.tick_params(axis='x', labelsize=6)
        ax.tick_params(axis='y', labelsize=6)

        self.do_common_args(ax, **kargs)

        return self.savefigure(fig, filename)

    def unified_scatter(self,
        labels,
        xdata,
        ydata,
        x,
        y,
        mode='PC',
        filename=None,
        spots=True,
        label=False,
        alpha=0.8,
        perc_weights=None,
        spot_cols='grey',
        overplot=None,
        spot_size=40,
        label_font_size=7,
        label_style=None,
        cut=None,
        squish_scales=False,
        only_plot_if_x_in_label=None,
        adjust_labels=False,
        cmap=None,
        cluster_data=None,
        draw_clusters=None,
        cluster_labels=None,
        cluster_centroids=None,
        **kargs):
        '''
        Unified for less bugs, more fun!
        '''
        ret_data = None

        if not "aspect" in kargs:
            kargs["aspect"] = "square"
        if 'figsize' not in kargs:
            kargs['figsize'] = (4,4)

        fig = self.getfigure(**kargs)
        ax = fig.add_subplot(111)

        if only_plot_if_x_in_label:
            newx = []
            newy = []
            newlab = []
            newcols = []
            for i, lab in enumerate(labels):
                if True in [l in lab for l in only_plot_if_x_in_label]:
                    if overplot and lab not in overplot: # Don't do twice
                        newx.append(xdata[i])
                        newy.append(ydata[i])
                        newlab.append(labels[i])
                        newcols.append(spot_cols[i])
            xdata = newx
            ydata = newy
            labels = newlab
            spot_cols = newcols
            #print zip(spot_cols, labels)

        if overplot: # Make sure some spots are on the top
            newx = []
            newy = []
            newcols = []
            for i, lab in enumerate(labels):
                if True in [l in lab for l in overplot]:
                    newx.append(xdata[i])
                    newy.append(ydata[i])
                    newcols.append(spot_cols[i])

        if draw_clusters:
            # unpack the cluster_data for convenience
            ax.set_prop_cycle(cycler(color=plot.get_cmap('tab20c').colors))
            n_clusters = cluster_data.n_clusters

            for labelk in range(n_clusters):
                cluster_center = cluster_centroids[labelk]
                this_x = [xdata[i] for i, l in enumerate(cluster_labels) if l == labelk]
                this_y = [ydata[i] for i, l in enumerate(cluster_labels) if l == labelk]
                ax.scatter(this_x, this_y, s=spot_size+1, alpha=1.0, edgecolors="none", zorder=5)

                #ax.plot(xdata[labelk], ydata[labelk], 'w', markerfacecolor=col, marker='.')

                ax.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor='none', markeredgecolor='black', alpha=0.4, markersize=6, zorder=6)
                ax.text(cluster_center[0], cluster_center[1], labelk, ha='center', zorder=7)

            ax.scatter(xdata, ydata, s=spot_size,
                        alpha=alpha, edgecolors="none",
                        c=spot_cols, cmap=cmap,
                        zorder=2)
        elif spots:
            ax.scatter(xdata, ydata,
                alpha=alpha, edgecolors="none",
                s=spot_size,
                c=spot_cols,
                cmap=cmap,
                zorder=2)
        else:
            # if spots is false then the axis limits are set to 0..1. I will have to send my
            # own semi-sensible limits:
            squish_scales = True

        if overplot:
            ax.scatter(newx, newy, s=spot_size+1, alpha=alpha, edgecolors="none", c=newcols, zorder=5)

        if label:
            texts = []
            for i, lab in enumerate(labels):
                if not spots and isinstance(spot_cols, list):
                    texts.append(ax.text(xdata[i], ydata[i], lab, size=label_font_size, color=spot_cols[i], style=label_style, ha='center'))
                else:
                    texts.append(ax.text(xdata[i], ydata[i], lab, size=label_font_size, style=label_style, color="black"))
            if adjust_labels:
                adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))

        # Tighten the axis
        if squish_scales:
            # do_common_args will override these, so don't worry
            dx = (max(xdata) - min(xdata)) * 0.05
            dy = (max(ydata) - min(ydata)) * 0.05
            ax.set_xlim([min(xdata)-dx, max(xdata)+dx])
            ax.set_ylim([min(ydata)-dy, max(ydata)+dy])

        if perc_weights is not None: # perc_weights is often a numpy array
            ax.set_xlabel("%s%s (%.1f%%)" % (mode, x, perc_weights[x-1])) # can be overridden via do_common_args()
            ax.set_ylabel("%s%s (%.1f%%)" % (mode, y, perc_weights[y-1]))
        else:
            ax.set_xlabel("%s%s" % (mode, x)) # can be overridden via do_common_args()
            ax.set_ylabel("%s%s" % (mode, y))

        if "logx" in kargs and kargs["logx"]:
            ax.set_xscale("log", basex=kargs["logx"])
        if "logy" in kargs and kargs["logy"]:
            ax.set_yscale("log", basey=kargs["logy"])

        if cut:
            rect = matplotlib.patches.Rectangle(cut[0:2], cut[2]-cut[0], cut[3]-cut[1], ec="none", alpha=0.2, fc="orange")
            ax.add_patch(rect)

            tdata = []
            for i in range(0, len(xdata)):
                if xdata[i] > cut[0] and xdata[i] < cut[2]:
                    if ydata[i] < cut[1] and ydata[i] > cut[3]:
                        if self.rowwise: # grab the full entry from the parent genelist
                            dat = {"pcx": xdata[i], "pcy": ydata[i]}
                            dat.update(self.parent.linearData[i])
                            tdata.append(dat)
                        else:
                            tdata.append({"name": lab[i], "pcx": xdata[i], "pcy": ydata[i]})
            if tdata:
                ret_data = genelist()
                ret_data.load_list(tdata)

        self.do_common_args(ax, **kargs)

        real_filename = self.savefigure(fig, filename)
        config.log.info("scatter: Saved '%s%s' vs '%s%s' scatter to '%s'" % (mode, x, mode, y, real_filename))
        return ret_data

    def dotbarplot(self, data, filename, yticktitle='Number', **kargs):
        """
        **Purpose**
            A drawing wrapper for the new style dot-mean-stderr plots.
            Each sample is a circle, the mean is a red line and the std err of the mean
            is shown where available.

            This version also supports heterogenous lengths of data.

        **Arguments**
            data (Required)
                A dictionary of label: [0, 1, 2, ... n] values.

                The x category labels will be taken from the dict key.

            order (Optional, default=data.keys())
                order for the conditions to be plotted

            filename (Required)
                filename to save the image to

            yticktitle (Optional default='Number')
                Added here just to emphasise you should really set your own y ticklabel.

        **Returns**
            The filename the data was saved to.

        """
        assert filename, 'You must specify a filename'
        assert isinstance(data, dict), 'data must be a dictionary'

        # ax 2:
        fig = self.getfigure(**kargs)
        fig.subplots_adjust(left=None, bottom=0.35, right=None, top=None, wspace=0.3, hspace=None)
        ax = fig.add_subplot(111)
        barh = list(data.values())
        labs = list(data.keys()) # Dicts are now ordered, so lab order should be fine.

        # Get the means/stdevs:
        means = [numpy.mean(data[k]) for k in data]
        stds = [numpy.std(data[k]) / math.sqrt(len(data[k])) for k in data]
        lengths = [len(data[k]) for k in data] # for working out if it is valid to plot errs or mean

        # convert the data array into a linear list of x and y:
        xs = numpy.arange(0.5, len(labs)+0.5)
        xd = []
        yd = []
        for i, k in enumerate(data):
            for d in data[k]:
                xd.append(xs[i])
                yd.append(d)


        ax.scatter(xd, yd, edgecolors='black', lw=0.5, c='none', s=20)
        if False not in [i>=3 for i in lengths]:
            ax.errorbar(xs, means, yerr=stds, barsabove=True, fmt='none', capsize=4, capthick=0.5, ls='-', color='black', lw=0.5)

        ax.set_ylim([0, max(yd)+10])
        ax.set_xlim([-0.2, len(labs)+0.2])
        ax.set_xticks(xs)
        ax.set_xticklabels(list(labs))
        ax.set_xticklabels(labs, rotation=45, rotation_mode="anchor", ha="right")
        ax.set_ylabel(yticktitle)

        if False not in [i>=2 for i in lengths]:
            ls = []
            for i in zip(xs, means):
                l = mlines.Line2D([i[0]-0.3, i[0]+0.3], [i[1], i[1]], c='red', lw=1)
                ax.add_line(l) #ax.scatter(xs, ms, color='red', marker='_')

        self.do_common_args(ax, **kargs)

        real_filename = self.savefigure(fig, filename)
        config.log.info("dotbarplot: Saved dotbarplot to '%s'" % (real_filename))
        return(real_filename)

    def proportional_bar(self,
        filename,
        data_dict,
        key_order=None,
        title='',
        cols=None,
        **kargs):
        '''
        **Purpose**
            Draw a bar plot, but with proporional bars.

        **Arguments**
            filename (Required)
                filename to save the figure to.

            data_dict (Required)
                {
                'row_name1': {'class1': 0, 'class2': 0},
                'row_name2': {'class1': 0, 'class2': 0},
                }

            key_order (Optional)
                order for the row_names;

            ...

        '''
        assert filename, 'A filename to save the image to is required'
        assert data_dict, 'data_dict is required'
        assert isinstance(data_dict, dict), 'data_dict is not a dict'

        if not cols:
            cols = plot.rcParams['axes.prop_cycle'].by_key()['color']

        # get all of the classes:
        if not key_order:
            all_keys = [] # preserve order
            for k in data_dict:
                for kk in data_dict[k]:
                    if kk not in all_keys:
                        all_keys.append(kk)
            print('Found {0} keys'.format(all_keys))
        else:
            all_keys = key_order

        vals = {k: [] for k in all_keys}

        labs = []
        for k in data_dict:
            labs.append(k)
            for kk in all_keys:
                vals[kk].append(float(data_dict[k][kk]))

        scaled = {k: [] for k in all_keys}
        sums = None
        for k in all_keys:
            if sums is None:
                sums = numpy.zeros(len(vals[k]))
            sums += vals[k]

        for k in all_keys:
            vals[k] = numpy.array(vals[k])
            scaled[k] = numpy.array(vals[k])
            scaled[k] /= sums
            scaled[k] *= 100

        plot_hei = (0.8) - (0.04*len(labs))

        plot.rcParams['pdf.fonttype'] = 42
        fig = plot.figure(figsize=[4,3])
        fig.subplots_adjust(left=0.35, right=0.95, bottom=plot_hei,)
        ax = fig.add_subplot(111)
        ax.set_prop_cycle('color', cols)

        ypos = numpy.arange(len(data_dict))

        # data_dict = {'bar_row': {'class': 0, class2': 0}}

        bots = numpy.zeros(len(labs))
        for k in vals:
            ax.barh(ypos, scaled[k], 0.7, label=k, left=bots)
            for y, v, s, b in zip(ypos, vals[k], scaled[k], bots):
                ax.text(b+(s//2), y, '{0:,.0f} ({1:.0f}%)'.format(v, s), ha='center', va='center', fontsize=6)
            bots += scaled[k]

        ax.set_yticks(ypos)
        ax.set_yticklabels(labs)

        ax.set_xlim([-2, 102])
        ax.set_xticks([0, 50, 100])
        ax.set_xticklabels(['0%', '50%', '100%'])
        ax.set_title(title, size=6)
        ax.legend()
        plot.legend(loc='upper left', bbox_to_anchor=(0.0, -0.4), prop={'size': 6})
        [t.set_fontsize(6) for t in ax.get_yticklabels()]
        [t.set_fontsize(6) for t in ax.get_xticklabels()]
        fig.savefig(filename)
        fig.savefig(filename.replace('.png', '.pdf'))
        self.do_common_args(ax, **kargs)
        fig.savefig(filename)

        real_filename = self.savefigure(fig, filename)
        config.log.info("proportional_bar: Saved '{0}'".format(real_filename))
        return real_filename

    def boxplots_vertical(self,
        filename,
        data,
        qs=None,
        title=None,
        xlims=None,
        sizer=0.022,
        vert_height=4,
        cols='lightgrey',
        bot_pad=0.1,
        showmeans=False,
        **kargs):

        assert filename, 'A filename to save the image to is required'

        plot.rcParams['pdf.fonttype'] = 42

        mmheat_hei = 0.1+(sizer*len(data))

        fig = self.getfigure(**kargs)
        fig.subplots_adjust(left=0.4, right=0.8, top=mmheat_hei, bottom=bot_pad)
        ax = fig.add_subplot(111)
        ax.tick_params(right=True)

        dats = list(data.values())
        r = ax.boxplot(
            dats,
            showfliers=False,
            whis=True,
            patch_artist=True,
            widths=0.5,
            vert=False,
            showmeans=showmeans)

        #print([i.get_data() for i in r['medians']])

        plot.setp(r['medians'], color='black', lw=2) # set nicer colours
        plot.setp(r['boxes'], color='black', lw=0.5)
        plot.setp(r['caps'], color="grey", lw=0.5)
        plot.setp(r['whiskers'], color="grey", lw=0.5)

        ax.set_yticks(numpy.arange(len(data.keys()))+1)
        ax.set_yticklabels(data.keys())

        xlim = ax.get_xlim()[1]
        if xlims:
            ax.set_xlim(xlims)
            xlim = xlims[1]

        if qs:
            for i, k, p in zip(range(0, len(data)), data, r['boxes']):
                ax.text(xlim+(xlim/8), i+1, '{:.1f}'.format(qs[k]), ha='left', va='center', fontsize=6,)

        if isinstance(cols, list):
            for i, k, p in zip(range(0, len(data)), data, r['boxes']):
                p.set_facecolor(cols[i])
        else:
            for i, k, p in zip(range(0, len(data)), data, r['boxes']):
                p.set_facecolor(cols)

        if title:
            ax.set_title(title, fontsize=6)

        [t.set_fontsize(6) for t in ax.get_yticklabels()]
        [t.set_fontsize(6) for t in ax.get_xticklabels()]

        self.do_common_args(ax, **kargs)

        fig.savefig(filename)

        real_filename = self.savefigure(fig, filename)
        config.log.info(f"boxplots_vertical: Saved '{real_filename}'")
        return real_filename

"""

PCA analysis for glbase expression objects.

"""

from operator import itemgetter

import numpy
import matplotlib.pyplot as plot
import matplotlib.patches
from mpl_toolkits.mplot3d import Axes3D, art3d
import scipy.cluster.vq
from sklearn.decomposition import PCA

from . import config
from .draw import draw
from .genelist import genelist

class manifold_pca:
    def __init__(self, parent=None, rowwise=False, feature_key_name=None, whiten=False, **kargs):
        """
        **Purpose**
            A custom class for PCA analysis of expression object data.

        **Arguments**
            parent (Required)
                The parent expression/genelist object.

            feature_key_name (ORequired)
                The key to use in the genelist to label.

            whiten (Optional, default=False)
                [From sklearn]
                When True (False by default) the components_ vectors are divided by
                n_samples times singular values to ensure uncorrelated outputs with
                unit component-wise variances.

                Whitening will remove some information from the
                transformed signal (the relative variance scales of the components) but
                can sometime improve the predictive accuracy of the downstream estimators
                by making there data respect some hard-wired assumptions.

            rowwise (Optional, default=False)
                Perform PCA on the columns or the rows

        """
        assert feature_key_name, "You must specifiy 'feature_key_name' from the original list to label the features"

        self.parent = parent

        self.matrix = self.parent.getExpressionTable().T # get a copy

        self.__draw = draw()
        self.cols = "black"
        self.rowwise = rowwise
        self.__model = None
        self.whiten = whiten
        self.feature_labels = parent[feature_key_name]
        self.labels = parent.getConditionNames() # It makes more sense to get a copy incase someone does something that reorders the list in between
        self.valid = False # Just check it's all calc'ed.

    def train(self, number_of_components, **kargs):
        '''
        **Purpose**
            Train the PCA on some array

        **Arguments**
            number_of_components (Required)
                the number of PC to collect

            whiten (Optional, default=False)
                [From sklearn]
                When True (False by default) the components_ vectors are divided by
                n_samples times singular values to ensure uncorrelated outputs with
                unit component-wise variances.

                Whitening will remove some information from the
                transformed signal (the relative variance scales of the components) but
                can sometime improve the predictive accuracy of the downstream estimators
                by making there data respect some hard-wired assumptions.

        **Returns**
            None

        '''
        if 'whiten' in kargs:
            self.whiten = kargs['whiten']
        self.__model = PCA(n_components=number_of_components, whiten=self.whiten)
        self.__transform = self.__model.fit_transform(self.matrix) # U, sample loading
        self.__components = self.__model.components_.T # V, The feature loading
        #self.__transform = self.__model.transform(self.matrix) # project the data into the PCA
        #self.__inverse_transform = self.__model.inverse_transform(self.matrix)
        self.valid = True

    def __repr__(self):
        return("<glbase.pca>")

    def __str__(self):
        ret = ["PCA object",
            "    Expression: %s" % self.parent.name,
            "    Trained   : %s" % self.valid,
            "    Whiten    : %s" % self.whiten,
            ]
        return("\n".join(ret))

    '''
    def __len__(self):
        """
        (Override)
        return the number of dimensions.
        """
        return(len(self.__d))
    '''

    def project(self, new_expn):
        """
        **Purpose**
            Give me a new expression object containing a new set of conditions and
            project those new conditions into the previously generated PCA.

            Note that it must have the same rows and row_names as the original data.

        **Arguments**
            new_expn (Required)
                a new expression object with at least one (new) condition.

        **returns**
            True if successful, and all subsequent pca methods will use the new projected data.
        """
        raise AsserionError('Not implemented')

    def explained_variance(self, filename=None, percent_variance=True, **kargs):
        """
        **Purpose**
            plot a graph of PC loading, percent variance for each componenet

        **Arguments**
            filename (Required)
                The filename to save the image to.

            <common figure arguments are supported>

        **Returns**
            None.
        """
        assert filename, "explained_variance: Must provide a filename"
        assert self.valid, 'model has not been trained, use pca.train()'

        if "aspect" not in kargs:
            kargs["aspect"] = "wide"

        expn_var = numpy.array(self.__model.explained_variance_ratio_) * 100.0

        fig = self.__draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        x = numpy.arange(len(expn_var))
        ax.bar(x, expn_var, ec="none", color="grey")
        ax.set_xlabel("Principal components")
        if percent_variance:
            ax.set_ylabel('Percent Variance')
        else:
            ax.set_ylabel("Loading")
        ax.set_xticklabels(x+1)
        ax.set_xticks(x)
        ax.set_xlim([-0.5, len(expn_var)-0.5])
        self.__draw.do_common_args(ax, **kargs)
        real_filename = self.__draw.savefigure(fig, filename)

        config.log.info("explained_variance: Saved PC loading '%s'" % real_filename)

    def get_loading_percents(self, **kargs):
        """
        **Purpose**
            Returns the percent of variance

        **Arguments**
            None

        **Returns**
            Returns an array for each PC and it's percent variance
        """
        return(numpy.array(self.__model.explained_variance_ratio_) * 100.0)

    def scatter(self, x, y, filename=None, spot_cols='grey', spots=True, label=False, alpha=0.8, overplot=None,
        spot_size=40, label_font_size=7, label_style='normal', cut=None, squish_scales=False, only_plot_if_x_in_label=None, **kargs):
        """
        **Purpose**
            plot a scatter plot of PCx against PCy.

        **Arguments**
            x, y (Required)
                PC dimension to plot as scatter
                Note that PC begin at 1 (and not at zero, as might be expected)

            filename (Required)

            spot_cols (Optional, default="black" or self.set_cols())
                list of colours for the samples, should be the same length as
                the number of conditions.

                if labels == True and spots == False and spot_cols is not None then
                    spot_cols will be used to colour the labels.

            label (Optional, default=False)
                label each spot with the name of the condition

            only_plot_if_x_in_label (Optional, default=None)
                Only plot an individual scatter if X is in the label name.

                This must be a list or tuple of names

                Allows you to effectively remove points from the PCA plot.

            overplot (Optional, default=False)
                send a list of condition names and these spots will be plotted twice, with
                spot_size +1 on the top layer.

            spots (Optional, default=True)
                Draw the spots

            alpha (Optional, default=0.8)
                alpha value to use to blend the individual points

            spot_size (Optional, default=40)
                Size of the spots on the scatter

            label_font_size (Optional, default=7)
                Size of the spot label text, only valid if label=True

            label_style (Optional, default='normal')
                add a 'style' for the text labels, one of:
                'normal', 'italic', 'oblique'

            cut (Optional, default=None)
                Send a rectangle of the form [topleftx, toplefty, bottomrightx, bottomrighty], cut out all of the items within that
                area and return their label and PC score

            squish_scales (Optional, default=False)
                set the limits very aggressively to [minmin(x), minmax(y)]

        **Returns**
            None
        """
        assert filename, "scatter: Must provide a filename"
        assert self.valid, 'model has not been trained, use pca.train()'

        labels = self.labels
        xdata = self.__transform[:,x-1]
        ydata = self.__transform[:,y-1]

        return self.__draw.unified_scatter(
            labels,
            xdata,
            ydata,
            x=x,
            y=y,
            filename=filename,
            spot_cols=spot_cols,
            spots=spots,
            label=label,
            alpha=alpha,
            overplot=overplot,
            perc_weights=self.get_loading_percents(),
            spot_size=spot_size,
            label_font_size=label_font_size,
            cut=cut,
            squish_scales=squish_scales,
            only_plot_if_x_in_label=only_plot_if_x_in_label,
            **kargs
        )

    def feature_scatter(self, x, y, filename=None, spot_cols='grey', spots=True, label=False, alpha=0.8,
        topbots=False, spot_size=40, label_font_size=7, cut=None, squish_scales=False,
        label_style='normal', **kargs):
        """
        **Purpose**
            plot a scatter plot of the loading for the features for PCx and PCy

        **Arguments**
            x, y (Required)
                PC dimension to plot as scatter
                Note that PC begin at 1 (and not at zero, as might be expected)

            label_key (Required)
                a label key in the original genelist to use as a label for the spots

            filename (Required)

            topbots (Optional, default=False)
                take on the top and bottom N on each PC axis

            spot_cols (Optional, default="black" or self.set_cols())
                list of colours for the samples, should be the same length as
                the number of conditions.

                if labels == True and spots == False and spot_cols is not None then
                    spot_cols will be used to colour the labels.

            label (Optional, default=False)
                label each spot with the name of the gene

            spots (Optional, default=True)
                Draw the spots

            alpha (Optional, default=0.8)
                alpha value to use to blend the individual points

            spot_size (Optional, default=40)
                Size of the spots on the scatter

            label_font_size (Optional, default=7)
                Size of the spot label text, only valid if label=True

            label_style (Optional, default='normal')
                add a 'style' for the text labels, one of:
                'normal', 'italic', 'oblique'

            cut (Optional, default=None)
                Send a rectangle of the form [leftx, topy, rightx, bottomy], cut out all of the items within that
                area and return their label and PC score

            squish_scales (Optional, default=False)
                set the limits very aggressively to [minmin(x), minmax(y)]

        **Returns**
            None
            You can get PC data from pca.get_uvd()
        """
        assert filename, "feature_scatter: Must provide a filename"

        labels = self.feature_labels
        xdata = self.__components[:,x-1]
        ydata = self.__components[:,y-1]

        return self.__draw.unified_scatter(
            labels,
            xdata,
            ydata,
            x=x,
            y=y,
            filename=filename,
            spot_cols=spot_cols,
            spots=spots,
            label=label,
            alpha=alpha,
            perc_weights=self.get_loading_percents(),
            label_style=label_style,
            spot_size=spot_size,
            label_font_size=label_font_size,
            cut=cut,
            squish_scales=squish_scales,
            only_plot_if_x_in_label=False,
            **kargs
        )

    def scatter3d(self, x, y, z, filename=None, spot_cols=None, label=False, stem=True,
        label_font_size=6, rotation=134, elevation=48, squish_scales=False,
        spot_size=40, depthshade=True, **kargs):
        """
        **Purpose**
            plot a scatter plot of PCx against PCy against PCz
            This is the 3D variant

        **Arguments**
            x, y, z (Required)
                PC dimension to plot as scatter
                Note that PC begin at 1 (and not at zero, as might be expected)

            spot_cols (Optional, default="black" or self.set_cols())
                list of colours for the samples, should be the same length as
                the number of conditions.

            spot_size (Optional, default=40)
                The relative spot size.

            label (Optional, default=False)
                label each spot with the name of the condition

            label_font_size (Optional, default=6)
                label font size.

            stem (Optional, default=False)
                Draw stems from the spot down to the base of the graph. (i.e. z=0)

            rotation (Optional, default=134)
                The rotation along the X, Y plane

            elevation (Optional, default=48
                The rotation along the Z plane.

            depthshade (Optional, default=True)
                turn on or off matplotlib depth shading of the points in the 3D acise

        **Returns**
            None
        """
        assert filename, "scatter: Must provide a filename"

        xdata = self.__transform[:,x-1]
        ydata = self.__transform[:,y-1]
        zdata = self.__transform[:,z-1]
        perc_weights = self.get_loading_percents()

        fig = self.__draw.getfigure(**kargs)
        ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=elevation, azim=rotation)

        cols = self.cols
        if spot_cols:
            cols = spot_cols

        sc = ax.scatter(xdata, ydata, zdata,
            #edgecolors="none",
            c=cols,
            s=spot_size,
            depthshade=depthshade
            )

        sc.set_edgecolor('none')

        if label:
            for i, lab in enumerate(self.labels):
                ax.text(xdata[i], ydata[i], zdata[i], lab, size=label_font_size, ha="center", va="bottom")

        if stem: # stem must go after scatter for sorting. Actually, not true right? matplotlib uses zorder for that...
            z_min = min(zdata)
            for x_, y_, z_ in zip(xdata, ydata, zdata):
                line = art3d.Line3D(*list(zip((x_, y_, z_min), (x_, y_, z_))), marker=None, c="grey", alpha=0.1)
                ax.add_line(line)

        ax.set_xlabel("PC%s (%.1f%%)" % (x, perc_weights[x-1])) # can be overridden via do_common_args()
        ax.set_ylabel("PC%s (%.1f%%)" % (y, perc_weights[y-1]))
        ax.set_zlabel("PC%s (%.1f%%)" % (z, perc_weights[z-1]))

        if "logx" in kargs and kargs["logx"]:
            ax.set_xscale("log", basex=kargs["logx"])
        if "logy" in kargs and kargs["logy"]:
            ax.set_yscale("log", basey=kargs["logy"])

        # squish_scales must be True as there seems to be a bug in 3D plots in guessing the scales
        # Don't worry about kargs, do_common_args will overwrite.
        ax.set_xlim([min(xdata), max(xdata)])
        ax.set_ylim([min(ydata), max(ydata)])
        ax.set_zlim([min(zdata), max(zdata)])

        #self.__draw.do_common_args(ax, **kargs)

        real_filename = self.__draw.savefigure(fig, filename)

        config.log.info("scatter3d(): Saved 'PC%s' vs 'PC%s' vs 'PC%s' scatter to '%s'" % (x, y, z, real_filename))

    def loading(self, filename=None, PC=-1, top=50, bot=50, label_key=None, all=None,
        position_override=None, selected_only=None, **kargs):
        """
        **Purpose**
            Get the loading for the items for a particular PC

        **Arguments**
            filename(Optional)
                filename to save the loading barchart to. This is optional, so you can use this function
                just top return top and bot

            PC (Required)
                The Principal Component to use

            label_key (Required)
                The key in the expression object to use to label the bar chart

            top & bot (Optional, default=50)
                the number of top and bottom genes to plot on the bar graph and
                also to return as a genelist. Set both to None to get all componenets

            all (Optional, default=False)
                if all is True, return all of the items. Make sure selected_only=None (or False)
                if you use all

            position_override (Optional, defualt=[0.3,0.03,0.3,0.96])
                specify your own positions for the barchart

            selected_only (Optional, default=None)
                Only show these selected items loading, row_loading only. Make sure all=False

        **Returns**
            topbot of the loading in a new genelist with an extra key "loadingPC<PC#>"

            The object will be based on the original expression object, and will be sorted
            based on the PC loading. This means you can do some neat stuff like:

            new_expn = pca.gene_loading(..., top=50, bot=50)
            new_expn.heatmap()
            new_expn.boxplot()

        """
        assert label_key, 'label_key is a required argument'
        assert PC>0, 'PC must be >0'
        PC -= 1 # Pad the PC so that the expected PC is returned rather than the zero-based PC.

        if not self.rowwise:
            return(self.__row_loading(filename=filename, PC=PC, top=top, bot=bot, label_key=label_key, all=all,
                position_override=position_override, selected_only=selected_only, **kargs))
        else:
            return(self.__condition_loading(filename=filename, PC=PC, top=top, bot=bot, label_key=label_key, all=all,
                position_override=position_override, **kargs))

    def __row_loading(self, filename=None, PC=-1, top=50, bot=50, label_key=None, all=False,
        position_override=None, selected_only=None, **kargs):
        """
        Internal handler for loading()
        """
        assert not self.rowwise, "loading(): You probably don't mean to use this when rowwise=True"
        assert PC >= 0, "loading(): PC of <1 specified"

        if "aspect" not in kargs:
            kargs["aspect"] = "long"

        if not position_override:
            position_override = [0.3,0.03,0.3,0.96]

        if not all and not selected_only and 'hlines' not in kargs: # put a horizontal line in at the halfway point
            kargs['hlines'] = [top-0.5,]

        if 'vlines' not in kargs:
            kargs['vlines'] = [0,]

        data = self.__components[:,PC]
        labs = self.parent[label_key]
        packed_data = [{label_key: i[0], "l": i[1]} for i in zip(labs, data)]

        sorted_data = sorted(packed_data, key=itemgetter("l"))
        data = [i["l"] for i in sorted_data]
        labs = [i[label_key] for i in sorted_data]

        if all:
            data = data
            labs = labs
        elif selected_only:
            data = []
            labs = []
            for item in sorted_data:
                if item[label_key] in set(selected_only):
                    data.append(item['l'])
                    labs.append(item[label_key])
        else:
            if bot > 0 and top > 0: # data[-0:] returns the entire list and data[0:0] returns [] !
                data = data[0:top] + data[-bot:]
                labs = labs[0:top] + labs[-bot:]
            elif top > 0:
                data = data[0:top]
                labs = labs[0:top]
            elif bot > 0:
                data = data[-bot:]
                labs = labs[-bot:]

        if filename:
            if 'size' not in kargs:
                kargs['size'] = (3,8)
            fig = self.__draw.getfigure(**kargs)
            ax = fig.add_subplot(111)
            ax.set_position(position_override)

            x = numpy.arange(len(data))
            ax.barh(x, data, ec="none", color="grey")
            ax.set_ylabel("Rows")
            ax.set_xlabel("Loading")
            ax.set_yticklabels(labs)
            ax.set_yticks(x)
            ax.set_ylim([-0.8, len(data)-0.2])
            [t.set_fontsize(6) for t in ax.get_yticklabels()]

            self.__draw.do_common_args(ax, **kargs)
            real_filename = self.__draw.savefigure(fig, filename)

            config.log.info("loading(): Saved PC loading '%s'" % real_filename)

        # work out the list to return
        newgl = genelist()
        newgl.load_list([{label_key: i[0], "pc_loading": i[1]} for i in zip(labs, data)]) # relist it so that top bot are used
        newexpn = newgl.map(genelist=self.parent, key=label_key, greedy=False)
        newexpn.sort("pc_loading")
        return(newexpn)

    def __condition_loading(self, filename=None, PC=-1, top=50, bot=50, label_key=None, all=False, **kargs):
        """
        internal alias for loading()
        """
        assert self.rowwise, "loading(): You probably don't mean to use this when rowwise=False"
        assert PC >= 0, "loading(): PC of <1 specified"

        if "aspect" not in kargs:
            kargs["aspect"] = "long"

        data = self.__components[:,PC]
        labs = self.parent._conditions
        packed_data = [{label_key: i[0], "l": i[1]} for i in zip(labs, data)]

        sorted_data = sorted(packed_data, key=itemgetter("l"))
        data = [i["l"] for i in sorted_data]
        labs = [i[label_key] for i in sorted_data]

        if all:
            data = data
            labs = labs
        else:
            if bot > 0 and top > 0: # data[-0:] returns the entire list and data[0:0] returns [] !
                data = data[0:top] + data[-bot:]
                labs = labs[0:top] + labs[-bot:]
            elif top > 0:
                data = data[0:top]
                labs = labs[0:top]
            elif bot > 0:
                data = data[-bot:]
                labs = labs[-bot:]

        if filename:
            fig = self.__draw.getfigure(**kargs)
            ax = fig.add_subplot(111)
            ax.set_position([0.3,0.03,0.6,0.96])

            x = numpy.arange(len(data))
            ax.barh(x-0.4, data, ec="black", color="grey")
            ax.set_ylabel("Columns")
            ax.set_xlabel("Loading")
            ax.set_yticklabels(labs)
            ax.set_yticks(x)
            ax.set_ylim([-0.5, len(data)-0.5])
            [t.set_fontsize(6) for t in ax.get_yticklabels()]

            self.__draw.do_common_args(ax, **kargs)
            real_filename = self.__draw.savefigure(fig, filename)

            config.log.info("loading(): Saved PC loading '%s'" % real_filename)

        # work out the list to return
        newgl = genelist()
        newgl.load_list([{"name": i[0], "pc_loading": i[1]} for i in zip(labs, data)]) # relist it so that top bot are used
        newgl.sort("pc_loading")
        return(newgl)

"""

PCA(SVD) analysis for glbase expression objects.

"""

from operator import itemgetter

import numpy
import numpy.linalg as LA
import matplotlib.pyplot as plot
import matplotlib.patches
from mpl_toolkits.mplot3d import Axes3D, art3d # not always available?
import scipy.cluster.vq 

from . import config
from .draw import draw
from .genelist import genelist

class manifold_svd:
    def __init__(self, parent=None, rowwise=False, label_key=None, whiten=False, mean_subtraction=False, **kargs):
        """
        **Purpose**
            A custom class for PCA(SVD) analysis of expression object data.
            
            Not recommended to use directly, but if you must then:
            
            expn = expression(...)
            expn.svd = pca(expn)
            
            expn.svd.loading(...)
            ...
            
        **Arguments**
            parent (Required)
                The parent expression/genelist object.
                
            whiten (Optional, default=False)
                'Normalise' the data. Each feature is divided by its standard deviation 
                across all observations to give it unit variance
                    
            mean_subtraction (Optional, default=False)
                subtract the mean from the matrix before PCA. Performed BEFORE Whitening
                
            rowwise (Optional, default=False)
                Perform PCA on the columns or the rows
                
            label_key (Optional, Required if rowwise=True)
                The key to use in the genelist to label. Required if rowwise=True
        """
        if rowwise: 
            assert label_key, "You must specifiy 'label_key' if rowwise=True"
        
        self.parent = parent
        
        matrix = numpy.array(parent.serialisedArrayDataList)
               
        self.__draw = draw()
        self.cols = "black"
        self.rowwise = rowwise
                               
        if rowwise:
            if mean_subtraction:
                config.log.info("pca(): mean subtract")
                matrix -= numpy.mean(matrix, axis=0)
            if whiten:
                config.log.info("pca(): whiten")
                matrix = scipy.cluster.vq.whiten(matrix)
            #print matrix
            self.__u, self.__d, self.__v = LA.svd(matrix, full_matrices=False)
            self.labels = parent[label_key]
            self.label_key = label_key
        else:
            matrix = matrix.T
            if mean_subtraction:
                config.log.info("pca(): mean subtract")
                matrix -= numpy.mean(matrix, axis=0)
            if whiten:
                config.log.info("pca(): whiten")
                matrix = scipy.cluster.vq.whiten(matrix)
            #print matrix
            self.__u, self.__d, self.__v = LA.svd(matrix, full_matrices=False)
            self.labels = parent.getConditionNames() # It makes more sense to get a copy incase someone 
            # does something that reorders the list in between 
            self.label_key = label_key
            
        self.__u = numpy.array(self.__u)
        self.valid = True # Just check it's all calc'ed.

    def __repr__(self):
        return("<glbase.pca>")
        
    def __str__(self):
        ret = ["SVD object",
            "    Expression: %s" % self.parent.name,
            "    Dimensions: %s" % len(self.__d),
            "    u: %s" % str(self.__u.shape),
            "    v: %s" % str(self.__v.shape)
            ]
        return("\n".join(ret))

    def __len__(self):
        """
        (Override)
        return the number of dimensions.
        """
        return(len(self.__d))

    def project(self, new_expn):
        """
        THIS IS CURRENTLY NOT WORKING
        
        **Purpose**
            Give me a new expression object containing a new set of conditions and
            project those new conditions into the previously generated PCA. 
            
        **Arguments**
            new_expn (Required)
                a new expression object with at least one condition.
                
        **returns**
            True if successful, and all subsequent pca methods will use the new projected data.
        """
        """
        data = numpy.array(self.parent.serialisedArrayDataList)
        import sklearn
        skpca = sklearn.decomposition.PCA()
        X_r = skpca.fit(data).transform(data)
        
        self.__v = X_r
        """
        # old martrisx
        matrix = numpy.array(self.parent.serialisedArrayDataList)
        U, S, V = numpy.linalg.svd(matrix.T, full_matrices=False)
        
        print("matrix", matrix.shape)
        
        # set-ups
        self.parent = new_expn
        if self.rowwise:
            self.labels = new_expn[self.label_key]
        else:
            self.labels = new_expn.getConditionNames()
        
        matrix = numpy.array(self.parent.serialisedArrayDataList)
        S = numpy.diag(S)
        print("U", U.shape)
        print("V", V.shape)
        print("S", S.shape)
        print("matrix", matrix.shape)
        
        #data = np.dot(U, np.dot(S, V))
        #X_transformed = np.dot(X_transformed, self.V.T)
        print(numpy.dot(S, V).shape)

        pr = numpy.dot(matrix, S)
        print("pr", pr.shape)
        #y = x*W;
        #y0 = Y(1,:);
        #sum(abs(y0 - y)) %
        
        # I want a new v. U and D are the same.
        
        self.__v = pr
        
        print(U)
        print()
        print(pr)
        
        print(numpy.allclose(U, pr))    
        print(numpy.allclose(matrix.T, numpy.dot(U, numpy.dot(S, V))))
        return(True)           
        
    def get_uvd(self):
        """
        **Purpose**
            Get the u, v, d matrices.
            
        **Arguments**
            None
            
        **Returns**
            A dict, containing:
            {"u": u,
            "v": v,
            "d": d}
        """
        ret = {"u": self.__u,
            "v": self.__v,
            "d": self.__d}
        return(ret)
    
    def max(self):
        """
        **Purpose**
            Return the maximum number of PC for this experiment.
        
        **Arguments**
            None
            
        **Returns**
            The number of PCs. Strictly, len(d)
        """
        return(len(self.__d))

    def set_cols(self, sample_colours):
        """
        **Purpose**
            Set the colours for the individual samples.
            
        **Arguments**
            sample_colours (Required)
                A list of sample colours (or None) for scatter plots and the like.
                Must be the same length as the number of conditions.
                
        **Returns**
            None
        """
        self.cols = sample_colours
    
    def pcloading(self, filename=None, percent_variance=True, **kargs):
        """
        **Purpose**
            plot a graph of PC loading.
            
        **Arguments**
            filename (Required)
                The filename to save the image to.
                
            percent_variance (Optional, default=True)    
                Report as the percent variance explained
                
            <common figure arguments are supported>
        
        **Returns**
            None.
        """
        assert filename, "loading(): Must provide a filename"
        
        if "aspect" not in kargs:
            kargs["aspect"] = "wide"
        
        newd = self.__d
        if percent_variance:
            newd2 = numpy.array(self.__d) **2
            newd = (newd2 / sum(newd2)) * 100.0
        
        fig = self.__draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        x = numpy.arange(len(newd))
        ax.bar(x-0.4, newd, ec="black", color="grey")
        ax.set_xlabel("Principal components")
        if percent_variance:
            ax.set_ylabel('Percent Variance')
        else:
            ax.set_ylabel("Loading")
        ax.set_xticklabels(x+1)
        ax.set_xticks(x)
        ax.set_xlim([-0.5, len(newd)-0.5])
        self.__draw.do_common_args(ax, **kargs)
        real_filename = self.__draw.savefigure(fig, filename)
        config.log.info("loading(): Saved PC loading '%s'" % real_filename)

    def get_loading_percents(self, exclude_first_pc=False, **kargs):
        """
        **Purpose**
            Returns the percent of variance:
            
            d^2 / sum(d^2)                 
            
        **Arguments**
            exclude_first_pc (Optional, default=False)
                exclude the first PC. 
                In a lot of my data I use unnormalised data, this usually means the first
                PC contains the 'shape' of the data and can dominate the overall 
                percent variance. Set this to True to ignore PC1.
            
        **Returns**
            Returns a dict in the form:
            {
            PC1: float(percent), 
            PC2: float(percent), ... 
            }
        """
        res = {}
        
        if exclude_first_pc:
            newd2 = numpy.array(self.__d[1:]) **2
        else:
            newd2 = numpy.array(self.__d) **2
        newd = (newd2 / sum(newd2)) * 100.0
        
        for i, p in enumerate(newd):
            if exclude_first_pc:
                res[i+2] = p
            else:
                res[i+1] = p
        
        return(res)
        
    def scatter(self, x, y, filename=None, spot_cols=None, spots=True, label=False, alpha=0.8, 
        spot_size=40, label_font_size=7, cut=None, squish_scales=False, only_plot_if_x_in_label=None, **kargs): 
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
                
                Allows you to effectively remove points from the PCA plot.
            
            spots (Optional, default=True)
                Draw the spots
            
            alpha (Optional, default=0.8)
                alpha value to use to blend the individual points
                
            spot_size (Optional, default=40)
                Size of the spots on the scatter
                
            label_font_size (Optional, default=7)
                Size of the spot label text, only valid if label=True
        
            cut (Optional, default=None)
                Send a rectangle of the form [topleftx, toplefty, bottomrightx, bottomrighty], cut out all of the items within that
                area and return their label and PC score    
                
            squish_scales (Optional, default=False)
                set the limits very aggressively to [minmin(x), minmax(y)]
        
        **Returns**
            None
            You can get PC data from pca.get_uvd()
        """
        assert filename, "scatter(): Must provide a filename"     

        labels = self.labels 
        xdata = self.__v[x-1]
        ydata = self.__v[y-1]
        
        return_data = self.__unified_scatter(labels, xdata, ydata, x=x, y=y, filename=filename, spot_cols=spot_cols, spots=spots, label=label, alpha=alpha, 
        spot_size=spot_size, label_font_size=label_font_size, cut=cut, squish_scales=squish_scales, only_plot_if_x_in_label=only_plot_if_x_in_label, **kargs)
        
        return(return_data)
        
    def __unified_scatter(self, labels, xdata, ydata, x, y, filename=None, spot_cols=None, spots=True, label=False, alpha=0.8, 
        spot_size=40, label_font_size=7, cut=None, squish_scales=False, only_plot_if_x_in_label=None, **kargs):
        '''
        Unified for less bugs, more fun!        
        '''
        perc_weights = self.get_loading_percents(exclude_first_pc=True)
        
        ret_data = None  
        
        if not "aspect" in kargs:
            kargs["aspect"] = "square"
        
        fig = self.__draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        
        cols = self.cols
        if spot_cols:
            cols = spot_cols            
        
        if only_plot_if_x_in_label:
            newx = []
            newy = []
            newlab = []
            newcols = []
            for i, lab in enumerate(labels):
                if only_plot_if_x_in_label in lab:
                    newx.append(xdata[i])
                    newy.append(ydata[i])
                    newlab.append(labels[i])
                    newcols.append(spot_cols[i])
            xdata = newx
            ydata = newy
            labels = newlab
            cols = newcols
            
        if spots:
            ax.scatter(xdata, ydata, s=spot_size, alpha=alpha, edgecolors="none", c=cols)
        else:
            # if spots is false then the axis limits are set to 0..1. I will have to send my
            # own semi-sensible limits:
            ax.set_xlim([min(xdata), max(xdata)])
            ax.set_ylim([min(ydata), max(ydata)])
            
        if label:
            for i, lab in enumerate(labels):
                if not spots and isinstance(spot_cols, list):
                    ax.text(xdata[i], ydata[i], lab, size=label_font_size, ha="center", va="top", color=spot_cols[i])
                else:
                    ax.text(xdata[i], ydata[i], lab, size=label_font_size, ha="center", va="top", color="black")
        
        # Tighten the axis
        if squish_scales:
            if not "xlims" in kargs:
                ax.set_xlim([min(xdata), max(xdata)])
        
            if not "ylims" in kargs:
                ax.set_ylim([min(ydata), max(ydata)])
        
        ax.set_xlabel("PC%s (%.1f%%)" % (x, perc_weights[x+1])) # can be overridden via do_common_args()
        ax.set_ylabel("PC%s (%.1f%%)" % (y, perc_weights[y+1]))
        
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
            
        self.__draw.do_common_args(ax, **kargs)
        
        real_filename = self.__draw.savefigure(fig, filename)
        config.log.info("scatter: Saved 'PC%s' vs 'PC%s' scatter to '%s'" % (x, y, real_filename)) 
        return(ret_data)

    def loading_scatter(self, x, y, label_key, filename=None, spot_cols=None, spots=True, label=False, alpha=0.8, 
        topbots=False,
        spot_size=40, label_font_size=7, cut=None, squish_scales=False, **kargs): 
        """
        **Purpose**
            plot a scatter plot of the loading values for PCx and PCy
        
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
        
            cut (Optional, default=None)
                Send a rectangle of the form [leftx, topy, rightx, bottomy], cut out all of the items within that
                area and return their label and PC score    
                
            squish_scales (Optional, default=False)
                set the limits very aggressively to [minmin(x), minmax(y)]
        
        **Returns**
            None
            You can get PC data from pca.get_uvd()
        """
        assert filename, "loading_scatter: Must provide a filename"     
        assert label_key, "loading_scatter: Must provide a label_key for the label names"
        assert label_key in self.parent, "loading_scatter(): I can't find '%s' label_key in the original genelist" % label_key

        ret_data = None      
        xdata = self.__u[:,x-1]
        ydata = self.__u[:,y-1]
        perc_weights = self.get_loading_percents(exclude_first_pc=True)
        
        labs = self.parent[label_key]
        if topbots:
            # Get the top and bot from the X and Y sorted PCs:
            sortable_data = list(zip(xdata, ydata, self.parent[label_key]))
            sorted_by_x = sorted(sortable_data, key=lambda sortable_data: sortable_data[0])
            x_tbs = list(sorted_by_x[0:topbots]) + list(sorted_by_x[-topbots:])
            sorted_by_y = sorted(sortable_data, key=lambda sortable_data: sortable_data[1])
            y_tbs = list(sorted_by_y[0:topbots]) + list(sorted_by_y[-topbots:])
            
            # Merge duplicates:
            all_items = list(set(x_tbs + y_tbs))
            
            xdata = [i[0] for i in all_items]
            ydata = [i[1] for i in all_items]
            labs = [i[2] for i in all_items]
            
        #print xdata, ydata
        
        if not "aspect" in kargs:
            kargs["aspect"] = "square"
        
        fig = self.__draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        
        cols = self.cols
        if spot_cols:
            cols = spot_cols            
        
        if spots:
            ax.scatter(xdata, ydata, s=spot_size, alpha=alpha, edgecolors="none", c=cols)
        else:
            # if spots is false then the axis limits are set to 0..1. I will have to send my
            # own semi-sensible limits:
            ax.set_xlim([min(xdata), max(xdata)])
            ax.set_ylim([min(ydata), max(ydata)])
        
        if label:
            for i, lab in enumerate(labs):
                if not spots and isinstance(spot_cols, list):
                    ax.text(xdata[i], ydata[i], lab, size=label_font_size, ha="center", va="top", color=spot_cols[i])
                else:
                    ax.text(xdata[i], ydata[i], lab, size=label_font_size, ha="center", va="top", color="black")
        
        # Tighten the axis
        if squish_scales:
            if not "xlims" in kargs:
                ax.set_xlim([min(xdata), max(xdata)])
        
            if not "ylims" in kargs:
                ax.set_ylim([min(ydata), max(ydata)])
        
        ax.set_xlabel("PC%s (%.1f%%)" % (x, perc_weights[x])) # can be overridden via do_common_args()
        ax.set_ylabel("PC%s (%.1f%%)" % (y, perc_weights[y]))
                
        if cut:
            rect = matplotlib.patches.Rectangle(cut[0:2], cut[2]-cut[0], cut[3]-cut[1], ec="none", alpha=0.2, fc="orange")
            ax.add_patch(rect)

            tdata = []
            labels = self.parent[label_key] # Just get once or big hit!
            for i in range(0, len(xdata)):
                if xdata[i] > cut[0] and xdata[i] < cut[2]:
                    if ydata[i] < cut[1] and ydata[i] > cut[3]:
                        tdata.append({"name": labels[i], "pcx": xdata[i], "pcy": ydata[i]})
            if tdata:
                ret_data = genelist()
                ret_data.load_list(tdata)
            
        self.__draw.do_common_args(ax, **kargs)
        
        real_filename = self.__draw.savefigure(fig, filename)
        config.log.info("loading_scatter: Saved 'PC%s' vs 'PC%s' scatter to '%s'" % (x, y, real_filename)) 
        return(ret_data)

    def scatter3d(self, x, y, z, filename=None, spot_cols=None, label=False, stem=False, 
        label_font_size=6, rotation=134, elevation=48, interactive=False, squish_scales=False, 
        spot_size=40, **kargs): 
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
                
            interactive (Optional, default=False)
                if True then spawn the matplotlib show() view. Note that this will pause
                execution of your script.
                Note that by default glbase uses a non-GUI matplotlib setting.
                
                You will need to fiddle around with matplotlib.use() before importing glbase         
            
            squish_scales (Optional, default=False)
                set the limits very aggressively to [minmax(x), minmax(y), minmax(z)]
        
        **Returns**
            None
            You can get PC data from pca.get_uvd()
        """
        assert filename, "scatter(): Must provide a filename"     
      
        xdata = self.__v[x-1]
        ydata = self.__v[y-1]
        zdata = self.__v[z-1]
            
        fig = self.__draw.getfigure(**kargs)
        ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=elevation, azim=rotation)
        
        cols = self.cols
        if spot_cols:
            cols = spot_cols            
        
        ax.scatter(xdata, ydata, zdata, edgecolors="none", c=cols, s=spot_size)
        if label:
            for i, lab in enumerate(self.labels):
                ax.text(xdata[i], ydata[i], zdata[i], lab, size=label_font_size, ha="center", va="bottom")
        
        if stem: # stem must go after scatter for sorting. Actually, not true right? matplotlib uses zorder for that...
            z_min = min(zdata)
            for x_, y_, z_ in zip(xdata, ydata, zdata):        
                line = art3d.Line3D(*list(zip((x_, y_, z_min), (x_, y_, z_))), marker=None, c="grey", alpha=0.1)
                ax.add_line(line)
        
        ax.set_xlabel("PC%s" % (x,)) # can be overridden via do_common_args()
        ax.set_ylabel("PC%s" % (y,))
        ax.set_zlabel("PC%s" % (z,))
        
        if "logx" in kargs and kargs["logx"]:
            ax.set_xscale("log", basex=kargs["logx"])
        if "logy" in kargs and kargs["logy"]:
            ax.set_yscale("log", basey=kargs["logy"])
        
        if squish_scales:   
            # Don't worry about kargs, do_common_args will overwrite.
            ax.set_xlim([min(xdata), max(xdata)])
            ax.set_ylim([min(ydata), max(ydata)])
            ax.set_zlim([min(zdata), max(zdata)])
        
        self.__draw.do_common_args(ax, **kargs)
        if "zlims" in kargs:
            ax.set_zlim([kargs["zlim"][0], kargs["zlim"][1]])
        
        if interactive:
            fig.show() # hope you are not on a cluster!
                
        real_filename = self.__draw.savefigure(fig, filename)
        
        config.log.info("scatter3d(): Saved 'PC%s' vs 'PC%s' vs 'PC%s' scatter to '%s'" % (x, y, z, real_filename))     

    def condition_loading(self, filename=None, PC=-1, top=50, bot=50, label_key='PC-score', all=False, **kargs):
        '''
        **Purpose**
            Get the column loading
            
        **Arguments**
            Filename (Required)
            
            PC (Required)
                The PC to get the condition loading for
        **Returns**
            A list of conditions, sorted in the same orientation as the 
        '''
        assert filename, 'condition_loading: You must specifiy a filename'
        assert PC >= 0, "condition_loading: PC of <1 specified"
        
        if "aspect" not in kargs:
            kargs["aspect"] = "long"
        
        data = self.__v[PC-1]
        labs = self.parent._conditions
        packed_data = [{label_key: i[0], "l": i[1]} for i in zip(labs, data)]
        #print packed_data
        
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

    def gene_loading(self, filename=None, PC=-1, top=50, bot=50, label_key=None, **kargs):
        """
        Deprecated.
        
        NOTE: This is an alias of loading(), please use loading() and NOT gene_loading()
        in future this method will be deleted.
        """
        return(self.loading(filename=filename, PC=PC, top=top, bot=bot, label_key=label_key, **kargs))

    def loading(self, filename=None, PC=-1, top=50, bot=50, label_key=None, all=None, **kargs):
        """
        **Purpose**
            Get the loading for the items for a particular PC
            
            (technically, get u for a particular PC)
            
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
                if all is True, return all of the items.
        
        **Returns**
            topbot of the loading in a new genelist with an extra key "loadingPC<PC#>"
            
            The object will be based on the original expression object, and will be sorted
            based on the PC loading. This means you can do some neat stuff like:
            
            new_expn = pca.gene_loading(..., top=50, bot=50)
            new_expn.heatmap()
            new_expn.boxplot()
            
        """
        PC -= 1 # Pad the PC so that the expected PC is returned rather than the zero-based PC.
        
        if not self.rowwise:
            return(self.__row_loading(filename=filename, PC=PC, top=top, bot=bot, label_key=label_key, all=all, **kargs))
        else:
            return(self.__condition_loading(filename=filename, PC=PC, top=top, bot=bot, label_key=label_key, all=all, **kargs))
        
    def __row_loading(self, filename=None, PC=-1, top=50, bot=50, label_key=None, all=False, **kargs):
        """
        Internal handler for loading()
        """
        assert not self.rowwise, "loading(): You probably don't mean to use this when rowwise=True"
        assert PC >= 0, "loading(): PC of <1 specified"
        
        if "aspect" not in kargs:
            kargs["aspect"] = "long"
        
        data = self.__u[:,PC]
        labs = self.parent[label_key]
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
            ax.set_ylabel("Rows")
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
        
        data = self.__u[:,PC]
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
        
"""

tSNE analysis for glbase expression objects.

This should really be merged with MDS and inherited...

"""

from operator import itemgetter

import numpy, random
import matplotlib.pyplot as plot
import matplotlib.patches
from mpl_toolkits.mplot3d import Axes3D, art3d
import scipy.cluster.vq
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import MiniBatchKMeans, AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph
from sklearn.neighbors import NearestCentroid
from scipy.cluster.hierarchy import dendrogram

from . import config
from .draw import draw
from .genelist import genelist

class base_manifold:
    def __init__(self, parent=None, name='none', manifold_type='base_manifold'):
        self.manifold_type = manifold_type
        self.parent = parent
        self.name = name
        self.configured = False
        self.trained = False
        self.clusters = False
        self.cluster_labels = None
        self.centroids = None
        self.__draw = draw()

    def __repr__(self):
        return "<glbase.{0}>".format(self.manifold_type)

    def __str__(self):
        ret = ["{0}} object".format(self.manifold_type),
            "\tExpression: %s" % self.parent.name,
            "\tConfigured: %s" % self.configured,
            "\tTrained   : %s" % self.trained,
            ]
        return "\n".join(ret)

    def configure(self,
        rowwise: str = False,
        feature_key_name: str = None,
        whiten: bool = False,
        random_state = None,
        verbose: int = 2,
        **kargs):
        """
        **Purpose**
            Configure the {0} Manifold

        **Arguments**
            rowwise (Optional, default=False)
                perform manifold on the rows, rather than the columns

            feature_key_name (Optional, default=False)
                if rowwise=True then this must be set to a key name
                in the expression object ot extract the row labels from.

            random_state (Optional, default=None)
                tSNE is non-determinisic
                set this to the seed you wish to use, otherwise a random number used

            whiten (Optional, default=False)
                set the data to unit variance

        """.format(self.manifold_type)
        if rowwise:
            # rowwise here is not needed
            assert feature_key_name, 'If rowwise=True then feature_key_name must also be valid'
            assert feature_key_name in list(self.parent.keys()), 'feature_key_name "%s" not found in this expression object' % feature_key_name
            self.labels = self.parent[feature_key_name]
            self.data_table = self.parent.getExpressionTable()
        else:
            self.labels = self.parent.getConditionNames()
            self.data_table = self.parent.getExpressionTable().T

        self.random_state = random_state
        random.seed(self.random_state)

        self.verbose = verbose

        self.whiten = whiten
        self.configured = True

    def scatter(self, filename=None, spot_cols='grey', spots=True, label=False, alpha=0.8,
        spot_size=40, label_font_size=7, cut=None, squish_scales=False,
        only_plot_if_x_in_label=None, draw_clusters=True, **kargs):
        """
        **Purpose**
            plot a scatter plot of the {0}.

        **Arguments**
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

                Allows you to effectively remove points from the tSNE plot.

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

            draw_clusters (Optional, default=True)
                colour the spots and label by clusters if True

        **Returns**
            None
        """.format(self.manifold_type)
        assert filename, "scatter: Must provide a filename"

        labels = self.labels
        xdata = self.npos[:, 0]
        ydata = self.npos[:, 1]

        if not self.clusters: draw_clusters=False # set to false if no clusters available;

        return_data = self.__draw.unified_scatter(labels, xdata, ydata, x=1, y=2, filename=filename,
            mode='{0} '.format(self.manifold_type), perc_weights=None,
            spot_cols=spot_cols, spots=spots, label=label, alpha=alpha,
            spot_size=spot_size, label_font_size=label_font_size, cut=cut, squish_scales=squish_scales,
            only_plot_if_x_in_label=only_plot_if_x_in_label,
            cluster_data=self.clusters, cluster_labels=self.cluster_labels, cluster_centroids=self.centroids, draw_clusters=draw_clusters,
            **kargs)

        return return_data

    def cluster(self, method=None, num_clusters=None, filename=None):
        '''
        **Purpose**
            Report louvain or leiden clusters for a trained 2D {0}

        **Arguments**
            method (Required)
                Sklearn method:

                https://scikit-learn.org/stable/modules/clustering.html#k-means

                Implemented:
                    'KMeans': The k-means (MiniBatchKMeans) algorithm. Requires a 'num_clusters' argument
                    'AgglomerativeClustering': AgglomerativeClustering. Requires a 'num_clusters' argument

            num_clusters (Required)
                The expected number of clusters.

        **Returns**
            THe cluster model and cluster labels

        '''.format(self.manifold_type)
        assert self.trained, '{0} not trained'.format(self.manifold_type)

        if self.clusters:
            config.log.warning('Overwriting exisitng cluster data')
        self.clusters = None
        self.__cluster_mode = method

        valid_methods = {'KMeans', 'AgglomerativeClustering'}
        assert method in valid_methods, 'method {0} not found'.format(method)

        xdata = self.npos[:, 0]
        ydata = self.npos[:, 1]

        if method == 'KMeans':
            assert num_clusters, 'if method is KMeans then you need a num_clusters'

            mbk = MiniBatchKMeans(init='k-means++',
                n_clusters=num_clusters,
                batch_size=100,
                n_init=50,
                max_no_improvement=10,
                verbose=1,
                random_state=self.random_state)
            labels = mbk.fit_predict(self.npos)
            self.clusters = mbk
            self.cluster_labels = labels
            self.centroids = mbk.cluster_centers_

        elif method == 'AgglomerativeClustering':
            knn_graph = kneighbors_graph(self.npos, num_clusters, include_self=False)

            self.__model = AgglomerativeClustering(
                linkage='ward',
                connectivity=knn_graph,
                n_clusters=num_clusters,
                affinity='euclidean')

            self.__full_model = AgglomerativeClustering( # For the tree;
                distance_threshold=0,
                n_clusters=None,
                linkage='ward',
                connectivity=knn_graph,
                affinity='euclidean'
                )

            self.__full_model_fp = self.__full_model.fit(self.npos)

            labels = self.__model.fit_predict(self.npos)
            self.clusters = self.__model
            self.cluster_labels = labels

            clf = NearestCentroid()
            clf.fit(self.npos, labels)
            self.centroids = clf.centroids_

        config.log.info('tsne.cluster: {0} clustered'.format(method))
        return self.clusters, self.cluster_labels, self.centroids

    def cluster_tree(self, filename, **kargs):
        """
        **Purpose**
            Draw the relationship between clusters as a tree.
            Only valid if clusering mode was 'AgglomerativeClustering'

        **Arguments**
            filename (Required)
                filename to save the image to
        """
        assert filename, 'You must specify a filename'
        assert self.__cluster_mode == 'AgglomerativeClustering', 'cluster_tree can only be used if the cluster method was AgglomerativeClustering'

        fig = self.__draw.getfigure()

        # Create linkage matrix and then plot the dendrogram

        # create the counts of samples under each node
        counts = numpy.zeros(self.__full_model_fp.children_.shape[0])
        n_samples = len(self.__full_model_fp.labels_)
        for i, merge in enumerate(self.__full_model_fp.children_):
            current_count = 0
            for child_idx in merge:
                if child_idx < n_samples:
                    current_count += 1  # leaf node
                else:
                    current_count += counts[child_idx - n_samples]
            counts[i] = current_count

        linkage_matrix = numpy.column_stack([
            self.__full_model_fp.children_,
            self.__full_model_fp.distances_,
            counts]).astype(float)

        ax = fig.add_subplot(111)
        # Plot the corresponding dendrogram
        dendrogram(linkage_matrix, ax=ax,
            truncate_mode='level', p=self.__model.n_clusters,
            **kargs)

        self.__draw.savefigure(fig, filename)

"""

UMAP analysis for glbase expression objects.

Requires umap-learn

"""

from operator import itemgetter

import numpy, random
from sklearn.decomposition import PCA

from . import config
from .genelist import genelist
from .base_manifold import base_manifold

if config.UMAP_LEARN_AVAIL:
    from umap import UMAP

class umap(base_manifold):
    def __init__(self, parent=None, name='none'):
        base_manifold.__init__(self, parent=parent, name=name, manifold_type='tsne')

    def train(self, num_pc, n_neighbors=None, min_dist=0.3):
        """
        **Purpose**
            Train the UMAP on the first <num_pc> components of a PCA

            UMAP is generally too computationally heavy to do on a full dataset, so you
            should choose the first few PCs to train the tSNE. Check the pca module
            for a PCA interface you can use to select the best PCs

        **Arguments**
            n_neighbors (Required)
                Estimated number of neighbours

            min_dist (Optional, default=0.3)
                minimum distance between points

        **Returns**
            None
        """
        assert self.configured, 'umap is not configured, run configure() first'
        assert n_neighbors, 'You must specify an estimate for n_neighbors'

        if isinstance(num_pc, int):
            self.__model = PCA(n_components=num_pc, whiten=self.whiten)
            self.__transform = self.__model.fit_transform(self.data_table)
            self.__pcas = self.__transform

        elif isinstance(num_pc, list):
            self.__model = PCA(n_components=max(num_pc)+1, whiten=self.whiten)
            self.__transform = self.__model.fit_transform(self.data_table)
            # get only the specific PCs
            self.__pcas = numpy.array([self.__transform[:,c-1] for c in num_pc]).T
        else:
            raise AssertionError('num_pcs must be either an integer or a list')

        self.__model = UMAP(
            n_components=2,
            n_neighbors=2,
            metric='correlation',
            random_state=self.random_state,
            verbose=self.verbose)

        self.npos = self.__model.fit_transform(self.__pcas)

        self.trained = True

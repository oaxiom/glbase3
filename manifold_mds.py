"""

MDS analysis for glbase expression objects.

"""

from operator import itemgetter

import numpy, random
from sklearn.decomposition import PCA
from sklearn.manifold import MDS

from . import config
from .genelist import genelist

from .base_manifold import base_manifold

class manifold_mds(base_manifold):
    def __init__(self, parent=None, name='none'):
        base_manifold.__init__(self, parent=parent, name=name, manifold_type='MDS')

    def train(self, num_pc):
        """
        **Purpose**
            Train the MDS on the first <num_pc> components of a PCA

            MDS is generally too computationally heavy to do on a full dataset, so you
            should choose the first few PCs to train the MDS. Check the pca module
            for a PCA interface you can use to select the best PCs

        **Arguments**
            num_pc (Required)
                The number of PCs of a PCA to use for MDS

                If it is an integer, MDS will use [1:num_pc]

                If it is a list MDS will only use those specific PCs.

        **Returns**
            None
        """
        assert self.configured, 'mds is not configured, run configure() first'

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

        self.__mds = MDS(n_components=2, n_jobs=1, n_init=20, verbose=self.verbose, random_state=self.random_state)
        self.npos = self.__mds.fit_transform(self.__pcas)

        self.trained = True

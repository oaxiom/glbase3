"""

tSNE analysis for glbase expression objects.

This should really be merged with MDS and inherited...

"""

from operator import itemgetter

import numpy, random
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from . import config
from .draw import draw
from .genelist import genelist
from .base_manifold import base_manifold

class tsne(base_manifold):
    def __init__(self, parent=None, name='none'):
        base_manifold.__init__(self, parent=parent, name=name, manifold_type='tsne')

    def configure(self, rowwise=False, feature_key_name=None, whiten=False,
        random_state=None, **kargs):
        """
        **Purpose**
            Configure the tSNE

        **Arguments**
            rowwise (Optional, default=False)
                perform PCA/tSNE on the rows, rather than the columns

            feature_key_name (Optional, default=False)
                if rowwise=True then this must be set to a key name
                in the expression object ot extract the row labels from.

            random_state (Optional, default=None)
                tSNE is non-determinisic
                set this to the seed you wish to use, otherwise a random number used

            whiten (Optional, default=False)
                set the data to unit variance

        """
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

        self.whiten = whiten
        self.configured = True

    def train(self, num_pc, perplexity=30):
        """
        **Purpose**
            Train the tSNE on the first <num_pc> components of a PCA

            tSNE is generally too computationally heavy to do on a full dataset, so you
            should choose the first few PCs to train the tSNE. Check the pca module
            for a PCA interface you can use to select the best PCs

        **Arguments**
            num_pc (Required)
                The number of PCs of a PCA to use for tSNE

                If it is an integer, tSNE will use [1:num_pc]

                If it is a list tSNE will only use those specific PCs.

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

        self.__model = TSNE(n_components=2,
            perplexity=perplexity,
            init='pca',
            random_state=self.random_state) # I make this deterministic
        self.npos = self.__model.fit_transform(self.__pcas)

        self.trained = True

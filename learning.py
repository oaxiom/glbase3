"""

learning for various classification tasks and techniques on expression objects

"""

from operator import itemgetter
import numpy, random

import matplotlib.pyplot as plot
import matplotlib.cm as cm
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import MiniBatchKMeans

from . import config
from .draw import draw
from .genelist import genelist

def make_meshgrid(x, y, h=.02):
    """Create a mesh of points to plot in

    Parameters
    ----------
    x: data to base x-axis meshgrid on
    y: data to base y-axis meshgrid on
    h: stepsize for meshgrid, optional

    Returns
    -------
    xx, yy : ndarray
    """
    x_min, x_max = x.min() - 1, x.max() + 1
    y_min, y_max = y.min() - 1, y.max() + 1
    xx, yy = numpy.meshgrid(numpy.arange(x_min, x_max, h),
                         numpy.arange(y_min, y_max, h))
    return xx, yy

def plot_contours(ax, clf, xx, yy, **params):
    """
    Plot the decision boundaries for a classifier.

    Parameters
    ----------
    ax: matplotlib axes object
    clf: a classifier
    xx: meshgrid ndarray
    yy: meshgrid ndarray
    params: dictionary of params to pass to contourf, optional

    """

    Z = Z.reshape(xx.shape)
    return ax.contourf(xx, yy, Z, **params)

class learning:
    def __init__(self,
        parent=None,
        name='none'):

        self.parent = parent
        self.name = name
        self.__draw = draw()

    def __repr__(self):
        return "<glbase.learning>"

    def __str__(self):
        ret = ["learning object",
            "\t"
            ]
        return "\n".join(ret)

    def configure(self,
        rowwise=False,
        feature_key_name=None,
        whiten=False,
        random_state=None,
        **kargs):
        """
        **Purpose**
            Configure the learning

        **Arguments**
            None

        """
        pass

    def svm_genes(self,
        features=None,
        fig_filename=None,
        **kargs):
        '''
        **Purpose**
            Perfrom SVM on genes for some features

        '''
        assert fig_filename, 'fig_filename=None'
        from sklearn import svm

        X = self.parent.numpy_array_all_data.T # (n_samples, n_features)
        Y = features# (n_samples)

        assert len(features) == X.shape[0], 'features ({0}) is not the same length as the condition names ({1})'.format(len(features), X.shape[1])

        clf = svm.SVC(kernel='linear', C=1.0)

        clf.fit(X, Y)

        #Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])

        config.log.info('learning.svm_genes: SVM fitted')

        fig = self.__draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        X0, X1 = X[:, 0], X[:, 4]
        #xx, yy = make_meshgrid(X0, X1)

        #plot_contours(ax, clf, xx, yy, cmap=cm.coolwarm, alpha=0.8)
        ax.scatter(X0, X1, c=Y, cmap=cm.coolwarm, s=20, edgecolors='k')

        self.__draw.savefigure(fig, fig_filename)
        config.log.info('learning.svm_genes: Saved {0}'.format(fig_filename))

    def nbayes(self,
        **kargs):
        """
        **Purpose **

        Naive Bayes (Gaussian) estimator

        **Arguments**

        """
        from sklearn.naive_bayes import GaussianNB

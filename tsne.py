"""

tSNE analysis for glbase expression objects.

This should really be merged with MDS

"""

from operator import itemgetter

import numpy, random
import matplotlib.pyplot as plot
import matplotlib.patches
from mpl_toolkits.mplot3d import Axes3D, art3d 
import scipy.cluster.vq 
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import config
from draw import draw
from genelist import genelist

class tsne:
    def __init__(self, parent=None, name='none'):
        self.parent = parent
        self.name = name
        self.configured = False
        self.trained = False
        self.__draw = draw()

    def __repr__(self):
        return("<glbase.mds>")
        
    def __str__(self):
        ret = ["tSNE object",
            "\tExpression: %s" % self.parent.name,
            "\tConfigured: %s" % self.configured,
            "\tTrained   : %s" % self.trained,
            ]
        return("\n".join(ret))
    
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
            
            random_state (Optional, default=1234567)
                MDS is non-determinisic 
                set this to the seed you wish to use.
            
            whiten (Optional, default=False)
                set the data to unit variance
                
            
        """
        if rowwise:
            # rowwise here is not needed
            assert feature_key_name, 'If rowwise=True then feature_key_name must also be valid'
            assert feature_key_name in self.parent.keys(), 'feature_key_name "%s" not found in this expression object' % feature_key_name
            self.labels = self.parent[feature_key_name]
            self.data_table = self.parent.getExpressionTable()
        else:
            self.labels = self.parent.getConditionNames()
            self.data_table = self.parent.getExpressionTable().T
        
        self.random_state = random_state
        random.seed(self.random_state)
        
        self.whiten = whiten        
        self.configured = True
        
    def train(self, num_pc):
        """
        **Purpose**
            Train the MDS on the first <num_pc> components of a PCA
        
            MDS is generally too computationally heavy to do on a full dataset, so you 
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
            raise AssertionError, 'num_pcs must be either an integer or a list'
        self.__model = TSNE(n_components=2, random_state=self.random_state) # I make this deterministic
        self.npos = self.__model.fit_transform(self.__pcas)
        
        self.trained = True
        
    def scatter(self, filename=None, spot_cols='grey', spots=True, label=False, alpha=0.8, 
        spot_size=40, label_font_size=7, cut=None, squish_scales=False, 
        only_plot_if_x_in_label=None, **kargs): 
        """
        **Purpose**
            plot a scatter plot of the tSNE.
        
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
        
        **Returns**
            None
        """
        assert filename, "scatter: Must provide a filename"     

        labels = self.labels 
        xdata = self.npos[:, 0]
        ydata = self.npos[:, 1]
        
        return_data = self.__draw.unified_scatter(labels, xdata, ydata, x=1, y=2, filename=filename, 
            mode='tSNE ', perc_weights=None, 
            spot_cols=spot_cols, spots=spots, label=label, alpha=alpha, 
            spot_size=spot_size, label_font_size=label_font_size, cut=cut, squish_scales=squish_scales, 
            only_plot_if_x_in_label=only_plot_if_x_in_label, **kargs)
        
        return(return_data)
    
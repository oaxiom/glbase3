"""
SOMPY was released under an MIT license:

The MIT License (MIT)

Copyright (c) 2014 Vahid Moosavi, www.vahidmoosavi.com

Permission is hereby granted, free of charge, to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in the Software
without restriction, including without limitation the rights to use, copy, modify,
merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

-------------

Modifications were added by Andrew P. Hutchins and are included as part of glbase,
also released under the MIT license

"""

import sys, os, timeit, math, textwrap, random, itertools
from time import time

import numpy as np
import scipy.spatial as spdist
from scipy.sparse import csr_matrix
#from sklearn.decomposition import RandomizedPCA
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn import neighbors
from sklearn.neighbors import NearestNeighbors
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.cluster.hierarchy import distance, linkage, dendrogram
from scipy.spatial.distance import pdist

try:
    from joblib import Parallel, delayed
except ImportError:
    pass

# glbase parts:
from . import config
from .draw import draw
from .genelist import genelist

def chunk_based_bmu_find(x, y, Y2):
    # This must be a pure function for Parallel
    dim = x.shape[1]
    dlen = x.shape[0]
    nnodes = y.shape[0]
    bmu = np.empty((dlen,2))

    # it seems that small batches for large dlen is really faster:
    # that is because of ddata in loops and n_jobs. for large data it slows down due to memory needs in parallel
    blen = min(50,dlen)
    i0 = 0
    d = None
    t = time()
    npd = np.dot # Speedness
    while i0+1 <= dlen:
        Low = (i0)
        High = min(dlen,i0+blen)+1

        i0 += blen
        ddata = x[Low:High]

        d = (npd(y, ddata.T) * -2) + Y2.reshape(nnodes, 1) # inline for hopeful numpy speedups

        bmu[Low:High,0] = np.argmin(d, axis=0)
        bmu[Low:High,1] = np.min(d, axis=0)

        del ddata
        d = None

    return bmu

class manifold_SOM(object):
    def __init__(self, parent, name):
        """
        **Purpose**
            self organising maps
        """
        self.name = name
        self.parent = parent
        self.codebook = None
        self.draw = draw()

        # Temp stores for gene/cluster recovery
        self.nearest_genes = None
        self.last_k_neighbours = 1

        # Stores for run order:
        self.config_run = False
        self.dlabel = False

    def config(self, threshold_value=False, digitize=False, nodenames=None, compnames=None, mapsize=None,
        CVrange=None, CVplot=None, initmethod='fullpca', seed=12345678, components=2,
        init_whiten=True,
        image_debug=False, **kargs):
        """
        **Purpose**
            configure the SOM, change options, etc.

        **Arguments**
            threshold_value (Optional, default=False)
                values in the grid below this threshold are set to 0.


            digitize (Optional, default=False)
                digitize the gene expression table. Only useable if threshold_value is used first.

                expects a number for the number of levels to digitize the data

            nodenames (Optional)
                A key name in the parent genelist to extract the SOM node names from.
                Typically gene names or accessions

            compnames (Optional)
                Change the individual SOM names (generally the conditions).

                By defualt uses parent.getConditionNames()

            mapsize (Optional)
                change the size of the map (in pixels/units/nodes)

            initmethod (Optional, default='fullpca')
                select the initialisation method for the SOM.

                fullpca uses the complete PCA model

                If the number of componenets (see below) is >2 then MDS is invoked to
                determine the final model.

                initmethod = ['fullpca']

            init_whiten (Optional, default=True)
                whiten (convert to unit variance) the data before doing the initialization specified above

            components (Optional, default=2)
                The number of PCA/MDS/etc components to use for the initialisation

            CVplot (Optional, default=False)
                Provide a filename for the CV plot, showing the Coefficient of variance (y axis)
                against the probe intensity (x axis).
                Generally a measure of heterogeneity, you will want to cut the lowest variance
                genes, using CVrange

                Assumes data is log2 transformed

            CVrange (Optional, default=False)
                Only keep genes within the CV ranges specified using [low, high]

            seed (Optional, default=12345678)
                a seed for random number usage.

            image_debug (Optional, default=False)
                Output a series of debug images showing the SOM and net training at each iteration (SLOW!)

                Set to a string for the path to put the images into.

            # Other (currently private) settings:
            algos = ['seq','batch']
            all_neigh = ['gaussian','manhatan','bubble','cut_gaussian','epanechicov' ]
            alfa_types = ['linear','inv','power']
        """
        random.seed(seed)
        self.seed = seed
        self.initmethod = initmethod
        self.algtype = 'batch'
        self.alfaini = .5
        self.alfafinal = .005
        self.neigh = 'gaussian'
        self.mapsize = mapsize
        self.image_debug = image_debug
        self.components = components
        self.init_whiten = init_whiten

        self.compnames = compnames or self.parent.getConditionNames()
        if nodenames:
            self.dlabel = self.parent[nodenames]
            self.node_names_key_name = nodenames

        if digitize:
            assert threshold_value, "If using digitize you must also set threshold_value"

        self.data_raw = self.parent.getExpressionTable() # funtion now returns a copy

        # CV plots and cutting
        cvs = None
        if CVplot:
            cvs = []
            yData = []
            xData = []
            for gene in self.parent:
                het = 2**np.array(gene["conditions"]) # This genes expression, must be 'unlogged'
                CV = het.std()/het.mean()
                yData.append(CV)
                xData.append(het.mean())
                cvs.append({"CV": CV, "tcount": het.mean(), "conditions": gene["conditions"]})

            fig = self.draw.getfigure(**kargs)
            ax = fig.add_subplot(111)
            ax.scatter(np.log2(xData), yData, c="grey", s=3, alpha=0.3, edgecolors="none")
            if CVrange:
                ax.axhline(CVrange[0], ls=":", color="red")
                ax.axhline(CVrange[1], ls=":", color="red")
            ax.set_xlabel("log2(expression value)")
            ax.set_ylabel("CV")
            actual_filename = self.draw.savefigure(fig, CVplot)
            config.log.info("som.config: saved '%s'" % actual_filename)

        if CVrange: # Do slicing here, not inside CVplot
            cvs = [] # Just generate again, fast enough
            for gene in self.parent:
                het = 2**np.array(gene["conditions"]) # This genes expression, must be 'unlogged'
                CV = het.std()/het.mean()
                yData.append(CV)
                xData.append(het.mean())
                cvs.append({"CV": CV, "tcount": het.mean(), "conditions": gene["conditions"]})

            res = []
            if nodenames:
                self.node_names_key_name = nodenames # I'm going to reload nodenames, this is safe as only used for labelling purposes, not to refer to previous list
                self.dlabel = []

            for index, cv in enumerate(cvs):
                if cv["CV"] >= CVrange[0] and cv["CV"] <= CVrange[1]:
                    res.append(self.data_raw[index])
                    if nodenames:
                        self.dlabel.append(self.parent[nodenames][index])
            self.data_raw = np.array(res)

            config.log.info("som.config: CV range sliced %s rows, %s remaining" % (len(self.parent)-self.data_raw.shape[0], self.data_raw.shape[0]))

            #print self.data_raw.shape, len(self.dlabel)

        self.data = np.copy(self.data_raw)       # get after CV slicing
        # I'm still not sure abourt the exact point raw_data should be stored.

        if threshold_value:
            assert digitize, "If using threshold_values you must also use digitize"
            self.data -= threshold_value
            self.data[self.data<0.0] = 0

            if digitize == 1:
                self.data[self.data>0.0] = 1
            elif digitize:
                self.data /= self.data.max()
                self.data *= digitize
                self.data = np.ceil(self.data) # always move up as the zeros have already been set

        self.data = normalize(self.data, method='var') # Already unit normalized?!

        self.finalise_config()

    def finalise_config(self):
        '''
        **Purpose**
            Finish the config
            Semi-internal method
        '''
        # Finalize setup:
        self.dim = self.data.shape[1]
        self.dlen = self.data.shape[0]

        self.set_topology(mapsize=self.mapsize, compname=self.compnames) # Can take None for mapsize
        self.calc_map_dist()

        if self.dlabel:
            assert len(self.dlabel) == self.dlen, "data labels length does not match the data length (%s != %s)" % (len(self.dlabel), self.dlen)

        self.config_run = True

    def set_topology(self, mapsize=None, mapshape='planar', lattice='rect', mask=None, compname=None):
        """
        **Purpose**
            set SOM topology

        **Arguments**
            mapsize
            all_mapshapes = ['planar', 'toroid', 'cylinder']
            all_lattices = ['hexa', 'rect']
        """
        self.mapshape = mapshape
        self.lattice = lattice

        #to set mask
        self.mask = np.ones([1,self.dim]) if not mask else mask
        #to set map size
        if mapsize is None: # It seems that this code does not guess correctly if mapsize is None...
            sq = math.ceil(math.sqrt(len(self.parent)))
            self.mapsize = (int(sq),int(sq))
        else:
            if len(mapsize) == 2:
                self.mapsize = [1, np.max(mapsize)] if np.min(mapsize) == 1 else mapsize
            elif len(mapsize) == 1:
                #s =  int (mapsize[0]/2)
                self.mapsize = [1 ,mapsize[0]]
                print('input was considered as node numbers')
                print('map size is [{0},{1}]'.format(s,s))
        self.nnodes = self.mapsize[0]*self.mapsize[1]

        # to set component names
        if not compname:
            try:
                cc = []
                for i in range(self.dim):
                    cc.append ('Variable-'+ str(i+1))
                    self.compnames = np.asarray(cc)[np.newaxis,:]
            except:
                print('no data yet: please first set training data to the SOM')
        else:
            assert len(compname) == self.dim, 'compname should have the same size'
            self.compname = np.asarray(compname)[np.newaxis,:]

    def calc_map_dist(self):
        """
        #calculating the grid distance, which will be called during the training steps
        #currently just works for planar grids
        """
        cd = self.nnodes
        UD2 = np.zeros((cd, cd))
        for i in range(cd):
            UD2[i,:] = self.grid_dist(i).reshape(1,cd)
        self.UD2 =  UD2

    def view_map(self, which_dim='all', pack=False, text_size=2.8, filename=None, grid=False, text=True):
        """
        **Purpose**
            ?

        """
        if np.min(self.mapsize) > 1:
            if not pack:
                self.view_2d(text_size, which_dim=which_dim)
            else:
                self.view_2d_pack(text_size, which_dim=which_dim, filename=filename, grid=grid, text=text)

        elif np.min(self.mapsize) == 1:
            self.view_1d(text_size, which_dim=which_dim)

    def init_map(self):
        """
        **Purpose**
            Initialise the map

        """
        valid_methods = ('random', 'pca', 'fullpca', 'mds')
        assert self.initmethod in valid_methods, "%s' is not a valid init method, shoose one of '%s'" % (self.initmethod, ', '.join(valid_methods))

        dim = 0
        n_nod = 0
        if self.initmethod == 'random':
            #It produces random values in the range of min-max of each dimension based on a uniform distribution
            mn = np.tile(np.min(self.data, axis=0), (self.nnodes,1))
            mx = np.tile(np.max(self.data, axis=0), (self.nnodes,1))
            self.codebook = mn + (mx-mn)*(np.random.rand(self.nnodes, self.dim))
        elif self.initmethod in ('pca', 'fullpca') or self.initmethod == 'mds':
            codebooktmp = self.lininit() #it is based on two largest eigenvalues of correlation matrix
            self.codebook = codebooktmp

    def train(self, trainlen=None, n_job=1, verbose=True):
        """
        **Purpose**
            Main training loop

        """
        assert self.config_run, "som.train: Please config() before train()"

        t0 = time()
        data = self.data
        nnodes = self.nnodes
        dlen = self.dlen
        dim = self.dim
        mapsize = self.mapsize

        #initialization
        if verbose:
            print('initialization method = %s' % self.initmethod)
            t0 = time()

        self.init_map()
        if verbose:
            print('initialization done in %f seconds' % round(time()-t0 , 3))

        self.batchtrain(njob=n_job, phase='rough', verbose=verbose)
        self.batchtrain(njob=n_job, phase='finetune', verbose=verbose)

        err = np.mean(self.bmu[1])
        if verbose:
            ts = round(time() - t0, 3)
            print()
            print("Total time elapsed: %.1f seconds" % ts)
            print("final quantization error: %.2f" % err)

    def project_data(self, data, k=1):
        """
        **Purpose**
            Project new data into an exisiting SOM
            to project a data set to a trained SOM and find the index of bmu
            It is based on nearest neighborhood search module of scikitlearn, but it is not that fast.
        """
        codebook = self.codebook
        data_raw = self.data_raw
        clf = neighbors.KNeighborsClassifier(n_neighbors=k)
        labels = np.arange(0, codebook.shape[0])
        clf.fit(codebook, labels)

        # the codebook values are all normalized
        # we can normalize the input data based on mean and std of original data
        data = normalize_by(data_raw, data, method='var')
        return clf.predict(data)

    def project_data_k_neigbours_for_index(self, data, index, k=1):
        """
        **Purpose**
            Project new data into an exisiting SOM

        """
        clf = neighbors.KNeighborsClassifier(n_neighbors=k)
        labels = np.arange(0, self.codebook.shape[0])

        print("codebook", self.codebook.shape)
        print("codebook slice", self.codebook[:,index].shape)
        print("labels", labels.shape)

        clf.fit(self.codebook, labels)
        print(self.codebook, labels)

        # the codebook values are all normalized
        # we can normalize the input data based on mean and std of original data
        data = normalize_by(self.data_raw, data, method='var')
        data = data.T
        print(data[:,index].shape)
        return clf.kneighbors(data[:,index], n_neighbors=k, return_distance=True)

    def __get_all_genes(self, codebook=None):
        """
        (Internal)

        Wrapper for the SOM projection. returns a bunch of stuff used in several functions.
        """
        if codebook is None:
            codebook = self.codebook

        data_tr = self.data_raw
        proj = self.project_data(data_tr)
        msz = self.mapsize
        coord = self.ind_to_xy(proj)

        aa = [[0]] * len(codebook)
        res = {}
        for n, i in enumerate(coord):
            if (i[0], i[1]) not in res:
                res[(i[0], i[1])] = []
            res[(i[0], i[1])].append(self.dlabel[n]) # Why not i[2]?
        return(data_tr, proj, coord, res)

    def get_all_genes(self):
        """
        **Purpose**
            Extract the locations of all of the genes

        **Arguments**
            None

        **Returns**
            A dictionary in the form:
            {
            (x0, y0): [gene1, gene2, .., genen],
            (x1, y1): [gene1, gene2, .., genen]
            }
        """
        _, _, _, res = self.__get_all_genes()

        ks = list(res.keys())
        ks.sort()
        return(res)

    def plot_gene_density(self, filename=None, topological=False, **kargs):
        """
        **Purpose**
            plot a map of the gene density in the SOM

        **Arguments**
            filename (Required)
                The filename to save the image to

            topological (Optional, default=False)
                By default a 2D image is drawn. Setting this to True will draw a
                3D topological image of the gene_density.
        """
        data_tr, proj, coord, res = self.__get_all_genes()

        hist_grid = np.zeros(self.mapsize)

        for n, i in enumerate(coord):
            hist_grid[i[0], i[1]] += 1

        fig = self.draw.getfigure()
        if not topological:
            ax = fig.add_subplot(111)
            hm = ax.imshow(hist_grid, extent=[0, self.mapsize[0], 0, self.mapsize[1]],
                aspect="auto", origin='lower', interpolation='none')
        else:
            X, Y = np.mgrid[0:self.mapsize[0], 0:self.mapsize[1]]

            ax = Axes3D(fig, elev=75, azim=95) #rect=[0, 0, .95, 1], )

            hm = ax.plot_surface(X, Y, hist_grid, cmap=cm.PuOr_r, shade=True,
                cstride=1, rstride=1, linewidth=0, antialiased=False) # antialiased stops transparent edges
            #surf = ax.contourf(X, Y, mp, cmap=cm.PuOr_r)
            ax.set_zlim(hist_grid.min(), hist_grid.max())

        ax.set_xticklabels('')
        ax.set_yticklabels('')
        ax.set_xlim(0, self.mapsize[0]-1)
        ax.set_ylim(0, self.mapsize[1]-1)
        cb = fig.colorbar(hm, orientation="vertical")
        cb.set_label("Genes in node")

        self.draw.do_common_args(fig, **kargs)
        real_filename = self.draw.savefigure(fig, filename)
        config.log.info("som.hit_map: Saved '%s'" % real_filename)

    def threshold_SOM_nodes(self, normalised_threshold=0.7, text_size=7, filename=None,
        alternate_soms=None, som_names=None, img_number_of_cols=5,return_masks=False,
        **kargs):
        """
        **Purpose**
            Threshold out nodes from a SOM, extract the genes and save them into a file.

            Additionally save an image file indicating the parts that have been thresholded

        **Arguments**
            normalised_threshold (Optional, default=0.7)
                A value between 0 and 1.0 to slice out. Essentially a % cut of the SOM.

            filename (Optional)
                save the SOM map with thresholding to a file

            img_number_of_cols (Optional, default=5)
                The number of columns of SOMs to plot.

            text_size (Optional, default=7)
                tet size for the SOM title. Only used if filename is also used.

            som_names (Optional, default=None)
                The name of the SOM(s) in a list in the data set (generally in a compname)

            alternate_soms (Optional, default=None)
                alternatively you can provide your own SOM(s) in a dictionary to threshold out

                the dict should look like: {name: SOM, ...}

            return_masks (Optional, default=False)
                By default threshold_SOM_nodes returns a set of genelists containing all of the data in
                the thresholds. If return_masks=True then instead threshold_SOM_nodes
                returns the mask arrays in a dictionary:
                {'som_name': numpy.array(mask), ...}
                Each mask is in the same format and orientation as the SOMs, with a 1 for
                thresholded and 0 for not thresholded.
        """
        soms_to_do = None
        if alternate_soms:
            assert isinstance(alternate_soms, dict), "som.threshold_SOM_nodes: expected a dict for alternate_soms"
            soms_to_do = alternate_soms
            som_order = list(soms_to_do.keys())
            som_order.sort() # Try to be helpful
        elif som_names:
            codebook = denormalize_by(self.data_raw, self.codebook) # This is wrong?
            if not isinstance(som_names, list):
                som_names = [som_names]
            som_indeces = [self.compnames.index(i) for i in som_names]
            soms_to_do = {name: codebook[:,som] for name, som in zip(som_names, som_indeces)}
            som_order = som_names
            # Check all som_anems are in the compnames...
        assert soms_to_do, 'som.threshold_SOM_nodes: You must use one of alternate_soms or som_names'

        data_tr, proj, coords, res = self.__get_all_genes()

        index_to_xy = np.arange(0, self.nnodes)
        index_to_xy = index_to_xy.reshape(self.mapsize[0], self.mapsize[1])[::-1] # The map/images are inverted for some reason
        #print "index_to_xy:", index_to_xy

        if filename:
            # Work out columns, etc
            no_rows = len(soms_to_do)/img_number_of_cols + 1
            no_cols = len(soms_to_do) if no_rows <=1 else img_number_of_cols
            h = .2
            w = .001
            fig = self.draw.getfigure(size=(no_cols*3.0*(1+w), no_rows*1.5*(1+h)))

        res = {}
        masks = {} # In case user wants the masks

        for som_index, som in enumerate(som_order):
            # Work out the actual_thresold for each SOM
            actual_threshold = ((soms_to_do[som].max()-soms_to_do[som].min()) * normalised_threshold)# + soms_to_do[som].min()

            mask = np.copy(soms_to_do[som])
            mask[mask<=actual_threshold] = 0
            mask[mask>actual_threshold] = 1

            #print "Mask:", mask

            # coord is a list of coords. Convert into t (x,y,): index dict:
            xy_coords = {}
            for n, i in enumerate(coords):
                if (i[0], i[1]) not in xy_coords:
                    xy_coords[(i[0], i[1])] = []
                xy_coords[(i[0], i[1])].append(self.dlabel[n])

            res_this_som = []
            for i in xy_coords:
                if mask[index_to_xy[i[0], i[1]]] == 1:
                    res_this_som += xy_coords[i]

            if filename:
                ax = fig.add_subplot(no_rows, no_cols*2, (som_index*2)+1) # Top is nromal SOM, bottom is the threshold
                ax.imshow(soms_to_do[som].reshape(self.mapsize[0], self.mapsize[1])[::-1],
                    extent=[0, self.mapsize[0], 0, self.mapsize[1]],
                    aspect="auto", origin='lower',
                    interpolation=config.get_interpolation_mode(filename))
                ax.set_xticklabels('')
                ax.set_yticklabels('')
                ax.set_xlim(0, self.mapsize[0]-1)
                ax.set_ylim(0, self.mapsize[1]-1)
                ax.tick_params(top="off", bottom="off", left="off", right="off")
                ax.set_title("\n".join(textwrap.wrap(som, 22)), fontsize=text_size)

                ax = fig.add_subplot(no_rows, no_cols*2, (som_index*2)+2) # Top is nromal SOM, bottom is the threshold
                ax.imshow(mask.reshape(self.mapsize[0], self.mapsize[1])[::-1],
                    extent=[0, self.mapsize[0], 0, self.mapsize[1]],
                    cmap=cm.binary_r,
                    aspect="auto", origin='lower',
                    interpolation=config.get_interpolation_mode(filename))
                ax.set_xticklabels('')
                ax.set_yticklabels('')
                ax.set_xlim(0, self.mapsize[0]-1)
                ax.set_ylim(0, self.mapsize[1]-1)
                ax.tick_params(top="off", bottom="off", left="off", right="off")
            res[som] = res_this_som
            masks[som] = mask.reshape(self.mapsize[0], self.mapsize[1])[::-1] # Be helpful

        # repack into genelists
        if not return_masks:
            newres = {}
            for k in res:
                newres[k] = genelist()
                if res[k]: # sometimes res can be empty?
                    newres[k].load_list([{self.node_names_key_name: i} for i in res[k]])

        if filename:
            #fig.tight_layout()
            plt.subplots_adjust(hspace=h,wspace=w)
            actual_filename = self.draw.savefigure(fig, filename, bbox_inches='tight') # bbox is not the same as tight_layout
            config.log.info('som.threshold_SOM_nodes: Saved "%s"' % actual_filename)
        if return_masks:
            return(masks)
        return(newres)

    def predict_by(self, data, Target, K=5, wt='distance'):
        """
        uniform
        """
        # here it is assumed that Target is the last column in the codebook
        #and data has dim-1 columns
        codebook = self.codebook
        data_raw = self.data_raw

        dim = codebook.shape[1]
        ind = np.arange(0,dim)
        indX = ind[ind != Target]
        X = codebook[:,indX]
        Y = codebook[:,Target]
        n_neighbors = K
        clf = neighbors.KNeighborsRegressor(n_neighbors, weights=wt)
        clf.fit(X, Y)
        # the codebook values are all normalized
        #we can normalize the input data based on mean and std of original data
        dimdata = data.shape[1]
        if dimdata == dim:
            data[:,Target] == 0
            data = normalize_by(data_raw, data, method='var')
            data = data[:,indX]
        elif dimdata == dim -1:
            data = normalize_by(data_raw[:,indX], data, method='var')
            #data = normalize(data, method='var')
        Predicted_values = clf.predict(data)
        Predicted_values = denormalize_by(data_raw[:,Target], Predicted_values)
        return Predicted_values

    def predict(self, X_test, K=5, wt='distance'):
        """
        uniform
        """
        # Similar to SKlearn we assume that we have X_tr, Y_tr and X_test
        # here it is assumed that Target is the last column in the codebook
        # and data has dim-1 columns
        data_raw = self.data_raw

        dim = self.codebook.shape[1]
        Target = data_raw.shape[1]-1
        X_train = self.codebook[:,:Target]
        Y_train= self.codebook[:,Target]
        n_neighbors = K
        clf = neighbors.KNeighborsRegressor(n_neighbors, weights = wt)
        clf.fit(X_train, Y_train)
        # the codebook values are all normalized
        #we can normalize the input data based on mean and std of original data
        X_test = normalize_by(data_raw[:,:Target], X_test, method='var')
        Predicted_values = clf.predict(X_test)
        Predicted_values = denormalize_by(data_raw[:,Target], Predicted_values)
        return Predicted_values

    def find_K_nodes(self, data, K=5):
        """
        (SOM Internal)
        **Purpose**
            ?

        """
        # we find the k most similar nodes to the input vector
        codebook = self.codebook
        neigh = NearestNeighbors(n_neighbors = K)
        neigh.fit(codebook)
        # the codebook values are all normalized
        #we can normalize the input data based on mean and std of original data
        data = normalize_by(self.data_raw, data, method='var')
        return neigh.kneighbors(data)

    def ind_to_xy(self, bm_ind):
        """
        (SOM Internal)
        **Purpose**
            ?

        """
        rows = self.mapsize[0]
        cols = self.mapsize[1]
        #bmu should be an integer between 0 to no_nodes
        out = np.zeros((bm_ind.shape[0], 3))
        out[:,2] = bm_ind
        out[:,0] = rows-1-bm_ind/cols
        out[:,1] = bm_ind%cols
        return out.astype(int)

    def predict_Probability(self, data, Target, K=5):
        """
        (SOM Internal)
        **Purpose**
            ?

        """
        # here it is assumed that Target is the last column in the codebook
        #and data has dim-1 columns
        codebook = self.codebook
        data_raw = self.data_raw
        dim = codebook.shape[1]
        ind = np.arange(0, dim)
        indX = ind[ind != Target]
        X = codebook[:,indX]
        Y = codebook[:,Target]
        n_neighbors = K
        clf = neighbors.KNeighborsRegressor(n_neighbors, weights='distance')
        clf.fit(X, Y)
        # the codebook values are all normalized
        # we can normalize the input data based on mean and std of original data
        dimdata = data.shape[1]
        if dimdata == dim:
            data[:,Target] == 0
            data = normalize_by(data_raw, data, method='var')
            data = data[:,indX]
        elif dimdata == dim -1:
            data = normalize_by(data_raw[:,indX], data, method='var')
            #data = normalize(data, method='var')
        weights,ind= clf.kneighbors(data, n_neighbors=K, return_distance=True)
        weights = 1./weights
        sum_ = np.sum(weights,axis=1)
        weights = weights/sum_[:,np.newaxis]
        labels = np.sign(codebook[ind,Target])
        labels[labels>=0] = 1

        #for positives
        pos_prob = labels.copy()
        pos_prob[pos_prob<0] = 0
        pos_prob = pos_prob*weights
        pos_prob = np.sum(pos_prob, axis=1)[:,np.newaxis]

        #for negatives
        neg_prob = labels.copy()
        neg_prob[neg_prob>0] = 0
        neg_prob = neg_prob*weights * -1
        neg_prob = np.sum(neg_prob, axis=1)[:,np.newaxis]

        #Predicted_values = clf.predict(data)
        #Predicted_values = denormalize_by(data_raw[:,Target], Predicted_values)
        return np.concatenate((pos_prob,neg_prob),axis=1)

    def node_Activation(self, data, wt='distance', Target=None):
        """
        (SOM Internal)
        **Purpose**
            ?

        """
        if not Target:
            codebook = self.codebook
            data_raw = self.data_raw
            clf = neighbors.KNeighborsClassifier(n_neighbors=self.nnodes)
            labels = np.arange(0, codebook.shape[0])
            clf.fit(codebook, labels)
            # the codebook values are all normalized
            # we can normalize the input data based on mean and std of original data
            data = normalize_by(data_raw, data, method='var')
            weights,ind = clf.kneighbors(data)

            weights = 1.0 / weights

            ##Softmax function
            S_  = np.sum(np.exp(weights),axis=1)[:,np.newaxis]
            weights = np.exp(weights)/S_

            return weights , ind

    def para_bmu_find(self, x, y, njb=1):
        dlen = x.shape[0]
        Y2 = None
        Y2 = np.einsum('ij,ij->i', y, y)
        bmu = None
        b = None
        #here it finds BMUs for chunk of data in parallel

        b = Parallel(n_jobs=njb, pre_dispatch='3*n_jobs')(delayed(chunk_based_bmu_find)\
            (x[i*dlen // njb:min((i+1)*dlen // njb, dlen)],y, Y2) \
            for i in range(njb))

        bmu = np.asarray(list(itertools.chain(*b))).T
        del b
        return bmu

    def update_codebook_voronoi(self, training_data, bmu, H, radius):
        """
        (SOM Internal)
        **Purpose**
        #First finds the Voronoi set of each node. It needs to calculate a smaller matrix.
        Super fast comparing to classic batch training algorithm
        # it is based on the implemented algorithm in som toolbox for Matlab by
        Helsinki university
        """
        #bmu has shape of 2,dlen, where first row has bmuinds
        # we construct ud2 from precomputed UD2 : ud2 = UD2[bmu[0,:]]
        nnodes = self.nnodes
        dlen = self.dlen
        dim = self.dim
        New_Codebook = np.empty((nnodes, dim))
        inds = bmu[0].astype(int)
        row = inds
        col = np.arange(dlen)
        val = np.tile(1,dlen)
        P = csr_matrix( (val,(row,col)), shape=(nnodes,dlen) )
        S  = np.empty((nnodes, dim))
        S = P.dot(training_data)

        # H has nnodes*nnodes and S has nnodes*dim  ---> Nominator has nnodes*dim
        Nom = np.empty((nnodes,nnodes))
        Nom = H.T.dot(S)
        #assert( Nom.shape == (nnodes, dim))
        nV = np.empty((1,nnodes))
        nV = P.sum(axis = 1).reshape(1, nnodes)
        #assert(nV.shape == (1, nnodes))
        Denom = np.empty((nnodes,1))
        Denom = nV.dot(H.T).reshape(nnodes, 1)
        #assert( Denom.shape == (nnodes, 1))
        New_Codebook = np.divide(Nom, Denom)
        Nom = None
        Denom = None
        #assert (New_Codebook.shape == (nnodes,dim))
        #setattr(som, 'codebook', New_Codebook)
        return np.around(New_Codebook, decimals = 6)

    def batchtrain(self, njob=1, phase=None, verbose=True):
        """
        (SOM Internal)
        **Purpose**
            ?

            Batch training which is called for rought training as well as finetuning

        """
        t0 = time()
        nnodes = self.nnodes
        dlen = self.dlen
        dim = self.dim
        mapsize = self.mapsize

        #############################################
        # seting the parameters
        initmethod = self.initmethod
        mn = np.min(mapsize)
        if mn == 1:
            mpd = float(nnodes*10)/float(dlen)
        else:
            mpd = float(nnodes)/float(dlen)

        ms = max(mapsize[0],mapsize[1])
        if mn == 1:
            ms = ms/5.0

        if phase == 'rough':
            #training length
            trainlen = int(np.ceil(40*mpd))
            #radius for updating
            if initmethod == 'random':
                radiusin = max(1, np.ceil(ms/2.))
                radiusfin = max(1, radiusin/8.)
            elif initmethod in ('pca', 'fullpca', 'mds'):
                radiusin = max(1, np.ceil(ms/8.))
                radiusfin = max(1, radiusin/4.)

        elif phase == 'finetune':
            #train lening length
            trainlen = int(np.ceil(40*mpd))
            #radius for updating
            if initmethod == 'random':
                radiusin = max(1, ms/8.) #from radius fin in rough training
                radiusfin = max(1, radiusin/16.)
            elif initmethod in ('pca', 'fullpca', 'mds'):
                radiusin = max(1, np.ceil(ms/8.)/2)
                radiusfin = 1

        radius = np.linspace(radiusin, radiusfin, trainlen)
        ##################################################

        UD2 = self.UD2
        #New_Codebook_V = np.empty((nnodes, dim))
        New_Codebook_V = self.codebook # ?!?!

        #X2 is part of euclidean distance (x-y)^2 = x^2 +y^2 - 2xy that we use for each data row in bmu finding.
        #Since it is a fixed value we can skip it during bmu finding for each data point, but later we need it calculate quantification error
        X2 = np.einsum('ij,ij->i', self.data, self.data)

        if verbose:
            print('%s training...' % phase)
            print('radius_ini: %.2f, radius_final: %.2f, trainlen: %d' %(radiusin, radiusfin, trainlen))

        for i in range(trainlen):
            #in case of Guassian neighborhood
            H = np.exp(-1.0*UD2/(2.0*radius[i]**2)).reshape(nnodes, nnodes)
            bmu = None
            bmu = self.para_bmu_find(self.data, New_Codebook_V, njb=njob)

            #updating the codebook
            New_Codebook_V = self.update_codebook_voronoi(self.data, bmu, H, radius)

            if verbose:
                print("Iteration: %d, quantization error: %.2f" % (i+1, np.mean(np.sqrt(bmu[1] + X2))))

        self.codebook = New_Codebook_V
        bmu[1] = np.sqrt(bmu[1] + X2)
        self.bmu = bmu

    def grid_dist(self, bmu_ind):
        """
        som and bmu_ind
        depending on the lattice "hexa" or "rect" we have different grid distance
        functions.
        bmu_ind is a number between 0 and number of nodes-1. depending on the map size
        bmu_coord will be calculated and then distance matrix in the map will be returned
        """
        try:
            lattice = self.lattice
        except:
            lattice = 'hexa'
            print('lattice not found! Lattice as hexa was set')

        if lattice == 'rect':
            return(self.rect_dist(bmu_ind))
        elif lattice == 'hexa':
            try:
                msize = self.mapsize
                rows = msize[0]
                cols = msize[1]
            except:
                rows = 0.
                cols = 0.

            #needs to be implemented
            print('to be implemented' , rows , cols)
            return np.zeros((rows,cols))

    def rect_dist(self,bmu):
        """
        (SOM Internal)
        **Purpose**
            ?
        """
        #the way we consider the list of nodes in a planar grid is that node0 is on top left corner,
        #nodemapsz[1]-1 is top right corner and then it goes to the second row.
        #no. of rows is map_size[0] and no. of cols is map_size[1]
        try:
            msize = self.mapsize
            rows = msize[0]
            cols = msize[1]
        except:
            pass

        #bmu should be an integer between 0 to no_nodes
        if 0 <= bmu <= (rows*cols):
            c_bmu = int(bmu%cols)
            r_bmu = int(bmu/cols)
        else:
            print('wrong bmu')

        #calculating the grid distance
        if np.logical_and(rows>0, cols>0):
            r,c = np.arange(0, rows, 1)[:,np.newaxis] , np.arange(0,cols, 1)
            dist2 = (r-r_bmu)**2 + (c-c_bmu)**2
            return dist2.ravel()
        else:
            print('please consider the above mentioned errors')
            return np.zeros((rows,cols)).ravel()

    def view_2d(self, text_size, which_dim='all', what='codebook'):
        """
        **Purpose**
            ?
        """
        msz0, msz1 = getattr(self, 'mapsize')
        if what == 'codebook':
            if hasattr(self, 'codebook'):
                codebook = self.codebook
                data_raw = self.data_raw
                codebook = denormalize_by(data_raw, codebook)
            else:
                print('First initialize codebook')

            if which_dim == 'all':
                dim = self.dim
                indtoshow = np.arange(0, dim).T
                ratio = float(dim)/float(dim)
                ratio = np.max((.35,ratio))
                sH, sV = 16,16*ratio*1
                plt.figure(figsize=(sH,sV))

            elif type(which_dim) == int:
                dim = 1
                indtoshow = np.zeros(1)
                indtoshow[0] = int(which_dim)
                sH, sV = 6,6
                plt.figure(figsize=(sH,sV))

            elif type(which_dim) == list:
                max_dim = codebook.shape[1]
                dim = len(which_dim)
                ratio = float(dim)/float(max_dim)
                #print max_dim, dim, ratio
                ratio = np.max((.35,ratio))
                indtoshow = np.asarray(which_dim).T
                sH, sV = 16,16*ratio*1
                plt.figure(figsize=(sH,sV))

            no_row_in_plot = dim/6 + 1 #6 is arbitrarily selected
            no_col_in_plot = dim if no_row_in_plot <= 1 else 6
            print(indtoshow)
            print(self.compnames)

            for axisNum in range(1, dim + 1):
                ax = plt.subplot(no_row_in_plot, no_col_in_plot, axisNum)
                ind = int(indtoshow[axisNum-1])
                mp = codebook[:,ind].reshape(msz0, msz1)
                pl = plt.pcolor(mp[::-1])
                print(ind, self.compnames[ind])
                plt.title(self.compnames[ind])
                font = {'size': text_size*sH/no_col_in_plot}
                plt.rc('font', **font)
                plt.axis('off')
                plt.axis([0, msz0, 0, msz1])
                ax.set_yticklabels([])
                ax.set_xticklabels([])
                plt.colorbar(pl)
            plt.show()

    def tree(self, filename=None, label_size=6, **kargs):
        """
        **Purpose**
            cluster the samples together based on the SOM

        **Arguments**
            filename (Required)
                filename to save an image to.

            label_size (Optional, default=7)
                The size of the text attached to the leaf labels.
        """
        assert filename, "You must specify a filename"
        assert self.codebook is not None, "You must initialize the codebook first"

        codebook = denormalize_by(self.data_raw, self.codebook).T

        #print codebook.shape, self.dim

        if "size" not in kargs: # resize if not specified
            kargs["size"] = (3,8)

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        ax.set_position([0.01, 0.01, 0.4, 0.98])

        dist = pdist(codebook, metric="euclidean")
        link = linkage(dist, 'complete', metric="euclidean")
        a = dendrogram(link, orientation='left', labels=self.compnames)

        ax.set_frame_on(False)
        ax.set_xticklabels("")
        ax.tick_params(top="off", bottom="off", left="off", right="off")
        [t.set_fontsize(label_size) for t in ax.get_yticklabels()]

        real_filename = self.draw.savefigure(fig, filename)
        config.log.info("som.tree(): Saved '%s'" % real_filename)
        return(None)

    def view_2d_pack(self, text_size, which_dim='all', filename=False, grid=True, text=True):
        """
        **Purpose**
            Draw a grid of the SOM maps, each map in 2D with a pseudo-color

        **Arguments**
            filename (Required)
                filename to save the image to.

            text_size
                font size for the title on each SOM

            which_dim (Optional, default='all')
                ?


        """
        assert self.codebook is not None, "You must initialize the codebook first"

        msz0, msz1 = self.mapsize
        codebook = denormalize_by(self.data_raw, self.codebook)

        if which_dim == 'all':
            dim = self.dim
            indtoshow = np.arange(0,dim).T
            ratio = float(dim)/float(dim)
            ratio = np.max((.35,ratio))
            sH, sV = 16,16*ratio*1

        elif type(which_dim) == int:
            dim = 1
            indtoshow = np.zeros(1)
            indtoshow[0] = int(which_dim)
            sH, sV = 6,6

        elif type(which_dim) == list:
            max_dim = codebook.shape[1]
            dim = len(which_dim)
            ratio = float(dim)/float(max_dim)
            #print max_dim, dim, ratio
            ratio = np.max((.35,ratio))
            indtoshow = np.asarray(which_dim).T
            sH, sV = 16,16*ratio*1

        no_row_in_plot = dim/20 + 1
        no_col_in_plot = dim if no_row_in_plot <=1 else 20
        compname = self.compnames
        h = .2
        w = .001
        fig = plt.figure(figsize=(no_col_in_plot*1.5*(1+w),no_row_in_plot*1.5*(1+h)))

        for axisNum in range(1, dim + 1):
            ax = fig.add_subplot(no_row_in_plot, no_col_in_plot, axisNum)
            ind = int(indtoshow[axisNum-1])
            mp = codebook[:,ind].reshape(msz0, msz1)

            if grid:
                pl = plt.pcolor(mp[::-1])
            else:
                plt.imshow(mp[::-1], interpolation=config.get_interpolation_mode(filename))
                plt.axis('off')

            if text:
                print(ind, self.compnames[ind])
                plt.title("\n".join(textwrap.wrap(self.compnames[ind], 19)))
                font = {'size': text_size}
                plt.rc('font', **font)

            plt.axis([0, msz0, 0, msz1])
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.xaxis.set_ticks([i for i in range(msz1)])
            ax.yaxis.set_ticks([i for i in range(msz0)])
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
                #ax.grid(True,linestyle='-', linewidth=0.5,color='k')
        plt.subplots_adjust(hspace=h,wspace=w)

        if filename:
            fig.savefig(filename, bbox_inches='tight', transparent=False, dpi=200)

        plt.close(fig)

    def lininit(self):
        """
        **Purpose**
            Initialise the seed network using PCA, RanomizedPCA or MDS when PC>2.
        """
        # X = UsigmaWT
        # XTX = Wsigma^2WT
        # T = XW = Usigma #Transformed by W EigenVector, can be calculated by
        # multiplication PC matrix by eigenval too
        # Further, we can get lower ranks by using just few of the eigen vevtors
        # T(2) = U(2)sigma(2) = XW(2) ---> 2 is the number of selected eigenvectors
        # This is how we initialize the map, just by using the first two first eigen vals and eigenvectors
        # Further, we create a linear combination of them in the new map by giving values from -1 to 1 in each
        # Direction of SOM map
        # it shoud be noted that here, X is the covariance matrix of original data

        msize = self.mapsize
        rows = msize[0]
        cols = msize[1]
        nnodes = self.nnodes

        if np.min(msize) > 1:
            coord = np.zeros((nnodes, 2))
            for i in range(0, nnodes):
                coord[i,0] = int(i/cols) #x
                coord[i,1] = int(i%cols) #y
            mx = np.max(coord, axis=0)
            mn = np.min(coord, axis=0)
            coord = (coord - mn)/(mx-mn)
            coord = (coord - .5)*2

            #data = np.copy(self.data)
            me = np.mean(self.data, 0)
            #data = (data - me)
            codebook = np.tile(me, (nnodes,1))

            if isinstance(self.components, int):
                max_comp = self.components
            elif isinstance(self.components, list):
                max_comp = max(self.components)
            else:
                raise AssertionError('components must be either an integer or a list')

            if self.initmethod in ('fullpca'):
                if self.initmethod == 'fullpca':
                    pca = PCA(n_components=max_comp, whiten=self.init_whiten)

                # Note that the init is done on the untransformed as the scipy implementation whitens and centers.
                eigvec = pca.fit_transform(self.data_raw.T) # see mds.py for transpose orders
                if isinstance(self.components, list): # get specific PCs
                    eigvec = np.array([eigvec[:,c-1] for c in self.components]).T

                eigval = pca.explained_variance_

            if self.components >2:
                config.log.info('Invoked MDS as PCs>2')
                # new code
                mds = MDS(n_components=2, n_jobs=1, n_init=20, verbose=0, random_state=self.seed) # I make this deterministic
                eigvec = mds.fit_transform(eigvec) # Not eigen values or vectors, but names maintained for code clarity
                eigval = mds.stress_ # For the normalization

                if self.image_debug:
                    self.draw.unified_scatter(self.compnames, eigvec[:, 0], eigvec[:, 1], x=1, y=2, filename="%s_md.png" % (self.image_debug,),
                        mode='MDS ', squish_scales=True)
                eigvec = eigvec.T # transpose after the draw call above

            if self.image_debug:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.scatter(eigvec[0], eigvec[1])
                fig.savefig('%s_init_int.png' % self.image_debug)
                self.eigvec = eigvec

            for j in range(nnodes):
                for i in range(eigvec.shape[0]):
                    codebook[j,:] = codebook[j, :] + coord[j,i] * eigvec[i,:]

            return np.around(codebook, decimals=6)

        elif np.min(msize) == 1:
            raise AssertionError('mapsize too small')
            coord = np.zeros((nnodes, 1))
            for i in range(0,nnodes):
                #coord[i,0] = int(i/cols) #x
                coord[i,0] = int(i%cols) #y
            mx = np.max(coord, axis = 0)
            mn = np.min(coord, axis = 0)
            #print coord

            coord = (coord - mn)/(mx-mn)
            coord = (coord - .5)*2
            #print coord
            data = getattr(self, 'data')
            me = np.mean(data, 0)
            data = (data - me)
            codebook = np.tile(me, (nnodes,1))
            pca = RandomizedPCA(n_components=1) #Randomized PCA is scalable
            #pca = PCA(n_components=2)
            pca.fit(data)
            eigvec = pca.components_
            eigval = pca.explained_variance_
            norms = np.sqrt(np.einsum('ij,ij->i', eigvec, eigvec))
            eigvec = ((eigvec.T/norms)*eigval).T; eigvec.shape

            for j in range(nnodes):
                for i in range(eigvec.shape[0]):
                    codebook[j,:] = codebook[j, :] + coord[j,i]*eigvec[i,:]
            return np.around(codebook, decimals = 6)

def normalize(data, method='var'):
    #methods  = ['var','range','log','logistic','histD','histC']
    if method == 'var':
        me = np.mean(data, axis=0)
        st = np.std(data, axis=0)
        return (data-me)/st

def normalize_by(data_raw, data, method='var'):
    #methods  = ['var','range','log','logistic','histD','histC']
    # to have the mean and std of the original data, by which SOM is trained
    me = np.mean(data_raw, axis=0)
    st = np.std(data_raw, axis=0)
    if method == 'var':
        return (data-me)/st

def denormalize_by(data_by, n_vect, n_method='var'):
    #based on the normalization
    if n_method == 'var':
        me = np.mean(data_by, axis=0)
        st = np.std(data_by, axis=0)
        return n_vect * st + me
    else:
        print('data was not normalized before')
        return n_vect

if __name__ == '__main__':
    from . import config
    config.SKLEARN_AVAIL = True # SKlearn import strangeness
    from .expression import expression
    from .helpers import glload
    import cProfile, pstats

    expn = glload('example/shared_raw_data/som_test_data_set.glb')
    expn += expn
    expn += expn # Dataset is too small

    # Can't call directly with expn.som
    s = SOM(parent=expn, name=expn.name)
    s.config(nodenames="name", threshold_value=math.log(10**1.6, 2), digitize=10)

    cProfile.run("s.train()", "profile.pro")
    p = pstats.Stats("profile.pro")
    p.strip_dirs().sort_stats("time").print_stats()

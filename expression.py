"""
**Purpose**
    An all-purpose container for transcript/gene expression data.

**To do**

"""

import sys, os, csv, string, math, copy, heapq, itertools, functools

from operator import itemgetter

import numpy, scipy
from numpy import array, arange, meshgrid, zeros, linspace, mean, object_, std # this should be deprecated
from scipy.cluster.hierarchy import distance, linkage, dendrogram
from scipy.cluster.vq import vq, kmeans, whiten, kmeans2
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plot
import matplotlib.cm as cm

from . import config, utils
from .flags import *
from .base_expression import base_expression
from .draw import draw
from .progress import progressbar
from .errors import AssertionError, ArgumentError
from .genelist import genelist
from .location import location
from .svd import svd
from .stats import stats

if config.NETWORKX_AVAIL and config.PYGRAPHVIZ_AVAIL:
    from .network import network
    # 2to3 problem:
    #if config.PYDOT_AVAIL:
    #    from bayes import bayes # Requires graphviz, but will probably port it later.

if config.SKLEARN_AVAIL:
    #from .learning import learning
    from .pca import pca
    from .mds import mds
    from .tsne import tsne
    from .som import SOM
    from .somde import somde
    from .umap import umap

if config.NETWORKX_AVAIL and config.PYGRAPHVIZ_AVAIL and config.SKLEARN_AVAIL:
    from .mdsquish import mdsquish

class expression(base_expression):
    def __init__(self, loadable_list=None, filename=None, format=None, expn=None, gzip=False, **kargs):
        """
        **Purpose**
            The base container for expression data.

            Please note that:
            Expression analysis in glbase requires easily manipulatable data.
            Examples for normalisation are
            genespring output, or output from R and PMA, Cufflinks, RSEM,
            EDASeq, DESeq, etc...

        **Arguments**
            filename (Required, one of loadable_list or filename)

                the filename of the microarray to load.

            loadable_list (Required, one of loadable_list or filename)
                a genelist-like object I can use to construct the list from. If you use this then
                expn should be a list of keys to extract the expression data from.

            expn (Required)
                If filename:

                Some sort of descriptor telling me where the expression data actually is.
                For example:
                    "column[3:]" (Take each column from column3 onwards, until the end column)
                Or:
                    "column[4::2]" (Take every second item, from column 4 until the end)
                Or:
                    "column[5:8]" (Take column 5 through 7 - 7 is not a typo. The lists are
                        zero-ordered and closed)

                In fact you can even give compound statements:
                    "[column[7], column[17]]" (Take columns 7 and 17)

                The only rules are it must be a valid piece of python code.

                If a loadable_list:

                This should be the name of the keys used to extract the expresion data

            err (Optional)
                Some sort of descriptor for where to get the error data from.
                NOT IMPLEMENTED

            cv_err (Optional)
                Some sort of descriptor for where to get confidence intervals from.
                This should be a tuple, as confidence intervals are 'lo' and 'hi'.
                NOT IMPLEMENTED

            cond_names (Optional)
                A list of the condition names (in order) if glbase is not working them
                out itself.

            name (Optional)
                By default expression will use the filename (removing any .txt, .tsv, etc) from
                the ends of the names

            silent (Optional, default=False)
                Do not output any reports (primarily this is for internal functions)

            gzip (Optional)
                is the filename gzipped?
        """
        if loadable_list:
            base_expression.__init__(self, loadable_list=loadable_list, expn=expn, **kargs)
        else:
            base_expression.__init__(self, filename=filename, expn=expn, format=format, **kargs)

    def __repr__(self):
        return("glbase.expression")

    def __getitem__(self, index):
        """
        Confers:

        a = expn["condition_name"]

        and inherits normal genelist slicing behaviour
        """
        if index in self._conditions:
            return(self.getDataForCondition(index))
        return(base_expression.__getitem__(self, index)) # otherwise inherit

    def __getattr__(self, name):
        """
        Confers this idiom:

        expn.network.genes(...)

        # These work and do not reinit as it works like this:
        # expn.bayes.learn()
        #       First call of bayes, calls __getattr__() and self.bayes = bayes()
        #       Subsequent calls will get the self.bayes slot instead of __getattr__()
        # expn.bayes.run()
        #       second call of bayes, self.bayes is already an attrib, so it will do self.bayes.run()
        #       with not new init.

        """

        if name == "network":
            assert config.NETWORKX_AVAIL, "Asking for a network object but networkx is not available"
            assert config.PYDOT_AVAIL, "Asking for a network object but pydot is not available"
            assert config.PYGRAPHVIZ_AVAIL, "Asking for a network object but pygraphviz is not available"
            self.network = network(self)
            return self.network

        elif name == "svd":
            return self.get_svd()

        elif name == "pca":
            return self.get_pca()

        elif name == "stats":
            self.stats = stats(self)
            return self.stats

        elif name == "mdsquish":
            assert config.NETWORKX_AVAIL, "Asking for mdsquish but networkx is not available"
            assert config.PYDOT_AVAIL, "Asking for a mdsquish object but pydot is not available"
            assert config.PYGRAPHVIZ_AVAIL, "Asking for a network object but pygraphviz is not available"
            self.mdsquish = mdsquish(self)
            return self.mdsquish

        elif name == "som":
            assert config.SKLEARN_AVAIL, "Asking for som but sklearn not available"
            sq = math.ceil(math.sqrt(len(self)))
            self.som = SOM(parent=self, name=self.name)
            return self.som

        elif name == 'somde':
            sq = math.ceil(math.sqrt(len(self)))
            self.somde = somde(parent=self, name=self.name)
            return self.somde

        elif name == 'mds':
            self.mds = mds(parent=self, name=self.name)
            return self.mds

        elif name == 'tsne':
            self.tsne = tsne(parent=self, name=self.name)
            return self.tsne

        elif name == 'umap':
            assert config.NETWORKX_AVAIL, "Asking for a UMAP object but umap-learn is not available"
            self.umap = umap(parent=self, name=self.name)
            return self.umap

        elif name == "bayes":
            assert config.NETWORKX_AVAIL, "Asking for a bayes object but networkx/graphviz is not available"
            assert config.PYDOT_AVAIL, "Asking for a bayes object but pydot is not available"
            assert config.PYGRAPHVIZ_AVAIL, "Asking for a network object but pygraphviz is not available"
            self.bayes = bayes(self)
            return self.bayes

        elif name == 'learning':
            assert config.SKLEARN_AVAIL, "Asking for som but sklearn not available"
            self.learning = learning(self)
            return self.learning

        raise AttributeError("'%s' object has no attribute '%s'" % (self.__repr__(), name))

    def get_pca(self, rowwise=False, label_key=None, **kargs):
        """
        **Purpose**
            Get the pca object for this expression data. For performing principal
            component analysis (PCA) on the expression data set.

            A typical workflow for PCA of expression data would go something like:

            expn = expression(filename="...", ...)
            pca = expn.get_pca()
            print pca.max() # maximum number of PCs.
            pca.loading(filename="...") # PC loading bar chart.
            pca.plot(1,2, filename="...") # PC

        **Arguments**
            whiten (Optional, default=False)
                'Normalise' the data. Each feature is divided by its standard deviation
                across all observations to give it unit variance

            mean_subtraction (Optional, default=False)
                subtract the mean from the matrix.

            rowwise (Optional, default=False)
                Perform PCA on the rows (probably genes), instead of on the conditions (the default)

            label_key (Optional, Required if rowwise=True)
                The key to use in the genelist to use as a label on plots, etc.

                Ignored if rowwise=False

        **Returns**
            A 'pca' object.
            See the documentation for pca for more details.

            You can get the pca object again in expn.pca
        """
        # This is currently hidden:
        self.__pca = pca(self, rowwise=rowwise, label_key=label_key, **kargs)
        return(self.__pca)

    def get_svd(self, rowwise=False, label_key=None, **kargs):
        """
        **Purpose**
            Get the svd object for this expression data. For performing principal
            component analysis (PCA/SVD) on the expression data set.

            A typical workflow for PCA of expression data would go something like:

            expn = expression(filename="...", ...)
            svd = expn.get_svd()
            print svd.max() # maximum number of PCs.
            svd.loading(filename="...") # PC loading bar chart.
            svd.plot(1,2, filename="...") # PC

        **Arguments**
            whiten (Optional, default=False)
                'Normalise' the data. Each feature is divided by its standard deviation
                across all observations to give it unit variance

            mean_subtraction (Optional, default=False)
                subtract the mean from the matrix.

            rowwise (Optional, default=False)
                Perform svd on the rows (probably genes), instead of on the conditions (the default)

            label_key (Optional, Required if rowwise=True)
                The key to use in the genelist to use as a label on plots, etc.

                Ignored if rowwise=False

        **Returns**
            A 'svd' object.
            See the documentation for svd for more details.

            You can get the svd object again in expn.svd
        """
        self.svd = svd(self, rowwise=rowwise, label_key=label_key, **kargs)
        return(self.svd)

    def sort_sum_expression(self, selected_conditions=None):
        """
        sort by the sum of conditions

        **Arguments**
            selected_conditions (Optional, default=None)
                You can send a list of condition names to use for the sum if you want.

                If this is none, then it uses all conditions

        """
        if selected_conditions:
            selected_condition_indeces = [self._conditions.index(i) for i in selected_conditions]
            comparator = lambda x: sum([x["conditions"][i] for i in selected_condition_indeces])
            self.linearData = sorted(self.linearData, key=comparator)
        else:
            self.linearData = sorted(self.linearData, key=lambda x: sum(x["conditions"]))
        self._optimiseData()
        return(True)

    def sort_column_sum_expression(self):
        """
        sort the conditions, left to right based on the sum of the column expression

        Like all sort_* methods this function is IN PLACE

        **Arguments**
            None

        **Returns**
            None
        """
        col_sums = numpy.sum(self.numpy_array_all_data, axis=0)

        new_condition_order = list(zip(self._conditions, col_sums))
        new_condition_order = sorted(new_condition_order, key=itemgetter(1))
        new_condition_order = [i[0] for i in new_condition_order]

        # Use an in-place variant of sliceConditions:
        newtab = []
        for name in new_condition_order:
            newtab.append(self.serialisedArrayDataDict[name])

        self._conditions = new_condition_order # Must update here for err rearrangement

        if "err" in self.linearData[0]:
            # Make an err version of serialisedArrayDataDict
            errs = numpy.copy([i["err"] for i in self.linearData])

            err_tab = None
            for name in new_condition_order:
                idx = self._conditions.index(name)
                err_col = errs[:,idx]
                if err_tab is None:
                    err_tab = err_col
                else:
                    err_tab = numpy.vstack((err_tab, err_col))

            for index, item in enumerate(self.linearData):
                item["err"] = list(err_tab[:,index])

        self.numpy_array_all_data = numpy.array(newtab).T
        self._load_numpy_back_into_linearData() # _conditions must be up to date
        self._optimiseData()
        return(True)

    def multi_sort(self, keys):
        """
        **Purpose**
            Sort a genelist using multiple keys.

            This version is tailored for expression objects, and will accept the name of a condition

        **Arguments**
            keys (Required)
                A list of key names to sort. Sorting will first sort keys[0] then key[1] through key[n]

        **Returns**
            returns True if it completes.
            sorts the list IN PLACE.
        """
        #assert key, "No such key '%s'" % key
        #assert key in self.linearData[0], "Data does not have key '%s'" % key

        comparers = []
        for k in keys:
            if k in self._conditions:
                kind = self._conditions.index(k)
                comparers.append(lambda x: x["conditions"][kind])
            else:
                comparers.append(itemgetter(k))

        def comparer(left, right):
            for fn in comparers:
                result = (fn(left)>fn(right))-(fn(left)<fn(right)) # py2.7: cmp(fn(left), fn(right)) and py3.6 bodge: (a>b)-(a<b)
                if result:
                    return result
            else:
                return 0

        self.linearData = sorted(self.linearData, key=functools.cmp_to_key(comparer))
        self._optimiseData()
        return(True)

    def sort_conditions(self, reverse=False):
        """
        **Purpose**
            Sort the Condition Names

            NOTE: In place method

        **Arguments**
            reverse (Optional, default=False)
                reverse the order of the
        """
        neworder = self.getConditionNames() # gives me a copy now.
        neworder = sorted(neworder, key=lambda s: s.lower())
        if reverse:
            neworder.reverse()

        copymask = [] # work out a mask to extract the correct array columns.
        for name in neworder:
            for i, c in enumerate(self._conditions):
                if c == name:
                    copymask.append(i)
                    break

        for index, item in enumerate(self.linearData):
            item["conditions"] = [item["conditions"][c] for c in copymask]
            if "err" in item:
                item["err"] = [item["err"][c] for c in copymask]

        self._conditions = neworder
        self._optimiseData()
        return(None)

    def findGene(self, **kargs):
        """
        Deprecated method
        """
        raise AssertioError("findGene() is deprecated, see getRowsByKey()")

    # ----------- overrides/extensions ---------------------------------

    def getGenomeName(self):
        if not self.genome:
            return("Genome not bound")
        return(self.genome.getName())

    def getColumns(self, return_keys=None, strip_expn=False):
        """
        **Purpose**
            return a new expression only containing the columns specified in return_keys (a list)

            This version for expression objects will preserve the conditions and expresion data.

        **Arguments**
            return_keys (Required)
                A list of keys to keep

            strip_expn (Optional, default=False)
                If True then remove the expression and err keys (if present).

                This will return a genelist

        **Returns**
            A new expression object or if strip_expn=True then a genelist
        """
        assert isinstance(return_keys, list), "getColumns: return_keys must be a list"
        not_found = []
        for k in return_keys:
            if k not in list(self.keys()):
                not_found.append(k)
        assert False not in [k in list(self.keys()) for k in return_keys], "key(s): '%s' not found" % (', '.join(not_found),)
        assert len(return_keys) == len(set(return_keys)), 'return_keys list is not unique'

        if strip_expn:
            newl = genelist()
            newl.name = str(self.name)
        else:
            newl = self.shallowcopy()
            newl.linearData = []
            if not "conditions" in return_keys and "conditions" in self.linearData[0]:
                return_keys.append("conditions")
            if "err" not in return_keys and "err" in self.linearData[0]:
                return_keys.append("err")

        for item in self.linearData:
            newd = {} # Could be done with dict comprehension.
            for key in return_keys:
                newd[key] = item[key] # This is wrong? It will give a view?

            newl.linearData.append(newd)
        newl._optimiseData()

        config.log.info("getColumns: got only the columns: %s" % (", ".join(return_keys),))
        return(newl)

    def strip_errs(self):
        """
        **Purpose**
            Remove the any err keys if present.
            Sometimes the err keys can be preserved inappropriately (Mostly due to running it through
            a method which does not support rescaling errors). This function removes the keys from the list

            This is an IN PLACE method.

        **Arguments**
            None

        **Returns**
            None
        """
        if "err" not in self.linearData[0]: # No errors to strip, silently fail
            config.log.warning("strip_errs: Tried to remove error data, but no error data found")
            return(None)

        for item in self.linearData:
            del item["err"]
        self._optimiseData() # I think this does nothing at the moment, but just in case I ever fix err key handling
        return(None)

    def merge(self, key=None, *tables):
        """
        **Purpose**
            Merge a bunch of expression tables by key

            This is basically a hstack for the expression data. This is required because map() can have some
            undesired effects in its treatment of expression objects.

            This method though has the disadvantage that your lists must be identical to begin with.

        **Arguments**
            tables (Required)
                A list of expression-objects to merge

            key (Required)
                The key to base the merge on

        **Returns**
            A new expression object. Only 'self' key's are maintained
        """
        assert key, "merge: You must specify a key"
        lls = [len(i) for i in tables]
        assert len(set(lls)) == 1, "merge: the expression objects must be identically sized"

        newgl = self.deepcopy()
        newgl._conditions = sum([gl._conditions for gl in tables], self._conditions)

        for item in newgl:
            others = [i._findDataByKeyLazy(key=key, value=item[key]) for i in tables]
            if None in others:
                raise AssertionError("merge: %s:%s not found in table" % (key, item[key]))
            others = [i["conditions"] for i in others]
            item["conditions"] = sum(others, item["conditions"])
            if "err" in item:
                item["err"] = sum(others, item["err"])

        newgl._optimiseData()
        return(newgl)

    def sliceConditions(self, conditions=None, **kargs):
        """
        **Purpose**

            return a copy of the expression-data, but only containing
            the condition names specified in conditions

            Note that you can also use this method to change the order of the conditions.
            Just slice all of the keys in the order you want them to occur in.

            Additionally, you can use this to replicate a key, e.g.

            gl = gl.sliceConditions(["cond1", "cond2", "cond1"])

            will now give you an expression set with two 'cond1' conditions

        **Arguments**
            conditions (Required)
                A list, or other iterable of condition names to extract
                from the expression-data. Every condition name must be present
                on the expression-data.

        **Result**

            A new expression-data object with the same settings as the original,
            but containing only the expression-data conditions specified in
            the 'conditions' argument.
        """
        assert conditions, "sliceConditions: You must specify a list of conditions to keep"
        assert isinstance(conditions, list) or isinstance(conditions, set), "sliceConditions: You must specify a list or set of conditions to keep"
        for item in conditions:
            assert item in self._conditions, "sliceConditions: '%s' condition not found in this expression data" % item

        newgl = self.deepcopy()

        newtab = [newgl.serialisedArrayDataDict[name] for name in conditions]

        # err is not stored as a serialisedArrayDataDict, so have to make one here:
        if "err" in self.keys():
            err_table = numpy.array([i["err"] for i in newgl.linearData])

            err_serialisedArrayDataDict = {}
            for index, name in enumerate(self._conditions):
                if name in conditions: # only load those we are going to slice in
                    err_serialisedArrayDataDict[name] = err_table[:,index]

            new_err_tab = numpy.array([err_serialisedArrayDataDict[name] for name in conditions]).T

            # unpack it back into the err key:
            for i, row in enumerate(new_err_tab):
                newgl.linearData[i]["err"] = list(row)

        newgl._conditions = conditions
        newgl.numpy_array_all_data = numpy.array(newtab).T
        newgl._load_numpy_back_into_linearData() # _conditions must be up to date

        newgl._optimiseData()

        config.log.info("sliceConditions: sliced for %s conditions" % (len(newgl[0]["conditions"]),))
        return(newgl)

    def getDataForCondition(self, condition_name):
        """
        **Purposse**
            get all of the expression-data data for a particular condition
            name, returns a list of all the values.
            The list remains in the same order as the overall list,

            This method returns a view of the data (not a copy)
        """
        #print self.serialisedArrayDataDict.keys()
        assert condition_name in self.getConditionNames(), "getDataForCondition: No condition named '%s' in this expression object" % condition_name

        return(self.serialisedArrayDataDict[condition_name])

    def getExpressionTable(self):
        """
        **Purpose**
            Return the entire expression table as a numpy array.
            Note that rows and columns are not labelled.

        **Arguments**
            None
        """
        return(numpy.copy(self.numpy_array_all_data))

    def subtract_mean(self):
        """
        **Purpose**
            Moves the expression values to their (column) mean. Now they will vary around zero.

        """
        self.numpy_array_all_data -= numpy.mean(self.numpy_array_all_data, axis=0)
        self._load_numpy_back_into_linearData()

    def whiten(self):
        """
        **Purpose**
            whiten the expression values.

            i.e. move them all to unit variance (i.e. 0..1)
        """
        self.numpy_array_all_data = scipy.cluster.vq.whiten(self.numpy_array_all_data)
        self._load_numpy_back_into_linearData()

    def convert_to_Z_score(self, row_wise_variance=True):
        """
        **Purpose**
            Convert the expression to a row-wise Z-score

            This is an IN PLACE function

        **Arguments**
            row_wise_variance (Optional, default=True)
                use a row_wise Z-score. If False then use the variance from all genes/rows
                on the expression object
        """
        config.log.warning('Deprecation: Please use row_Z() method')
        self.row_Z(row_wise_variance=row_wise_variance)

    def row_Z(self, row_wise_variance=True):
        """
        **Purpose**
            Convert the expression to a row-wise Z-score

            This is an IN PLACE function

        **Arguments**
            row_wise_variance (Optional, default=True)
                use a row_wise Z-score. If False then use the variance from all genes/rows
                on the expression object
        """
        expn = self.numpy_array_all_data.T
        m = numpy.mean(expn, axis=0)
        if row_wise_variance:
            s = numpy.std(expn, axis=0)
        else:
            s = numpy.std(expn)

        z = (expn - m) / s
        z[numpy.isnan(z)] = 0
        self.numpy_array_all_data = z.T

        self._load_numpy_back_into_linearData()

    def column_Z(self, col_wise_variance=True):
        '''
        **Purpose**
            Perform column Z-score conversion

            This is an IN PLACE function

        **Arguments**
            col_wise_variance (Optional, default=True)
                For calculation of the Z-score, if set to True use only the variance in each
                column.

                If set to False use the entire array

        '''
        expn = self.numpy_array_all_data
        m = numpy.mean(expn, axis=0)
        if col_wise_variance:
            s = numpy.std(expn, axis=0)
        else:
            s = numpy.std(expn)

        z = (expn - m) / s
        self.numpy_array_all_data = z
        self._load_numpy_back_into_linearData()

    def normalize(self):
        """
        **Purpose**
            Normalise the numerical data between 0 .. 1

        **Returns**
            None
            THIS IS AN IN-PLACE CONVERSION
        """
        mins = numpy.min(self.numpy_array_all_data)
        maxs = numpy.max(self.numpy_array_all_data)
        rng = maxs - mins

        self.numpy_array_all_data = (self.numpy_array_all_data - mins) / rng
        self._load_numpy_back_into_linearData()

    def normalize_columns(self):
        """
        **Purpose**
            column based normalization of data

        **Returns**
            None
            THIS IS AN IN-PLACE CONVERSION
        """
        mins = numpy.min(self.numpy_array_all_data, axis=0)
        maxs = numpy.max(self.numpy_array_all_data, axis=0)
        rng = maxs - mins

        self.numpy_array_all_data = 1.0 - (((1.0 - -1.0) * (maxs - self.numpy_array_all_data)) / rng)
        self._load_numpy_back_into_linearData()

    def normalize_rows(self):
        """
        **Purpose**
            Perform row normalisation such that max = 1.0 and min = 0.0

        **Purpose**
            row based normalization of data

        **Returns**
            None
            THIS IS AN IN-PLACE CONVERSION
        """
        #for row in xrange(data.shape[0]): # from heatmap()
        #    mi = min(data[row,:])
        #    ma = max(data[row,:])
        #    data[row,:] = (data[row,:]-mi) / (ma-mi)

        ma = numpy.max(numpy.abs(self.numpy_array_all_data), axis=1)
        mi = numpy.min(numpy.abs(self.numpy_array_all_data), axis=1)
        d = ma - mi

        self.numpy_array_all_data = numpy.apply_along_axis(lambda a: (a-mi)/d, 0, self.numpy_array_all_data)#, self.numpy_array_all_data-mi
        #self.numpy_array_all_data = numpy.apply_along_axis(lambda a: a / , 0, self.numpy_array_all_data

        #self.numpy_array_all_data = (self.numpy_array_all_data-mi) / (ma-mi)

        self._load_numpy_back_into_linearData()
        return(None)

    def digitize(self, number_of_steps, min=None, max=None):
        """
        **purpose**
            Digitize the data into the supplied <number_of_steps>.

            This will convert all of the expression values from 0 .. number_of_steps bins.

            If min and max are specified, these values will be used in place of guessing from the min and max of the data.

        **Returns**
            None
            THIS IS AN IN-PLACE CONVERSION
        """
        if min:
            mins = min
        else:
            mins = self.numpy_array_all_data.min()

        if max:
            maxs = max
        else:
            maxs = self.numpy_array_all_data.max()

        rng = maxs - mins

        self.numpy_array_all_data -= mins
        self.numpy_array_all_data /= rng
        self.numpy_array_all_data *= number_of_steps
        self.numpy_array_all_data = numpy.ceil(self.numpy_array_all_data) # always move up as the zeros have already been set

        self._load_numpy_back_into_linearData()

    def coerce(self, new_type):
        """
        **Purpose**
            Semi-internal/obscure function. Coerces the data in condition into
            the type specified by new type. Primarily this is to convert the
            expression data from/to integers or floats for downstream R problems.

        **Arguments**
            new_type (Required)
                generally int or float

        **Returns**
            None
            THIS IS AN IN-PLACE CONVERSION
        """
        if new_type == int:
            self.numpy_array_all_data = self.numpy_array_all_data.astype(numpy.int32)
            self._load_numpy_back_into_linearData()
        else:
            for item in self.linearData:
                item["conditions"] = [new_type(i) for i in item["conditions"]]
        return(None)

    def heatmap(self, filename=None, row_label_key="name", row_color_threshold=None,
        optimal_ordering=True, **kargs):
        """
        **Purpose**

            draw a simple heatmap of the current expression-data data.

        **Arguments**
            filename (Required)
                the filename of the image to save. depending upon the current
                drawing settings it will save either a png (default) svg or eps.

            bracket (Optional, default = no bracketing performed)
                bracket the data within a certain range of values.
                For example to bracket to 0 .. 1 you would use the syntax::

                    result = array.heatmap(filename="ma.png", bracket=[0,1])

                Or for something like log2 normalised array data::

                    result = array.heatmap(filename="ma.png", bracket=[-2,2])

                "Bracket' chops off the edges of the data, using this logic::

                    if value > high_bracket then value := high_bracket
                    if value < low_bracket then value := low_bracket

                See normal for a method that modifies the data.

            row_label_key (Optional, default="name")
                A key in your genelist to use to label the rows. Examples would be gene names accesion
                numbers or something else.

            normal (Optional, default = no normalising)
                Unimplemented

            row_cluster (Optional, default = True)
                cluster the rows? True or False

            row_color_threshold (Optional, default=None)
                color_threshold to color the rows clustering dendrogram.

                See also scipy.hierarchy.dendrogram

            col_cluster (Optional, default = True)
                cluster the column conditions, True or False

            log (Optional, defualt=False, True|False of 2..n for log2, log10)
                log the y axis (defaults to e)
                send an integer for the base, e.g. for log10

                log=10

                for log2

                log=2

                for mathematical constant e

                log=True
                log="e"

            row_tree (Optional, default=False)
                provide your own tree for drawing. Should be a valid dendrogram tree.
                Probably the output from tree(), although you could role your own with Scipy.
                row_labels and the data
                will be rearranged based on the tree, so don't rearrange the data yourself.

            col_tree (Optional, default=False)
                provide your own tree for ordering the data by. See row_tree for details.
                This one is applied to the columns.

            highlights (Optional, default=None)
                sometimes the row_labels will be suppressed as there is too many labels on the plot.
                But you still want to highlight a few specific genes/rows on the plot.
                Send a list to highlights that matches entries in the row_names.

            digitize (Optional, default=False)
                change the colourmap (either supplied in cmap or the default) into a 'digitized' version
                that has large blocks of colours, defined by the number you send to discretize.
                In other words, place the expression values into the number of 'digitized' bins

            cmap (Optional, default=matplotlib.cm.RdBu)
                colour map for the heatmaps. Use something like this:

                import matplotlib.cm as cm

                gl.heatmap(..., cmap=cm.afmhot)

            col_norm (Optional, default=False)
                normalise each column of data between 0 .. max => 0.0 .. 1.0

            row_norm (Optional, default=False)
                similar to the defauly output of heatmap.2 in R, rows are normalised 0 .. 1

            row_font_size (Optional, default=guess suitable size)
                the size of the row labels (in points). If set this will also override the hiding of
                labels if there are too many elements.

            col_font_size (Optional, default=8)
                the size of the column labels (in points)

            heat_wid (Optional, default=0.25)
                The width of the heatmap panel. The image goes from 0..1 and the left most
                side of the heatmap begins at 0.3 (making the heatmap span from 0.3 -> 0.55).
                You can expand or shrink this value depending wether you want it a bit larger
                or smaller.

            heat_hei (Optional, default=0.85)
                The height of the heatmap. Heatmap runs from 0.1 to heat_hei, with a maximum of 0.9 (i.e. a total of 1.0)
                value is a fraction of the entire figure size.

            colbar_label (Optional, default="expression")
                the label to place beneath the colour scale bar

            grid (Optional, default=False)
                draw a grid around each cell in the heatmap.

            draw_numbers (Optional, default=False)
                draw the values of the heatmaps in each cell see also draw_numbers_threshold

            draw_numbers_threshold (Optional, default=-9e14)
                draw the values in the cell if > draw_numbers_threshold

            draw_numbers_fmt (Optional, default= '%.1f')
                string formatting for the displayed values

                You can also send arbitrary text here, (for example, if you wanted to
                mark significane with a '*' then you could set draw_numbers_fmt='*').

            draw_numbers_font_size (Optional, default=7)
                the font size for the numbers in each cell

            imshow (Optional, default=False)
                Embed the heatmap as an image inside a vector file. (Uses matplotlib imshow
                to draw the heatmap part of the figure. Allows very large matrices to
                be saved as an svg, with the heatmap part as a raster image and all other elements
                as vectors).

            sample_label_colbar (Optional, default=None)
                add a colourbar for the samples names. This is designed for when you have too many
                conditions, and just want to show the different samples as colours

                Should be a list of colours in the same order as the condition names

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
                An improved ordering for the tree, at some computational and memory cost.
                Can be trouble on very large heatmaps

                See https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html

        **Result**
            saves an image to the 'filename' location and

            returns a dictionary containing the real_filename and re-ordered labels after clustering (if any).

            returns the 'actual filename' that really gets saved (glbase
            will modify e.g. a '.png' ending to '.svg' or '.eps' etc. depending
            upon the current setings).
        """
        # checks for option here please.
        assert filename, "heatmap: you must specify a filename"
        assert row_label_key in list(self.keys()), 'row_label_key "%s" not found in this genelist' % row_label_key

        data = self.serialisedArrayDataList

        if "log" in kargs:
            data = self.__log_transform_data(self.serialisedArrayDataList, log=kargs["log"])

        # convert it into the serialisedArrayDataDict that _heatmap() expects.
        newdata = {}
        for index, name in enumerate(self._conditions):
            newdata[name] = data[index] # get the particular column

        res = self.draw.heatmap(data=newdata,
            row_names=self[row_label_key],
            col_names=self.getConditionNames(),
            filename=filename,
            row_color_threshold=row_color_threshold,
            optimal_ordering=optimal_ordering,
            **kargs)

        config.log.info("heatmap: Saved %s" % res["real_filename"])
        return(res)

    def __fold_change(self, c1v, c2v, log=2):
        """
        Calculate the fold-change of two values
        assumes values have already been padded.
        """
        try:
            if c2v > c1v:
                if log:
                    return(math.log((c2v/c1v), log))
                else:
                    return((c2v/c1v))
            else:
                if log:
                    return(-math.log(c1v/c2v, log))
                else:
                    return(-(c1v/c2v))

        except (OverflowError, ZeroDivisionError, ValueError):
            if c2v > c1v:
                config.log.error("(%.2f/%.2f) failed" % (c2v, c1v))
            else:
                config.log.error("(%.2f/%.2f) failed" % (c1v, c2v))
            raise Exception("__fold_change() encountered an error, possibly the pad value is too small, or you are trying to apply fold-change to log transformed data")

    def normaliseToCondition(self, condition_name, use_fold_change=False, keep_normed=False, pad=1e-6,
        log=2, **kargs):
        """
        **Purpose**
            normalise all other conditions to condition_name and delete
            condition name from the expression-data list

        **Arguments**
            keep_normed (boolean, Optional)
                keep the original data, but set it to 0

            use_fold_change (Optional, default=False)
                by default glbase performs::
                    2^ -[(control_sample) / (experimental_sample)]

                this will result in values bound between 0...1 with no change being 0.5.

                If use_fold_change is set to True a value will be calulated such that
                positive values are the fold-change upwards and -ve values are fold-change down
                wards from the control_sample

            pad (Optional, default=1.0E-06)
                The amount to pad the divisor so as to avoid divide by zero errors and
                overflow errors.

            log (Optional, default=2)
                log transform the 'fold-change' by the base defined in log (usually 2, so that logFC=1 = FC=2)
                Only used if use_fold_change = True. Set to False if you don't want to log transform FC.

        **Returns**
            returns the newly normalised list
        """
        names = self.getConditionNames()

        assert condition_name in names, "normaliseToCondition: condition name '%s' is not in this expression-data" % condition_name

        name_index = names.index(condition_name)
        #print name_index

        newl = self.shallowcopy()
        newl.linearData = []
        newl._conditions = []

        p = progressbar(len(self.linearData))
        for index, item in enumerate(self.linearData):
            #print item
            old_array_data = item["conditions"]
            #print old_array_data
            new_array_data = []
            toNormal = old_array_data[name_index]+pad # stop divByZero errors.
            for i, datum in enumerate(old_array_data):
                if names[i] != condition_name:
                    if use_fold_change:
                        c1v = toNormal
                        c2v = datum+pad
                        #print c2v, c1v, (c2v/c1v), math.pow(2, -(float(toNormal) / (float(datum)+pad)))
                        new_array_data.append(self.__fold_change(c1v, c2v, log=log))
                    else:
                        new_array_data.append(math.pow(2, -(float(toNormal) / (float(datum)+pad))))
                elif keep_normed:
                    if use_fold_change:
                        new_array_data.append(0.0)
                    else:
                        new_array_data.append(0.5)

            data_copy = utils.qdeepcopy(item)
            newl.linearData.append(data_copy)
            data_copy["conditions"] = new_array_data # load the new_array_data over the old one.

            p.update(index)

        # rebuild the condition names (optimiseData can't handle this)
        if keep_normed:
            newl._conditions = names
        else:
            # delete the old label
            newl._conditions = []
            for name in names:
                if name != condition_name:
                    newl._conditions.append(name)
                elif keep_normed:
                    newl._conditions.append(name)

        newl._optimiseData()
        config.log.info("normaliseToCondition: '%s'" % condition_name)
        return(newl)

    def norm_multi_fc(self, conds=None, pad=1e-6, log=2, **kargs):
        """
        **Purpose**
            transform the data, but do it relative to multiple base lines.

            [Commonly, but incorrectly referred to as 'normalisation']

            For example suppose you had the following samples:

            sampleA, sampleB, sampleC, sampleD, sampleE

            And say you wanted to transform the expression values into fold-change, but
            sampleB should be transformed relative to sampleA whilst sampleD and E are relative to sampleC.

            You'd do this:

            newexpn = expn.norm_multi_fc({"sampleA": ["sampleB"], "sampleC": ["sampleD", "sampleE"]})

            mean_replicates is now designed to estimate the error stored in "err" key if available.

            The estimate is not guaranteed to be accurate though so use with care. Unfortunately the expression class
            cannot store assymetric errors, so the max() of the error is instead taken.

        **Arguments**
            conds (Required)
                A dictionary, in the form::
                    {"norm_sample": [A, B, list of samples to normalise to norm_sample],
                    "another_norm_sample": [C, D, other samples],
                    ...
                    }

                Each sample in the lists [] will be normalised to the sample named as the key in the dict.

            pad (Optional, default=1e-6)
                A pad value to avoid Infinity errors and DivideByZero errors.

            log (Optional, default=2)
                Base to log transform the fold-change by. Set to None/False  if you don't want to
                log transform the data.
        **Returns**
            A new expression object.
        """
        assert list(conds.keys()), "norm_multi_fc: 'conds' must be a dictionary-like object"
        assert list(conds.values()), "norm_multi_fc: 'conds' must be a dictionary-like object"
        #if "err" in self.linearData[0]:
        #    config.log.warning("'err' key in this data, norm_multi_fc() will not correct the errors, I will delete this key")

        newl = []
        newe = None
        for item in self.linearData:
            newc = [0 for n in item["conditions"]]
            for k in conds:
                kin = self._conditions.index(k)
                newc[kin] = 0.0

                for n in conds[k]:
                    nin = self._conditions.index(n)
                    newc[nin] = self.__fold_change(item["conditions"][kin]+pad, item["conditions"][nin]+pad, log=log)
                    if "err" in newc:
                        del newc["err"]

            if "err" in item:
                newe = [0 for n in item["err"]]
                for k in conds:
                    kin = self._conditions.index(k)
                    newe[kin] = 0.0 # by definition.

                    for n in conds[k]:
                        nin = self._conditions.index(n)
                        expn_lo_err = item["conditions"][nin] - item["err"][nin] + pad
                        expn_base = item["conditions"][kin] + pad # I ignore any potential error here.
                        expn_hi_err = item["conditions"][nin] + item["err"][nin] + pad
                        up = self.__fold_change(expn_base, expn_lo_err, log=log)
                        dn = abs(self.__fold_change(expn_base, expn_hi_err, log=log))
                        newe[nin] = max(up, dn)
                        #print (up,dn)

            newi = item.copy()
            newi["conditions"] = newc
            if newe and "err" in newe:
                newi["err"] = newe
            newl.append(newi)

        newgl = self.deepcopy()
        newgl.linearData = newl
        newgl._optimiseData()
        return(newgl)

    def mean_replicates(self, *reps, **kargs):
        """
        **Purpose**
            replace replicates with a single column representing the mean of the replicates.

            Don't send the same condition name as part of two other replicates. this function
            will fail silently.

            The standard error will be stored in a new key, 'err' in the same order as in each
            condition. This key is used by some downstream drawing tools for error bars.

        **Arguments**
            Accepts a single (unnamed) argument and keywords (see below)
                The condition names of replicates to merge, as a set of lists of replicates
                specifying the name of the conditions to merge as replicates. e.g.:

                newexpn = expn.mean_replicates(["condA_rp1", "condA_rp2"], ["condB_rp1", "condB_rp2", "condB_rp3"])

            threshold (Optional, default=0.8)
                by default mean_replicates will measure the Pearson correlation between samples and
                if the score is below threshold will emit a warning

            output_pears (Optional, default=False)
                If a filename, mean_replicates will output a square table for all of the Pearson
                correlation values for each pair of RNA samples, where replicates are available.

            pearson_hist (Optional, default=False)
                Save a histogram of the distribution of Pearson correlations.

        **Returns**
            A new expression object containing mean-ed values for each set of replicates.
        """
        newe = []
        all_conds = self._conditions
        all_reps = set([x for sublist in reps for x in sublist]) # Flatten the 2D list. I still don't know how this works.
        done = []
        pearson_vals = []
        # Test that all condition names are available:
        #print [c in self._conditions for c in all_reps], self._conditions
        if False in [c in self._conditions for c in all_reps]:
            missing_conds = [c for c in all_reps if c not in self._conditions]
            raise AssertionError("mean_replicates: '%s' condition names not found" % (", ".join(sorted(missing_conds)),))

        threshold = 0.8
        if "threshold" in kargs and kargs["threshold"]:
            threshold = kargs["threshold"]

        output_pears = False
        pearson_hist_filename = False

        if "output_pears" in kargs and kargs["output_pears"]:
            output_pears = kargs["output_pears"]
        if "pearson_hist" in kargs and kargs["pearson_hist"]:
            pearson_hist_filename = kargs["pearson_hist"]

        new_serialisedArrayDataDict = {}
        errors = {}
        new_condition_name_list = []

        for cind, cond in enumerate(self._conditions):
            if cond in all_reps: # Its one of the reps to merge
                if cond not in done: # check it's not done already:

                    # get the p it is in:
                    p = [i for i in reps if cond in i][0] # the indeces of the replicates to merge, the [0] is so that if the rep is in two sets I don't end up with multiple sets
                    expn_vals = numpy.array([self.serialisedArrayDataDict[i] for i in p])
                    mean = sum(expn_vals) / float(len(p))
                    err = numpy.std(expn_vals, axis=0) / math.sqrt(len(expn_vals))
                    new_serialisedArrayDataDict[p[0]] = mean # merge into the 0th replicate key
                    errors[p[0]] = err
                    new_condition_name_list.append(p[0])

                    # add all reps to the done list:
                    [done.append(i) for i in p]
            else: # not to be merged or modified, so just add it to conditions.
                new_serialisedArrayDataDict[cond] = self.serialisedArrayDataDict[cond] # merge into the 0th replicate key
                errors[cond] = numpy.zeros(len(self.serialisedArrayDataDict[cond]))
                new_condition_name_list.append(cond)

        # reload self.serialisedArrayDataDict back into numpy_array_all_data
        newgl = self.deepcopy()
        newgl._conditions = new_condition_name_list

        # Get an numpy array of all of the data:
        newgl.serialisedArrayDataDict = new_serialisedArrayDataDict
        newgl.numpy_array_all_data = numpy.vstack([new_serialisedArrayDataDict[i] for i in new_condition_name_list])
        newgl.numpy_array_all_data = newgl.numpy_array_all_data.T

        # Do the same for errors:
        # Get an array of all errors:
        error_stack = numpy.vstack([errors[i] for i in new_condition_name_list]).T
        for i, row in enumerate(error_stack):
            newgl.linearData[i]["err"] = list(row)

        # Load back
        newgl._load_numpy_back_into_linearData() # Already calls _optimiseData()

        pear_out = numpy.zeros([len(self._conditions), len(self._conditions)])
        # check pairs for pearson correlation
        for r in reps:
            # r can be a list n entries long. I need to compare all vs all
            for i1, p1 in enumerate(r):
                for i2, p2 in enumerate(r):
                    if i1 != i2 and i1 < i2:
                        p1d = self.getDataForCondition(p1)
                        p2d = self.getDataForCondition(p2)
                        corr = scipy.stats.pearsonr(p1d, p2d)
                        if corr[0] < threshold:
                            config.log.warning("Samples '%s' vs '%s', pearson=%.2f" % (p1, p2, corr[0]))
                        if output_pears:
                            p1ind = self._conditions.index(p1)
                            p2ind = self._conditions.index(p2)
                            pear_out[p1ind, p2ind] = corr[0]
                            pear_out[p2ind, p1ind] = corr[0]
                        pearson_vals.append(corr[0])
                    elif i1 == i2 and output_pears:
                        p1ind = self._conditions.index(p1)
                        pear_out[p1ind, p1ind] = 1.0

        if output_pears:
            oh = open(output_pears, "w")
            oh.write("\t%s\n" % "\t".join(self._conditions))
            for i, row in enumerate(pear_out):
                d_out = []
                for d in row:
                    if d == 0:
                        d_out.append("")
                    else:
                        d_out.append(str(d))
                oh.write("%s\t%s\n" % (self._conditions[i], "\t".join(d_out)))
            oh.close()
            config.log.info("mean_replicates: Saved Pearson correlation table '%s'" % (output_pears,))

        if pearson_hist_filename:
            fig = self.draw.getfigure(**kargs)
            axis = fig.add_subplot(111)
            axis.hist(pearson_vals, bins=50, range=[0, 1], facecolor='grey', ec="none", alpha=1.0)
            axis.set_xlim([0,1])

            config.log.info("mean_replicates: Saved Pearson histogram '%s'" % self.draw.savefigure(fig, pearson_hist_filename))

        config.log.info("mean_replicates: Started with %s conditions, ended with %s" % (len(self._conditions), len(newgl[0]["conditions"])))
        return(newgl)

    def add_fc_key(self, key="fc", cond1=None, cond2=None, log=2, pad=1.0E-06, and_err=False, **kargs):
        """
        **Purpose**
            Add in a fold-change key for the fold change from cond1 to cond2.

            Note that it will pad the values by 1.0E-08 to avoid divide by zero errors.

            You can change this value with the 'pad' argument if it is too large/small

        **Arguments**
            key (Optional, default="fc")
                The key name to load the fold-change value into

            cond1 (Required)
                condition name 1

            cond2 (Required)
                condition name 2

            pad (Optional, default=1.0E-06)
                The amount to pad the divisor so as to avoid divide by zero errors and
                overflow errors.

            log (Optional, default=2)
                By default fold-changes are log2 transformed. This means a fold-change of
                2.0 becomes 1.0 and no change is no 0.0. Set this to None to disable this
                behaviour

            and_err (Optional, default=False)
                and estimate the err and load into a key name as specified by and_err.
                Note that your expression object must have an 'err' key to estimate the error.
                Also, expresion objects can't hold assymetric error bars. Hence the minimum value
                will be taken. Now, that may sound a bit odd, but as fold-changes are commonly
                log transformed taking the maximum value will overestimate the error, whilst the min
                value will show an accurate 'upward' error bar and a somewhat inaccurate downward error
                bar.

        **Returns**
            A new expression object containing the <fc> key.
        """
        assert cond1, "add_fc_key: you must specify a name for condition1"
        assert cond1 in self._conditions, "add_fc_key: appears '%s' not in this expression-object" % cond1
        assert cond2, "add_fc_key: you must specify a name for condition1"
        assert cond2 in self._conditions, "add_fc_key: appears '%s' not in this expression-object" % cond2
        if and_err:
            assert "err" in self.linearData[0], "add_fc_key: 'err' key not found in this genelist"

        newl = self.deepcopy() # full copy

        c1i = self._conditions.index(cond1)
        c2i = self._conditions.index(cond2)

        for item in newl:
            c1v = item["conditions"][c1i]+pad
            c2v = item["conditions"][c2i]+pad
            try:
                item[key] = self.__fold_change(item["conditions"][c1i]+pad, item["conditions"][c2i]+pad, log=log)
            except OverflowError as DivByZeroError:
                if c2v > c1v:
                    config.log.error("(%.2f/%.2f) failed" % (c2v, c1v))
                else:
                    config.log.error("(%.2f/%.2f) failed" % (c1v, c2v))
                raise Exception("add_fc_key: encountered an error, possibly the pad value '%s' is too small, or you need to log transform the data as you are getting an Infinity result" % pad)

            if and_err: # If it got here then the above try probably passed.
                expn_lo_err = item["conditions"][c2i] - item["err"][c2i] + pad
                expn_base = item["conditions"][c1i] + pad # I ignore any potential error here.
                expn_hi_err = item["conditions"][c2i] + item["err"][c2i] + pad
                up = abs(item[key] - (self.__fold_change(expn_base, expn_lo_err, log=log)))
                dn = abs(item[key] - self.__fold_change(expn_base, expn_hi_err, log=log))
                item[and_err] = min(up,dn)
                #print item["name"]
                #print item["conditions"]
                #print item["err"]
                #print expn_lo_err, expn_base, expn_hi_err, up, dn, item[key]

        newl._optimiseData()
        config.log.info("add_fc_key: Added fold-change key '%s'" % key)
        return(newl)

    def filter_by_fc(self, fckey=None, direction="any", value=2.0, **kargs):
        """
        **Purpose**
            Filter data which passes some sort of fold-change, defined in a previous key generated
            by add_fc_key()

        **Arguments**
            fckey (Required)
                The name of the fc_key to use, previously generated by add_fc_key()

            direction (Optional, default="any", values=["up", "down", "dn", "any"])
                The direction of change. "up" is +value, "down"/"dn" is -value and "any" is
                either +value or -value.

            value (Optional, default=2.0)
                The value of change required to pass the test.
                comparisons are evaluated as >= (i.e. greater than or equal to)

        **Returns**
            A new expression-object containing only the items that pass.
        """
        assert fckey, "filter_by_fc: 'fckey' argument is required"
        assert fckey in self.linearData[0], "filter_by_fc: '%s' not found in this expression object" % fckey
        assert direction in ("up", "down", "dn", "any"), "filter_by_fc(): direction argument '%s' not recognised" % direction

        new_expn = []

        for item in self.linearData:
            if direction == "up":
                if item[fckey] >= value:
                    new_expn.append(item)

            elif direction in ("down", "dn"):
                if item[fckey] <= -value:
                    new_expn.append(item)

            elif direction == "any":
                if item[fckey] >= value or item[fckey] <= -value:
                    new_expn.append(item)

        ret = expression(loadable_list=new_expn, cond_names=self._conditions)

        rep_d = {"up": "+", "dn": "-", "down": "-", "any": "+\-"}

        config.log.info("filter_by_fc: Filtered expression by fold-change '%s' %s%s, found: %s" % (fckey, rep_d[direction], value, len(ret)))
        return(ret)

    def filter_low_expressed(self, min_expression, number_of_conditions):
        """
        **Purpose**
            filter genes by a minimum_expression value in at least number_of_conditions

            Basically the command:

            sum([i > min_expression for i in <expression_data>]) >= number_of_conditions

        **Arguments**
            min_expression (Required)
                The minimum expression value required to pass the test

            number_of_conditions (Required)
                The number of conditions that must be greater than min_expression

        **Results**
            Returns a new genelist
        """

        newl = self.shallowcopy()
        newl.linearData = []

        for item in self.linearData:
            if sum([int(i > min_expression) for i in item["conditions"]]) >= number_of_conditions: # passed
                newl.linearData.append(item.copy())

        assert len(newl) > 0, "filter_low_expressed: The number of genes passing the filter was zero!"
        newl._optimiseData()
        config.log.info("filter_low_expression: removed %s items, list now %s items long" % (len(self) - len(newl), len(newl)))
        return(newl)

    def filter_by_mean_expression(self, min_mean_expression, max_mean_expression):
        """
        **Purpose**
            filter genes by a minimum_expression value in at least number_of_conditions

            Basically the command:

            mean(item["conditions"]) >= min_mean_expression
            and
            mean(item["conditions"]) <= max_mean_expression

        **Arguments**
            min_mean_expression (Required)
                The minimum mean expression value required to pass the test

            max_mean_expression (Required)
                The maximum mean expression value required to pass the test

        **Results**
            Returns a new genelist
        """

        newl = self.shallowcopy()
        newl.linearData = []

        for item in self.linearData:
            if numpy.mean(item["conditions"]) >= min_mean_expression and numpy.mean(item["conditions"]) <= max_mean_expression: # passed
                newl.linearData.append(item.copy())

        assert len(newl) > 0, "filter_by_mean_expression: The number of genes passing the filter was zero!"
        newl._optimiseData()
        config.log.info("filter_by_mean_expression: removed %s items, list now %s items long" % (len(self) - len(newl), len(newl)))
        return(newl)

    def filter_high_expressed(self, max_expression, number_of_conditions):
        """
        **Purpose**
            remove genes with a maximum expression value in at least number_of_conditions

            Basically the command:

            sum([int(i > max_expression) for i in item["conditions"]]) <= number_of_conditions

        **Arguments**
            max_expression (Required)
                The maximum expression value required to pass the test

            number_of_conditions (Required)
                The number of conditions that must be greater than min_expression

        **Results**
            Returns a new genelist
        """

        newl = self.shallowcopy()
        newl.linearData = []

        for item in self.linearData:
            if sum([int(i > max_expression) for i in item["conditions"]]) <= number_of_conditions: # passed
                newl.linearData.append(item.copy())

        assert len(newl) > 0, "filter_high_expressed: The number of genes passing the filter was zero!"
        newl._optimiseData()
        config.log.info("filter_high_expressed: removed %s items, list now %s items long" % (len(self) - len(newl), len(newl)))
        return(newl)

    def filter_by_expression(self, condition_name, minimum_expression, **kargs):
        """
        **Purpose**
            Keep only items in <condition_name> with >= minimum_expression

        **Arguments**
            condition_name (Required)
                The name of the condition to test

            minimum_expression (Required)
                the minimum expression value equal to or greater required

        **Returns**
            A new expression object with probes that fail to pass removed.
        """
        # kargs tidier to go here.

        newl = self.deepcopy()
        newl.linearData = []

        cindex = self._conditions.index(condition_name)

        for row in self.linearData:
            if row["conditions"][cindex] >= minimum_expression:
                newl.linearData.append(row)
        newl._optimiseData()
        config.log.info("filter_by_expression: %s items >= %s in '%s'" % (len(newl), minimum_expression, condition_name))
        return(newl)

    def filter_conditions_by_threshold(self, genes, key, threshold, **kargs):
        """
        **Purpose**
            filter out conditions based on a threshold test for
            a list of genes.

            One useage might be to remove all conditions not expressing some set of
            housekeeping genes:

            newe = expn.filter_conditions_by_threshold(["GAPDH", "ACTB", "EEF1A1"], "name", 1)

        **Arguments**
            genes (Required)
                a List of genes (or rows) to test,

            key (Required)
                key in the genelist to use to find the genes

            threshold (Required)
                all genes specified in genes must be > than this threshold.

        **Returns**
            A new expression object, with the conditions that failed removed.

        """
        cond_names = self.getConditionNames()

        cell_expressing = []
        genes_to_test = genes

        # I should get all the genes here, then just iterate through the genes, in case there are duplicate.s

        for g in genes_to_test:
            this_gene = e.getRowsByKey(key=key, values=g)[0] # what if more?
            for i, n in enumerate(cond_names):
                if this_gene["conditions"][i] >= threshold[g]:
                    cell_expressing.append(n)

        cell_expressing = list(set(cell_expressing))

        e = e.sliceConditions(cell_expressing)

        return(newe)

    def filter_by_value(self, value, absolute=False, **kargs):
        """
        **Purpose**
            Keep only items in <condition_name> with >= minimum_expression

        **Arguments**
            value (Required)
                Keep only rows with > value

            absolute (Optional, default=False)
                ignore the sign? (e.g. abs(value))

        **Returns**
            A new expression object with rows that fail to pass removed.
        """
        newl = self.deepcopy()
        newl.linearData = []

        for row in self.linearData:
            if absolute:
                vals = [abs(v) > value for v in row["conditions"]]
            else:
                vals = [v > value for v in row["conditions"]]
            if True in vals:
                newl.linearData.append(row)
        newl._optimiseData()
        config.log.info("filter_by_value: removed %s items, list now %s items long" % (len(self) - len(newl), len(newl)))
        return(newl)

    def filter_by_CV(self, minCV, maxCV, pad=0.1):
        """
        **Purpose**
            Filter each gene (row) to be within the range minCV, maxCV
            coefficient of variations. For good measure it will add in a
            'CV' key to the genelist

        **Arguments**
            minCV, maxCV (Required)
                The minimum and maximum CVs thresolds required

            pad (Optional, default=0.1)
                pad for the mean divisor to stop division by zero

        **Returns**
            A new expression object with rows that fail to pass removed.
        """
        newl = self.deepcopy()
        newl.linearData = []

        for row in self.linearData:
            CV = numpy.std(row['conditions']) / (numpy.mean(row['conditions'])+pad)
            if CV > minCV and CV < maxCV:
                row['CV'] = CV
                newl.linearData.append(row)
        newl._optimiseData()
        config.log.info("filter_by_CV: removed %s items, list now %s items long" % (len(self) - len(newl), len(newl)))
        return(newl)

    def bracket(self, min, max):
        """
        **Purpose**
            force the data to span only min and max:

            for each_value:
                if each_value > max:
                    each_value = max
                elif each_value <min:
                    each_value = min

        **Arguments**
            min, max (Required)
                minimum and maximum values

        **Returns**
            A new expn object.

        """
        newl = self.deepcopy()

        for item in newl:
            newc = []
            for c in item['conditions']:
                if c > max:
                    newc.append(max)
                elif c < min:
                    newc.append(min)
                else:
                    newc.append(c)
            item['conditions'] = newc
        newl._optimiseData()
        return(newl)

    def draw_cumulative_CV(self, filename, pad=0.1, label_genes=None, label_genes_key=None,
        label_fontsize=7, **kargs):
        '''
        **Purpose**
            Draw a cumulative CV plot for all genes (rows) in the expression object.

        **Arguments**
            filename (Required)
                The filename to save the image to.

            pad (Optional, default=0.1)
                pad for the mean divisor to stop division by zero

            label_genes (Optional, default=None)
                a list of names for rows (genes?) you want to mark on the plot

            label_genes_key (Optional, Required if label_genes is used)
                a key name to use to match the label_genes

            label_fontsize (Optional, default=7)
                Font size for the labels
        '''
        assert filename, "no filename specified"

        if label_genes_key:
            CV = [(row[label_genes_key], numpy.std(row['conditions']) / (numpy.mean(row['conditions'])+pad)) for row in self.linearData]
            CV = sorted(CV, key=itemgetter(1))
            labels = [i[0] for i in CV]
            CV = [i[1] for i in CV]
        else:
            CV = [numpy.std(row['conditions']) / (numpy.mean(row['conditions'])+pad) for row in self.linearData]
            CV.sort()

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        ax.plot(CV, color="red")
        for t in label_genes:
            try:
                x = labels.index(t)
                y = CV[x]
                ax.text(x, y, t, size=label_fontsize, ha='center', va='center') # Trust me, keep it va='center' it reduces ambiguity
            except ValueError:
                config.log.warning("label '%s' not found" % t)

        self.draw.do_common_args(ax, **kargs)
        realfilename = self.draw.savefigure(fig, filename)
        config.log.info("draw_cumulative_CV: Saved '%s'" % (realfilename))
        return(None)

    def draw_scatter_CV(self, filename, pad=0.1, label_genes=None, label_genes_key=None,
        label_fontsize=7, **kargs):
        '''
        **Purpose**
            Draw a scatter CV plot for all genes in the expression object.

            X axis is mean expression
            Y axis is CV

        **Arguments**
            filename (Required)
                The filename to save the image to.

            pad (Optional, default=0.1)
                pad for the mean divisor to stop division by zero

            label_genes (Optional, default=None)
                a list of names for rows (genes?) you want to mark on the plot

            label_genes_key (Optional, Required if label_genes is used)
                a key name to use to match the label_genes

            label_fontsize (Optional, default=7)
                Font size for the labels

            logx (Optional, default=2)
                You usually calculate CV on the raw (un transofrmed data)
                Generally it is easier to comprehend if the X axis (mean expression)
                is log2 transformed
        '''
        assert filename, "no filename specified"
        if 'logx' not in kargs:
            kargs['logx'] = 2 # Use the do_common_args system for clean overrides

        if label_genes_key:
            data = [(row[label_genes_key], numpy.std(row['conditions']) / (numpy.mean(row['conditions'])+pad), numpy.mean(row['conditions'])+pad) for row in self.linearData]
            labels = [i[0] for i in data]
            X = [i[2] for i in data]
            Y = [i[1] for i in data]
        else:
            data = [(numpy.std(row['conditions']) / (numpy.mean(row['conditions'])+pad), numpy.mean(row['conditions'])+pad) for row in self.linearData]
            X = [i[1] for i in data]
            Y = [i[0] for i in data]

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        ax.scatter(X, Y, color="grey", s=5, alpha=0.5, edgecolor='none')
        if label_genes:
            for t in label_genes:
                try:
                    i = labels.index(t)
                    x = X[i]
                    y = Y[i]
                    ax.text(x, y, t, size=label_fontsize, ha='center', va='center') # Trust me, keep it va='center' it reduces ambiguity
                except ValueError:
                    config.log.warning("label '%s' not found" % t)

        ax.set_xlabel('Mean')
        ax.set_ylabel('CV')

        self.draw.do_common_args(ax, **kargs)
        realfilename = self.draw.savefigure(fig, filename)
        config.log.info("draw_scatter_CV: Saved '%s'" % (realfilename))
        return(None)

    def scatter(self, x_condition_name, y_condition_name, filename=None, genelist=None, key=None,
        label=False, label_fontsize=12, **kargs):
        """
        **Purpose**
            draw an X/Y dot plot or scatter plot, get R^2 correlation etc.

        **Arguments**
            x_condition_name = the name of the er... X condition
            y_condition_name = the name of the er... Y condition

            genelist (Optional)
                If you send a genelist and a key then these items will be emphasised on the dotplot.

            key (Optional, Required if 'genelist' used)
                The key to match between the expression data and the genelist.

            label (Optional, default=False)
                If genelist and key are set, then you can label these spots with the label from
                genelist[key].
                If you want to label all of the items in the list, then just send the original genelist
                and a key to use and all of the spots will be labelled.

            label_fontsize (Optional, default=14)
                labels fontsize

            do_best_fit_line (Optional, default=False)
                draw a best fit line on the scatter

            print_correlation (Optional, default=None)
                You have to spectify the type of correlation to print on the graph.
                valid are:

                r = R (Correlation coefficient)
                r2 = R^2.
                pearson = Pearson

                You need to also set do_best_fit_line=True for this to work

            spot_size (Optional, default=5)
                The size of each dot.

            plot_diag_slope (Optional, default=False)
                Plot a diagonal line across the scatter plot.

            available key-word arguments:
            xlabel, ylabel, title, log (set this to the base to log the data by),
            xlims, ylims, spot_size,

        **Returns**
            the actual filename saved as and a new image in filename.
        """

        assert filename, "no filename specified"
        assert x_condition_name in self.serialisedArrayDataDict, "%s x-axis condition not found" % x_condition_name
        assert y_condition_name in self.serialisedArrayDataDict, "%s y-axis condition not found" % y_condition_name

        x_data = self.getDataForCondition(x_condition_name)
        y_data = self.getDataForCondition(y_condition_name)

        if not "xlabel" in kargs:
            kargs["xlabel"] = x_condition_name
        if not "ylabel" in kargs:
            kargs["ylabel"] = y_condition_name

        if genelist and key:
            matches = genelist.map(genelist=self, key=key) # make sure resulting object is array
            tx = matches.getDataForCondition(x_condition_name)
            ty = matches.getDataForCondition(y_condition_name)
            if label:
                kargs["spot_labels"] = matches[key]

            real_filename = self.draw.nice_scatter(x_data, y_data, filename, spots=(tx, ty), label_fontsize=label_fontsize, **kargs)
        else:
            real_filename = self.draw.nice_scatter(x_data, y_data, filename, **kargs)


        config.log.info("scatter: Saved '%s'" % real_filename)
        return(True)

    def boxplot(self, filename=None, showfliers=True, whis=1.5, **kargs):
        """
        **Purpose**

        Draw a boxplot of all conditions.

        **Arguments**

            filename (Required)
                filename to save as. The file extension may be modified
                depending the setting of the current config.DEFAULT_DRAWER

            log (True|False of 2..n for log2, log10)
                log the y axis (defaults to e)
                send an integer for the base, e.g. for log10

                log=10

                for log2

                log=2

                for mathematical constant e

                log=True
                log="e"

            showfliers (Optional, default=True)
                draw the 9/95 % outlier ticks on the plot

            whis (Optional, default=1.5)
                The location of the whiskers for the boxplot, see
                matplotlib for details of the settings

        **Results**

        saves an image with the correct filetype extension for the current
        config.DEFAULT_DRAWER.
        returns the actual filename used to save the file.
        """
        assert filename, "must provide a filename"

        data = self.serialisedArrayDataList

        if "log" in kargs and kargs["log"]:
            data = self.__log_transform_data(data, log=kargs["log"])

        # do plot
        actual_filename = self.draw.boxplot(data=data, filename=filename,
            labels=self.getConditionNames(), showfliers=showfliers, **kargs)

        config.log.info("boxplot: Saved %s" % actual_filename)
        return(actual_filename)

    def violinplot(self, filename=None, beans=False, **kargs):
        """
        **Purpose**

        Draw a beanplot/violinplot of all conditions.

        **Arguments**
            filename (Required)
                filename to save as. The file extension may be modified
                depending the setting of the current config.DEFAULT_DRAWER

            beans (Optional, default=False)
                Draw the beans on the beanplot, each dot corresponds to a sample.
                The algorithm for arranging the dits is a bit whacky, so may not look that nice.
                So by default it is off.

        **Results**
            saves an image with the correct filetype extension for the current
            config.DEFAULT_DRAWER
            returns the actual filename used to save the file.
        """
        assert filename, "must provide a filename"

        data = self.serialisedArrayDataDict

        # do plot
        actual_filename = self.draw.beanplot(data=data, filename=filename,
            order=self.getConditionNames(), beans=beans, **kargs)

        config.log.info("beanplot: Saved %s" % actual_filename)
        return(actual_filename)

    def multi_line(self, filename=None, alpha=None, **kargs):
        """
        **Purpose**
            Draw each gene as an alpha blended line. Ideal for visualising continuous data, for example time courses.

        **Arguments**
            filename (Required)
                The filename of the image to save.

            alpha (Optional, default=guess a value)
                The alpha blending value. By default glbase will attempt to guess a suitable
                alpha value. But this may be too light/dark for your tastes. You can set the value here, from
                0.0 (completely transparent) to 1.0 (completely opaque)

        **Return**
            None
        """
        assert filename, "You must specify a filename"

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        # Guess an alpha value:
        if not alpha:
            alpha = 8.0 / len(self)

        for item in self:
            ax.plot(item["conditions"], color="grey", alpha=alpha)

        self.draw.do_common_args(ax, **kargs)
        realfilename = self.draw.savefigure(fig, filename)
        config.log.info("multi_line: Saved '%s' with alpha=%.5f" % (realfilename, alpha))
        return(None)

    def log(self, base=math.e, pad=0.00001):
        """
        **Purpose**
            log transform the data

            NOTE: THis is one of the few IN-PLACE glbase commands.

        **Arguments**
            base (Optional, default=math.e)
                the base for the log transform.

            pad (Optional, default=1e-6)
                value to pad all values by to log(0) errors.

        **Returns**
            None
        """
        do_log = False

        if base == math.e:
            do_log = math.e
        elif isinstance(base, bool):
            do_log = math.e
        elif isinstance(base, int):
            do_log = base
        else:
            do_log = False

        for item in self:
            item["conditions"] = [math.log(v+pad, do_log) for v in item["conditions"]]
        self._optimiseData()
        return(None)

    def unlog(self, base=None, adjuster=0.00001):
        """
        **Purpose**
            return the raw data to the unlogged form. YOU MUST PROVIDE THE CORRECT BASE

            NOTE: THis is one of the few IN-PLACE glbase commands.

            Also note that continual repeated log() unlog() will cause the data to 'drift' away from its
            actual values.

        **Arguments**
            base (Required)
                the base for the log transform.

        **Returns**
            None
        """
        for item in self:
            item["conditions"] = [base**(v+adjuster) for v in item["conditions"]]
        self._optimiseData()
        return(None)

    def mult(self, number=None):
        """
        **Purpose**
            multiply the data by some number

            This method came about as I had a reason to multiply RNA-seq data to get it above zero
            Particularly concerning the TPM measure.

            NOTE: This is one of the few IN-PLACE glbase commands.

        **Arguments**
            number (Required)
                the number to multiply by

        **Returns**
            None
        """
        assert number, "must provide a number"

        for item in self.linearData:
            item["conditions"] = [v*number for v in item["conditions"]]
        self._optimiseData()
        return(None)

    def add(self, number=None):
        """
        **Purpose**
            add number to the expression data

            NOTE: This is one of the few IN-PLACE glbase commands.

        **Arguments**
            number (Required)
                the number to add to the expression data

        **Returns**
            None
        """
        assert number, "must provide a number"

        for item in self.linearData:
            item["conditions"] = [v+number for v in item["conditions"]]
        self._optimiseData()
        return(None)

    def abssub(self, number=None):
        """
        **Purpose**
            subtract number from all of the expression data, moving the value towards 0 with a limit of 0.

            This will perform the following logic::

                if expression_value > 0:
                    expression_value -= number
                    if expression_value < 0:
                        expression_value = 0

                elif expression_value < 0:
                    expression_value += number
                    if expression_value > 0:
                        expression_value = 0


            NOTE: This is an IN-PLACE glbase commands.

        **Arguments**
            number (Required)
                the number to add to the expression data

        **Returns**
            None
        """
        assert number, "must provide a number"

        for item in self.linearData:
            newcond = []
            for expression_value in item["conditions"]:
                if expression_value >= 0:
                    expression_value -= number
                    if expression_value < 0:
                        expression_value = 0

                elif expression_value <= 0:
                    expression_value += number
                    if expression_value > 0:
                        expression_value = 0
                newcond.append(expression_value)
            item["conditions"] = newcond
        self._optimiseData()
        return(None)

    def sub(self, number=None):
        """
        **Purpose**
            subtract number from all of the expression data.

            This is an IN-PLACE operation

        **Arguments**
            number (Required)
                the number to add to the expression data

        **Returns**
            None
        """
        assert number, "must provide a number"

        for item in self.linearData:
            newcond = [i-number for i in item["conditions"]]
            item["conditions"] = newcond
        self._optimiseData()
        return(None)

    def __log_transform_data(self, serialisedArrayDataList=None, log=math.e):
        """
        (Internal)

        transforms the data based on base
        helper for drawBoxPlot() and draw drawCurves()

        Zeros are trimmed from the data. This is important because now the output
        may not be synchronised to the input. Care should be taken with this transform.
        """
        assert serialisedArrayDataList, "[Internal] __log_transform_data() - no data provided"

        do_log = False

        if log == math.e:
            do_log = math.e
        elif isinstance(log, bool):
            do_log = math.e
        elif isinstance(log, int):
            do_log = log
        else:
            do_log = False

        if do_log:
            data = []
            for set in serialisedArrayDataList:
                row = []
                for index, item in enumerate(set):
                    if item != 0.0:
                        row.append(math.log(item, do_log))
                    else:
                        row.append(math.log(0.000001, do_log)) # Append a very small value
                data.append(row)
            return(data)
        else:
            return(serialisedArrayDataList)

    def drawCurves(self, filename=None, **kargs):
        """
        **Purpose**

        draw a bell-curve diagram of the expression-data expression.

        **Arguments**

            filename (Required)
                filename of the resulting

            window
                size of window for moving average

            modifier
                undocumented fudge for float based arrays.

            xlimits
                a tuple of the form: (minimum_x, maximum_y)

            log (True|False of 2..n for log2, log10)
                log the y axis (defaults to e)
                send an integer for the base, e.g. for log10

                log=10

                for log2

                log=2

                for mathematical constant e

                log=True
                log="e"

                Data points that are 0.0 will be trimmed from the data.
                This means the number of samples plotted may not
                be the same as the number of points in the microarray
                data.

            cumulative (True|False, default False)

                draw cumulative curves.

            verbose (True|False), default=False
                print out the means and standard deviations.

        **Result**

        saves an image to 'filename'
        returns the actual filename (the filename may be modified
        depending upon the current display driver)

        """
        assert filename, "no filename given"

        window_size = 200
        simple_args = ["filename", "window", "modifier"]
        xlimits = None
        do_log = False
        extra_args = {}
        data = self.serialisedArrayDataList

        for k in kargs:
            if k == "log":
                data = self.__log_transform_data(self.serialisedArrayDataList, log=kargs["log"])
            if k == "filename":
                filename = kargs[k]
            if k == "window":
                window_size = kargs[k]
            if k == "modifier": # undocumented fudge for flat arrays.
                binner_modifier = kargs[k]
            if k == "xlimits": # requires a tuple.
                xlimits = kargs[k]
            if k == "cumulative":
                extra_args["cumulative"] = kargs["cumulative"]

        # normalise data, for each condition.
        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        if "verbose" in kargs and kargs["verbose"]:
            print("name\tmean\tstd")

        for i, c in enumerate(self.getConditionNames()):
            ax.hist(data[i], bins=window_size, histtype="step", label=c, **extra_args)
            m = mean(data[i])
            d = std(data[i])
            ax.axvline(x=m, color="red")
            ax.axvline(x=m-d, color='grey', ls=":")
            ax.axvline(x=m+d, color='grey', ls=":")
            if "verbose" in kargs and kargs["verbose"]:
                print("%s\t%.2f\t%.2f" % (c, m, d))

        if xlimits:
            ax.xlim(xlimits)
        ax.legend()

        real_filename = self.draw.savefigure(fig, filename)
        config.log.info("curves: Saved '%s'" % filename)
        return(True)

    def getDataByCriteria(self, **kargs):
        """
        **Purpose**
            used to test for fold_up, sig_up, etc...

            function can be any helper_function, which uses data[] as it's set for each item.
            a set of already defined functions exist in flags.py
            function must accept these arguments: (data[], conditionNames)
            data is a dict of the form {"con_name1": value, "con_name2": value ...}
        """
        function = None
        normal = None
        for k in kargs:
            if k == "function":
                function = kargs[k]
            if k == "normal":
                normal = kargs[k]

        if not function:
            config.log.error("Criteria function unavailable.")
            return(False)

        newl = self.deepcopy()
        newl.linearData = []

        conNames = self.getConditionNames()
        for item in self:
            data = item["conditions"]
            # package as a dict:
            dd = {}
            for index, name in enumerate(conNames):
                dd[name] = data[index]
            if function(dd, conNames, normal, **kargs): # pass on other kargs
                newl.linearData.append(deepcopy(item))

        newl._optimiseData()
        return(newl)

    def _insertCondition(self, condition_name, condition_data, range_bind=None, **kargs):
        """
        (Internal)
        candidate for reveal?

        This should be done using the numpy array
        """
        self._conditions.append(condition_name)
        max_data = max(condition_data)
        min_data = min(condition_data)
        if len(condition_data) != len(self):
            config.log.error("Insertion of array data, wrongly sized")
            return(False)
        for index, item in enumerate(self):
            if range_bind:
                toAdd = (float((condition_data[index] - min_data)) / max_data)
                toAdd = (toAdd * (range_bind[0] + range_bind[1])) - range_bind[0]
            else:
                toAdd = condition_data[index]
            item["conditions"].append(toAdd)
        self._optimiseData()
        return(True)

    def cumulative_distributions(self, genelists=None, filename=None, key=None, **kargs):
        """
        **Purpose**

            draw a set of cumulative distributions, based on a selection
            of genelist-like objects that can be mapped to
            this expression object using 'key'

        **Arguments**

            genelists (Required)
                a list or other iterable of genelists

            filename (Required)
                the filename to save the image to.

            key (Required)
                the key to use to match the microarray to the genelist.

        **Returns**

            An image, saved to 'filename'
        """
        valig_args = ["genelists", "filename", "key"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError(self.cumulative_distributions, k)

        assert filename, "you must specify a valid filename"
        assert key, "you must specify a mapping key"
        assert genelists[0], "you must specify a valid list of genelists"
        assert key in genelists[0], "key '%s' not found in the genelists" % key
        assert key in self, "key '%s' not found in the expression-data" % key

        mapped_scores = []

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        axis = fig.add_subplot(111)

        for g in genelists:
            mapped = self.map(genelist=g, key=key)

            nmap = []
            # this will sum all items in array.
            for a in mapped:
                print(a["conditions"])
                s = sum(a["conditions"])
                nmap.append(s)

            # cumulate the nmap
            for i, v in enumerate(nmap):
                try:
                    nmap[i] = nmap[i] + nmap[i+1]
                except:
                    break

            print(nmap)

        # matplotlib junk is inappropriately here: to go later.

        axis.plot(nmap, label=g.name)

        axis.set_title("")
        #axis.show_legend()
        fig.savefig(filename)

    def drawBarChart(self, gene_symbols=None, filename=None, key=None, labels=None,
        errs_are_absolute=False, error_keys=None, fake_entries=True, **kargs):
        """
        **Purpose**

            Draw barcharts of the individual genes or a (more likely) list of genes.
            I will use these symbols and collect every matching entry in the expression-data
            and will draw a horizontal bar-chart

            You probably want barh_single_item() in preference over this (old) method

        **Arguments**

            gene_symbols (Required)
                a list or single item that I can find in the expression data.

            filename (Required)
                the filename to save the image to.

            key (Optional, default="name"|"symbol")
                The key to match the gene_symbols in.
                Ususally you will mean either name, gene symbol, or enst or equivalent
                annotation

            labels (Optional, default=<the same as 'key' argument, above>)
                specify an optional key to use for the actual labels on the chart
                If left unspecified then the 'key' argument is used for the labels.
                This allows you to search using 'key' but actually label the items with 'label'

            error_keys (Optional, default=None)
                If the array data has error values then you can send the key names.
                The input comes in two forms::

                    error_keys="stderr" # this will use one key that should match the
                                        # expression data

                    error_keys=["conf_lo", "conf_hi"] # If you have paired keys for
                    # error values. In this instace conf_lo holds the lower bound confidence
                    # interval and conf_hi the Higher bound.

            errs_are_absolute (Optional, default=False)
                If your error bars are absolute then set this to True.
                What it means is if your expression is (say) 10, and your error is 2.
                then the error bars will be drawn from 8 -> 12 (ie. the errors are relative
                to the expression value. If however you are using

            fake_entries (Optional, default=True)
                If the value you search for (usually a gene name) does not exist in the data
                then bar_chart will fake an empty entry, and the gene will still show up on the
                graph even though there is no data in your list. It will be given expression values of
                0 (actually, a very small float to avoid log errors).
                Set this to False if you don't want empty values to be placed on the bar_chart

            Other arguments understood by the figure
                title, xlabel, ylabel

                and 'cols' - where you can send a list of colours for the bars. The list
                must be as long as the data, or it will throw an error.

        **Example**

            Okay, this is pretty complicated (It's because the expression system has only a rudimentary
            understanding of error values, I may change this in the future)::

                x = expression(...)

                # I don't have any error values
                x.bar_chart(filename="meh1.png", gene_symbols=["Nanog", "Sox2", "Sox17", "Stat3"],
                    key="name")

                # My errors are confidence intervals relative to the expression value
                x.bar_chart(filename="meh2.png", gene_symbols=["Nanog", "Sox2", "Sox17", "Stat3"],
                    key="name", error_keys=["conf_lo", "conf_hi"], errors_are_absolute=True)

                # My errors are standard errors (or equivalent) and are symmetric around the mean
                x.bar_chart(filename="meh3.png", gene_symbols=["Nanog", "Sox2", "Sox17", "Stat3"],
                    key="name", error_keys="stderr")

        **Returns**

            A saved image in filename and None

        """
        assert filename, "no filename specified"
        assert gene_symbols, "no gene_symbols specified"
        if "cols" in kargs:
            assert len(kargs["cols"]) == len(self.linearData[0]["conditions"]), "the colour array is not the same length as the array data!"

        if not key:
            if "name" in self:
                key = "name"
            else:
                if not "symbol" in self:
                    key = "symbol"
                else:
                    raise AssertionError("No suitable key found")

        if not labels:
            labels = key

        if not isinstance(gene_symbols, list):
            gene_symbols = [gene_symbols]

        empty_array = [0.0001 for x in self._conditions]
        fake_entry = {"conditions": empty_array, "loc": location(loc="chr1:1000-1000")}
        # work out the error keys:
        if error_keys:
            if isinstance(error_keys, list): # load kargs for draw.bar_chart()
                kargs["err_up"] = error_keys[0]
                kargs["err_dn"] = error_keys[1]
                fake_entry = {"conditions": empty_array, error_keys[0]: empty_array, error_keys[1]: empty_array, "loc": location(loc="chr1:1000-1000")}
            else: # assume string
                kargs["err"] = error_keys
                fake_entry = {"conditions": empty_array, error_keys: empty_array, "loc": location(loc="chr1:1000-1000")}

        subset = []
        for g in gene_symbols:
            nl = self._findDataByKeyGreedy(key, g)
            if nl:
                subset = subset + nl
            else:
                new = fake_entry.copy()
                new.update({key: "%s (not detected)" % g, labels: "%s (not detected)" % g}) # I want to signify that these are empty
                subset.append(new)
                config.log.warning("%s not found. Spoofing an empty entry" % g)

        if subset:
            nl = genelist()
            nl.load_list(subset)

            # This is a kind of peculiar way of doing it...
            self.draw.bar_chart(filename=filename, genelist=nl, data="conditions",
                errs_are_absolute=errs_are_absolute,
                labels=labels, cond_names=self._conditions,
                **kargs)

    def hist(self, filename=None, range=None, suppress_zeros=True, log=None, kde=False,
        covariance=0.2, bins=50, color='grey', legend=True, **kargs):
        """
        **Purpose**
            Draw a normal histogram of the expression values

        **Arguments**
            filename (Required)
                the filename to save the resulting image to.

            covariance (Optional, default=0.2)
                undocumented wierdness (Sorry)

            range (Optional, default=(min(expression), max(expression)))
                This is the maximum value to build the KDE over. (ie. 0 ... mmax).
                You probably really want to chage this value, as small outliers will distort the
                distribution massively.
                Play around with the value until it gets you a nice normal-like distribution.

            suppress_zeros (Optional, default=False)
                microarray expression-data don't generally include zero values for expression.
                However, RNA-seq does include zero values for some transcripts.
                In some instances you may want to suppress these values.

                If this is set to True expression values of zero are not included in the
                histogram.

            log (Optional, default=False)
                log transform the data.
                At the moment only log2 is supported.

            color (Optional, default='grey')
                color in the histogram

            legend (Optional, default=True)
                Draw the legend.

            kde (Optional, default=False)
                use kernel density estimation to smooth the data

            covariance (Optional, default=0.2)
                Value for KDE smoothing.

            bins (Optional, default=50)
                number of bins to use. (Also affects the resolution of the kde smoothing)

        **Returns**

        """
        assert filename, "Need filename, leh!"

        binned_data = {}

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        all_data = {}

        for k in self._conditions:
            expn_values = self.getDataForCondition(k)
            if not range: # sample a range if none specified
                range = (min(expn_values), max(expn_values))

            if suppress_zeros:
                expn_values = [v for v in expn_values if int(v*1000000000) != 0]

            expn_values = numpy.array(expn_values)

            if log:
                expn_values = numpy.log2(expn_values)

            if kde:
                expn_values = utils.kde(expn_values, range=range, covariance=covariance, bins=bins)
                ax.plot(expn_values, label=k)
            else:
                ax.hist(expn_values, color=color, bins=bins, range=range, density=True, histtype='stepfilled', label=k)
            all_data[k] = expn_values

        if legend:
            ax.legend(ncol=1)

        self.draw.do_common_args(ax, **kargs)
        real_filename = self.draw.savefigure(fig, filename)

        config.log.info("hist: Saved '%s'" % real_filename)
        return({"data": all_data, "labels": self._conditions})

    def tree(self, mode="conditions", filename=None, row_name_key=None,
        cluster_mode="euclidean", color_threshold=None, label_size=6, cut=False,
        radial=False, optimal_ordering=True,
        **kargs):
        """
        **Purpose**
            Draw a hierarchical clustered tree of either the 'conditions' or 'rows'

        **Arguments**
            filename (Required)
                filename to save an image to.

            mode (Optional, default="conditions")
                cluster either the "conditions" or the "rows" or "genes". (rows are usually genes)

            row_name_key (Optional, default=None)
                A key to use for the row_names, if mode == "row"

            cluster_mode (Optional, default="euclidean")
                A metric to cluster the data by.

            color_threshold (Optional, default=None)
                By default tree() uses the Scipy/MATLAB default colours to signify
                links with similar distances. Set to -1 to change the tree to all blue.

                see scipy.cluster.hierarch.dendrogram

            label_size (Optional, default=7)
                The size of the text attached to the leaf labels.

            cut (Optional, default=False)
                the cut threshold on the tree, to select groups and return the cut.
                Will draw a line on the dendrgram tree indicating its position.
                The value should range from 0..1.0

                Note that if cut is used then row_name_key must also be valid

                Also note that the order of the clusters returned is random.

            bracket (Optional, default=False)
                The min and max to bracket the data. This is mainly here for compatibility with heatmap()
                So that the row/col trees in heatmap and here exactly match

            optimal_ordering (Optional, default=True)
                See https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html

        **Returns**
            if cut is False:
                The Tree in a dictionary {"dendrogram": ..., "linkage": ..., "distance": ...}
            if cut is True:
                A list of genelists, containing all of the items from the cut.

        """
        valid_modes = ("conditions", "rows", "genes")

        assert mode in valid_modes, "'%s' not a valid mode" % mode
        if cut:
            assert cut >= 0.0 and cut <= 1.0, "cut '%.2f' must be between 0 and 1" % cut
            assert row_name_key, "'row_name_key' must also be valid"

        if not "size" in kargs: # resize if not specified
            kargs["size"] = (3,6)

        fig = self.draw.getfigure(**kargs)

        data = numpy.array(self.serialisedArrayDataList)

        if "log" in kargs:
            data = self.__log_transform_data(self.serialisedArrayDataList, log=kargs["log"])

        if "bracket" in kargs: # done here so clustering is performed on bracketed data
            data = self.draw.bracket_data(data, kargs["bracket"][0], kargs["bracket"][1])
            vmin = kargs["bracket"][0]
            vmax = kargs["bracket"][1]

        if mode == "conditions": # Use the condition names for rows:
            row_names = self._conditions
        elif mode == "rows" or mode == "genes":
            data = data.T
            if row_name_key:
                row_names = self[row_name_key]
            else:
                row_names = None

        ax = fig.add_subplot(111)
        ax.set_position([0.01, 0.01, 0.4, 0.98])
        # from scipy;
        # generate the dendrogram
        dist = pdist(data, metric=cluster_mode)
        if cut or color_threshold: # override color_threshold with cut
            color_threshold = cut
            color_threshold = color_threshold*((dist.max()-dist.min())+dist.min()) # convert to local measure

        link = linkage(dist, 'complete', metric=cluster_mode, optimal_ordering=optimal_ordering)
        a = dendrogram(link, orientation='left', labels=row_names,
            ax=ax,
            color_threshold=color_threshold, get_leaves=True)

        ax.set_frame_on(False)
        ax.set_xticklabels("")
        ax.tick_params(top=False, bottom=False, left=False, right=False)

        if cut:
            # I believe that the old code is actually correct
            ret = []
            config.log.info("tree: Using local threshold '%.2f'" % color_threshold)
            clus = scipy.cluster.hierarchy.fcluster(link, color_threshold, 'distance')
            for i, net in enumerate(sorted(set(clus))):
                # get all of the gene names back out.
                idxs = [idx for idx, x in enumerate(clus) if x == net]
                newl = [{row_name_key: row_names[idx]} for idx in idxs]
                newgl = genelist()
                newgl.load_list(newl)
                ret.append(newgl)

            '''
            ret = []
            color_threshold
            config.log.info("tree: Using local cut threshold '%.3f'" % color_threshold)
            clus = scipy.cluster.hierarchy.fcluster(link, color_threshold, 'distance')

            print clus

            # Complicated way to do, but required to maintain the set order:
            uniq_clus_in_order_bottom_to_top = []
            for i in clus:
                if i not in uniq_clus_in_order_bottom_to_top:
                    uniq_clus_in_order_bottom_to_top.append(i)

            for net in uniq_clus_in_order_bottom_to_top:
                # get all of the gene names back out.
                idxs = [idx for idx, x in enumerate(clus) if x == net]
                newl = [{row_name_key: row_names[idx]} for idx in idxs]
                newgl = genelist()
                newgl.load_list(newl)
                ret.append(newgl)
            '''

            ax.axvline(color_threshold, color="grey", ls=":")
            config.log.info("tree: Found %s clusters" % len(ret))

        # Use the tree to reorder the data.
        row_names = a["ivl"]
        [t.set_fontsize(label_size) for t in ax.get_yticklabels()]

        #self.draw.do_common_args(ax, **kargs) # broken!?
        if filename:
            real_filename = self.draw.savefigure(fig, filename)
            config.log.info("tree: Saved '%s'" % real_filename)

        if not cut:
            ret = {"dendrogram": a, "linkage": link, "distance": dist}
        return(ret)

    def radial_tree(self, mode="conditions", filename=None, row_name_key=None,
        cluster_mode="euclidean", color_threshold=None, label_size=7, cut=False,
        radial=False,
        **kargs):
        """
        **Purpose**
            Draw a hierarchical clustered tree of either the 'conditions' or 'rows'

            In this version the output is a radial tree, with the root at the center.

        **Arguments**
            filename (Required)
                filename to save an image to.

            mode (Optional, default="conditions")
                cluster either the "conditions" or the "rows" or "genes". (rows are usually genes)

            row_name_key (Optional, default=None)
                A key to use for the row_names, if mode == "row"

            cluster_mode (Optional, default="euclidean")
                A metric to cluster the data by.

            color_threshold (Optional, default=None)
                By default tree() uses the Scipy/MATLAB default colours to signify
                links with similar distances. Set to -1 to change the tree to all blue.

                see scipy.cluster.hierarch.dendrogram

            label_size (Optional, default=7)
                The size of the text attached to the leaf labels.

            cut (Optional, default=False)
                the cut threshold on the tree, to select groups and return the cut.
                Will draw a line on the dendrgram tree indicating its position.
                The value should range from 0..1.0

                Note that if cut is used then row_name_key must also be valid

                Also note that the order of the clusters returned is random.

            bracket (Optional, default=False)
                The min and max to bracket the data. This is mainly here for compatibility with heatmap()
                So that the row/col trees in heatmap and here exactly match

        **Returns**
            if cut is False:
                The Tree in a dictionary {"dendrogram": ..., "linkage": ..., "distance": ...}
            if cut is True:
                A dict of genelists, containing all of the items from the cut.

        """
        config.log.warning('radial_tree: is currently experimental!')

        from . import radial_tree

        valid_modes = ("conditions", "rows", "genes")
        assert mode in valid_modes, "'%s' not a valid mode" % mode
        if cut:
            assert cut >= 0.0 and cut <= 1.0, "cut '%.2f' must be between 0 and 1.0" % cut
            assert row_name_key, "'row_name_key' must also be valid"

        if not "size" in kargs: # resize if not specified
            kargs["size"] = (6,6)

        fig = self.draw.getfigure(**kargs)
        data = numpy.array(self.serialisedArrayDataList) # get a copy

        if "log" in kargs:
            data = self.__log_transform_data(self.serialisedArrayDataList, log=kargs["log"])

        if "bracket" in kargs: # done here so clustering is performed on bracketed data
            # And will result in the same clustering as in heatmaps.
            data = self.draw.bracket_data(data, kargs["bracket"][0], kargs["bracket"][1])
            vmin = kargs["bracket"][0]
            vmax = kargs["bracket"][1]

        if mode == "conditions": # Use the condition names for rows:
            row_names = self._conditions
        elif mode == "rows" or mode == "genes":
            data = data.T
            if row_name_key:
                row_names = self[row_name_key]
            else:
                row_names = None


        # from scipy;
        # generate the dendrogram
        dist = pdist(data, metric=cluster_mode)
        if cut or color_threshold: # override color_threshold with cut
            color_threshold = cut
            color_threshold = color_threshold*((dist.max()-dist.min())+dist.min()) # convert to local measure

        link = linkage(dist, 'complete', metric=cluster_mode)

        rtree = radial_tree.tree(link)

        ax = fig.add_subplot(111)
        rtree = radial_tree.tree(link)
        rtree.draw()

        #ax.set_frame_on(False)
        #ax.set_xticklabels("")
        #ax.tick_params(top="off", bottom="off", left="off", right="off")

        if cut:
            # Probably too tough to support?
            # I believe that the old code is actually correct
            ret = []
            config.log.info("radial_tree: Using local threshold '%.2f'" % color_threshold)
            clus = scipy.cluster.hierarchy.fcluster(link, color_threshold, 'distance')
            for i, net in enumerate(sorted(set(clus))):
                # get all of the gene names back out.
                idxs = [idx for idx, x in enumerate(clus) if x == net]
                newl = [{row_name_key: row_names[idx]} for idx in idxs]
                newgl = genelist()
                newgl.load_list(newl)
                ret.append(newgl)

            ax.axvline(color_threshold, color="grey", ls=":")
            config.log.info("radial_tree: Found %s clusters" % len(ret))

        # Use the tree to reorder the data.
        #row_names = a["ivl"]
        #[t.set_fontsize(label_size) for t in ax.get_yticklabels()]

        #self.draw.do_common_args(ax, **kargs) # broken!?
        if filename:
            real_filename = self.draw.savefigure(fig, filename)
            config.log.info("radial_tree: Saved '%s'" % real_filename)

        if not cut:
            ret = {"dendrogram": a, "linkage": link, "distance": dist}
        return(ret)

    def gene_curve(self, key=None, values=None, filename=None, moving_average=None,
        **kargs):
        """
        **Purpose**
            Draw all of the genes specified on a single axis.

        **Arguments**
            key (Required)
                The key to use to match to items found in values

            values (Required)
                The values to match, should be a list or other iterable

            filename (Required)
                the filename to save to

            moving_average (Optional, default=False)
                use a moving average to smooth the gene curves.
                If set then should be an integer specifying the number of neighboring bins to use

        **Returns**
            A file in filename and None
        """
        assert key, "Must specify key"
        assert values, "must specify values"
        assert filename, "must specify filename to save to"

        if not "aspect" in kargs:
            kargs["aspect"] = "wide"

        keeps = frozenset(values)

        x_data = numpy.arange(len(self._conditions))

        data = {}
        for i in self.linearData:
            if i[key] in keeps:
                data[i[key]] = i["conditions"]

        fig = self.draw.getfigure(**kargs)

        ax = fig.add_subplot(111)

        for k in data:
            if moving_average:
                x_data, y_data = utils.movingAverage(data[k], window=moving_average)
                ax.plot(x_data, y_data, label=k)
            else:
                ax.plot(x_data, data[k], label=k)

        ax.legend()
        self.draw.do_common_args(ax, **kargs)
        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("gene_curve: Saved '%s'" % actual_filename)

    def correlation_heatmap(self, axis="conditions", filename=None, label_key=None, mode="r2", aspect="square", bracket=(0,1),
        optimal_ordering=True, **kargs):
        """
        **Purpose**
            Plot a heatmap of the (R, R^2, Pearson or Spearman) correlation for all pairs of
            samples in the expression object.

            Note that R and R^2 the calculation is very fast but pearson and spearman are quite a bit slower.

        **Arguments**
            filename (Required)
                the filename to save the image to.

            axis (Optional, default="conditions")
                perform the clustering based on conditions or rows, genes (actually just rows=genes)

                Valid axes are:
                    'conditions' # the default
                    'rows' or 'genes' # cluster by rows (or genes for typical expression data)

            label_key (Optional, but required if axis=genes or rows)
                You must specify a label key if the axis is rows or genes

            mode (Optional, default="r2")
                by default the R (Coefficient of determination) is squared. Set to 'r' for
                Coefficient of determination value.
                use 'pearson' for a Pearson correlation score and 'spearman' for a
                Spearman-ranked score.

            bracket (optional, default=(0,1))
                bracket the heatmap by these values.

            Other heatmap options that will work:
                heat_hei
                heat_wid
                ...

        **Returns**
            A dict containing::
                "data": 2D numpy array containing the correlation scores (depending upon the mode)
                "labels": the labels along the top and bottom of the array, sorted according to the
                clustering (if any)

            Note that the arrays in 'data' are not clustered and will not be in the same order as the figure.
        """
        assert filename, "You must specify a filename"
        assert axis in ("conditions", "genes", "rows"), "'%s' axis not found" % axis
        assert mode in ("r", "r2", "pearson", "spearman"), "'%s' mode not found" % mode

        if axis == "conditions":
            data_table = self.getExpressionTable().T
            labels = self._conditions
        elif axis == "rows" or axis == "genes":
            assert label_key, "You must specify a label_key if axis = rows/genes"
            data_table = self.getExpressionTable()
            labels = self[label_key]

        if mode in ("r", "r2"):
            arr = numpy.corrcoef(data_table) # PMCC or little r,
            if mode == "r2":
                arr *= arr

        else:
            arr = numpy.zeros((len(labels), len(labels)))
            ps = numpy.zeros((len(labels), len(labels)))
            p = progressbar(len(labels))
            for ic1, c1 in enumerate(labels):
                for ic2, c2 in enumerate(labels):
                    if ic1 == ic2: # stop identical labels mucking it up.
                        arr[ic1, ic2] = 1.0
                    else:
                        x = data_table[ic1,:]
                        y = data_table[ic2,:]

                        if mode == "pearson":
                            arr[ic1,ic2], ps[ic1, ic2] = scipy.stats.pearsonr(x, y)
                        elif mode == "spearman":
                            arr[ic1,ic2], ps[ic1, ic2] = scipy.stats.spearmanr(x, y)
                p.update(ic1)

        square = True
        if "heat_hei" in kargs or "heat_wid" in kargs:
            square=False

        results = self.draw.heatmap(filename=filename, data=arr, square=square,
            bracket=bracket, aspect=aspect, row_names=labels, col_names=labels,
            colbar_label="Correlation (%s)" % mode, optimal_ordering=optimal_ordering, **kargs)
        config.log.info("correlation_heatmap: Saved '%s'" % results["real_filename"])
        return({"data": results["reordered_data"], "labels": results["reordered_cols"]})

    def closest_correlate(self, target_condition, number=5, cut_off=0.7, method="r2",
        pretty_print=False, **kargs):
        """

        **Purpose**
            Find the closest correlating samples to <target_condition>

        **Arguments**
            target_condition (Required)
                The name of the target condition to find the closest correlates for

            number (Optional, default=5)
                report <number> closest correlates

            cut_off (Optional, default=0.7)
                Do not report a correlate if it falls below the <cut_off> value

            method (Optional, default="r2")
                The method to use for the correlation, must be one of:
                "r", "r2", "spearman", "pearson"

            pretty_print (Optional, default=False)
                By default closest_correlate dumps it all in a fiddly dict.
                With this enabled instead closest_correlate() will print it to the screen
                nicely formatted.

        **Returns**
            A list of dictionaries containing the closest correlations, their rank
            and their correlation score

        """
        assert target_condition in self._conditions, "closest_correlate: '%s' condition was not found" % target_condition
        assert method in ("r", "r2", "pearson", "spearman"), "closest_correlate: method '%s' not recognised" % method

        # get all the correlations:

        res = []
        x = self.getDataForCondition(target_condition)

        for ic1, c1 in enumerate(self._conditions):
            if c1 != target_condition:
                y = self.getDataForCondition(c1)
                if method in ("r", "r2"):
                    arbr = scipy.polyfit(x, y, 1)
                    xr = scipy.polyval(arbr, x)
                    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x,y)
                    if method == "r":
                        correlation = r_value
                    elif method == "r2":
                        correlation = r_value**2

                elif method == "pearson":
                    correlation, p_value = scipy.stats.pearsonr(x, y)
                elif method == "spearman":
                    correlation, p_value = scipy.stats.spearmanr(x, y)

                if correlation > cut_off:
                    res.append({"name": c1, "correlation": correlation})

        # sort the dict by correlation:
        res = sorted(res, key=itemgetter("correlation"))
        res.reverse()

        if pretty_print:
            for rank, item in enumerate(res):
                print("%s:\t%s (%.3f)" % (rank+1, item["name"], item["correlation"]))

        return(res)

    def cut(self, function):
        """
        **Purpose**
            Perform an action (any function(), often a lambda) on each item
            and return the resulting expression object.

            This function does (essentially) this::

            newlist = []
            for item in self:
                if function(item):
                    newlist.append(item)
            return(newlist)

        **Arguments**
            function (Required)
                A Python method or lambda, that takes one argument, which is each row
                of the item. Must return True or False

        **Returns**
            A new expression object.
        """

        assert function, "cut: you must supply a falid function"
        assert function({"name": None}), "cut: function appears to be invalid"

        newl = []
        for item in self.linearData:
            #print item
            #print function(item)
            if function(item):
                newl.append(item)

        if not newl:
            config.log.info("cut: result of cut() is empty!")
            return(None)

        cc = self.shallowcopy() # going to replace linearData.
        cc.linearData = newl
        cc._optimiseData()

        return(cc)

    def barh_single_item(self, key=None, value=None, filename=None, tree=None,
        plot_mean=True, plot_stdev=False, fold_change=False, tight_layout=False,
        bar_cols=None, vert_space=0.75, hori_space=0.5, **kargs):
        """
        **Purpose**
            Plot a horizontal bar for a single item (gene?), for all conditions.

            Uses an 'err' key for error bars if present

            Note that this will find only the first item in the list that matches.

        **Arguments**
            key (Required)
                The key to search for the item

            value (Required)
                the value to search for

            filename (Required)
                the filename to save the plot to

            tree (Optional, default=None)
                A tree used to order the conditions by. Should be the output from tree()
                and should only be clustering the conditions

            plot_mean (Optional, default=True)
                plot a grey line showing the mean of all of the expression values

            plot_stdev (Optional, default=False)
                plot a blue line at 1x stdev and a red line at 2x stdev

            bar_cols (Optional, default=None)
                a list of colours to use to colour the bars

            fold_change (Optional, default=False)
                by default barh_single_itme expects absolute levels of expression, but if you want to
                provide fold-change data (which centres around 0 then setting this to True will
                modify the plot and make it more friendly for fold-change plost
                (re-bracket, a line at 0, and 2 and 4 fold-change (assuming data is
                log2(fold-change)), plot_mean and plot_stdev=False).

            Typical arguments for plots are supported:
                xlabel - x-axis label
                ylabel - y-axis label
                title  - title
                xlims - x axis limits (Note that barh_single_item will clamp on [0, max(data)+10%] by default.
                ylims - y-axis limits
                xticklabel_fontsize - x tick labels fontsizes
                yticklabel_fontsize - y tick labels fontsizes

                hori_space (default=0.5)
                vert_space (default=0.75) - a special arg to help pad the barchart up when using very long plots with a ot of samples

            tight_layout (Optional, default=False)
                wether to use matplotlib tight_layout() on the plot

        **Returns**
            The item used to draw the barchart
        """
        assert key, "barh_single_item: you must specify a key"
        assert value, "barh_single_item: you must specify a value"
        assert filename, "barh_single_item: you must specify a filename"

        if tree and bar_cols:
            config.log.warning('Using tree and bar_cols, make sure your bar_cols are in the same order as the tree!')

        item = self._findDataByKeyLazy(key, value)

        if not item:
            config.log.warning("barh_single_item: '%s:%s' not found in this list, not saving" % (key, value))
            return(None)

        if not bar_cols:
            bar_cols = 'grey'

        err = None
        # re-order "conditions" by tree, if present:
        if tree:
            data = []
            conds = []
            if "err" in item:
                err = []
            if "dendrogram" in tree:
                order = tree["dendrogram"]["ivl"]
            else:
                # This may die if not sending the correct data
                order = tree["ivl"]

            # resort the data by order;
            for label in order:
                indx = self._conditions.index(label)
                data.append(item["conditions"][indx])
                conds.append(label)
                if "err" in item:
                    err.append(item["err"][indx])

        else:
            data = list(item["conditions"])
            conds = list(self._conditions)
            if "err" in item:
                err = list(item["err"])

        fig = self.draw.getfigure(aspect="long", **kargs)
        ax = fig.add_subplot(111)

        # Work out vert spacing
        hori_remaining = 1.0-hori_space-0.03
        vert_remaining = 1.0-vert_space
        ax.set_position([hori_remaining, vert_remaining/2, hori_space, vert_space]) # x,y,width,height

        y = numpy.arange(len(data))
        if err:
            # error bars'd graphs look better with black on grey
            ax.barh(y, data, xerr=err, ec="none", color=bar_cols, height=0.5, error_kw={'linewidth': 0.3, 'capthick': 0.3}, capsize=2.0)
        else:
            # no error bars, solid black is better
            ax.barh(y, data, ec="none", color=bar_cols, height=0.5)

        ax.set_yticklabels(conds)
        ax.set_yticks(y)#+0.25) # 0.25 should be half of the height of the bars so that text aligns with the bar
        ax.set_ylim([-0.5, len(data)-0.5])
        [item.set_markeredgewidth(0.0) for item in ax.yaxis.get_ticklines()]

        if fold_change:
            ax.axvline(-1, color="red", ls=":")
            ax.axvline(1, color="red", ls=":")
            ax.axvline(-2, color="orange", ls=":")
            ax.axvline(2, color="orange", ls=":")
            ax.axvline(-3, color="yellow", ls=":")
            ax.axvline(3, color="yellow", ls=":")
        else:

            if plot_mean:
                ax.axvline(numpy.average(data), color="grey", ls=":", lw=0.5)
                #ax.text(numpy.average(data)+0.1, len(data)+0.02, "m", size=6)

            if plot_stdev:
                mean = numpy.average(data)
                std = numpy.std(data)
                ax.axvline(mean + std, color="blue", ls=":", lw=0.5)
                ax.axvline(mean - std, color="blue", ls=":", lw=0.5)
                ax.axvline(mean + (std*2.0), color="red", ls=":", lw=0.5)
                #ax.text(mean + std+0.1, len(data)+0.02, "1x", size=6)
                #ax.text(mean + (std*2.0)+0.1, len(data)+0.02, "2x", size=6)

        if not "xlims" in kargs:
            if fold_change:
                mx = max(data)+(max(data)*0.1)
                if mx < 2.0:
                    mx = 2.0
                ax.set_xlim([-mx, mx])
            else:
                ax.set_xlim([0, max(data)+(max(data)*0.1)])

        if fold_change:
            ax.axvline(0, color="grey", ls="-")

        if not "yticklabel_fontsize" in kargs:
            kargs["yticklabel_fontsize"] = 10

        if not "title" in kargs:
            kargs["title"] = value

        self.draw.do_common_args(ax, **kargs)
        if tight_layout:
            fig.tight_layout()
        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("barh_single_item: Saved '%s'" % actual_filename)
        return(item)

    def time_course_plot(self, key=None, value=None, filename=None, timepoints=None, plotmean=True, **kargs):
        """
        **Purpose**
            Plot a nice graph from timecourse expression data (i.e. a line graph)

        **Arguments**
            key (Required)
                The key to search for the item

            value (Required)
                the value to search for

            filename (Required)
                The filename to save to

            timepoints (Optional)
                integer time points for your time course. If no timepoints are given glbase
                assumes your time points are equally spaced and in the same order from
                getConditionNames()

                however, if you supply a series of integers (or floats) then the plot will be arranged
                so that the x-axis

            plotmean (Optional, default=True)
                plot a line on the graph indicating the mean expression across the time course.

            Other standard draw() arguments are also accepted:
                xlabel - x-axis label
                ylabel - y-axis label
                title  - title
                xlims - x axis limits
                ylims - y-axis limits
                logx - set the x scale to a log scale argument should equal the base
                logy - set the y scale to a log scale
                legend_size - size of the legend, small, normal, medium
                xticklabel_fontsize - x tick labels fontsizes
                yticklabel_fontsize - y tick labels fontsizes
                vlines - A list of X points to draw a vertical line at
                hlines - A list of Y points to draw a vertical line at

        """
        assert filename, "time_course_plot: You must provide a filename"

        if not timepoints:
            timepoints = numpy.arange(len(self.getConditionNames())) + 1

        assert key, "time_course_plot: you must specify a key"
        assert value, "time_course_plot: you must specify a value"
        assert filename, "time_course_plot:you must specify a filename"

        data = self._findDataByKeyLazy(key, value)

        if not data:
            config.log.warning("time_course_plot: '%s:%s' not found in this list, not saving" % (key, value))
            return(None)

        data = list(data["conditions"])

        if not "aspect" in kargs:
            kargs["aspect"] = "wide"
        if not "size" in kargs:
            kargs["size"] = "small"

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        ax.plot(timepoints, data, '-o', color="black", mec="black", mfc="none")

        if plotmean:
            m = numpy.average(data)
            ax.axhline(m, ls=":", color="grey")

        if not "ylims" in kargs:
            kargs["ylims"] = (0, max(data)+(max(data)/10.0))

        ax.set_title("%s:%s" % (key, value))

        self.draw.do_common_args(ax, **kargs)
        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("time_course_plot: Saved '%s'" % actual_filename)
        return(actual_filename)

    def fc_scatter(self, cond1, cond2, filename=None, plot_diagonals=True, zoom_bracket=[-3,3],
        label_key="name", text_threshold=1.0, alpha=0.2, **kargs):
        """
        **Purpose**
            Draw a scatter plot of the fold-changes between two conditions

            There is not so much difference between this and scatter()
            The major advantage is it is set up to label particular genes, genes moving away
            from the diagonals and has more sensible defaults for a scatter plot based on fold-change.

            You could certainly emulate this using the normal scatter(). But think of this as a nice helper
            function.

            NOTE: This does NOT fold-change transform your data and assumes
            you have already done this. It also does not check if the data is fold-change data.

            use something like norm_multi_fc() or normaliseToCondition() before using this function.

            By default it will plot two plots, one zoomed in and one zoomed out to show all of the data

        **Arguments**
            cond1 (Required)
                the name of condition1

            cond2 (Required)
                the name of condition2

            filename (Required)
                the filename to save the scatter plot to

            plot_diagonals (Optional, default=True)
                plot guidelines showing the diagonals

            zoom_bracket (Optional, default=-[3,3])
                the bracket to use for the zoomed in plot

            label_key (Optional, default="name")

            text_threshold (Optional, default=1.0)
                Only draw the gene label if > abs(text_threshold) in any direction

            alpha (Optional, default=0.2)
                Blening fraction for the alpha (opacity) channel for the dots.

        **Returns**
            None
        """
        assert filename, "fc_scatter: you must specify a filename"

        # get the data
        c1d = self.getDataForCondition(cond1)
        c2d = self.getDataForCondition(cond2)
        names = self[label_key]

        # load the data for plotting
        pt_x = []
        pt_y = []
        labs = []
        cols = []
        for index, x in enumerate(c1d):
            pt_x.append(c1d[index])
            pt_y.append(c2d[index])
            cols.append("blue")

            # do the colouring here if away from the diagonal
            if label_key:
                labs.append(names[index])
            """
                    # get the distance from the diagonal:
            dx = (ex[l1] - ex[l2]) # i.e. x - y
            dy = (ex[l2] - ex[l1]) # i.e. y - x

            d = dx * dx # actually the squared distance
            print dx, dy, d
            if d > 0.05 or d < -0.05:
                cols.append("red")
                labs.append(g["name"])
                print "%s\t%.2f" % (g["name"], d)
                both.append({"ensg": g["ensg"]})
            else:
                cols.append("blue")
                labs.append("")
            """

        # plot
        if not "size" in kargs:
            kargs["size"] = (13,6)
        fig = self.draw.getfigure(**kargs)

        ax = fig.add_subplot(121)
        ax.scatter(pt_x, pt_y, alpha=alpha, edgecolor='none', color=cols)
        for i in range(len(labs)):
            if abs(pt_x[i]) > text_threshold or abs(pt_y[i]) > text_threshold:
                ax.text(pt_x[i], pt_y[i], labs[i], size=5, ha="center")

        # Diagonal slopes:
        ax.plot([5, -5], [5, -5], ":", color="grey")
        ax.plot([-5, 5], [5,-5], ":", color="grey")

        ax.axvline(0, color="grey", ls="-.")
        ax.axhline(0, color="grey", ls="-.")
        ax.set_xlim(zoom_bracket)
        ax.set_ylim(zoom_bracket)
        ax.set_xlabel(cond1)
        ax.set_ylabel(cond2)

        ax = fig.add_subplot(122)
        ax.scatter(pt_x, pt_y, alpha=alpha, edgecolor='none', color=cols)
        #for i in xrange(len(labs)):
        #    ax.text(pt_x[i], pt_y[i], labs[i], size=6, ha="center")

        # Diagonal slopes:
        if plot_diagonals:
            ax.plot([5, -5], [5, -5], ":", color="grey")
            ax.plot([-5, 5], [5,-5], ":", color="grey")

        minmax = max([abs(min(c1d)), max(c1d), abs(min(c2d)), max(c2d)])

        ax.axvline(0, color="grey", ls=":")
        ax.axhline(0, color="grey", ls=":")
        ax.set_xlim([-minmax,minmax])
        ax.set_ylim([-minmax,minmax])
        ax.set_xlabel(cond1)
        ax.set_ylabel(cond2)

        # best fit:
        """
        (ar, br) = polyfit(pt_x, pt_y, 1)
        print ar, br
        #xr = polyval([ar,br], pt_x)
        slope, intercept, r_value, p_value, std_err = linregress(pt_x,pt_y)

        print "m, x, r2:", slope, intercept, r_value*r_value

        mx = [min(pt_x), max(pt_x)]
        my = [(slope * min(pt_x)) + intercept, (slope * max(pt_x)) + intercept]

        ax.plot(mx, my, "r.-")
        """

        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("fc_scatter: Saved '%s'" % actual_filename)
        return(None)

    def kmeans(self, filename=None, key=None, seeds=None, dowhiten=True, plotseeds=True, **kargs):
        """
        **Purpose**
            perform k-means clustering on the expression data.

            There are two modes to k-means: Either randomly gather
            patterns from the data or the user can provide the initial patterns to
            cluster around.

            Output a stack of plots and return the groups

        **Arguments**
            filename (Optional)
                filename to save the profiles to

            key (Optional)
                A key to use to extract the initial values from

            seeds (Optional)
                If you send a number:
                    kmeans() will select 0..seeds randomly from the data
                If seeds is a genelist
                    For each item in the genelist, k-means cluster around that
                    value.

            dowhiten (Optional, default=True)
                'whiten' the data, i.e. standardise the variation.
                This can be a bad idea if your data is already heavily normalised (i.e.
                variation has already been controlled).

            plotseeds (Optional, default=False)
                don't plot the seeds

        **Returns**
            A dictionary containing keys either from the names of the seeds or the random
            group from the seed
        """

        # expression table needs to be transformed.
        data = numpy.array(self.serialisedArrayDataList).T

        if dowhiten:
            config.log.info("kmeans: whiten...")
            wt = whiten(data) # features are columns
        else:
            wt = data

        if isinstance(seeds, int):
            seeds = seeds
            seed_names = None
        else: # assume genelist
            assert key, "You must specify a key"
            cents = []
            seed_names = []
            for i in seeds:
                col = self.index(key, i[key])# I need to know the column numbers for each item in seeds:
                cents.append(self[col]["conditions"])
                seed_names.append(i[key])
            seeds = numpy.array(cents)

        config.log.info("kmeans: kmeans...")
        centroids, variance = kmeans(wt, seeds)


        config.log.info("kmeans: vq...")
        code, distance = vq(wt, centroids)

        clusters = {}
        for feature, cluster in enumerate(code):
            if not cluster in clusters:
                clusters[cluster] = []
            clusters[cluster].append(data[feature])

        # Guess a suitable arrangement
        sq = math.ceil(math.sqrt(len(clusters)))
        if not "size" in kargs:
            kargs["size"] = (sq*2, sq*2)

        fig = self.draw.getfigure(**kargs)

        fig.subplots_adjust(0.01, 0.01, 0.98, 0.98, wspace=0.15, hspace=0.15)

        for i, k in enumerate(clusters):
            ax = fig.add_subplot(sq, sq, i+1)
            ta = numpy.array(clusters[k])
            ax.plot(ta.T, alpha=0.1, color="grey")

            if seed_names:
                # plot he original centroid:
                if plotseeds:
                    ax.plot(cents[i], color="black")
                ax.set_title("%s (%s)" % (seed_names[i], len(clusters[k])), size=6)
            else:
                ax.set_title("Cluster %s (%s)" % (k+1, len(clusters[k])), size=6)

            ax.set_xlim([-1,len(self._conditions)])
            #ax.set_xticklabels("")
            #ax.set_yticklabels("")

        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("kmeans: Saved '%s'" % actual_filename)

    def bundle(self, filename=None, key=None, seeds=None, plotseeds=True, Z=0.9, mode="pearsonr",
        negative_correlations=False, **kargs):
        """
        **Purpose**
            perform 'gene-bundling' on the expression data.

            I still can't get k-means to correctly bundle together similar genes. So instead I will roll my own
            gene-similarity bundler.

            Output a stack of plots and return the groups

        **Arguments**
            filename (Optional)
                filename to save the profiles to

            seeds (Required)
                For each item in the list, bundle together genes based on similarity.

            key (Required)
                A key to use to extract the initial seeds from

            mode (Optional, default="pearsonr")
                Valid modes are:
                    pearsonr - Pearson correlation
                    spearmanr - Rank order correlation

            Z (Optional, default=0.9)
                Pearson correlation score required for a match (1.0 = perfect match, higher is better).
                Note that the means are unscaled, so this will recover patterns that are similar but may have
                much higher or lower mean signals.

            plotseeds (Optional, default=True)
                Plot the seeds on the mini-plots with a thick black line.

            negative_correlations (Optional, default=False)
                Also include negative correlations. Uses the same threshold as specified in 'Z', but
                now -Z is used to cut off correlations.

        **Returns**
            A dictionary in the form:

            {"<seed>": genelist(), "<seed2>": genelist(), ...}
        """
        assert mode in ("pearsonr", "spearmanr"), "bundle: '%s' mode not found!" % mode

        # Get the expression data.
        data = numpy.array(self.serialisedArrayDataList).T

        stds = numpy.std(data, axis=0)

        centroids = []
        for seed_bundle in seeds:
            expns = []
            if isinstance(seed_bundle, tuple): # bundle of bundles.
                # get all the expression for this bundle:
                for i in seed_bundle:
                    col = self.index(key, i) # I need to know the column numbers for each item in seeds:
                    expns.append(self[col]["conditions"])
                expn_data = numpy.average(expns, axis=0)
                assert expn_data.std() > 0.0, 'Variance is too low for gene "%s"' % seed_bundle
                centroids.append({"name": ",".join(seed_bundle), "avg": expn_data, "std": numpy.std(expn_data, axis=0), "members": []})
            else: # singlet
                col = self.index(key, seed_bundle) # I need to know the column numbers for each item in seeds:
                expn_data = numpy.array(self[col]["conditions"])
                assert expn_data.std() > 0.0, 'Variance is too low for gene "%s"' % seed_bundle
                centroids.append({"name": seed_bundle, "avg": expn_data, "std": stds, "members": []}) # Notice uses the entire array std

        # put all of the genes into a bundle:
        for gene in self:
            scores = []
            for cen in centroids:

                #variation = numpy.abs(gene["conditions"] - cen["avg"]) / cen["std"] # This would be better done with CV.
                #score = numpy.sum(variation) / len(cen["avg"])
                if mode == "pearsonr":
                    score = scipy.stats.pearsonr(cen["avg"], gene["conditions"])[0]
                elif mode == "spearmanr":
                    score = scipy.stats.spearmanr(cen["avg"], gene["conditions"])[0]
                scores.append(score)
            #print scores

            if max(scores) > Z: # It has at least one good score;
                best = scores.index(max(scores))
                centroids[best]["members"].append(gene) # put it in the best scoring group
                gene["score"] = max(scores) # return score
            elif negative_correlations and min(scores) < -Z:
                best = scores.index(min(scores))
                centroids[best]["members"].append(gene) # put it in the best scoring group
                gene["score"] = min(scores) # return score

                # I need to do a subsequent pass to add all other genes into a bin.

        # Guess a suitable arrangement for the figure
        sq = math.ceil(math.sqrt(len(centroids)))
        if not "size" in kargs:
            kargs["size"] = (sq*2, sq*2)

        fig = self.draw.getfigure(**kargs)

        fig.subplots_adjust(0.01, 0.01, 0.98, 0.98, wspace=0.15, hspace=0.15)

        xx = numpy.arange(len(self._conditions))

        for i, cen in enumerate(centroids):
            ax = fig.add_subplot(sq, sq, i+1)

            if cen["members"]:
                dd = numpy.array([g["conditions"] for g in cen["members"]])

                ax.plot(xx, dd.T, alpha=0.1, color="grey")

                if plotseeds:
                    ax.plot(cen["avg"], color="black") # err?
                ax.set_title("%s (%s)" % (cen["name"], len(cen["members"])), size=6)

            [t.set_fontsize(6) for t in ax.get_xticklabels()]
            [t.set_fontsize(6) for t in ax.get_yticklabels()]

            ax.set_xlim([-1,len(self._conditions)])

        res = {}
        for c in centroids:
            gl = genelist(name=c["name"])
            gl.load_list(c["members"])
            res[c["name"]] = gl

        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("bundle: Saved '%s'" % actual_filename)
        return(res)

    def volcano(self, condition_name, p_value_key, label_key=None, filename=None,
        label_fontsize=6, label_significant=0.01, **kargs):
        """
        **Purpose**
            draw a Volcano plot (fold change versus P/Q-value

            Assumes that your data is already in fold-change, or Z-score or some other
            value centered on 0.

        **Arguments**
            condition_name (Required)
                the name of the condition to use.

            key (Optional, Required if 'genelist' used)
                The key to match between the expression data and the genelist.

            label_key (Optional, default=None)
                Put a text label on top of those points that are <label_significant

                None = draw no labels.

            label_significant (Optional, default=0.01)
                Put a text label on top of those points that are <label_significant

            label_fontsize (Optional, default=6)
                labels fontsize

            spot_size (Optional, default=5)
                The size of each dot.

            available key-word arguments:
            xlabel, ylabel, title, log (set this to the base to log the data by),
            xlims, ylims, spot_size,

        **Returns**
            the actual filename saved as and a new image in filename.
        """
        assert filename, "no filename specified"
        assert condition_name in self.serialisedArrayDataDict, "%s condition not found" % x_condition_name
        assert p_value_key in self.keys(), '"%s" p_value_key not found in this list' % p_value_key

        x_data = self.getDataForCondition(condition_name)
        y_data = self[p_value_key]

        if not "xlabel" in kargs:
            kargs["xlabel"] = 'Fold-change'
        if not "ylabel" in kargs:
            kargs["ylabel"] = 'P-value'

        xlim = max(x_data)

        if label_key:
            assert label_key in self.keys(), 'label_key "%s" not found in this list' % label_key
            tx = []
            ty = []
            matches = []

            for x, p, label in zip(x_data, y_data, self['name']):
                if p < label_significant:
                    tx.append(x)
                    ty.append(p)
                    matches.append(label)
                kargs["spot_labels"] = matches

            real_filename = self.draw.nice_scatter(x=x_data, y=y_data,
                filename=filename,
                #spots=(tx, ty),
                xlims=[-xlim, xlim],
                plot_diag_slope=False,
                label_fontsize=label_fontsize,
                **kargs)
        #else:
        #    real_filename = self.draw.nice_scatter(x_data, y_data, filename, xlims=[-xlim, xlim],
        #        plot_diag_slope=False, **kargs)


        config.log.info("scatter: Saved '%s'" % real_filename)
        return(True)

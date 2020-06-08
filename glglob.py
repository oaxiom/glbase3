"""

glglobs are almagamated lists, useful for drawing comparisons between lists.

renamed to glglob as it clashes with a matplotlib and python module name

"""

import sys, os, csv, string, math, numpy, pickle

from numpy import array, zeros, object_, arange
from copy import deepcopy
from operator import itemgetter
from statistics import pstdev, mean

from . import config, utils
from .flags import *
from .base_genelist import _base_genelist
from .draw import draw
from .errors import AssertionError, NotImplementedError, GlglobDuplicateNameError
from .location import location
from .progress import progressbar
from .genelist import genelist
from .expression import expression

import matplotlib.pyplot as plot
import matplotlib.cm as cm
from scipy.stats import spearmanr, pearsonr

if config.SKLEARN_AVAIL: # These are optional
    from sklearn.cluster import DBSCAN

class glglob(_base_genelist): # cannot be a genelist, as it has no keys...
    def __init__(self, *args, **kargs):
        """
        **Purpose**
            An bunch of genelists to enable across-genelist meta comparisons

        **Arguments**
            list of genelists
                The first argument must be a list of genelists.

            Optional arguments:
                None specified.
        """
        _base_genelist.__init__(self)

        # args should be a list of lists.
        # we then store them in the linearData set

        if args: # So we can have empty glglobs.
            self.linearData = args
            self._optimiseData()
        else:
            self.linearData = []
        self.draw = draw(self)

    def __repr__(self):
        return("glbase.glglob")

    def __str__(self):
        # work out all of the types
        types = []
        for item in self.linearData:
            if item.__repr__() not in types:
                types.append(item.__repr__())
        return("glglob contains: %s items of type(s): %s" % (len(self.linearData), ", ".join(types)))

    def _optimiseData(self): # no keys, so would die.
        """
        (Override)

        actually does have keys: the list names
        """
        self.__list_name_lookback = {} # the index location.
        for index, item in enumerate(self.linearData):
            if item.name in self.__list_name_lookback:
                raise GlglobDuplicateNameError("self._optimiseData", item.name)
            else:
                self.__list_name_lookback[item.name] = index

    def loadCSV(self):
        config.log.error("glglobs cannot be represented as a CSV/TSV file, use glload() to load binary files")
        return(False)

    def saveCSV(self):
        config.log.error("glglobs cannot be represented as a CSV/TSV file, use .save() to save binary files")
        return(False)

    def __getitem__(self, value):
        """
        (Override)

        glglobs should be accesible by name.
        """
        if value in self.__list_name_lookback:
            return(self.linearData[self.__list_name_lookback[value]]) # get the actual item
        return(None)

    def __setitem__(self, value):
        """
        (Override)
        """
        config.log.error("glglobs cannot be written to")

    def compare(self, key=None, filename=None, method=None, delta=200, matrix_tsv=None,
        mode='fast',
        row_cluster=True, col_cluster=True, bracket=None,
        pearson_tsv=None,
        jaccard=False,
        **kargs):
        """
        **Purpose**
            perform a square comparison between the genelists according
            to some sort of criteria.

            Performs Pearson correlation between the patterns of overlap

        **Arguments**
            key (string, required)
                key to use when performing the comparison.

            filename (string, required)
                filename to save heatmap image as.

            method (string, "overlap|collide|map", required)
                method to use to compare the genelists

            delta (integer, optional, default=200)
                an optional command for delta expansion of the reads. This feeds into
                "overlap" and "collide" methods. And is ignored by "map"
                See the documentation for overlap and collide in genelist for
                exactly what these mean

            distance_score (string, "euclidean", optional defaults to "euclidean")
                Scoring method for distance caluclations.
                This is not implemented at the moment, and only uses a
                euclidean distance.

            output_pair_wise_correlation_plots (True|False)
                This is a debug option to output graphs showing the
                actual correlation matrices of the pair-wise correlations.
                Probably best only to do this if you
                know what you are doing!

            matrix_tsv (filename)
                If true save a tsv containing the overlap scores

            row_cluster (optional, True|False, default=True)
                cluster the rows of the heatmap

            col_cluster (optional, True|False, default=True)
                cluster the columns of the heatmap

            row_font_size (Optional, default=guess suitable size)
                the size of the row labels (in points). If set this will also override the hiding of
                labels if there are too many elements.

            col_font_size (Optional, default=8)
                the size of the column labels (in points)

            jaccard (Optional, default=False)
                Use the Jaccard index https://en.wikipedia.org/wiki/Jaccard_index
                instead of the (Classic) number of overlaps.

        **Result**
            returns the distance matrix if succesful or False|None if not.
            Saves an image to 'filename' containing a grid heatmap
        """
        valid_modes = ['fast', 'slow']
        valid_methods = ["overlap", "collide", "map"]
        valid_dist_score = ["euclidean"]
        distance_score = "euclidean" # override for the moment
        assert mode in valid_modes, "must use a valid method for comparison ({0})".format(", ".join(valid_modes))
        assert method in valid_methods, "must use a valid method for comparison (%s)" % ", ".join(valid_methods)
        assert filename, "Filename must be valid"
        assert key in self.linearData[0].linearData[0], "key '%s' not found" % (key,) # just check one of the lists
        assert distance_score in valid_dist_score, "%s is not a valid distance metric" % (distance_score,)

        if mode == 'slow':
            return self.__compare_slow(key=key, filename=filename, method=method, delta=delta,
                matrix_tsv=matrix_tsv, pearson_tsv=pearson_tsv,
                row_cluster=row_cluster, col_cluster=col_cluster, bracket=bracket,
                jaccard=jaccard,
                **kargs)
        elif mode == 'fast':
            return self.__compare_fast(key=key, filename=filename, method=method, delta=delta,
                matrix_tsv=matrix_tsv, pearson_tsv=pearson_tsv,
                row_cluster=row_cluster, col_cluster=col_cluster, bracket=bracket,
                jaccard=jaccard,
                **kargs)

    def __compare_fast(self, key=None, filename=None, method=None, delta=200, matrix_tsv=None,
        row_cluster=True, col_cluster=True, bracket=None, pearson_tsv=None, bin_size=5000,
        jaccard=False, **kargs):
        """
        **Purpose**
            The new style faster O(lin*lin/23) version

        """
        assert not jaccard, 'Jaccard not implemented for compare mode=fast'
        assert method != 'map', 'method=map not implemented when mode=fast'

        mat = {}
        num_samples = len(self.linearData)

        if jaccard:
            if not bracket:
                bracket = [0.0, 0.4]
        else: # Integers will do fine to store the overlaps
            if not bracket:
                bracket = [-0.2, 1]

        matrix = zeros( (len(self), len(self)), dtype=numpy.float64 ) # Must be float;

        config.log.info('Stage 1: Overlaps')
        p = progressbar(num_samples)
        # Prep the overlap table;
        for ia, this_peak_list in enumerate(self.linearData):
            for la in this_peak_list.linearData:
                chrom = la['loc']['chr']
                if chrom not in mat:
                    mat[chrom] = {}

                # check for an overlap;
                left = la['loc']['left'] - delta # overlap;
                rite = la['loc']['left'] + delta
                ctr = (left + rite) // 2
                if method == 'collide':
                    left = ctr - delta
                    rite = ctr + delta

                hit = None

                # super-fast mode (genome is binned)
                bin = int(ctr / bin_size)*bin_size# still keep proper locations if you want an accurate matrix_tsv
                bin = (bin, bin+bin_size)
                if bin not in mat[chrom]:
                    mat[chrom][bin] = [0] * num_samples # no hit found
                mat[chrom][bin][ia] = 1

                '''
                # fast mode (jiggle algortihm)
                for l in mat[chrom]:


                    if rite >= l[0] and left <= l[1]:
                        newctr = ((l[0] + left) // 2 + (l[1] + rite) // 2 ) // 2

                        hit = l
                        #hit = (newctr-delta, newctr+delta)
                        #if hit != l:
                        #    mat[chrom][hit] = mat[chrom][l] # jiggle together the peak overlaps;
                        #    del mat[chrom][l] # delete the old one, and expand the coords to jiggle the peak
                        break
                if hit:
                    mat[chrom][hit][ia] = 1
                else:
                    mat[chrom][(left,rite)] = [0] * num_samples # no hit found
                '''
            p.update(ia)

        # output the full location matrix;
        if matrix_tsv:
            names = [i.name for i in self.linearData]
            oh = open(matrix_tsv, "w")
            oh.write("%s\n" % "\t".join(['loc',] + names))

            for chrom in mat:
                for bin in sorted(mat[chrom]):
                    line = ['chr{0}:{1}-{2}'.format(chrom, bin[0], bin[1]), ]

                    line += [str(i) for i in mat[chrom][bin]]

                    oh.write('{0}\n'.format('\t'.join(line)))

            oh.close()
            config.log.info('compare: save {0} matrix_tsv'.format(matrix_tsv))

        # Go through the table once more and merge overlapping peaks?
        # Only if the jiggle=True

        # You can now dispose of the location data and convert the matrices to numpy arrays
        config.log.info('Stage 2: Clean matrix')
        total_num_peaks = 0
        for chrom in mat:
            mat[chrom] = numpy.array([mat[chrom][loc] for loc in mat[chrom]])
            total_num_peaks += mat[chrom].shape[0]
        config.log.info('Total number of peaks = {0:,}'.format(total_num_peaks))

        # Now it's simpler, take the column sums;
        config.log.info('Stage 3: Collect overlaps')

        peak_lengths = [len(a) for a in self.linearData]

        # convert to the triangular matrix:
        p = progressbar(len(peak_lengths))
        for ia, la in enumerate(peak_lengths):
            for ib, lb in enumerate(peak_lengths):
                if ia == ib:
                    matrix[ia, ib] = peak_lengths[ia] # should be filled in with the maximum possible overlap.
                    continue
                elif ia < ib: # make triangular
                    continue

                # the overlap is each row sums to 2
                res = 1 # If two lists collide to produce 0 hits it eventually ends up with nan, so put in a pseudo overlap;
                for chrom in mat:
                    #print(mat[chrom][ia,:], mat[chrom][ib:,])
                    s = mat[chrom][:,ia] + mat[chrom][:,ib] # down, across
                    res += len(s[s>=2])

                if jaccard:
                    res = (res / float(len(la) + len(lb) - res))

                matrix[ia, ib] = res

            p.update(ia)

        # fill in the gaps in the triangle
        for ia, la in enumerate(self.linearData):
            for ib, lb in enumerate(self.linearData):
                if ia < ib:
                    matrix[ia,ib] = matrix[ib,ia]

        config.log.info('Stage 4: Correlations')
        # Normalise
        if jaccard:
            result_table = matrix # normalised already;
        else:
            # data must be normalised to the maximum possible overlap.
            for ia, la in enumerate(peak_lengths):
                for ib, lb in enumerate(peak_lengths):
                    matrix[ia,ib] = (matrix[ia,ib] / min([la, lb]))

            print(matrix)

            corr_result_table = zeros( (len(self), len(self)) ) # square matrix to store the data.
            # convert the data to a pearson score.
            for ia, this_col in enumerate(matrix):
                for ib, other_col in enumerate(matrix):
                    if ia != ib:
                        corr_result_table[ia,ib] = pearsonr(this_col, other_col)[0] # [0] = r score, [1] = p-value
                    else:
                        corr_result_table[ia,ib] = 1.0

            result_table = corr_result_table

        if pearson_tsv:
            names = [i.name for i in self.linearData]
            oh = open(pearson_tsv, "w")
            oh.write("%s\n" % "\t".join([] + names))

            for ia, la in enumerate(names):
                oh.write("%s" % la)
                for ib, lb in enumerate(names):
                    oh.write("\t%s" % corr_result_table[ia,ib])
                oh.write("\n")
            oh.close()
            config.log.info('compare: save {0} pearson_tsv'.format(pearson_tsv))

        # need to add the labels and serialise into a dict of lists.
        dict_of_lists = {}
        row_names = []
        for index, item in enumerate(self.linearData):
            dict_of_lists[item.name] = result_table[index]
            row_names.append(item.name) # preserve order of row names.

        if "output_pair_wise_correlation_plots" in kargs and kargs["output_pair_wise_correlation_plots"]:
            # output the plot matrices.
            for ia, la in enumerate(self.linearData):
                for ib, lb in enumerate(self.linearData):
                    plot_data = []
                    if ia != ib:
                        x = matrix[ia,]
                        y = matrix[ib,]
                        self.draw._scatter(x, y, xlabel=row_names[ia], ylabel=row_names[ib],
                            filename="dpwc_plot_%s_%s.png" % (row_names[ia], row_names[ib]))

        if "aspect" in kargs:
            aspect = kargs["aspect"]
        else:
            aspect = "normal"

        # respect heat_wid, hei if present
        square = True
        if "heat_hei" in kargs or "heat_wid" in kargs:
            square=False

        #print dict_of_lists
        if jaccard:
            colbar_label = 'Jaccard index'
        else:
            colbar_label = 'Pearson correlation'

        # draw the heatmap and save:
        ret = self.draw.heatmap(data=dict_of_lists, filename=filename,
            colbar_label=colbar_label, bracket=bracket,
            square=square, cmap=cm.hot, cluster_mode="euclidean", row_cluster=row_cluster, col_cluster=col_cluster,
            row_names=row_names, col_names=row_names, aspect=aspect, **kargs)

        config.log.info("compare: Saved Figure to '%s'" % ret["real_filename"])
        return(dict_of_lists)

    def __compare_slow(self, key=None, filename=None, method=None, delta=200, matrix_tsv=None,
        row_cluster=True, col_cluster=True, bracket=None, pearson_tsv=None,
        jaccard=False, **kargs):
        """
        **Purpose**
            The old-style all-vs all intersect based
        """

        config.log.info("This may take a while, all lists are intersected by '%s' with '%s' key" % (method, key))

        if jaccard:
            if not bracket:
                bracket = [0.0, 0.4]
            # I need a float matrix? I thought the default is a float64?
            matrix = zeros( (len(self), len(self)), dtype=numpy.float64) # 2D matrix.
        else: # Integers will do fine to store the overlaps
            if not bracket:
                bracket = [-0.2, 1]
            matrix = zeros( (len(self), len(self)) ) # 2D matrix.

        for ia, la in enumerate(self.linearData):
            for ib, lb in enumerate(self.linearData):
                if ia == ib:
                    matrix[ia, ib] = len(la) # should be filled in with the maximum possible overlap.
                elif ia < ib: # make search triangular
                    pass
                else:
                    res = 0
                    if method == "collide":
                        res = la.collide(genelist=lb, loc_key=key, delta=delta)
                    elif method == "overlap":
                        res = la.overlap(genelist=lb, loc_key=key, delta=delta)
                    elif method == "map":
                        res = la.map(genelist=lb, key=key)

                    if res:
                        res = len(res)
                    else:
                        res = 1 # If two lists collide to produce 0 hits it eventually ends up with nan
                        # in the table which then buggers up the clustering below.
                    if jaccard:
                        res = (res / float(len(la) + len(lb) - res))
                    matrix[ia, ib] = res

        #print matrix

        # fill in the gaps in the triangle
        for ia, la in enumerate(self.linearData):
            for ib, lb in enumerate(self.linearData):
                if ia < ib:
                    matrix[ia,ib] = matrix[ib,ia]

        if matrix_tsv:
            names = [i.name for i in self.linearData]
            oh = open(matrix_tsv, "w")
            oh.write("%s\n" % "\t".join([] + names))

            for ia, la in enumerate(names):
                oh.write("%s" % la)
                for ib, lb in enumerate(names):
                    oh.write("\t%s" % matrix[ia,ib])
                oh.write("\n")
            oh.close()

        if jaccard:
            result_table = matrix
        else:
            # data must be normalised to the maximum possible overlap.
            for ia, la in enumerate(self.linearData):
                for ib, lb in enumerate(self.linearData):
                    matrix[ia,ib] = (matrix[ia,ib] / min([len(la), len(lb)]))
            corr_result_table = zeros( (len(self), len(self)) ) # square matrix to store the data.
            # convert the data to a pearson score.
            for ia, this_col in enumerate(matrix):
                for ib, other_col in enumerate(matrix):
                    if ia != ib:
                        corr_result_table[ia,ib] = pearsonr(this_col, other_col)[0] # [0] = r score, [1] = p-value
                    else:
                        corr_result_table[ia,ib] = 1.0

            result_table = corr_result_table

        if pearson_tsv:
            names = [i.name for i in self.linearData]
            oh = open(pearson_tsv, "w")
            oh.write("%s\n" % "\t".join([] + names))

            for ia, la in enumerate(names):
                oh.write("%s" % la)
                for ib, lb in enumerate(names):
                    oh.write("\t%s" % corr_result_table[ia,ib])
                oh.write("\n")
            oh.close()

        # need to add the labels and serialise into a doct of lists.
        dict_of_lists = {}
        row_names = []
        for index, item in enumerate(self.linearData):
            dict_of_lists[item.name] = result_table[index]
            row_names.append(item.name) # preserve order of row names.

        if "output_pair_wise_correlation_plots" in kargs and kargs["output_pair_wise_correlation_plots"]:
            # output the plot matrices.
            for ia, la in enumerate(self.linearData):
                for ib, lb in enumerate(self.linearData):
                    plot_data = []
                    if ia != ib:
                        x = matrix[ia,]
                        y = matrix[ib,]
                        self.draw._scatter(x, y, xlabel=row_names[ia], ylabel=row_names[ib],
                            filename="dpwc_plot_%s_%s.png" % (row_names[ia], row_names[ib]))

        if "aspect" in kargs:
            aspect = kargs["aspect"]
        else:
            aspect = "normal"

        # respect heat_wid, hei if present
        square = True
        if "heat_hei" in kargs or "heat_wid" in kargs:
            square=False

        #print dict_of_lists
        if jaccard:
            colbar_label = 'Jaccard index'
        else:
            colbar_label = 'Pearson correlation'

        # draw the heatmap and save:
        ret = self.draw.heatmap(data=dict_of_lists, filename=filename,
            colbar_label=colbar_label, bracket=bracket,
            square=square, cmap=cm.hot, cluster_mode="euclidean", row_cluster=row_cluster, col_cluster=col_cluster,
            row_names=row_names, col_names=row_names, aspect=aspect, **kargs)

        config.log.info("compare: Saved Figure to '%s'" % ret["real_filename"])
        return(dict_of_lists)

    def venn(self, key=None, filename=None, mode='map', **kargs):
        """
        **Purpose**
            draw 2, 3 or 4 venn Diagrams
            currently only equally size venndiagrams are supported.
            (proportional venn diagrams are experimental only, enable them
            using experimental_proportional_venn = True as an argument).

            your glglob should be loaded with several genelist-like objects.

            Note that you can do simple 2 overlap venn_diagrams using any
            pair of genelists with this sort of code:

            genelist.map(genelist=other_genelist, <...>, image_filename="venndiagram.png")

            Note also that the data from each genelist will be converted to unique
            values. So the final numbers may not match your original list sizes

        **Arguments**
            key
                key to use to map() between the two lists.

            filename
                save the venn diagram to this filename

            mode (Optional, default='map')
                set to 'collide' if you prefer to use the glbase location
                collide() to map the Venn.

                Note that 'key' should point to a location key, and you can pass an optional
                'delta' command

                Also note that location Venn overlaps are only partly accurate, as it is
                possible to have multiple overlaps .

            title (Optional)
                title for the figures
                defaults to <list> vs <list> vs ...

        **Returns**
            A venn diagram saved in filename.
        """
        valid_args = ["filename", "key", "title", "experimental_proportional_venn"]
        for k in kargs:
            if not k in valid_args:
                raise ArgumentError(self.map, k)

        assert len(self.linearData) <= 5, "currently glglob venn diagrams only support at most 5 overlaps"
        assert len(self.linearData) >= 2, "you must send at least two lists"
        assert filename, "no filename specified for venn_diagrams to save to"
        assert key, "Must specify a 'key' to map the two lists"

        if mode == 'collide': # Deflect to the collide routine
            return(self.__venn_collide(key=key, filename=filename))

        proportional = False
        if "experimental_proportional_venn" in kargs and kargs["experimental_proportional_venn"]:
            proportional=True

        # The code below is needlessly verbose for clarity.
        if len(self.linearData) == 2:
            A = set(self.linearData[0][key])
            B = set(self.linearData[1][key])

            AB = A & B

            realfilename = self.draw.venn2(len(A), len(B), len(AB),
                self.linearData[0].name, self.linearData[1].name,
                filename, **kargs)
            return(None)

        elif len(self.linearData) == 3:
            A = set(self.linearData[0][key])
            B = set(self.linearData[1][key])
            C = set(self.linearData[2][key])

            AB = A & B
            AC = A & C
            BC = B & C

            ABC = A & B & C

            # check for none's:
            if AB:
                AB = len(AB)
            else:
                AB = 0

            if AC:
                AC = len(AC)
            else:
                AC = 0

            if BC:
                BC = len(BC)
            else:
                BC = 0

            if ABC:
                ABC = len(ABC)
            else:
                ABC = 0

            realfilename = self.draw.venn3(len(A), len(B), len(C), AB, AC, BC, ABC,
                self.linearData[0].name, self.linearData[1].name, self.linearData[2].name,
                filename, **kargs)

        elif len(self.linearData) == 4:
            A = set(self.linearData[0][key])
            B = set(self.linearData[1][key])
            C = set(self.linearData[2][key])
            D = set(self.linearData[3][key])
            print(len(A), len(B), len(C), len(D))

            # Use set logic to work out the actual values:
            # I'm pretty sure this is accurate at this point.

            ABCD = A & B & C & D

            ABC = (A & B & C) - ABCD
            ABD = (A & B & D) - ABCD
            ACD = (A & C & D) - ABCD
            BCD = (B & C & D) - ABCD

            AB = (((A & B) - ABC) - ABD) - ABCD
            AC = (((A & C) - ABC) - ACD) - ABCD
            AD = (((A & D) - ABD) - ACD) - ABCD
            BC = (((B & C) - ABC) - BCD) - ABCD
            BD = (((B & D) - ABD) - BCD) - ABCD
            CD = (((C & D) - ACD) - BCD) - ABCD

            A = A - ABCD - ABC - ABD - ACD - AB - AC - AD
            B = B - ABCD - ABC - ABD - BCD - AB - BC - BD
            C = C - ABCD - ABC - ACD - BCD - AC - BC - CD
            D = D - ABCD - ABD - ACD - BCD - AD - BD - CD

            # A, B, C, D, AB, AC, AD, BC, BD, CD, ABC, ABD, ACD, BCD, ABCD
            lists = [len(A), len(B), len(C), len(D),
                len(AB), len(AC), len(AD), len(BC), len(BD), len(CD),
                len(ABC), len(ABD), len(ACD), len(BCD), len(ABCD)]

            labs = [self.linearData[i].name for i in range(4)]

            realfilename = self.draw.venn4(lists, labs, filename, **kargs)

        elif len(self.linearData) == 5:
            raise NotImplementedError

            A = set(self.linearData[0][key])
            B = set(self.linearData[1][key])
            C = set(self.linearData[2][key])
            D = set(self.linearData[3][key])
            E = set(self.linearData[4][key])

            AB = A & B
            AC = A & C
            AD = A & D
            AE = A & E
            BC = B & C
            BD = B & D
            BE = B & E
            CD = C & D
            CE = C & E
            DE = D & E

            ABC = A & B & C
            ABD = A & B & D
            ABE = A & B & E
            ACD = A & C & D
            ACE = A & C & E
            ADE = A & D & E
            BCD = B & C & D
            BCE = B & C & E
            BDE = B & D & E
            CDE = C & D & E

            ABCD = A & B & C & D
            ABCE = A & B & C & E
            ACDE = A & C & D & E
            BCDE = B & C & D & E

            ABCDE = A & B & C & D & E

            # A, B, C, D, E,
            # AB, AC, AD, AE, BC, BD, BE, CD, CE, DE,
            # ABC, ABD, ABE, ACD, ACE, ADE, BCD, BCE, BDE, CDE
            # ABCD, ABCE, ACDE, BCDE,
            # ABCDE
            lists = [len(A), len(B), len(C), len(D), len(E),
                len(AB), len(AC), len(AD), len(AE), len(BC), len(BD), len(BE), len(CD), len(CE), len(DE),
                len(ABC), len(ABD), len(ABE), len(ACD), len(ACE), len(ADE), len(BCD), len(BCE), len(BDE), len(CDE),
                len(ABCD), len(ABCE), len(ACDE), len(BCDE), (ABCDE)]

            labs = [self.linearData[i].name for i in range(5)]

            realfilename = self.draw.venn5(lists, labs, filename, **kargs)

        config.log.info("venn: Saved Figure to '%s'" % realfilename)
        return(None)

    def __venn_collide(self, key, filename, delta=200, **kargs):
        #assert key in self.keys()
        assert len(self.linearData) >= 3, "currently glglob venn diagrams only support at least 3-way overlaps"
        assert len(self.linearData) <= 3, "currently glglob venn diagrams only support at most 3-way overlaps"

        if len(self.linearData) == 3:
            A = self.linearData[0]
            B = self.linearData[1]
            C = self.linearData[2]

            AB = A.collide(genelist=B, key=key, delta=delta)
            AC = A.collide(genelist=C, key=key, delta=delta)
            BC = B.collide(genelist=C, key=key, delta=delta)

            ABC = AB.collide(genelist=BC, key=key, delta=delta)

            # check for none's:
            if AB:
                AB = len(AB)
            else:
                AB = 0

            if AC:
                AC = len(AC)
            else:
                AC = 0

            if BC:
                BC = len(BC)
            else:
                BC = 0

            if ABC:
                ABC = len(ABC)
            else:
                ABC = 0

            # You only need to provide the lengths, the overlaps are calculated in venn3:
            realfilename = self.draw.venn3(len(A), len(B), len(C), AB, AC, BC, ABC,
                self.linearData[0].name, self.linearData[1].name, self.linearData[2].name,
                filename, **kargs)
        return(None)

    def moving_average_maps(self, mode="graph", compare_array=None, filename=None, key=None,
        normalise=True, **kargs):
        """
        **Purpose**
            Draw moving average maps in a variety of formats.

        **Arguments**
            mode (Optional, defaults to "graph")
                "graph"
                    draw a series of line graphs on the same graph.
                "heatmap"
                    draw as a heatmap
                "stacked_plots"
                    draw as a series of stacked line plots.

            compare_array (Required)
                The name of the array to use as a compare array. You have to draw the
                moving average by comparing against some array already in
                the glglob.

            filename
                filename to save the image as.

            key
                key to use to match between the compare_array and the rest of the
                data in this glglob

            normalise (Optional, default=True)
                normalise the data, True or False.

        **Returns**
            The actual filename used to save and an image saved
            to 'filename'
        """
        # get the compare array:
        compare_array = self[compare_array]
        assert compare_array, "the compare array was not found in this glglob"

        res = {}
        last_val = 0
        for pl in self.linearData:
            if pl.name != compare_array.name: # a simple == will fail here, so I use the names to compare
                res[pl.name] = []
                last_val = 0
                names = pl[key] # get the name column.
                for i, v in enumerate(compare_array):
                    if v[key] in names:
                        #last_val += 1
                        #res[pl.name].append(last_val)
                        res[pl.name].append(1)
                    else:
                        res[pl.name].append(0)

                res[pl.name] = utils.movingAverage(res[pl.name], int(len(res[pl.name]) * 0.2))[1] # keep only y,
                typical_length = len(res[pl.name]) # measure length to generate x axis later.

        # append the compare_array for comparison.

        res[compare_array.name] = [i[0] for i in compare_array["conditions"]][0:typical_length] # this will be a list of lists (with 1 entry) - flatten the list.
        # the arange fakes the x axis for plots.
        # I have to slice the array as it's slightly shorter than the movingAverage plots...
        # You will only notice this if the list is short and you can see the names

        if normalise: #this will normalise the y-axis
            for k in res:
                # normalise to 0-->100
                min_val = min(res[k])
                max_val = max(res[k]) - min_val

                for i, v in enumerate(res[k]):
                    res[k][i] = ((v-min_val) / max_val) * 100.0

        if mode == "graph":
            # fake the xdata:
            xdata = arange(0,typical_length)
            plot.cla()
            fig = plot.figure(figsize=(8,5))
            axis = fig.add_subplot(111)
            for pl in res:
                axis.plot(xdata, res[pl][1], label=pl)

            axis.set_title("")
            axis.legend(loc=2, markerscale=0.1)#, prop={"size": "xx-small"})
            fig.savefig(filename)
            real_filename = filename
        elif mode == "heatmap":
            # At the moment the data is in the form: [ [ (x,y), (x,y), ]:class ... []]
            # I need to strip out the x data.
            scalebar_name = "score"
            if normalise:
                scalebar_name = "normalised score"
            # use key to build an array ready for draw._heatmap
            colnames=[]
            for k in res:
                res[k] = res[k]
                colnames.append(k)

            real_filename = self.draw._heatmap(data=res, filename=filename, col_names=colnames, row_names=compare_array["name"],
                row_cluster=False, col_cluster=True, colour_map=cm.Blues, #vmax=1,
                colbar_label=scalebar_name, aspect="square")
        config.log.info("moving_average_maps: Saved image to '%s'" % real_filename)

    def overlap_heatmap(self, filename=None, score_key=None, resolution=1000, optics_cluster=True):
        """
        **Purpose**
            Draw a heatmap from ChIP-seq binding data (or lists of genomic coordinates),
            showing one row for each binding site, particularly where the sites overlap.

            This will compare (and overlap) each of the lists, then produce a large heatmap
            with each row for a unique genomic location.

            For the purposes of this tool, the genome is divided up into blocks (based on
            resolution) and binding is then tested against these blocks.

        **Arguments**
            filename (Required)
                filename to save the heatmap to.

            score_key (Optional, default=None)
                by default each row will be scored in a binary manner. i.e. Is the binding site present or not?
                However, if your ChIP-seq lists have some sort of intensity score then this
                can be used instead to score the overlap.

            resolution (Optional, default=1000)
                Number of base pairs to divide the genome up into.

            optics_cluster (Optional, default=True)
                By default overlap_heatmap() will cluster the data based on OPTICS:

                https://en.wikipedia.org/wiki/OPTICS_algorithm
                If optics_cluster=False


        **Returns**
            A heatmap and the heatmap data (as a Numpy array) for any further processing.
        """
        assert filename, "Must provide a filename"
        assert config.SKLEARN_AVAIL, "the Python package sklearn is required for overlap_heatmap()"

        # Iterate through each genelist and build a non-redundant table of genomic
        # blocks

        chr_blocks = {}
        total_rows = 0

        for index, gl in enumerate(self.linearData):
            for item in gl:
                # Add chromosome to the cache if not present.
                if item["loc"]["chr"] not in chr_blocks:
                    chr_blocks[item["loc"]["chr"]] = {}

                block_id = "bid:%s" % (math.floor(item["loc"]["left"] / resolution) * 1000, )
                if block_id not in chr_blocks[item["loc"]["chr"]]:
                    chr_blocks[item["loc"]["chr"]][block_id] = [0 for x in range(len(self.linearData))] # one for each gl
                    total_rows += 1

                if score_key:
                    chr_blocks[item["loc"]["chr"]][block_id][index] = item[score_key]
                else:
                    chr_blocks[item["loc"]["chr"]][block_id][index] = 1

        config.log.info("overlap_heatmap(): Found %s unique genomic regions" % total_rows)

        # Build the table for the heatmap
        tab = numpy.zeros([len(self.linearData), total_rows])

        crow = 0
        for c in chr_blocks:
            for bid in chr_blocks[c]:
                for i in range(len(chr_blocks[c][bid])): # or len(self.linearData)
                    tab[i, crow] = chr_blocks[c][bid][i]
                crow += 1

        tab = tab.T
        # dendrogram dies here, so need other ways to cluster
        # DBSCAN consumes too much unnecessary memory.
        """
        alg = DBSCAN(eps=0.2)
        print alg.fit(tab)
        print alg.labels_
        clusters = numpy.unique(alg.labels_)
        print clusters
        # reorder the list based on cluster membership
        newd = {}
        for index, c in enumerate(alg.labels_):
            if c not in newd:
                newd[c] = []
            newd[c].append(tab[index])

        # load it back into a numpy array
        tab = None
        for c in clusters:
            new = numpy.vstack(newd[c])
            if tab is None:
                tab = new
            else:
                tab = numpy.vstack([tab, new])
        """
        # Yay, roll my own clustering!

        # I already know how many possible clusters there will be.
        #num_clusters = math.factorial(len(self.linearData))

        # build a cluster table, containing all possible variants for this len(self.linearData)
        clusters = {}
        for row in tab:
            # Make an identifier for the cluster:
            id = tuple([bool(i) for i in row])
            if id not in clusters:
                clusters[id] = []
            clusters[id].append(row)

        # I want to sort the clusters first:
        sorted_clusters = []
        for c in clusters:
            sorted_clusters.append({"id": c, "score": sum(c)})
        sorted_clusters = sorted(sorted_clusters, key=itemgetter("score"))

        # Flattent the arrays and load it back into a numpy array
        tab = None
        for c in sorted_clusters:
            new = numpy.vstack(clusters[c["id"]])
            if tab is None:
                tab = new
            else:
                tab = numpy.vstack([tab, new])

        ret = self.draw.heatmap(data=tab, filename=filename, col_names=[gl.name for gl in self.linearData], row_names=None,
                row_cluster=False, col_cluster=True, colour_map=cm.Reds, heat_wid=0.7, heat_hei=0.7, bracket=[0,tab.max()])

        config.log.info("overlap_heatmap: Saved overlap heatmap to '%s'" % ret["real_filename"])
        return(tab)

    def __peak_cluster(self, list_of_peaks, merge_peaks_distance):
        # Merge overlapping peaks
        chr_blocks = {}
        total_rows = 0
        #merged_peaks = {}
        p = progressbar(len(list_of_peaks))
        for idx, gl in enumerate(list_of_peaks):
            for p1 in gl["loc"]:
                #p1 = p1.pointify().expand(merge_peaks_distance) # about 10% of the time is in __getitem__ from the loc, so unpack it;
                cpt = (p1.loc["left"] + p1.loc['right']) // 2
                p1_chr = p1['chr']
                p1_left = cpt - merge_peaks_distance
                p1_right = cpt + merge_peaks_distance
                if not p1_chr in chr_blocks:
                    chr_blocks[p1_chr] = {}

                binary = [0 for x in range(len(list_of_peaks))] # set-up here in case I need to modify it.

                for p2 in chr_blocks[p1_chr]: # p2 is now a block_id tuple
                    #if p1.qcollide(p2):
                    if p1_right >= p2[0] and p1_left <= p2[1]: # unfolded for speed.
                        binary = chr_blocks[p1_chr][p2]["binary"] # preserve the old membership

                        # remove the original entry
                        del chr_blocks[p1_chr][p2]
                        total_rows -= 1

                        # Add in a new merged peak:
                        cpt = (((p1_left+p2[0])//2) + ((p1_right+p2[1])//2)) // 2 # pointify()

                        p1_left=cpt-merge_peaks_distance
                        p1_right=cpt+merge_peaks_distance
                        # Don't get confused here, p1 is added onto the block heap below:
                        break

                # modify binary to signify membership for this peaklist
                binary[idx] = 1

                # Add p1 onto the blocklist
                block_id = (p1_left, p1_right)
                if block_id not in chr_blocks[p1_chr]:
                    chr_blocks[p1_chr][block_id] = {"binary": binary,
                        "pil": [0 for x in range(len(list_of_peaks))]} # one for each gl, load pil with dummy data.
                    total_rows += 1 # because the result is a dict of dicts {"<chrname>": {"bid": {data}}, so hard to keep track of the total size.

            p.update(idx)
        return total_rows, chr_blocks

    def chip_seq_cluster(self, list_of_peaks, merge_peaks_distance=400, sort_clusters=True,
        _get_chr_blocks=False, **kargs):
        """
        **Purpose**
            Combine and merge all peaks, extract the read pileups then categorize the peaks into
            similar groupings. Return a new list of genelists, one genelist for each grouping
            that contains the list of genomic locations in each group.

            Return a glbase expression object with each row a merged (unique) peak, each
            column is a peak

            Be careful, the resulting objects can get very huge!

            The order of the genomic locations and order of the groups must be maintained
            between the heatmap and the returned data.

            NOTE: I sort of named this function incorectly with the whole 'cluster' business.
            Although it's not wrong to label the returned groups as clusters it is certainly
            confusing and may imply that some sort of k-means or hierarchical clustering
            is performed. No clustering is performed, instead groups are made based on a binary
            determination from the list_of_peaks. So below, where I refer to 'cluster'
            I really mean group. Later I may add k-means clustering, which may make things even more
            confusing.

            Here is a detailed explanation of this function:

            1. Join all of the peaks into a redundant set of coordinates

            2. Merge all of the genomic regions to produce a single list of unique genomic regions
            (this is what it means by "chip_seq_cluster_heatmap(): Found <number> unique
            genomic regions")

            3. Build a table of all possible peak combinations:

            e.g. for two chip-seq lists, A and B:

            listA only: [True, False]
            listB only: [False, True]
            listA and listB: [True, True]

            It is these that are the 'clusters' (or groups). In this case there would
            be just 3 groups. The more lists the more possible groups.

            Note that groups with no members are culled.

        **Arguments**
            list_of_peaks (Required)
                A list of genelists of peaks from your ChIP-seq data to interrogate. The order of the libraries
                Genomic location data should be stored in a 'loc' key in the genelist.

            merge_peaks_distance (Optional, default=400)
                Maximum distance that the centers of any two peaks can be apart before the two peaks are merged into
                a single peak. (taking the mean of the peak centers)

            sort_clusters (Optional, default=True)
                sort the clusters from most complex to least complex.
                Note that chip_seq_cluster_heatmap cannot preserve the order of the peaks
                (it's impossible), so setting this to false will just randomise the order of the clusters
                which may not be particularly helpful.

        **Returns**
            Returns a glbase expression object, with rows as unique genomic peaks and
            columns as each peak list.
            The values will be filled with 0 or 1, if it was a peak or not a peak.
        """
        assert list_of_peaks, 'list_of_peaks is empty'
        assert len(list_of_peaks[0]) > 0, 'list_of_peaks lists appear to be empty'

        # get a non-redundant list of genomic regions based on resolution.
        chr_blocks = {} # stores a binary identifier
        pil_blocks = {}
        total_rows = 0

        peak_lengths = sum([len(p) for p in list_of_peaks])
        config.log.info("chip_seq_cluster_heatmap: Started with {0} redundant peaks".format(peak_lengths))
        total_rows, chr_blocks = self.__peak_cluster(list_of_peaks, merge_peaks_distance)
        config.log.info("chip_seq_cluster: Found %s unique genomic regions" % total_rows)

        if _get_chr_blocks:
            return chr_blocks

        # Convert the chr_blocks into a expression object
        tab = []
        for chrom in chr_blocks:
            for loc in chr_blocks[chrom]:
                l = location(chr=chrom, left=loc[0], right=loc[1])
                cid = int("".join([str(i) for i in chr_blocks[chrom][loc]["binary"]]), 2)
                #print cid
                tab.append({'loc': l, 'conditions': chr_blocks[chrom][loc]['binary'], 'cid': cid})

        e = expression(loadable_list=tab, cond_names=[p.name for p in list_of_peaks])
        if sort_clusters:
            e.sort('cid')

        return e

    def chip_seq_cluster_heatmap(self, list_of_peaks, list_of_trks, filename=None, normalise=False, bins=20,
        pileup_distance=1000, merge_peaks_distance=400, sort_clusters=True, cache_data=False, bracket=None,
        range_bracket=None, frames=False, titles=None, read_extend=200, imshow=True, cmap=cm.plasma,
        log_pad=None, log=2,
        size=None, **kargs):
        """
        **Purpose**
            Combine and merge all peaks, extract the read pileups then categorize the peaks into
            similar groupings. Return a new list of genelists, one genelist for each grouping
            that contains the list of genomic locations in each group. Finally, draw a nice
            heatmap to <filename>.

            The order of the genomic locations and order of the groups must be maintained
            between the heatmap and the returned data.

            NOTE: I sort of named this function incorectly with the whole 'cluster' business.
            Although it's not wrong to label the returned groups as clusters it is certainly
            confusing and may imply that some sort of k-means or hierarchical clustering
            is performed. No clustering is performed, instead groups are made based on a binary
            determination from the list_of_peaks. So below, where I refer to 'cluster'
            I really mean group. Later I may add k-means clustering, which may make things even more
            confusing.

            Here is a detailed explanation of this function:

            1. Join all of the peaks into a redundant set of coordinates

            2. Merge all of the genomic regions to produce a single list of unique genomic regions
            (this is what it means by "chip_seq_cluster_heatmap(): Found <number> unique
            genomic regions")

            3. Build a table of all possible peak combinations:

            e.g. for two chip-seq lists, A and B:

            listA only: [True, False]
            listB only: [False, True]
            listA and listB: [True, True]

            It is these that are the 'clusters' (or groups). In this case there would
            be just 3 groups. The more lists the more possible groups.

            Note that groups with no members are culled.

            4. for <each genome location> get the pileup from the approprate track. expand
            around the genomic location of the bin by <pileup_distance> and build a heat map.
            (This is the really slow bit). The option <bins> will bin the pileup
            (i.e. chr1:10001-10400 would have 400 base pairs, but would be divided into twenty bins.
            This makes drawing the heatmap feasable as too many squares for the heatmap
            will consume RAM and CPU.

            6. Order the heatmap from most complex bin (genome regions with peaks in all libraries)
            to least complex bin (chip-seq library specific peaks) and draw.

        **Arguments**
            list_of_peaks (Required)
                A list of genelists of peaks from your ChIP-seq data to interrogate. The order of the libraries
                MUST be the same as the order of the list_of_trks. genomic location data
                should be stored in a 'loc' key in the genelist.

            list_of_trks (Required)
                A list of trks to draw the sequence tag reads from to build the pileups.
                This list MUST be in the same order as the list_of_peaks.

            filename (Optional, default=None)
                If set to a string a heatmap will be saved to filename.

            normalise (Optional, default=False)
                Normalize the read pileup data within each library to the size of the
                library to assist in cross-comparison of ChIP-seq libraries.

            merge_peaks_distance (Optional, default=400)
                Maximum distance that the centers of any two peaks can be apart before the two peaks are merged into
                a single peak. (taking the mean of the peak centers)

            pileup_distance (Optional, default=1000)
                distance around the particular bin to draw in the pileup.

            read_extend (Optional, default=200)
                The size in base pairs to extend the read. If a strand is present it will expand
                from the 3' end of the read.

            bins (Optional, default=20)
                number of bins to use for the pileup. Best to use conservative numbers (10-50) as
                large numbers of bins can consume huge amounts of memory.

            sort_clusters (Optional, default=True)
                sort the clusters from most complex to least complex.
                Note that chip_seq_cluster_heatmap cannot preserve the order of the peaks
                (it's impossible), so setting this to false will just randomise the order of the clusters
                which may not be particularly helpful.

            log (Optional, default=2)
                Use logarithms for the heatmap. Possible options are 2 and 10.

            cmap (Optional, default=matplotlib.cm.YlOrRd)
                A colour map for the heatmap.

            titles (Optional, default=peaks.name)
                Supply your own titles for the top of the heatmap columns

            range_bracket (Optional, default=None, exclusive with range_bracket)
                Okay, hold your hats, this is complicated.
                range_bracket will bracket the range of values as a fraction between [min(data), max(data)]
                i.e. If range_bracket=0.5 (the default) then the data is bracketed as:
                [max(data)*range_bracket[0], max(data)*range_bracket[1]]. The practical upshot of this is it allows you to shift
                the colour bracketing on the heatmap around without having to spend a long time finding
                a suitable bracket value.

                Bracketing is performed AFTER log.

                Typical bracketing would be something like [0.4, 0.9]

                By default glbase attempts to guess the best range to draw based on the
                median and the stdev. It may not always succeed.

            bracket (Optional, default=None, exclusive with range_bracket)
                chip_seq_cluster_heatmap() will make a guess on the best bracket values and will output that
                as information. You can then use those values to set the bracket manually here.
                This is bracketing the data AFTER log transforming.

            cache_data (Optional, default=False)
                cache the pileup data into the file specified in cache_data. This speeds up analysis.
                Note that storage of data is AFTER normalisation, resolution, pileup_distance,
                bins, but before sort_clusters and before heatmap drawing.
                This allows you to store the very slow part of chip_seq_cluster_heatmap()
                and so iterate through different heatmap drawing options without having to
                do the whole pileup again.

                note that if cache_data file does not exist then it will be created and
                pileup data generated. If the file does exist, data will be read from that
                file and used for heatmap drawing.

            frames (Optional, default=False)
                Draw black frames around the heatmaps and category maps. I prefer without,
                so that is the default!

            imshow (Optional, default=False)
                Embed the heatmap as an image inside a vector file. (Uses matplotlib imshow
                to draw the heatmap part of the figure. Allows very large matrices to
                be saved as an svg, with the heatmap part as a raster image and all other elements
                as vectors).

        **Returns**
            Returns a list of genelists, one genelist for each major category. The genelist
            contains a list of locations belonging to that group. Note that the genomic locations
            may not exactly match with the original provided locations as chip_seq_cluster_heatmap()
            will merge redundant peaks by taking the mid point between two close peaks.

            The order of the genomic locations and order of the groups will be maintained
            between the heatmap and the returned data.

            The formal returned data is in a dict so that it can give information about each cluster grouping:

            {"<cluster_id>": {"genelist": <a genelist object>, "cluster_membership": (True, True, ..., False)}, ...}

            The "cluster_membership" key returns a tuple of the same length as the number
            list_of_peaks (and in the same order) indicating that this particular cluster
            represents binding (True) or not (False) in each original list_of_peaks.
        """
        assert not (range_bracket and bracket), "You can't use both bracket and range_bracket"
        assert len(list_of_peaks) == len(list_of_trks), 'len(list_of_peaks) != len(list_of_trks)'

        # get a non-redundant list of genomic regions based on resolution.
        chr_blocks = {} # stores a binary identifier
        pil_blocks = {}
        total_rows = 0
        resolution = merge_peaks_distance # laziness hack!

        # Confirm that all lists contain a 'loc' key
        assert False not in ['loc' in gl.keys() for gl in list_of_peaks], 'One of your peak data (list_of_peaks) does not contain a "loc" key'

        peak_lengths = sum([len(p) for p in list_of_peaks])
        config.log.info("chip_seq_cluster_heatmap: Started with {0} redundant peaks".format(peak_lengths))
        total_rows, chr_blocks = self.__peak_cluster(list_of_peaks, merge_peaks_distance)
        config.log.info("chip_seq_cluster: Found %s unique genomic regions" % total_rows)

        # Get the size of each library if we need to normalize the data.
        if normalise:
            # get and store the read_counts for each library to reduce an sqlite hit.
            read_totals = [trk.get_total_num_reads()/float(1e6) for trk in list_of_trks]

        # I will need to go back through the chr_blocks data and add in the pileup data:
        bin_size = int((resolution+resolution+pileup_distance) / bins)
        block_len = (resolution+resolution+pileup_distance+pileup_distance) # get the block size
        data = None

        # sort out cached_data
        if cache_data and os.path.isfile(cache_data): # reload previously cached data.
            oh = open(os.path.realpath(cache_data), "rb")
            chr_blocks = pickle.load(oh)
            oh.close()
            config.log.info("chip_seq_cluster_heatmap: Reloaded previously cached pileup data: '%s'" % cache_data)
            # this_loc will not be valid and I test it for length below, so I need to fake one.

        else:
            # No cached data, so we have to collect ourselves.
            config.log.info('chip_seq_cluster_heatmap: Collecting pileup data...')
            p = progressbar(len(list_of_trks))
            # New version that grabs all data and does the calcs in memory, uses more memory but ~2-3x faster
            for pindex, trk in enumerate(list_of_trks):
                for index, chrom in enumerate(chr_blocks):
                    # The chr_blocks iterates across all chromosomes, so this only hits the db once per chromosome:
                    del data
                    data = trk.get_array_chromosome(chrom, read_extend=read_extend) # This will use the fast cache version if available.

                    for block_id in chr_blocks[chrom]:
                        left = block_id[0] - pileup_distance
                        right = block_id[1] + pileup_distance
                        # It's possible to ask for data beyond the edge of the actual data. So trim the right vals
                        if right > len(data):
                            right = len(data)
                        if left > len(data):
                            left = len(data)

                        dd = data[left:right]
                        #dd = trk.get(loc=None, c=chrom, left=left, rite=right) # single gets are faster than the above messy stuff;

                        if len(dd) < block_len: # This should be a very rare case...
                            num_missing = block_len - len(dd)
                            ad = numpy.zeros(num_missing)
                            dd = numpy.append(dd, ad)
                            dd = list(dd)

                        if normalise:
                            # normalise before bin?
                            pil_data = [av/read_totals[pindex] for av in dd]

                        chr_blocks[chrom][block_id]["pil"][pindex] = [sum(dd[i:i+bin_size]) for i in range(0, len(dd), bin_size)] #pil_data = utils.bin_sum_data(dd, bin_size)
                p.update(pindex)

            if cache_data: # store the generated data for later.
                oh = open(cache_data, "wb")
                pickle.dump(chr_blocks, oh, -1)
                oh.close()
                config.log.info("chip_seq_cluster_heatmap: Saved pileup data to cache file: '{0}'".format(cache_data))

        # assign each item to a group and work out all of the possible groups
        cluster_ids = []
        for chrom in chr_blocks:
            for block_id in chr_blocks[chrom]:
                cluster_id = tuple([bool(i) for i in chr_blocks[chrom][block_id]["binary"]])
                if cluster_id not in cluster_ids:
                    cluster_ids.append(cluster_id)
                chr_blocks[chrom][block_id]["cluster_membership"] = cluster_id

        # I want to sort the groups from the most complex to the least complex.
        if sort_clusters:
            sorted_clusters = []
            for c in cluster_ids:
                sorted_clusters.append({"id": c, "score": sum(c)})
            sorted_clusters = sorted(sorted_clusters, key=itemgetter("score"))
            # This result is actually least to most, but as the heatmap is drawn bottom to top it makes sense to
            # preserve this order.
        else:
            pass
            #URK!

        # build the super big heatmap table
        tab_wid = block_len * len(list_of_peaks)
        tab_spa = len(list_of_peaks) # distance between each set of blocks.
        tab = None

        ret_data = {}
        list_of_tables = [None for i in list_of_peaks]
        groups = []
        pileup_data = {}

        # And arrange according to the groups.
        for cluster_index, cluster_id in enumerate(sorted_clusters):
            for chrom in chr_blocks:
                for block_id in chr_blocks[chrom]:
                    if chr_blocks[chrom][block_id]["cluster_membership"] == cluster_id["id"]:
                        for peaks in range(len(list_of_peaks)):
                            row = chr_blocks[chrom][block_id]["pil"][peaks]
                            if list_of_tables[peaks] is None:
                                # append together all pileup data in a long row and stick on the tab array.
                                list_of_tables[peaks] = [row,]
                            else:
                                list_of_tables[peaks].append(row)

                            # store the pileup_data for later.
                            if (cluster_index+1) not in pileup_data:
                                pileup_data[cluster_index+1] = [None for i in list_of_peaks]

                            if pileup_data[cluster_index+1][peaks] is None: # numpy testing.
                                pileup_data[cluster_index+1][peaks] = numpy.array(row, dtype=numpy.float64)
                            else:
                                pileup_data[cluster_index+1][peaks] += row

                        # Also add it into the return data.
                        if cluster_index+1 not in ret_data:
                            ret_data[cluster_index+1] = {"genelist": genelist(name="cluster_%s" % (cluster_index+1,)), "cluster_membership": cluster_id["id"]}
                        this_loc = location(loc="chr%s:%s-%s" % (chrom, int(block_id[0]), int(block_id[1]))) # does not include the pileup_distance
                        ret_data[cluster_index+1]["genelist"].linearData.append({"loc": this_loc})
                        groups.append(cluster_index+1)

        # finish off the pileup_data:
        for cid in pileup_data:
            for pid in range(len(pileup_data[cid])):
                pileup_data[cid][pid] /= len(ret_data[cid]["genelist"])# for i in pileup_data[cid][pid]]

        self.__pileup_data = pileup_data
        self.__pileup_names = [g.name for g in list_of_peaks] # names for each sample, taken from peaks.
        self.__pileup_groups_membership = sorted_clusters
        self.__pileup_group_sizes = [groups.count(i) for i in range(0, len(sorted_clusters)+1)]

        config.log.info("chip_seq_cluster_heatmap: There are %s groups" % len(sorted_clusters))

        # rebuild the genelist quickdata and make genelist valid:
        for cid in ret_data:
            ret_data[cid]["genelist"]._optimiseData()

        colbar_label = "Tag density"

        if log:
            if not log_pad:
                log_pad = 0.1

            for index in range(len(list_of_tables)):
                if log == 2:
                    list_of_tables[index] = numpy.log2(numpy.array(list_of_tables[index])+log_pad)
                    colbar_label = "Log2(Tag density)"
                elif log == 10:
                    list_of_tables[index] = numpy.log10(numpy.array(list_of_tables[index])+log_pad)
                    colbar_label = "Log10(Tag density)"
                else:
                    list_of_tables[index] = numpy.array(list_of_tables[index])

        if normalise:
            colbar_label = "Normalised %s" % colbar_label

        self.__pileup_y_label = colbar_label

        tab_max = max([tab.max() for tab in list_of_tables]) # need to get new tab_max for log'd values.
        tab_min = min([tab.min() for tab in list_of_tables])
        #tab_median = numpy.median([numpy.median(tab) for tab in list_of_tables])
        tab_mean = mean([numpy.average(tab) for tab in list_of_tables])
        tab_stdev = numpy.std(numpy.array([tab for tab in list_of_tables]))

        config.log.info("chip_seq_cluster_heatmap: min=%.2f, max=%.2f, mean=%.2f, stdev=%.2f" % (tab_min, tab_max, tab_mean, tab_stdev))
        if range_bracket:
            bracket = [tab_max*range_bracket[0], tab_max*range_bracket[1]]
        elif bracket:
            bracket = bracket # Fussyness for clarity.
        else: # guess a range:
            bracket = [tab_mean, tab_mean+(tab_stdev*2.0)]
            config.log.info("chip_seq_cluster_heatmap: suggested bracket = [%s, %s]" % (bracket[0], bracket[1]))

        #real_filename = self.draw.heatmap2(filename=filename, row_cluster=False, col_cluster=False,
        #    data=tab, colbar_label=colbar_label, bracket=bracket)
        if filename:
            if not titles:
                titles = [p.name for p in list_of_peaks]
            real_filename = self.draw.multi_heatmap(filename=filename, groups=groups, titles=titles, imshow=imshow, size=size,
                list_of_data=list_of_tables, colour_map=cmap, colbar_label=colbar_label, bracket=bracket, frames=frames)

        config.log.info("chip_seq_cluster_heatmap: Saved overlap heatmap to '%s'" % real_filename)
        return ret_data

    def chip_seq_cluster_pileup(self, filename=None, multi_plot=True, **kargs):
        """
        **Purpose**
            This is an addendum to chip_seq_cluster_heatmap(). You only need run this
            directly after chip_seq_cluster_heatmap() and it will draw aggregate pileup
            graphs either all on the same graph (with a legend) or will draw a multi_plot
            with many graphs (when multi_plot=True, the default).

            Note, you do not need to respecify the data for this, but you must run
            chip_seq_cluster_heatmap() first, before this function.

            By default the yscale is locked to the maximum value - so the plots are all to the same
            scale.

        **Arguments**
            filename (Required)
                A base file name to save the images to. This function will save multiple files, one
                for each cluster/group. The file name will be modified in this manner:

                if filename=="base.png", the modified versions will be: "base_cid[1..n].png"
                i.e. "_cid<num>" will be inserted before the final filetype (in this case a png
                file).

            multi_plot (Optional, default=True) ONLY True IS IMPLEMENTED
                If True, plot all pileups on separate graphs, plotted sequentially horizontally as part of the same
                figure.

                If False then plot them all on the same graph, and with a legend.

        **Returns**
            The pileup_data as a dict in the form:
            {<cluster_id>: [array_data1, array_data2 ... array_dataN],
            <cluster_id>: [...],
            ...
            }
        """
        assert filename, "chip_seq_cluster_pileup(): you must provide a filename"
        assert self.__pileup_names, "chip_seq_cluster_pileup(): You must run chip_seq_cluster_heatmap() first"

        #print self.__pileup_data

        base_filename = ".".join(filename.split(".")[0:-1])

        num_plots = len(self.__pileup_data[1])
        if not "size" in kargs:
            kargs["size"] = (4*num_plots, 7)

        maxx = self.__pileup_data[1][0].shape[0]
        for cid in self.__pileup_data:
            this_filename = "%s_cid%s.png" % (base_filename, cid) # savefigure will modify png if needed.

            fig = self.draw.getfigure(**kargs)
            fig.suptitle("Group: %s #members: %s Membership: %s" % (cid, self.__pileup_group_sizes[cid], self.__pileup_groups_membership[cid-1]["id"]))
            # get the max x and y axes:
            maxy = max([a.max() for a in self.__pileup_data[cid]])
            miny = min([a.min() for a in self.__pileup_data[cid]])

            for cfig, data in enumerate(self.__pileup_data[cid]):
                ax = fig.add_subplot(1, len(self.__pileup_data[cid]), cfig+1)
                #print data
                x = numpy.arange(len(data))
                ax.plot(x, data)

                ax.set_xlim([0, maxx-2]) # -2 to trim off the unsightly tail due to binning.
                [t.set_visible(False) for t in ax.get_xticklabels()]
                ax.set_ylim([miny, maxy])
                if cfig >= 1: # nice bodge to blank labels on subsequent graphs.
                    [t.set_visible(False) for t in ax.get_yticklabels()]
                else:
                    ax.set_ylabel(self.__pileup_y_label)
                ax.set_title("%s (%s)" % (self.__pileup_names[cfig], self.__pileup_groups_membership[cid-1]["id"][cfig]))

                self.draw.do_common_args(ax, **kargs)

            self.draw.savefigure(fig, this_filename)
        return(self.__pileup_data)

    def genome_dist_radial(self, genome, layout, filename=None, randoms=None, **kargs):
        """
        **Purpose**
            Measure genome distributions of a list of genome coordinates relative to a list of TSS's.

            As seen in Hutchins et al., 2013 NAR Figure 1D. The version here is slightly generalised compared
            to the version used in that paper.

            Also, one disadvantage here is the lack of a key...

        **Arguments**
            genome (Required)
                A genome with a tss_loc key or loc key to annotate against

            layout (Required)
                You need to specify a tuple describing the layout arrangement i.e. the number of radial plot rows and columns

            randoms (Optional)
                A list of random peaks to treat as the background binding or binding pattern expected by chance alone.

            filename (Required)
                filename to save the radial plots to

        """

        assert genome[0], "genome_dist_radial: genome appears to be empty"
        assert "tss_loc" in list(genome.keys()), "genome_dist_radial: genome does not have a 'tss_loc' key"
        # check layout == len(self.linearData)

        annotation = genome

        res = {}

        for p in self.linearData:
            data, back, back_err, cats = p.genome_distribution(annotation, randoms, filename=None)

            res[p.name] = {"data": data, "back": back, "err": back_err}

        fig = self.draw.getfigure(**kargs)
        fig.subplots_adjust(0.02, 0.02, 0.97, 0.97, wspace=0.1, hspace=0.1)

        # Work out values for the histograms
        erad = 0.69813170079773 # each segment gets 0.69813170079773 (or thereabouts) rads
        eradh = erad / 2.0
        eradq = eradh / 2.0
        theta = numpy.arange(0.0, 2*numpy.pi, 2*numpy.pi/len(data))
        width = (numpy.pi/4)*len(data) # in rads?
        width = 0.5

        # colour for each segment
        colors = ["#FFF800", # (255, 248, 0)
            "#000E7C", # (0, 14, 177)
            "#001EFF", # (0, 30, 255)
            "#6275FF", # (98, 117, 255)
            "#B1BAFF", # (177, 186, 255)
            "#FFB7B1", # (255, 183, 177)
            "#FF6E62", # (255, 110, 98)
            "#FF1300", # (255, 19, 0)
            "#7C0900"] # (124, 9, 0)

        for i, k in enumerate(res): # ugh. random order...
            ax = fig.add_subplot(layout[0], layout[1], i+1, polar=True)

            if res[k]["back"]:
                axes[k].bar(theta-0.10, res[k]["back"], width=erad, bottom=0.0, alpha=0.8, ec="none", color="grey")
            ax.bar(theta, res[k]["data"], width=erad-0.20, bottom=0.0, alpha=0.9, ec="none", color=colors)
            ax.set_title(k, size=7)
            ax.set_xticks(theta-0.10)
            ax.set_xticklabels("")
            l = ax.get_ylim()
            #print k, ["%s%%" % i for i in range(l[0], l[1]+5, l[1]//len(axes[k].get_yticklabels()))]
            #print [str(t) for t in axes[k].get_yticklabels()]
            [t.set_fontsize(10) for t in ax.get_yticklabels()]
            #print ["%s%%" % (i*10, ) for i, t in enumerate(axes[k].get_yticklabels())]
            #print [t.get_text() for t in axes[k].get_yticklabels()]
            ax.set_yticklabels(["%s%%" % i for i in range(int(l[0]), int(l[1]+5), int(l[1]//len(ax.get_yticklabels())))][1:])

        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("genome_dist_radial: Saved '%s'" % actual_filename)

    def GO_heatmap(self, filename, p_value_limit=0.01, num_top=5, pvalue_key='pvalue',
            size=[8, 6], bracket=[1.3,4], row_cluster=True, col_cluster=False, # heatmap args
            heat_wid=0.15, cmap=cm.Reds, border=True, row_font_size=7,
            heat_hei='proportional', grid=True, ontology=None, draw_numbers_fmt='%.1f',
            draw_numbers=True, draw_numbers_threshold=2.0, draw_numbers_font_size=5, do_negative_log10=True,
            **kargs):
        '''
        **Purpose**
            Produce a heatmap of GO categories from a glglob of GO genelists (glgo's)

        **Arguments**
            filename (Required)
                filename to save the resulting heatmap to

            p_value_limit (Optional, default=0.01)
                minimum p-value to include in list

            pvalue_key (Optional, default='p_value')
                The key in your GO lists that contains some sort of significance
                result.

            num_top (Optional, default=5)
                Generally these heatmaps do not have much space to contain that
                many categories, so you need to take the top N from each list.

                GO_heatmap will 'fill in' categories in other lists even if they do not
                fulfill the 'p_value_limit' and 'num_top' criteria.

                However, this only occurs if your GO lists contain all GO terms. If
                you have already truncated the lists for some p-value then glbase cannot
                fill in the missing data.

            ontology (Optional, default=False)
                DAVID will give you a table containing all GO categories. Use this to specify using only
                a single ontology to use. Assumes the genelists have a 'ontology' key.

            This function will also accept all glbase heatmap arguments (see expression.heatmap).
            A few args have altered defaults:

            heat_hei (Optional, default='proportional')
                Sets the heatmap to a fixed y-size for each row.
                Set to a normal heat_wid value if you prefer.

            bracket (Optional, default=[1.3, 4.0])
                the bracket for the min and max of the heatmap. This sort of bracket
                assumes your data is -log10 transformed and so the p-value would
                range from 0.05 to 0.0001

            do_negative_log10 (Optional, default=True)
                By default convert the value in pvalue into the -log10()
                Set this to False if you don't want to convert

        **Returns**
            The resorted row names (as a list) and a heatmap in filename

        '''

        format = {'force_tsv': True, 'pvalue': 1, 'name': 0}

        main_cluster_membership = {}
        number_of_clusters = len(self)
        go_store = {}
        main_cluster_membership = {}
        cond_names_idx = {}

        for idx, go in enumerate(self.linearData):
            cond_names_idx[go.name] = idx
            if not go:
                config.log.warning("GO_heatmap: GO list '%s' was empty, skipping" % go.name)
                continue

            go.sort(pvalue_key)
            #go.reverse() # huh?

            #print(go)

            if ontology:
                this_ont = go.getRowsByKey('ontology', ontology)
                topN = this_ont[0:num_top]
            else:
                topN = go[0:num_top]

            for item in topN:
                if do_negative_log10:
                    if float(item[pvalue_key]) < p_value_limit:
                        if item['name'] not in go_store:
                            go_store[item['name']] = [-1] * (number_of_clusters)
                        go_store[item['name']][idx] = -math.log10(item['pvalue'])
                else:
                    if float(item[pvalue_key]) > -math.log10(p_value_limit): # i.e. 0.01
                        if item['name'] not in go_store:
                            go_store[item['name']] = [-1] * (number_of_clusters)
                        go_store[item['name']][idx] = item['pvalue']


        # fill in the holes:
        for go in self.linearData:
            for k in go_store:
                this_k = go.get(key='name', value=k, mode='lazy') # by default
                if this_k:

                    if do_negative_log10:
                        if float(item[pvalue_key]) < p_value_limit:
                            if item['name'] not in go_store:
                                go_store[item['name']] = [-1] * (number_of_clusters)
                            go_store[k][cond_names_idx[go.name]] = -math.log10(float(this_k[0]['pvalue']))
                    else:
                        if float(item[pvalue_key]) > -math.log10(p_value_limit): # i.e. 0.01
                            if item['name'] not in go_store:
                                go_store[item['name']] = [-1] * (number_of_clusters)
                            go_store[k][cond_names_idx[go.name]] = float(this_k[0]['pvalue'])

        newe = []

        for k in go_store:
            newe.append({'name': k.replace("~", ":"), 'conditions': go_store[k]}) # REPAIR DAVID GO names

        cond_names = sorted(zip(list(cond_names_idx.keys()), list(cond_names_idx.values())), key=itemgetter(1))
        cond_names = [i[0] for i in cond_names]

        goex = expression(loadable_list=newe, cond_names=cond_names)
        if len(goex) == 0:
            config.log.warning('GO list was empty, skipping')
            return(False)

        if heat_hei == 'proportional':
            heat_hei=0.012*len(goex)

        res = goex.heatmap(filename=filename, size=size, bracket=bracket,
            row_cluster=row_cluster, col_cluster=col_cluster,
            heat_wid=heat_wid, cmap=cmap, border=border,
            row_font_size=row_font_size, heat_hei=heat_hei, grid=grid,
            draw_numbers=draw_numbers, colbar_label='-log10(%s)' % pvalue_key,
            draw_numbers_threshold = -math.log10(p_value_limit),
            draw_numbers_fmt=draw_numbers_fmt,
            draw_numbers_font_size=draw_numbers_font_size)
        config.log.warning("GO_heatmap: Saved heatmap '%s'" % filename)
        return(reversed(res["reordered_rows"]))

    def measure_density(self, trks, peaks, norm_by_library_size=True, log=False,
        read_extend=0, pointify=True, expand=1000,
        **kargs):
        """
        **Purpose**
            get the seq tag density from the trks, and return as an expression object

        **Arguments**
            trks (Required)
                a list of tracks/flats

            peaks (Required)
                a list of peaks, a genelist containing a 'loc' key

            read_extend (Optional, default=200)
                read extend the sequence tags in the tracks by xbp

            norm_by_library_size (Optional, default=True)
                normalise the result by the [library size/1e6]

            log (Optional, default=False)
                log transform the resulting matrix

            pointify (Optional, default=True)
                convert the genomic locations to the center of the peak

            expand (Optional, default=500)
                expand the left and right flanks of the genomic coordiantes by <expand> base pairs
                Performed AFTER pointify

        **Returns**
            an expression object, with the conditions as the tag density from the tracks

        """
        assert isinstance(trks, list), 'measure_density: trks must be a list'
        assert 'loc' in list(peaks.keys()), 'measure_density: no loc key found in peaks'
        all_trk_names = [t["name"] for t in trks]
        assert len(set(all_trk_names)) == len(all_trk_names), 'track names are not unique. Please change the track["name"] to unique names'

        peaks = peaks.deepcopy()
        if pointify:
            peaks = peaks.pointify()
        if expand:
            peaks = peaks.expand('loc', expand)

        peaks.sort('loc')

        newl = []
        curr_chrom = None
        curr_data = None
        curr_n = 0
        for p in peaks:
            p["conditions"] = [0.0 for t in trks]

        all_chroms = len(set([i['chr'] for i in peaks['loc']])) * len(trks)
        all_sizes = [t.get_total_num_reads() / 1e6 for t in trks]

        for it, t in enumerate(trks):
            pb = progressbar(all_chroms)
            curr_chrom = None
            for p in peaks:
                p_loc = p['loc']
                if p_loc['chr'] != curr_chrom:
                    del curr_data
                    curr_data = t.get_array_chromosome(p_loc['chr'], read_extend=read_extend) # this is a numpy array
                    curr_chrom = p_loc['chr']
                    pb.update(curr_n)
                    curr_n += 1

                d = curr_data[p_loc['left']:p_loc['right']]

                if len(d) == 0: # fell off edge of array
                    p["conditions"][it] = 0 # Need to put a value in here
                    continue

                if norm_by_library_size:
                    p["conditions"][it] = numpy.average(d) / all_sizes[it]
                else:
                    p["conditions"][it] = numpy.average(d)

        expn = expression(loadable_list=peaks.linearData, cond_names=[t["name"] for t in trks])
        if log:
            expn.log(2, .1)

        return(expn)

    def redefine_peaks(self, super_set_of_peaks, list_of_flats, filename=None, Z_threshold=1.2,
        peak_window=200, lambda_window=5000, **kargs):
        """
        **Purpose**
            A strategy for re-calling peaks from some arbitrary set of flat files.

            It's a bit similar to the strategy published in Li et al., Cell Stem Cell 2017.
            However, this version uses a local lambda model drawn from a 10 kb window surrounding the peak.

            Breifly, the strategy is:

            Peak calling is conservative on any single ChIP-seq library. To get better sensitivity
            I pool information from other libraries by making a superset of peaks
            (all possible peaks in some set of ChIP-seq) and then 're-calling' the peaks
            in each library by modeling the enrichment. This allows me to rescue weak peaks.

            It then builds a model, constructs a Z-score and then only keeps those peaks
            that are greater than the threshold.

            This allows you to rescue weak peaks, and also to clean the false-negatives that are very common
            when cross-comparing peaks.

        **Arguments**
            super_set_of_peaks (Required)
                A genelist, containing a 'loc' key.

                This can be any list of peaks. I suggest you create this using the
                glglob.chip_seq_cluster() function. However, you can also ue bedtools
                or your favourite grouping strategy.

                I will pointify and expand the genomic locations to make them all uniform in size.

            list_of_flats (Required)
                a list of flats, make sure each flat has a unique 'name' slot as that
                will be used to store the result.

                For each flat the peaks will be recalled, and returned as a key in the returned dictionary.

            filename (Optional)
                the basefilename to save the model images to, one file for each flat.

            Z_threshold (Optional, default=1.2)
                The Z score cut off to call a peak or non-peak.

            peak_window (Optional, default=200)
                The window around the peak center, or summit to measure the peak enrichment.

            lambda_window (Optional, default=5000)
                The window around the

                Note that the actual window will be from the peak_window to the extent:

                left flank = -lambda_window <- lambda_window-peak_window
                right flank = peak_window -> lambda_window-peak_window

                i.e. with the default settings:

                    left flank      peak      right flank
                |-----4900 bp-----|-200bp-|-----4900 bp-----|

        **Returns**
            A dictionary, in the form:
            {flat[0]['name']: genelist(),
            flat[0]['name']: genelist(),
            ...
            flat[0]['name']: genelist()}

            where each genelist contains a 'loc' key, and several new keys, one for each flat:
            *_lam10 = the local lambda score
            *_lam10std = the local lambda standard deviation
            *_peak_score = the maximum peak height for this peak.

        """
        assert isinstance(list_of_flats, list), 'list_of_flats must be a list'
        assert False not in [i['name'] for i in list_of_flats], 'list_of_flats seems to not contain flats or flat-like objects'
        assert len(set([i['name'] for i in list_of_flats])) == len(list_of_flats), 'the "name" slots of the flats are not unique'
        assert 'loc' in super_set_of_peaks.keys(), 'super_set_of_peaks does not contain a "loc" (genomic location) key'

        peak_window = peak_window // 2
        lambda_window = lambda_window - peak_window

        rets = {f['name']: [] for f in list_of_flats}

        super_set_of_peaks = super_set_of_peaks.pointify().expand('loc', peak_window) # Make peaks symmetric
        super_set_of_peaks = [p['loc'].loc for p in super_set_of_peaks]

        # First I estimate the local background
        for f in list_of_flats:
            sam_name  = f['name'].replace('.flat', '')
            config.log.info('Doing {0}'.format(sam_name))

            this_chrom = None
            this_data = None
            prog = progressbar(len(super_set_of_peaks))

            for i, p in enumerate(super_set_of_peaks):
                p_loc_chrom = p['chr']
                p_loc_left = p['left']
                p_loc_rite = p['right']

                # Chrom cache version
                if p_loc_chrom != this_chrom:
                    this_data = f.get_array_chromosome(p_loc_chrom)
                    this_chrom = p_loc_chrom

                # I guess this is possible to be longer than the chrom:
                lambd_left = p_loc_left-lambda_window
                lambd_left = (lambd_left if lambd_left>0 else 0)
                lambd_rite = p_loc_rite+lambda_window
                lambd_rite = (lambd_rite if lambd_rite<len(this_chrom) else len(this_chrom))

                left_flank = this_data[lambd_left:p_loc_left-peak_window]
                rite_flank = this_data[p_loc_rite+peak_window:lambd_rite]
                peak_data = this_data[p_loc_left:p_loc_rite]

                # The above can fail, as peaks can come from dense data, and then be tested against a sparse flat
                if len(peak_data) == 0:
                    p['peak_score'] = 0 # fill the entries in, with 0 due to missing data in the array.
                    p['lam10'] = 0
                    p['lam10std'] = 0
                    continue

                all_lambda = left_flank + rite_flank
                mean_lambda = sum(all_lambda) / len(all_lambda)
                p['lam10'] = mean_lambda
                p['lam10std'] = pstdev(all_lambda)
                p['peak_score'] = max(peak_data) # should this be the max?
                prog.update(i)

            lam10 = [p['lam10'] for p in super_set_of_peaks]
            avg = mean(lam10)
            std = pstdev(lam10)
            config.log.info('Average background: %.3f' % avg)
            config.log.info('Average STDev: %.3f' % std)

            thresh = avg + (std * Z_threshold)
            config.log.info('Guessed threshold value of {1:.2f} (For a Z of {0})'.format(Z_threshold, thresh))

            # Plot the histogram:
            if filename:
                fig = self.draw.getfigure(**kargs)
                ax = fig.add_subplot(111)
                ax.hist(lam10, bins=50, range=[0,50], histtype='step', label='Background')
                ax.hist([p['peak_score'] for p in super_set_of_peaks], bins=50, range=[0,50], histtype='step', label='Peaks')
                ax.axvline(avg, ls=':', color='red')
                ax.axvline(avg+std, ls=':', color='green')
                ax.legend()
                self.draw.savefigure(fig, '{0}_{1}.png'.format(filename, sam_name))

            # redefine peaks:
            prog = progressbar(len(super_set_of_peaks))
            for i, p in enumerate(super_set_of_peaks):
                # First, filter on the global Z:
                if p['peak_score'] > thresh:
                    # Then filter on the localz:
                    if p['peak_score'] > (p['lam10'] + (p['lam10std']*Z_threshold)):
                        p_add = {'loc': location(chr=p['chr'], left=p['left'], right=p['right'])}
                        p_add['peak_height'] = p['peak_score']
                        try:
                            p_add['Z-score'] = ((p['peak_score'] - p['lam10']) / p['lam10std'])
                        except ZeroDivisionError:
                            p_add['Z-score'] = 100
                        rets[f['name']].append(p_add)

                prog.update(i)

        for f in rets:
            bed = genelist()
            bed.load_list(rets[f])
            bed.name = f
            rets[f] = bed

        config.log.info('New peak lengths:')
        for f in rets:
            config.log.info('    %s: %s peaks' % (f, len(rets[f])))

        return(rets)

    def chip_seq_heatmap(self,
        list_of_peaks,
        list_of_trks,
        filename:str = None,
        normalise=False,
        bins:int = 100,
        pileup_distance:int = 1000,
        cache_data=False,
        bracket=None,
        range_bracket=None,
        frames=True,
        row_labels=None,
        col_labels=None,
        read_extend:int = 200,
        imshow:bool = True,
        cmap=cm.plasma,
        log_pad=None,
        log=2,
        per_column_bracket=False,
        size=None,
        **kargs):
        """
        **Purpose**
            Draw heatmaps for the indicated list of peaks.

            peaks will be piled up on top of each other in proportional blocks (separated by a line).

            pileups will be plotted above, and each

        **Arguments**
            list_of_peaks (Required)
                A list of genelists of peaks from your ChIP-seq data to interrogate.
                The peaks will be stacked from bottom to top

            list_of_trks (Required)
                A list of trks to draw the sequence tag reads from to build the pileups.

            filename (Optional, default=None)
                If set to a string a heatmap & pileup will be saved to filename.

            normalise (Optional, default=False)
                Normalize the read pileup data within each library to the size of the
                library to assist in cross-comparison of ChIP-seq libraries.

            pileup_distance (Optional, default=1000)
                distance around the particular bin to draw in the pileup.

            read_extend (Optional, default=200)
                The size in base pairs to extend the read. If a strand is present it will expand
                from the 3' end of the read.

            row_labels (Optional, default= from the peak.name)
                row labels for the heatmaps;

            col_labels (Optional, default=from the trk['name'])
                column labels

            bins (Optional, default=100)
                number of bins to use for the pileup. Best to use conservative numbers (30-200) as
                large numbers of bins can consume huge amounts of memory.

            log (Optional, default=2)
                Use logarithms for the heatmap. Possible options are 2 and 10.

            cmap (Optional, default=matplotlib.cm.plasma)
                A colour map for the heatmap.

            ######## Bracket system:

            There area bunch of args here.

            per_column_bracket sets a bracket for each column (track) in the heatmaps. If this is True,
            you should use the range_bracket system.

            If per_column_bracket is False, then you can use the range_bracket system (recommended),
            but you cna also set your own bracket with the bracket option.

            per_column_bracket (Optional, default=False)
                Have a bracket for each column, or (when False) a single bracket for all of the data

            range_bracket (Optional, default=[0.4, 0.9], exclusive with range_bracket)
                chip_seq_heatmap can have column-wise (track-wise) specific brackets. So, only range_bracket is valid

                Okay, hold your hats, this is complicated.

                range_bracket will bracket the range of values as a fraction between [min(data), max(data)]
                i.e. If range_bracket=0.5 (the default) then the data is bracketed as:
                [max(data)*range_bracket[0], max(data)*range_bracket[1]]. The practical upshot of this is it allows you to shift
                the colour bracketing on the heatmap around without having to spend a long time finding
                a suitable bracket value.

                Bracketing is performed AFTER log.

                Typical bracketing would be something like [0.4, 0.9]

                By default glbase attempts to guess the best range to draw based on the
                median and the stdev. It may not always succeed.

            cache_data (Optional, default=False)
                cache the pileup data into the file specified in cache_data. This speeds up reanalysis.
                Note that storage of data is AFTER normalisation, resolution, pileup_distance,
                bins, but before heatmap drawing.

                This allows you to store the slow part of chip_seq_heatmap()
                and so iterate through different heatmap drawing options without having to
                do the whole pileup again.

                note that if cache_data file does not exist then it will be created and
                pileup data generated. If the file does exist, data will be read from that
                file and used for heatmap drawing.

            frames (Optional, default=True)
                Draw black frames around the heatmaps and category maps.

            imshow (Optional, default=False)
                Embed the heatmap as an image inside a vector file. (Uses matplotlib imshow
                to draw the heatmap part of the figure. Allows very large matrices to
                be saved as a reasonably sized svg/pdf, with the heatmap part as a raster image
                and all other elements as vectors).

        **Returns**
            Returns None
        """
        assert not (range_bracket and bracket), "You can't use both bracket and range_bracket"
        assert False not in ['loc' in gl.keys() for gl in list_of_peaks], 'At least one of your peak data (list_of_peaks) does not contain a "loc" key'

        total_rows = 0

        # Get the size of each library if we need to normalize the data.
        if normalise:
            # get and store the read_counts for each library to reduce an sqlite hit.
            read_totals = [trk.get_total_num_reads()/float(1e6) for trk in list_of_trks]

        # I will need to go back through the chr_blocks data and add in the pileup data:
        bin_size = int((pileup_distance+pileup_distance) / bins)
        #block_len = pileup_distance+pileup_distance # get the block size
        data = None

        # Populate the datastores:
        matrix = {}
        pileup = {}
        for tindex, _ in enumerate(list_of_trks):
            matrix[tindex] = {} # Populate the final matrix
            pileup[tindex] = {} # Populate the final matrix
            for plidx, peaklist in enumerate(list_of_peaks):
                matrix[tindex][plidx] = {} # Populate the final matrix
                pileup[tindex][plidx] = {} # Populate the final matrix
                for pindex, peak in enumerate(peaklist):
                    matrix[tindex][plidx][pindex] = None
                    pileup[tindex][plidx][pindex] = None

        # Populate the order data so I can use the chromosome cache system;
        porder = {}
        for plidx, peaklist in enumerate(list_of_peaks):
            porder[plidx] = {} # make an index hitter so that order is preserved:
            for pindex, peak in enumerate(peaklist):
                p_loc_chrom = peak['loc']['chr']
                if p_loc_chrom not in porder[plidx]:
                    porder[plidx][p_loc_chrom] = []
                porder[plidx][p_loc_chrom].append(pindex)

        # sort out cached_data
        if cache_data and os.path.isfile(cache_data): # reload previously cached data.
            oh = open(os.path.realpath(cache_data), "rb")
            matrix = pickle.load(oh)
            oh.close()
            config.log.info("chip_seq_heatmap: Reloaded previously cached pileup data: '%s'" % cache_data)
            # sanity check the matrix data
            assert isinstance(matrix, dict), '{0} does not match the expected data, suggest you rebuild'.format(cache_data)
            assert isinstance(matrix[0], dict), '{0} does not match the expected data, suggest you rebuild'.format(cache_data)
            assert isinstance(matrix[0][0], (numpy.ndarray, numpy.generic)), '{0} does not match the expected data, suggest you rebuild'.format(cache_data)
            assert len(matrix) == len(list_of_trks), '{0} does not match the expected data, suggest you rebuild'.format(cache_data)
            assert len(matrix[0]) == len(list_of_peaks), '{0} does not match the expected data, suggest you rebuild'.format(cache_data)
            for it, t in enumerate(list_of_trks):
                for ip, p in enumerate(list_of_peaks):
                    #print(matrix[it][ip].shape, (len(p), bins))
                    assert matrix[it][ip].shape == (len(p), bins), '{0} does not match the expected data, suggest you rebuild'.format(cache_data)
        else:
            # No cached data, so we have to collect ourselves.
            config.log.info('chip_seq_heatmap: Collecting pileup data...')
            p = progressbar(len(list_of_trks))
            # New version that grabs all data and does the calcs in memory, uses more memory but ~2-3x faster
            for tindex, trk in enumerate(list_of_trks):
                for plidx, peaklist in enumerate(list_of_peaks):
                    for chrom in porder[plidx]:
                        # The chr_blocks iterates across all chromosomes, so this only hits the db once per chromosome:
                        data = trk.get_array_chromosome(chrom, read_extend=read_extend) # This will use the fast cache version if available.

                        for pidx, peak in enumerate(porder[plidx][chrom]): # peak is the index to look into
                            left = peaklist.linearData[peak]['loc']['left']
                            rite = peaklist.linearData[peak]['loc']['right']
                            cpt = (left + rite) // 2
                            left = cpt - pileup_distance
                            rite = cpt + pileup_distance

                            # It's possible to ask for data beyond the edge of the actual data. truncate...
                            if rite > len(data):
                                rite = len(data)
                            if left > len(data):
                                left = len(data)

                            dd = data[left:rite]

                            if normalise:
                                pil_data = [av/read_totals[tindex] for av in dd]

                            # Fill in the matrix table:
                            #pileup[tindex][plidx][pidx] += pil_data
                            matrix[tindex][plidx][peak] = [sum(dd[i:i+bin_size]) for i in range(0, len(dd), bin_size)]
                            #print(matrix[tindex][plidx][pidx])
                            #chr_blocks[chrom][block_id]["pil"][tindex] = [sum(dd[i:i+bin_size]) for i in range(0, len(dd), bin_size)] #pil_data = utils.bin_sum_data(dd, bin_size)
                p.update(tindex)

            # convert to numpy arrays;
            for tindex, _ in enumerate(list_of_trks):
                for plidx, peaklist in enumerate(list_of_peaks):
                    twoD_list = []
                    for pindex, _ in enumerate(peaklist): # preserve original order;
                        twoD_list.append(matrix[tindex][plidx][pindex])
                    matrix[tindex][plidx] = numpy.array(twoD_list)

            if cache_data: # store the generated data for later.
                oh = open(cache_data, "wb")
                pickle.dump(matrix, oh, -1)
                oh.close()
                config.log.info("chip_seq_heatmap: Saved pileup data to cache file: '{0}'".format(cache_data))

        colbar_label = "Tag density"


        if log:
            if not log_pad:
                log_pad = 0.1

            if not range_bracket: # suggest reasonable range;
                range_bracket = [0.6, 0.9]

            for tindex, _ in enumerate(list_of_trks):
                for plidx, peaklist in enumerate(list_of_peaks):

                    if log == 2:
                        matrix[tindex][plidx] = numpy.log2(matrix[tindex][plidx]+log_pad)
                        colbar_label = "Log2(Tag density)"
                    elif log == 10:
                        matrix[tindex][plidx] = numpy.log10(matrix[tindex][plidx]+log_pad)
                        colbar_label = "Log10(Tag density)"
                    else:
                        raise AssertionError('log={0} not found'.format(log))

        else:
            if not range_bracket: # suggest reasonable range;
                range_bracket = [0.0, 0.1]

        if normalise:
            colbar_label = "Normalised %s" % colbar_label

        if not row_labels:
            row_labels = [p.name for p in list_of_peaks]

        if not col_labels:
            col_labels = [t['name'] for t in list_of_trks]

        bracket = None
        brackets = None

        if per_column_bracket:
            # Suggest brackets:

            t_stats = []
            brackets = []
            for tindex, _ in enumerate(list_of_trks): # I can have track-wise brackets;
                for plidx, peaklist in enumerate(list_of_peaks):
                    tab_max = max([tab.max() for tab in matrix[tindex][plidx]]) # need to get new tab_max for log'd values.
                    tab_min = min([tab.min() for tab in matrix[tindex][plidx]])
                    #tab_median = numpy.median([numpy.median(tab) for tab in list_of_tables])
                    tab_mean = mean([numpy.average(tab) for tab in matrix[tindex][plidx]])
                    tab_std = numpy.std(numpy.array([tab for tab in matrix[tindex][plidx]]))
                    t_stats.append((tab_max, tab_min, tab_mean, tab_std))

                config.log.info('chip_seq_heatmap: trk={0} min={1:.2f}, max={2:.2f}, mean={3:.2f}, stdev={4:.2f}'.format(list_of_trks[tindex]['name'], tab_min, tab_max, tab_mean, tab_std))

                tab_range = tab_max - tab_min
                top = tab_min + tab_range*range_bracket[0]
                bot = tab_min + tab_range*range_bracket[1]
                brackets.append([top, bot])
                config.log.info("chip_seq_heatmap: trk={0}, suggested bracket=({1:.2f}, {2:.2f})".format(list_of_trks[tindex]['name'], brackets[tindex][0], brackets[tindex][1]))
        else:
            # Suggest brackets:
            if bracket:
                pass # USe the arg
            else: # Guess a bracket for all heatmaps;
                t_stats = []
                brackets = []
                for tindex, _ in enumerate(list_of_trks): # I can have track-wise brackets;
                    for plidx, peaklist in enumerate(list_of_peaks):
                        tab_max = max([tab.max() for tab in matrix[tindex][plidx]]) # need to get new tab_max for log'd values.
                        tab_min = min([tab.min() for tab in matrix[tindex][plidx]])
                        #tab_median = numpy.median([numpy.median(tab) for tab in list_of_tables])
                        tab_mean = mean([numpy.average(tab) for tab in matrix[tindex][plidx]])
                        tab_std = numpy.std(numpy.array([tab for tab in matrix[tindex][plidx]]))
                        t_stats.append((tab_max, tab_min, tab_mean, tab_std))

                    config.log.info('chip_seq_heatmap: trk={0} min={1:.2f}, max={2:.2f}, mean={3:.2f}, stdev={4:.2f}'.format(list_of_trks[tindex]['name'], tab_min, tab_max, tab_mean, tab_std))

                    tab_range = tab_max -tab_min
                    top = tab_min + tab_range*range_bracket[0]
                    bot = tab_min + tab_range*range_bracket[1]
                    brackets.append([top, bot])

                bracket = [max([b[0] for b in brackets]), max([b[1] for b in brackets])]
                brackets = None
                config.log.info("chip_seq_heatmap: suggested bracket=({0:.2f}, {1:.2f})".format(bracket[0], bracket[1]))

        if filename:
            real_filename = self.draw.grid_heatmap(
                data_dict_grid=matrix,
                filename=filename,
                row_labels=row_labels,
                col_labels=col_labels,
                colbar_label=colbar_label,
                imshow=imshow,
                size=size,
                colour_map=cmap,
                bracket=bracket,
                brackets=brackets,
                frames=frames
                )

        config.log.info("chip_seq_heatmap: Saved overlap heatmap to '{0}'".format(real_filename))
        return None

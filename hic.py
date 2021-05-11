'''

hic.py

Analysis for HiC data.

TODO:
-----
. merge_hiccys could be adapted to measure variance, and so extract DE?

'''

import pickle, numpy, math, gzip
from operator import itemgetter
from shutil import copyfile

import matplotlib.cm as cm
import scipy
from scipy import ndimage, interpolate
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from .genelist import genelist
from . import config
from .draw import draw
from .progress import progressbar
from .location import location
from .format import minimal_bed

if config.STATSMODELS_AVAIL:
    from statsmodels.nonparametric.smoothers_lowess import lowess
else:
    raise AssertionError('Asked for a hic, but statsmodels is not available')

import matplotlib.pyplot as plot

if config.H5PY_AVAIL:
    import h5py
else:
    raise AssertionError('Asked for a hic, but h5py is not available')


def reshap_mats2(mat, dimX, dimY):
    '''
    Sometimes the matrices are slightly off, reshape them to the same sizes.
    reshapes matB to the dimensions provided

    This deos not work!

    '''
    interpolator_func = interpolate.interp2d(range(mat.shape[0]), range(mat.shape[1]), mat, kind='linear')
    return interpolator_func(dimX, dimY)


def reshap_mats(mat, dimX, dimY):
    '''
    Sometimes the matrices are slightly off, reshape them to the same sizes.
    reshapes matB to the dimensions provided

    NOTE: All matrices MUST be square len(x) == len(y)

    '''
    X = numpy.linspace(0, mat.shape[0], mat.shape[0])
    f = interpolate.RectBivariateSpline(X, X, mat)
    Xnew = numpy.linspace(0, mat.shape[0], dimX)
    return f(Xnew, Xnew)

def reshap_mats_irregular(mat, dimX, dimY, t):
    '''
    Non-square matrix version
    '''
    X = numpy.linspace(0, mat.shape[0], mat.shape[0])
    Y = numpy.linspace(0, mat.shape[0], mat.shape[1])
    f = interpolate.RectBivariateSpline(X, Y, mat)
    Xnew = numpy.linspace(0, mat.shape[0], dimX)
    Ynew = numpy.linspace(0, mat.shape[1], dimY)

    return f(Xnew, Ynew)

def merge_hiccys(new_hic_filename, name, *hics):
    '''
    **Purpose**
        Take the mean values for two or more hic objects.
        Save the result as a new hiccy File

        This algorithm assumes you know what you are doing, and you are trying to merge
        hiccy objects with the same resolution

    **Arguments**
        new_hic_filename (Required)
            filename to save the new hiccy object into

        name (Required)
            name of the new hiccy object

        hics (Required)
            A list of filenames to merge

    **Returns**
        None
    '''
    # TODO: Merge OE and AB matrices (regenerate, surely?)
    assert new_hic_filename, 'You must specify a filename in new_hic_filename'
    assert name, 'You must specify a name'
    assert hics, 'a list of hics was not specified'
    assert isinstance(hics, tuple), 'hics must be a tuple of >=2'
    assert (len(hics) > 1), 'hics must be >=2 length'

    # For the first one to merge, do an OS copy for speed and setup
    #copyfile(hics[0], new_hic_filename)

    # I bind all the hics:
    h0 = hic(filename=hics[0], name=hics[0])
    hics = [hic(filename=f, name=f) for f in hics[1:]]

    # Do a setup;
    newhic = hic(filename=new_hic_filename, name=name, new=True, _readplus=False)

    newhic.hdf5_handle.attrs['name'] = name
    newhic.hdf5_handle.attrs['inter_chrom_only'] = h0.hdf5_handle.attrs['inter_chrom_only']
    newhic.hdf5_handle.attrs['OE'] = False
    newhic.hdf5_handle.attrs['AB'] = False
    newhic.hdf5_handle.attrs['version'] = h0.hdf5_handle.attrs['version']
    newhic.hdf5_handle.attrs['num_bins'] = h0.hdf5_handle.attrs['num_bins']
    newhic.hdf5_handle.attrs['bin_size'] = h0.hdf5_handle.attrs['bin_size']
    newhic.all_chrom_names = hics[0].all_chrom_names # Made into a set later
    newhic.draw = draw()

    for chrom in newhic.all_chrom_names:
        data = numpy.array(h0.mats[chrom])
        for h in hics:
            newdata = h.mats[chrom]

            #if newdata.shape != data.shape:
            #    newdata = reshap_mats(newdata, data.shape[0], data.shape[1])
            data += newdata

        data /= (len(hics)+1)

        grp = newhic.hdf5_handle.create_group('matrix_{}'.format(chrom))
        grp.create_dataset('mat', data.shape, dtype=numpy.float32, data=data)
        config.log.info('Added chrom=%s to table' % chrom)

    # save the bin data as an emulated dict: Just copy from h0
    grp = newhic.hdf5_handle.create_group('bin_lookup')
    for chrom in newhic.all_chrom_names:
        grp = newhic.hdf5_handle.create_group('bin_lookup/chrom_%s' % chrom)
        flat_bin = numpy.array(h0.hdf5_handle['bin_lookup/chrom_%s/bins' % chrom][()])
        grp.create_dataset('bins', (flat_bin.shape), dtype=int, data=flat_bin, chunks=True, compression='lzf')

    dat = [str(n).encode("ascii", "ignore") for n in newhic.all_chrom_names]
    newhic.hdf5_handle.create_dataset('all_chrom_names', (len(newhic.all_chrom_names), 1), 'S10', dat)

    newhic._OEmatrix()
    newhic._ABcompart()
    newhic.close()
    config.log.info('Merged {0} matrices'.format(len(hics)+1,))

# I'm not using the base_genelist class, so, you need to add in defs as needed.
class hic:
    def __init__(self, filename=None, name='', new=False, inter_chrom_only=True, _readplus=False):
        """
        **Purpose**
            store for HiC data.

        **Arguments**
            filename (Required)
                hic is an on-disk array. So, you need to provide
                a filename to save the data into.

                I suggest a '.hiccy' extension, so as not to confuse with JuiceBox's .hic

            new (Optional, default=False)
                set this to True if you want to make a new hiccy

            name (Optional)
                a name to use as a label, etc

            inter_chrom_only (Required, default=True)
                At the moment only inter-chromsomal contacts are recorded.

        """
        assert inter_chrom_only, 'inter_chrom_only != True, only inter-chromosomal links are supported'
        self.version = '1.1'
        self.readonly = True
        self.pca_valid = False
        self.tsne_trained = False
        self.filename = filename
        # These need to also be packaged into the hiccy
        if new:
            self.readonly = False
            #try:
            self.hdf5_handle = h5py.File(filename, 'w')  # This should be x to stop overwriting an exisiting file

            # h5py hierarchy (Made at load_time)
            # Things stored:
            # str name
            # matrix # A numpy style matrix
            # chrom_edges slice locations for all the chromsome names
            self.hdf5_handle.attrs['name'] = name
            self.hdf5_handle.attrs['inter_chrom_only'] = inter_chrom_only
            self.hdf5_handle.attrs['OE'] = False
            self.hdf5_handle.attrs['version'] = self.version
            self.all_chrom_names = [] # Made into a set later
            self.draw = draw()
        else: # old
            if _readplus:
                self.hdf5_handle = h5py.File(filename, 'r+')
                self.readonly = False
            else:
                self.hdf5_handle = h5py.File(filename, 'r')
            # TODO: fetch it back out as a set:
            self.tad_lookup = None

            dat = self.hdf5_handle['all_chrom_names'][()]
            dat = [i[0] for i in dat]# flatten
            self.all_chrom_names = [n.decode("ascii", "ignore") for n in dat]

            # Rescue bin_lookup from disk:
            # save the bin data as an emulated dict:
            self.bin_lookup_by_binID = {}
            self.bin_lookup_by_chrom = {}
            for chrom in self.all_chrom_names:
                flat_bin = self.hdf5_handle['bin_lookup/chrom_%s/bins' % chrom][()]
                self.bin_lookup_by_chrom[chrom] = [(row[0], row[1], row[2], row[3]) for row in flat_bin]

                for row in flat_bin:
                    self.bin_lookup_by_binID[row[0]] = (chrom, row[1], row[2], row[3]) # row[3] == newBinID

            #print(self.bin_data)

            self.mats = {}
            self.OE = {}
            self.AB = {}
            for chrom in self.all_chrom_names:
                self.mats[chrom] = self.hdf5_handle['matrix_%s/mat' % chrom]
                self.OE[chrom] = self.hdf5_handle['OE_{}/OE'.format(chrom)]
                self.AB[chrom] = self.hdf5_handle['AB_{}/AB'.format(chrom)]

            self.draw = draw()
            config.log.info('Bound "%s" Hic file' % filename)
        return None

    def visit(self):
        self.hdf5_handle.visit(print)

    def __len__(self):
        return self['num_bins']

    def __getitem__(self, index):
        """
        Confers:

        a = hic["condition_name"]

        """
        return self.hdf5_handle.attrs[index]

    def keys(self):
        return self.hdf5_handle.keys()

    def _optimiseData(self):
        pass

    def close(self):
        self.hdf5_handle.flush()
        self.hdf5_handle.close()

    def __find_binID_spans(self, loc):
        """
        (Internal)

        Fing the binIDs for the span loc.

        Returns a tuple containing the local binIDs for these coordiantes,
        and the new updated genomic coordinates taken from the binspans

        """
        assert self.readonly, 'To __find_binID_spans, this must be read-only, set new=False'

        # These are in global bins.
        binLeft = 1e50
        binRight = -1
        locLeft = 0
        locRight = 0
        mostLeft = 1e50
        mostRight = -1

        chrom = loc.loc['chr']
        if 'chr' not in chrom:
            chrom = 'chr{0}'.format(chrom)

        #print(self.bin_lookup_by_chrom[loc['chr']])
        for bin in self.bin_lookup_by_chrom[chrom]:
            # bin = (binID, left, right)
            #print(bin, loc)
            if bin[2] > loc['left'] and bin[0] < binLeft:
                binLeft = bin[0]
                locLeft = bin[1]
                #print('bL', binLeft, binRight, loc)
            if loc['right'] > bin[1] and bin[0] > binRight:
                binRight = bin[0]
                locRight = bin[2]
                #print('bR', binLeft, binRight, loc)

            # To work out the local bin positions:
            if bin[0] < mostLeft:
                mostLeft = bin[0]
            if bin[0] > mostRight:
                mostRight = bin[0]
                mostRightLoc = bin[2]

        # in case locations are asked for off the edge of the chromsome:
        if binLeft < 0:
            binLeft = 0
            locLeft = 0
        if binRight > mostRight:
            binRight = mostRight
            locLeft = mostRightLoc

        newloc = location(chr=chrom, left=locLeft+1, right=locRight-1) # the +1/-1 stops the bins from overlapping by 1bp, and solves a lot of problems

        binLeft = binLeft - mostLeft
        binRight = binRight - mostLeft
        mostRight = self.bin_lookup_by_binID[mostRight][3] # Convert to newBinID:
        mostLeft = self.bin_lookup_by_binID[mostLeft][3]
        #print (binLeft, binRight, newloc, mostLeft, mostRight)
        assert binRight-binLeft > 10, 'the genome view (loc) is too small, and contains < 10 bins'

        return (binLeft, binRight, newloc, mostLeft, mostRight)

    def __quick_find_binID_spans(self, loc_chrom=None, loc_left=None, loc_rite=None, do_assert_check=True):
        """
        (Internal)

        Fing the binIDs for the span loc.

        Returns a tuple containing the local binIDs for these coordinates,

        """
        # These are in global bins.
        binLeft = 1e50
        binRight = -1
        mostLeft = 1e50
        mostRight = -1

        #print(self.bin_lookup_by_chrom[loc['chr']])
        for bin in self.bin_lookup_by_chrom[loc_chrom]:
            # bin = (binID, left, right)
            #print(bin, loc)
            if bin[2] > loc_left and bin[0] < binLeft:
                binLeft = bin[0]
            if loc_rite > bin[1] and bin[0] > binRight:
                binRight = bin[0]

            # To work out the local bin positions:
            if bin[0] < mostLeft:
                mostLeft = bin[0]
            if bin[0] > mostRight:
                mostRight = bin[0]

        # in case locations are asked for off the edge of the chromsome:
        if binLeft < 0:
            binLeft = 0
        if binRight > mostRight:
            binRight = mostRight

        binLeft = binLeft - mostLeft
        binRight = binRight - mostLeft
        if do_assert_check: assert binRight-binLeft > 2, 'the genome view (loc) is too small, and contains < 3 bins'

        return (binLeft, binRight)

    def __find_binID_chromosome_span(self, chrom):
        """
        (Internal)

        Get all of the binID (mostLeft and mostRight) for the indicated chromosome
        """
        assert self.readonly, 'To __find_binID_chromosome_span, this must be read-only, set new=False'
        # These are in global bins.
        locLeft = 0
        locRight = 0
        mostLeft = 1e50
        mostRight = -1

        for bin in self.bin_lookup_by_chrom[chrom]:
            # To work out the local bin positions:
            if bin[0] < mostLeft:
                mostLeft = bin[0]
            if bin[0] > mostRight:
                mostRight = bin[0]
                mostRightLoc = bin[2]

        #newloc = location(chr=loc['chr'], left=mostLeft+1, right=mostRight-1) # the +1/-1 stops the bins from overlapping by 1bp, and solves a lot of problems
        mostRight = self.bin_lookup_by_binID[mostRight][3] # Convert to newBinID:
        mostLeft = self.bin_lookup_by_binID[mostLeft][3]

        return (mostLeft, mostRight)

    def _OEmatrix(self):
        '''

        Call after loading mat to build the OE matrices;

        '''
        fig = plot.figure()

        for axidx, chrom in enumerate(self.all_chrom_names):
            ax = fig.add_subplot(5, 5, axidx+1)
            ax.set_title(chrom, fontsize=6)
            ax.tick_params(axis='both', labelsize=6)

            # Simple enough, take the diagonal mean (E), then for each point (O) O/E

            mats = self.hdf5_handle['matrix_{}/mat'.format(chrom)]
            grp = self.hdf5_handle.create_group('OE_{}'.format(chrom))
            #with numpyp.errstate(divide='ignore', invalid='ignore'):

            # get diagonal means;
            means = []
            flipped = numpy.array(mats) # numpy.fliplr(mats)
            for d in range(flipped.shape[0]):
                s = flipped.diagonal(offset=d)
                s = s[s > 0] # filter no datas
                if s.shape[0] > 0:
                    m = numpy.sum(s) / s.shape[0]
                else:
                    m = 0
                means.append(m)
                #means.append(d)
            oldmeans = numpy.array(means) # [::-1]

            # smooth means;
            fit = lowess(oldmeans, numpy.arange(len(means)), frac=0.04, it=6, is_sorted=True)

            ax.plot(means, lw=0.3, alpha=0.5, c='red')
            # The first ~2 points are a bad fit, so just replace them with the raw:
            means = fit[:,1]
            means[0] = oldmeans[0]
            means[1] = oldmeans[1]
            #means[2] = oldmeans[2]
            #print(means[0:4], oldmeans[0:4])
            ax.plot(fit[:,0], means, lw=0.6, alpha=0.2, c='black')
            sliding_means = means #[::-1] # initial;
            full_window = numpy.concatenate((means[::-1], means[1:], means[::-1], means[1:]), axis=None) # for sliding window;

            # Now get O/E for each point;
            newmat = numpy.array(mats)
            for x in range(mats.shape[0]):
                for y in range(mats.shape[0]):
                    with numpy.errstate(divide='ignore', invalid='ignore'):
                        newmat[x,y] = newmat[x,y] / sliding_means[y] # actual O/E calc;
                        #newmat[x,y] = sliding_means[y]
                        #newmat[x,y] = y
                t = mats.shape[0] - x
                sliding_means = full_window[t-2:t+mats.shape[0]+1] # increment;
                # inc
            newmat[numpy.isnan(flipped)] = 0
            newmat[numpy.isinf(flipped)] = 0
            newmat[flipped == 0] = 0 # remove no data

            newmat = numpy.corrcoef(newmat)

            grp.create_dataset('OE', newmat.shape, dtype=numpy.float32, data=newmat)
        config.log.info('Calculated O/E data')

        fig.savefig('OE_decay_{}.pdf'.format(self['name']))

        self.hdf5_handle.attrs['OE'] = True
        return

    def _ABcompart(self):
        '''

        Build A/B compartments;

        '''

        # From: https://github.com/dekkerlab/cworld-dekker/blob/master/scripts/python/getEigenVectors.py

        numPCs = 3

        """
        performs eigen vector analysis, and returns 3 best principal components
        result[0] is the first PC, etc
        """
        '''
        for chrom in self.all_chrom_names:
            newmat = self.hdf5_handle['OE_{}/OE'.format(chrom)]
            with numpy.errstate(divide='ignore', invalid='ignore'):
                A = numpy.corrcoef(newmat) # Pearson Corr of normalised distance

            A = numpy.array(A)
            M = (A - numpy.mean(A.T, axis=1)).T
            covM = numpy.dot(M, M.T)
            latent, coeff = scipy.sparse.linalg.eigsh(covM, numPCs)
            #return (np.transpose(coeff[:,::-1]),latent[::-1])
            pc1 = numpy.transpose(coeff[:,::-1])

            grp = self.hdf5_handle.create_group('AB_{}'.format(chrom))
            grp.create_dataset('AB', pc1.shape, dtype=numpy.float32, data=pc1)

        config.log.info('Calculated A/B compartments')
        self.hdf5_handle.attrs['AB'] = True
        return
        '''

        for chrom in self.all_chrom_names:
            m = numpy.array(self.hdf5_handle['OE_{}/OE'.format(chrom)])

            grp = self.hdf5_handle.create_group('AB_{}'.format(chrom))
            with numpy.errstate(divide='ignore', invalid='ignore'):
                #m = numpy.corrcoef(m)
                #m = numpy.log2(m)
                m[numpy.isnan(m)] = 0
                m[numpy.isinf(m)] = 0

            w, v = numpy.linalg.eig(m)
            if hasattr(v, 'mask'):
                v.mask = False
            e = numpy.real(v[:, 0]) # first PC

            grp.create_dataset('AB', e.shape, dtype=numpy.float32, data=e)
        config.log.info('Calculated A/B compartments')

        self.hdf5_handle.attrs['AB'] = True
        return


    def load_hicpro_matrix(self, matrix_filename, bed_file):
        """
        **Purpose**
            Load a .matrix file output from HiC-Pro

            bin1, bin2, norm_score

        **Arguments**
            filename (Required)
                the Filename to load.
        """
        assert not self.readonly, 'To load, this must be read-only, set new=True'

        # First go through the bed file and get all the bins, chroms and sizes.
        bins = []
        bin_size = None
        self.all_chrom_names = []
        oh = open(bed_file, 'r')
        for lin in oh:
            # format = chr, left, right, bin#
            lin = lin.strip().split('\t')

            chr = lin[0].replace('chr', '')
            try:
                chr = int(chr)
            except ValueError:
                pass # It's chrX chrVII etc.

            bins.append((chr, int(lin[1]), int(lin[2]), int(lin[3])))
            if not bin_size:
                bin_size = int(lin[2]) - int(lin[1])
            if chr not in self.all_chrom_names:
                self.all_chrom_names.append(chr)

        oh.close()
        self.all_chrom_names = set(self.all_chrom_names)
        self.hdf5_handle.attrs['num_bins'] = len(bins)
        self.hdf5_handle.attrs['bin_size'] = bin_size

        config.log.info('Found {0} bins, each bin = {1} bp'.format(self['num_bins'], self['bin_size']))
        # the bins are not sorted,
        bin_order = sorted(bins, key=lambda v: (isinstance(v[0], str), v[0], v[1])) # Sort chroms by str, then int
        bin_lookup = {} # python 3: no need OrderedDict

        for newbinID, bin in enumerate(bin_order):
            # Basically, bin_lookup[oldbinID] = (chrom, left, right, newbinID)
            bin_lookup[bin[3]] = (bin[0], bin[1], bin[2], newbinID) # left, right, bin#

        # Get the edges of the chorms, and so the matrix sizes
        self.chrom_edges = {}
        for chrom in self.all_chrom_names:
            self.chrom_edges[chrom] = [1e20,-1]

            for bin in bin_lookup:
                if bin_lookup[bin][0] != chrom:
                    continue

                # fill in:
                if bin_lookup[bin][3] < self.chrom_edges[chrom][0]:
                    self.chrom_edges[chrom][0] = bin_lookup[bin][3]
                if bin_lookup[bin][3] > self.chrom_edges[chrom][1]:
                    self.chrom_edges[chrom][1] = bin_lookup[bin][3]
        # Have to then convert them all to the new bins:
        # Although not strictly necessary in the matrix-based, I preserve this fiddly step so that
        # when I implement the inter-chrom interactions it is easier.
        matrices = {}
        bins = {} # The bin -> matrix convertor
        for chrom in self.all_chrom_names:
            # Now I know how big the arrays need to be, and all of the chromosomes:
            size =  self.chrom_edges[chrom][1]-self.chrom_edges[chrom][0]
            matrices[chrom] = numpy.zeros((size+1, size+1), dtype='float32')
            #bins[k] =
        # TODO: Support inter-chromosomal contacts with a sparse array:

        oh = open(matrix_filename, 'r')
        p = progressbar(self['num_bins'] * self['num_bins']) # expected maximum
        for idx, lin in enumerate(oh):
            lin = lin.strip().split('\t')
            # First check the two bins are on the same chromosome:
            bin1 = bin_lookup[int(lin[0])]
            bin2 = bin_lookup[int(lin[1])]

            if bin1[0] == bin2[0]:
                chrom = bin1[0]
                chrom_edge = self.chrom_edges[chrom][0]
                x = bin1[3]-chrom_edge
                y = bin2[3]-chrom_edge

                matrices[chrom][x,y] = float(lin[2])
                matrices[chrom][y,x] = float(lin[2])
            else:
                # TODO: Support for interchrom with a sparse array:
                pass

            p.update(idx)
        p.update(self['num_bins'] * self['num_bins']) # the above is unlikely to make it to 100%, so fix the progressbar.
        oh.close()

        # save the bin data as an emulated dict:
        grp = self.hdf5_handle.create_group('bin_lookup')
        for chrom in self.all_chrom_names:
            grp = self.hdf5_handle.create_group('bin_lookup/chrom_%s' % chrom)
            flat_bin = []

            for oldbinID in bin_lookup:
                if bin_lookup[oldbinID][0] == chrom:
                    flat_bin.append([oldbinID, bin_lookup[oldbinID][1], bin_lookup[oldbinID][2], bin_lookup[oldbinID][3]]) # i.e. oldId, left, right, newID
            flat_bin = numpy.array(flat_bin)
            grp.create_dataset('bins', (flat_bin.shape), dtype=int, data=flat_bin, chunks=True, compression='lzf')

            #self.hdf5_handle.create_dataset('bin_lookup', (len(to_store), 1), dtype='S10', data=to_store)

        # A sparse array would be better. The matrix is only around ~40% full in 100k, and
        # likely this gets much lower as the resolution drops...
        config.log.info('Loaded matrix %s*%s' % (self['num_bins'], self['num_bins']))
        for chrom in self.all_chrom_names:
            config.log.info('Added chrom=%s to table' % chrom)
            grp = self.hdf5_handle.create_group('matrix_%s' % chrom)
            grp.create_dataset('mat', matrices[chrom].shape, dtype=numpy.float32, data=matrices[chrom])
        config.log.info('Saved all matrices to hdf5')

        dat = [str(n).encode("ascii", "ignore") for n in self.all_chrom_names]
        self.hdf5_handle.create_dataset('all_chrom_names', (len(self.all_chrom_names), 1), 'S10', dat)

        self._OEmatrix()
        self._ABcompart()

        return True

    def load_cooler_pixels_dump(self, matrix_filename):
        """
        **Purpose**
            Load a file output by cooler with the command:

            cooler dump --header --join --balanced ${out}.cooler | gzip > ${out}

            e.g:

            chrom1	start1	end1	chrom2	start2	end2	count	balanced
            chr1	0	150000	chr1	0	150000	22	0.00115439
            chr1	0	150000	chr1	150000	300000	129
            chr1	0	150000	chr1	300000	450000	6
            chr1	0	150000	chr1	450000	600000	10
            chr1	0	150000	chr1	600000	750000	37
            chr1	0	150000	chr1	750000	900000	34	0.00165211
            chr1	0	150000	chr1	900000	1050000	15
            chr1	0	150000	chr1	1050000	1200000	9
            chr1	0	150000	chr1	1200000	1350000	8

            Only rows with a balanced score are kept. Otherwise set to 0.

        **Arguments**
            filename (Required)
                the Filename to load.
        """
        assert not self.readonly, 'To load, this must be read-only, set new=True'

        # First go through the file and get all the bins, chroms and sizes.
        bins = []
        bin_size = None

        # Sample the files to get the chromosome sizes:
        max_chrom_size = {}
        config.log.info('Preparse {0}'.format(matrix_filename))
        oh = gzip.open(matrix_filename, 'rt')
        for lin in oh:
            if 'chrom1' in lin:
                continue

            lin = lin.strip().split('\t')

            if not bin_size: # sample the bin size;
                bin_size = int(lin[2]) - int(lin[1])

            chr1 = lin[0]
            if chr1 not in max_chrom_size:
                max_chrom_size[chr1] = 0
            if int(lin[2]) > max_chrom_size[chr1]:
                max_chrom_size[chr1] = int(lin[2])

            # same for other read, as you could imagine a situation where a right bin is not seen in the left bin?
            chr2 = lin[3]
            if chr2 not in max_chrom_size:
                max_chrom_size[chr2] = 0
            if int(lin[5]) > max_chrom_size[chr2]:
                max_chrom_size[chr2] = int(lin[5])

        for chrom in sorted(max_chrom_size):
            config.log.info('Observed chrom {0} sizes = {1:,} bp'.format(chrom, max_chrom_size[chrom]))

        # Build the bin lookups based on the maximum chromsome sizes and the sampled
        # binsizes;
        bin_id = 0
        bin_lookup_by_chr_left = {}
        bin_lookup_by_id = {}
        for chrom in max_chrom_size:
            for left in range(0, max_chrom_size[chrom]+bin_size, bin_size):
                bin_loc = (chrom, left)
                if bin_loc not in bin_lookup_by_chr_left:
                    bin_lookup_by_chr_left[bin_loc] = bin_id
                    bin_lookup_by_id[bin_id] = bin_loc
                    bin_id += 1

        oh.close()

        self.all_chrom_names = set(max_chrom_size)
        self.hdf5_handle.attrs['num_bins'] = bin_id
        self.hdf5_handle.attrs['bin_size'] = bin_size

        config.log.info('Found {0:,} bins, each bin = {1:,} bp'.format(self['num_bins'], self['bin_size']))
        # the bins are not sorted,
        bin_order = sorted(bins, key=lambda v: (v[0], v[1]))

        # Get the edges of the chroms, and so the matrix sizes
        self.chrom_edges = {}
        for chrom in self.all_chrom_names:
            self.chrom_edges[chrom] = [1e20,-1]

            for bin in bin_lookup_by_chr_left:
                if bin[0] != chrom:
                    continue

                # get the edges:
                bin_id = bin_lookup_by_chr_left[bin]
                if bin_id < self.chrom_edges[chrom][0]:
                    self.chrom_edges[chrom][0] = bin_id
                if bin_id > self.chrom_edges[chrom][1]:
                    self.chrom_edges[chrom][1] = bin_id

        # Have to then convert them all to the new bins:
        # Although not strictly necessary in the matrix-based, I preserve this fiddly step so that
        # when I implement the inter-chrom interactions it is easier.
        matrices = {}
        bins = {} # The bin -> matrix convertor
        for chrom in self.all_chrom_names:
            # Now I know how big the arrays need to be, and all of the chromosomes:
            size = self.chrom_edges[chrom][1]-self.chrom_edges[chrom][0]
            matrices[chrom] = numpy.zeros((size+1, size+1), dtype='float32')
            # TODO: Support for interchrom with a sparse array:

        oh = gzip.open(matrix_filename, 'rt')
        p = progressbar(self['num_bins'] * self['num_bins']) # expected maximum
        for idx, lin in enumerate(oh):
            if 'chrom1' in lin:
                continue
            lin = lin.strip().split('\t')

            if len(lin) < 8: # balanced is blank;
                continue

            # First check the two bins are on the same chromosome:
            if lin[0] == lin[3]:
                chrom = lin[0]

                bin1 = bin_lookup_by_chr_left[(chrom, int(lin[1]))]
                bin2 = bin_lookup_by_chr_left[(chrom, int(lin[4]))]

                x = bin1-self.chrom_edges[chrom][0]
                y = bin2-self.chrom_edges[chrom][0]

                matrices[chrom][x,y] = float(lin[7]) # i.e. weight
                matrices[chrom][y,x] = float(lin[7])
            else:
                # TODO: Support for interchrom with a sparse array:
                pass

            p.update(idx)
        p.update(self['num_bins'] * self['num_bins']) # the above is unlikely to make it to 100%, so fix the progressbar.
        oh.close()

        # save the bin data as an emulated dict:
        grp = self.hdf5_handle.create_group('bin_lookup')
        for chrom in self.all_chrom_names:
            grp = self.hdf5_handle.create_group('bin_lookup/chrom_%s' % chrom)
            flat_bin = []

            for bin in bin_lookup_by_chr_left:
                if bin[0] == chrom:
                    bin_id = bin_lookup_by_chr_left[bin]
                    flat_bin.append([bin_id-self.chrom_edges[chrom][0], bin[1], bin[1]+bin_size, bin_id]) # i.e. oldId, left, right, newID
            flat_bin = numpy.array(flat_bin)

            grp.create_dataset('bins', (flat_bin.shape), dtype=int, data=flat_bin, chunks=True, compression='lzf')

            #self.hdf5_handle.create_dataset('bin_lookup', (len(to_store), 1), dtype='S10', data=to_store)

        # A sparse array would be better. The matrix is only around ~40% full in 100k, and
        # likely this gets much lower as the resolution drops...
        config.log.info('Loaded matrix %s*%s' % (self['num_bins'], self['num_bins']))
        for chrom in self.all_chrom_names:
            config.log.info('Added chrom=%s to table' % chrom)
            grp = self.hdf5_handle.create_group('matrix_%s' % chrom)
            grp.create_dataset('mat', matrices[chrom].shape, dtype=numpy.float32, data=matrices[chrom])
        config.log.info('Saved matrices to hdf5')

        dat = [str(n).encode("ascii", "ignore") for n in self.all_chrom_names]
        self.hdf5_handle.create_dataset('all_chrom_names', (len(self.all_chrom_names), 1), 'S10', dat)

        self._OEmatrix()
        self._ABcompart()

        return True

    def save_np3_column_matrix(self, filename, nohead=False):
        """
        **Purpose**
            Some other tools ask for a n+3 chromosome matrix,

            In the form
            chrom   left    right   0   0   0   0   0   0    .... #bins

            if nohead=True, then save just the matrix

            These tables are chromsome local only.

            Additionally, it will save them as 'filename_chromX.matrix'
            with one chromosome in each matrix file

        **Arguments**
            filename (Required)
                filename to save the matrix file to.

                The main aim of this tool is input for TopDom

        """
        assert self.readonly, 'must be readonly to save_np3_column_matrix. set new=False'
        assert filename, 'You need to specify a filename'

        for chrom in self.all_chrom_names:
            chrom_name = chrom
            if 'chr' not in chrom:
                chrom_name = 'chr%s' % chrom_name

            actual_filename = '%s_chrom%s.matrix' % (filename, chrom)
            oh = open(actual_filename, 'w')

            mostLeft, mostRight = self.__find_binID_chromosome_span(chrom)
            mat = self.mats[chrom][()]
            bins = self.bin_lookup_by_chrom[chrom]
            for m, b in zip(mat, bins):
                if nohead:
                    lin = [str(i) for i in m]
                else:
                    lin = [chrom_name, b[1], b[2]] + list(m)
                    lin = [str(i) for i in lin]
                oh.write('{0}\n'.format('\t'.join(lin)))
            oh.close()

            config.log.info('Saved save_np3_column_matrix() "%s"' % actual_filename)

    def load_tad_calls(self, filename, format='bed'):
        """
        **Purpose**
            Load and bind a set of TAD calls from some other tool.

        **Arguments**
            filename (Required)
                filename to load the TAD calls from

            format (Optional, default=TopDom)
                One of:
                'bed':
                    Load a 3 column BED file with the calls
                'TopDom':
                    Load a TopDom BED file.
                    This includes the 'gap/domain/boundary'
                    information, which will also be plotted.

        **Returns**
            The TAD BED as a genelist
        """
        assert self.readonly, 'must be readonly to load_tad_calls. set new=False'
        assert filename, 'You need to specify a filename'

        if format == 'bed':
            format = minimal_bed
        elif format == 'TopDom':
            format = minimal_bed
            format.update({"tad_type": 3})

        self.tad_calls = genelist(filename, format=format)
        self.tad_calls.sort('loc')
        # make a quicklookup by chrom for TAD extraction
        self.tad_lookup = {}
        for chrom in self.all_chrom_names:
            self.tad_lookup[chrom] = []
        for item in self.tad_calls:
            self.tad_lookup[item['loc']['chr']].append(item)

        return self.tad_calls

    def heatmap(self,
        filename,
        chr=None,
        loc=None,
        key='matrix',
        bracket=None,
        colour_map=cm.inferno_r,
        log2=False,
        **kargs):
        """
        **Purpose**
            Draw an interaction heatmap

        **Arguments**
            filename (Required)
                filename to save the image to

            chr (Optional, default=None (i.e. all chromsomes))
            loc (Optional, default=None)
                Exclusive: You can use only one of chr and loc
                You can set the view to an entire chromosome with chr,
                or you can set a specific location with loc

            key (Optional, default='matrix')
                The data matrix name to use for plotting.
                Valid options are
                    matrix: normal contact probability
                    OE: generated with OEmatrix
                    AB: generated with ABmatrix

            aspect (Optional, default='square')
                image aspect for the heatmap

            colbar_label (Optional, default='log2(Density)')
                Label for the colour bar.

            colour_map (Optional, default=cm.inferno_r)
                matplotlib colour map to use to colour.

            bracket (Optional, default=None)
                clip the data within the ranges [low, high]

            log2 (Optional, default=False)
                transform the data to log2 before plotting;

        **Returns**
            None, and a file in filename.

        """
        assert self.readonly, 'must be readonly to draw a heatmap. set new=False'
        assert filename, "heatmap: you must specify a filename"
        if chr and loc:
            raise AssertionError('chr and loc both contain values. You can only use one')

        dataset_to_use = {
            'matrix': self.mats,
            'OE': self.OE,
            'AB': self.OE,
            }
        # TODO: sanity checking on matrix availability;
        cmap_to_use = {
            'matrix': cm.viridis,
            'OE': cm.RdBu_r,
            'AB': cm.BrBG_r,
            }

        assert key in dataset_to_use, '{} is not a valid dataset key'.format(key)

        if chr:
            data = dataset_to_use[key][str(chr)]
            #ABdata = self.AB[str(chr)]
        elif loc:
            if not isinstance(loc, location):
                loc = location(loc)

            chrom = loc.loc['chr']
            if 'chr' not in chrom:
                chrom = 'chr{}'.format(chrom)

            # Need to get the binIDs and the updated location span
            localLeft, localRight, loc, _, _ = self.__find_binID_spans(loc)

            data = dataset_to_use[key][chrom][localLeft:localRight, localLeft:localRight]
            #ABdata = self.AB[chrom][localLeft:localRight]
        else:
            raise NotImplementedError('chr=None not implemented')
            data = self.matrix # use the whole lot

        if not "aspect" in kargs:
            kargs["aspect"] = "square"

        colbar_label = "Density"
        if "colbar_label" in kargs:
            colbar_label = kargs['colbar_label']

        fig = self.draw.getfigure(**kargs)

        # positions of the items in the plot:
        heatmap_location  = [0.20,   0.01,   0.79,   0.79]
        ABtop =             [0.20,   0.81,   0.79,   0.07]
        ABlef =             [0.12,   0.01,   0.07,   0.79]
        scalebar_location = [0.20,   0.97,   0.79,   0.02]

        if bracket: # done here so clustering is performed on bracketed data
            #data = self.draw.bracket_data(numpy.log2(self.matrix+0.1), bracket[0], bracket[1])
            if log2:
                with numpy.errstate(divide='ignore'):
                    data = numpy.log2(data)
                    data[numpy.isneginf(data)] = 0
                colbar_label = 'Log2(Density)'
            # Faster numpy"
            data = numpy.clip(data, bracket[0], bracket[1])
            vmin = bracket[0]
            vmax = bracket[1]
        else:
            if log2:
                with numpy.errstate(divide='ignore'):
                    data = numpy.log2(data)
                colbar_label = 'Log2(Density)'
            #data[data == -numpy.inf] = 0
            data = numpy.array(data)
            vmin = data.min()
            vmax = data.max()

        # ---------------- (A/B plots) ---------------------
        if key == 'AB':
            ABdata = numpy.array(self.AB[str(chr)])

            ax1 = fig.add_subplot(142)
            ax1.set_position(ABtop)
            ax1.plot(ABdata)
            ax1.set_xlim([0, len(ABdata)])
            ax1.axhline(0.0, lw=1.0, c='black')
            ax1.tick_params(left=None, bottom=None)
            ax1.set_xticklabels('')
            ax1.set_yticklabels('')

            ax2 = fig.add_subplot(143)
            ax2.set_position(ABlef)
            ax2.plot(range(len(ABdata)), ABdata)
            ax2.set_ylim([0, len(ABdata)])
            ax2.axvline(0.0, lw=1.0, c='black')
            #ax2.tick_params(left=None, bottom=None)
            #ax2.set_xticklabels('')
            #ax2.set_yticklabels('')


        # ---------------- (heatmap) -----------------------
        ax3 = fig.add_subplot(141)

        ax3.set_position(heatmap_location) # must be done early for imshow
        hm = ax3.imshow(data, cmap=cmap_to_use[key], vmin=vmin, vmax=vmax, aspect="auto",
            origin='lower', extent=[0, data.shape[1], 0, data.shape[0]],
            interpolation=config.get_interpolation_mode(filename))

        #ax3.set_frame_on(True)
        ax3.set_position(heatmap_location)
        ax3.set_xlim([0,data.shape[1]])
        ax3.set_ylim([0,data.shape[0]])
        ax3.set_yticklabels("")
        ax3.set_xticklabels("")

        ax3.tick_params(top=False, bottom=False, left=False, right=False)
        [t.set_fontsize(5) for t in ax3.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(5) for t in ax3.get_xticklabels()]

        ax0 = fig.add_subplot(144)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)

        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0)
        cb.set_label(colbar_label)
        [label.set_fontsize(5) for label in ax0.get_xticklabels()]

        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("heatmap: Saved '%s'" % actual_filename)

    def tri_plot(self, filename, chr=None, loc=None, bracket=None, colour_map=cm.inferno_r,
        **kargs):
        """
        **Purpose**
            Draw the side-on half triangle plot.

        **Arguments**
            filename (Required)
                filename to save the image to.

            chr (Optional, default=None (i.e. all chromsomes))
            loc (Optional, default=None)
                Exclusive: You can use only one of chr and loc
                You can set the view to an entire chromosome with chr,
                or you can set a specific location with loc

            colour_map (Optional, default=cm.inferno_r)
                matplotlib colour map to use to colour.

            bracket (Optional, default=None)
                clip the data within the ranges [low, high]

        **Returns**
            None and a file in filename.
        """
        assert self.readonly, 'must be readonly to draw a heatmap. set new=False'
        assert filename, "heatmap: you must specify a filename"
        if chr and loc:
            raise AssertionError('chr and loc both contain values. You can only use one')

        if chr:
            data = self.mats[str(chr).replace('chr', '')]
        elif loc:
            if not isinstance(loc, location):
                loc = location(loc)

            # Need to get the binIDs and the updated location span
            localLeft, localRight, loc, _, _ = self.__find_binID_spans(loc)

            data = self.mats[loc['chr']][localLeft:localRight, localLeft:localRight]
        else:
            data = self.matrix # use the whole lot

        if bracket: # done here so clustering is performed on bracketed data
            #data = self.draw.bracket_data(numpy.log2(self.matrix+0.1), bracket[0], bracket[1])
            with numpy.errstate(divide='ignore'):
                data = numpy.log2(data)
            # Faster numpy"
            data = numpy.clip(data, bracket[0], bracket[1])
            vmin = bracket[0]
            vmax = bracket[1]
        else:
            with numpy.errstate(divide='ignore'):
                data = numpy.log2(data)
            vmin = data.min()
            vmax = data.max()

        if not "aspect" in kargs:
            kargs["aspect"] = "square"
        if not "colbar_label" in kargs:
            kargs["colbar_label"] = "log2(Density)"

        heatmap_location =  [0.05,   0.01,   0.90,   0.80]
        scalebar_location = [0.05,  0.97,   0.90,   0.02]

        if 'figsize' not in kargs:
            kargs['figsize'] = (8,3)
        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        dst = ndimage.rotate(data, 45, order=0, reshape=True, prefilter=False, cval=0)

        ax.set_position(heatmap_location)
        hm = ax.imshow(dst, cmap=colour_map, vmin=vmin, vmax=vmax, aspect="auto",
            origin='lower', extent=[0, data.shape[1], 0, data.shape[0]],
            interpolation=config.get_interpolation_mode(filename))

        ax.set_xlim([0,data.shape[1]])
        ax.set_ylim([data.shape[1]/2,(data.shape[1]/2)+data.shape[1]/4]) # get the middle start position, then about 1/X of the way up
        #print([dst.shape[1]/2,data.shape[1]])
        #ax.set_yticklabels("")
        #ax.set_xticklabels("")

        ax.yaxis.tick_right()
        ax.tick_params(top=False, bottom=False, left=False, right=False)
        [t.set_fontsize(5) for t in ax.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(5) for t in ax.get_xticklabels()]

        ax0 = fig.add_subplot(122)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)
        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0, cmap=colour_map)
        cb.set_label(kargs["colbar_label"])
        [label.set_fontsize(5) for label in ax0.get_xticklabels()]

        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("tri_plot: Saved '%s'" % actual_filename)

        return()

    def __plot_tad_calls(self, ax, ax_position, loc, tad_calls):
        """
        (Internal)

        Plot Tad calls from a flat bed with no 'tad_type' key

        """
        borders = [tad['loc']['left'] for tad in tad_calls]# This is easy, just take all of the left most positions:
        tad_loc = loc.expand(len(loc)//5) # I want to collect off the edges so the plots spread off the edge

        # TAD triangle plot
        tads = []
        for tad in tad_calls:
            # check its inside the mostLeft and mostRight:
            if tad_loc.qcollide(tad['loc']):
                tads.append((tad['loc']['left'], tad['loc']['right']))

        xx = []
        yy = []
        for b in tads:
            midH = (b[1] - b[0])/2
            mid = midH + b[0]

            xx.append(b[0])
            yy.append(0)

            xx.append(mid)
            yy.append(midH)
        ax.set_position(ax_position)
        ax.plot(xx, yy)
        ax.tick_params(top=False, bottom=False, left=False, right=False)
        ax.set_xlim([loc['left'], loc['right']]) # yes, to the loc
        ax.set_yticklabels("")
        #ax.set_xticklabels("")

    def __plot_tad_calls_topdom(self, ax, ax_position, loc, tad_calls):
        """
        (Internal)

        Plot Tad calls from a flat bed with a 'tad_type' key

        """
        borders = [tad['loc']['left'] for tad in tad_calls]# This is easy, just take all of the left most positions:
        tad_loc = loc.expand(len(loc)//5) # I want to collect off the edges so the plots spread off the edge

        # TAD triangle plot
        tads = []
        for tad in tad_calls:
            # check its inside the mostLeft and mostRight:
            if tad_loc.qcollide(tad['loc']):
                tads.append((tad['loc']['left'], tad['loc']['right'], tad['tad_type']))

        xx = []
        yy = []
        for b in tads:
            midH = (b[1] - b[0])/2
            mid = midH + b[0]

            xx.append(b[0])
            yy.append(0)

            xx.append(mid)
            yy.append(midH)
        ax.set_position(ax_position)
        ax.plot(xx, yy)
        ax.tick_params(top=False, bottom=False, left=False, right=False)
        ax.set_xlim([loc['left'], loc['right']]) # yes, to the loc
        ax.set_yticklabels("")
        #ax.set_xticklabels("")

    def tri_plot_and_freq_plot(self, filename, expn=None, expn_cond_name=None,
        chr=None, loc=None,
        bracket=None, colour_map=cm.inferno_r,
        **kargs):
        """
        **Purpose**
            Draw the side-on half triangle plot. and a frequency plot below

        **Arguments**
            filename (Required)
                filename to save the image to.

            expn (Required)
                an expression object, with conditions to use in the frequency plot.

            expn_cond_name

            chr (Optional, default=None (i.e. all chromsomes))
            loc (Optional, default=None)
                Exclusive: You can use only one of chr and loc
                You can set the view to an entire chromosome with chr,
                or you can set a specific location with loc

            colour_map (Optional, default=cm.inferno_r)
                matplotlib colour map to use to colour.

            bracket (Optional, default=None)
                clip the data within the ranges [low, high]

        **Returns**
            None and a file in filename.
        """
        assert self.readonly, 'must be readonly to draw a heatmap. set new=False'
        assert filename, "heatmap: you must specify a filename"
        if chr and loc:
            raise AssertionError('chr and loc both contain values. You can only use one')

        tad_calls = None

        if chr:
            data = self.mats[str(chr).replace('chr', '')]
            if self.tad_lookup:
                tad_calls = self.tad_lookup[loc['chrom']]

            mostLeft, mostRight = self.__find_binID_chromosome_span(str(chr).replace('chr', ''))

            # I'm a bit confused by this code :*(
            this_chrom = [0] * (mostRight-mostLeft+1)
            cindex = expn.getConditionNames().index(expn_cond_name)
            for i in expn.linearData:
                # Take the old BinID and convert it to the new binID:
                # If inside this part of the chromosome:
                if str(i['loc']['chr']) == str(chr.replace('chr', '')):
                    local_bin_num = (self.bin_lookup_by_binID[i['bin#']][3] - mostLeft)
                    this_chrom[local_bin_num] = i['conditions'][cindex]

            plot_y = this_chrom
            plot_x = numpy.arange(0, len(plot_y))

        elif loc:
            if not isinstance(loc, location):
                loc = location(loc)

            # Need to get the binIDs and the updated location span
            localLeft, localRight, loc, mostLeft, mostRight = self.__find_binID_spans(loc)
            data = self.mats[loc['chr']][localLeft:localRight, localLeft:localRight]

            # I just assume the bins match
            this_chrom = [0] * (localRight-localLeft+1)
            cindex = expn.getConditionNames().index(expn_cond_name)
            for i in expn.linearData:
                if i['loc']['chr'] == loc['chr']:
                    # Take the old BinID and convert it to the new binID:
                    # If inside this part of the chromosome:
                    if loc.qcollide(i['loc']):
                        local_bin_num = (self.bin_lookup_by_binID[i['bin#']][3] - mostLeft) - localLeft
                        this_chrom[local_bin_num] = i['conditions'][cindex]
            plot_y = this_chrom
            plot_x = numpy.arange(0, len(plot_y))


            if self.tad_lookup:
                tad_calls = self.tad_lookup[loc['chr']]
        else:
            raise(AssertionError, 'Not implemented here! :(')

        if bracket: # done here so clustering is performed on bracketed data
            #data = self.draw.bracket_data(numpy.log2(self.matrix+0.1), bracket[0], bracket[1])
            with numpy.errstate(divide='ignore'):
                data = numpy.log2(data)
            # Faster numpy"
            data = numpy.clip(data, bracket[0], bracket[1])
            vmin = bracket[0]
            vmax = bracket[1]
        else:
            with numpy.errstate(divide='ignore'):
                data = numpy.log2(data)
            vmin = data.min()
            vmax = data.max()

        if not "aspect" in kargs:
            kargs["aspect"] = "square"
        if not "colbar_label" in kargs:
            kargs["colbar_label"] = "log2(Density)"

        scalebar_location = [0.05,  0.97,   0.90,   0.02]
        heatmap_location  = [0.05,  0.50,   0.90,   0.40]
        tadtri_location   = [0.05,  0.40,   0.90,   0.10]
        tadborder_location= [0.05,  0.35,   0.90,   0.05] # Bashes into;
        freq_location     = [0.05,  0.01,   0.90,   0.28]

        if 'figsize' not in kargs:
            kargs['figsize'] = (8,5)

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(512)

        dst = ndimage.rotate(data, 45, order=0, reshape=True, prefilter=False, cval=0)

        ax.set_position(heatmap_location)
        hm = ax.imshow(dst, cmap=colour_map, vmin=vmin, vmax=vmax, aspect="auto",
            origin='lower', extent=[0, data.shape[1], 0, data.shape[0]],
            interpolation=config.get_interpolation_mode(filename))

        ax.set_xlim([0,data.shape[1]])
        ax.set_ylim([data.shape[1]/2,(data.shape[1]+data.shape[1]/10)]) # the /10 determines how much it gets stretched.
        ax.set_yticklabels("")
        ax.set_xticklabels("")

        ax.yaxis.tick_right()
        ax.tick_params(top=False, bottom=False, left=False, right=False)
        [t.set_fontsize(5) for t in ax.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(5) for t in ax.get_xticklabels()]

        ax0 = fig.add_subplot(511)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)
        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0, cmap=colour_map)
        t = cb.set_label(kargs["colbar_label"], fontsize=6)
        #t.set_fontsize(6)
        cb.ax.tick_params(labelsize=5)
        #[label.set_fontsize(5) for label in ax0.get_xticklabels()]

        if tad_calls:
            ax = fig.add_subplot(514)
            if 'tad_type' in self.tad_lookup[loc['chr']]:
                self.__plot_tad_calls(ax, tadtri_location, loc, tad_calls)
            else: # assume a flat BED:
                self.__plot_tad_calls_topdom(ax, tadtri_location, loc, tad_calls)
        # end TADS

        # freq plot
        ax1 = fig.add_subplot(515)
        ax1.set_position(freq_location)
        #ax1.set_frame_on(False)
        # get a chrom if required:
        ax1.plot(plot_x, plot_y)
        ax1.set_xlim([0, len(plot_y)])
        ax1.tick_params(top=False, bottom=False, left=True, right=False)
        ax1.set_xticklabels("")

        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("tri_plot: Saved '%s'" % actual_filename)

        return()

    def density_plot(self, filename, vmin=0, vmax=50000, **kargs):
        """
        **Purpose**
            Draw a density histogram of the matrix scores.

        **Arguments**
            filename (Required)
                filenmame to save the image to

        **Returns**
            None and image in filename
        """
        assert self.readonly, 'must be readonly to draw a density_plot. set new=False'
        assert filename, "density_plot: you must specify a filename"

        fig = self.draw.getfigure(**kargs)
        ax1 = fig.add_subplot(121) # All
        ax2 = fig.add_subplot(122) # inter-chromosome

        max = 0
        for i in self.matrix: # this will return each row
            if i.max() > max:
                max = i.max()

        data = [0] * (int(max)+1)
        for r in self.matrix:
            for i in r:
                if int(i) > 0: # There is a lot of 0
                    data[int(i)] += 1

        ax1.plot(numpy.arange(0, len(data)), data)

        l = numpy.log2(numpy.array(data))
        ax2.plot(numpy.arange(0, len(data)), l)

        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("density_plot: Saved '%s'" % actual_filename)

        return()

    def pca(self, number_of_components=10, chrom=None):
        """
        **Purpose**
            Train a PCA on the matrix
        """
        assert self.readonly, 'must be readonly to draw a density_plot. set new=False'
        assert chrom, 'You must specify a chromosome for pca()'

        matrix = self.mats[str(chrom).replace('chr', '')][:]
        matrix = numpy.log2(matrix+0.1)

        self.__pcalabels = self.bin_lookup_by_chrom[str(chrom).replace('chr', '')]

        self.__model = PCA(n_components=number_of_components, whiten=True)
        self.__transform = self.__model.fit_transform(matrix) # U, sample loading
        self.__components = self.__model.components_.T # V, The feature loading

        self.pca_valid = chrom
        config.log.info("pca: Trained PCA")

    def explained_variance(self, filename=None, percent_variance=True, **kargs):
        """
        **Purpose**
            plot a graph of PC loading, percent variance for each componenet

        **Arguments**
            filename (Required)
                The filename to save the image to.

            <common figure arguments are supported>

        **Returns**
            None.
        """
        assert filename, "explained_variance: Must provide a filename"
        assert self.pca_valid, 'model has not been trained, use pca()'

        if "aspect" not in kargs:
            kargs["aspect"] = "wide"

        expn_var = numpy.array(self.__model.explained_variance_ratio_) * 100.0

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        x = numpy.arange(len(expn_var))
        ax.bar(x, expn_var, ec="none", color="grey")
        ax.set_xlabel("Principal components")
        if percent_variance:
            ax.set_ylabel('Percent Variance')
        else:
            ax.set_ylabel("Loading")
        ax.set_xticklabels(x+1)
        ax.set_xticks(x)
        ax.set_xlim([-0.5, len(expn_var)-0.5])
        self.draw.do_common_args(ax, **kargs)
        real_filename = self.draw.savefigure(fig, filename)

        config.log.info("explained_variance: Saved PC loading '%s'" % real_filename)

    def get_loading_percents(self, **kargs):
        """
        **Purpose**
            Returns the percent of variance

        **Arguments**
            None

        **Returns**
            Returns an array for each PC and it's percent variance
        """
        return(numpy.array(self.__model.explained_variance_ratio_) * 100.0)

    def scatter(self, mode, x=None, y=None, filename=None, spot_cols=None, cmap=None, spots=True, label=False, alpha=0.5, overplot=None,
        spot_size=5, label_font_size=7, label_style='normal', cut=None, squish_scales=False, **kargs):
        """
        **Purpose**
            plot a scatter plot of PCx against PCy.

        **Arguments**
            mode (Required)
                One of 'pca' or 'tsne'
                The mode to take the scatter data from.

            x, y (Required when mode=='pca')
                PC dimension to plot as scatter
                Note that PC begin at 1 (and not at zero, as might be expected)

            filename (Required)

            spot_cols (Optional, default="black" or self.set_cols())
                list of colours for the samples, should be the same length as
                the number of conditions.

                if labels == True and spots == False and spot_cols is not None then
                    spot_cols will be used to colour the labels.

            only_plot_if_x_in_label (Optional, default=None)
                Only plot an individual scatter if X is in the label name.

                This must be a list or tuple of names

                Allows you to effectively remove points from the PCA plot.

            spots (Optional, default=True)
                Draw the spots

            spot_size (Optional, default=40)
                Size of the spots on the scatter

            label_font_size (Optional, default=7)
                Size of the spot label text, only valid if label=True

            label_style (Optional, default='normal')
                add a 'style' for the text labels, one of:
                'normal', 'italic', 'oblique'

            cut (Optional, default=None)
                Send a rectangle of the form [topleftx, toplefty, bottomrightx, bottomrighty], cut out all of the items within that
                area and return their label and PC score

            squish_scales (Optional, default=False)
                set the limits very aggressively to [minmin(x), minmax(y)]

        **Returns**
            None
        """
        assert filename, "scatter: Must provide a filename"
        assert self.pca_valid, 'PCA model has not been trained, or the model was trained on the wrong chromosome, use pca()'

        if mode == 'pca':
            assert x and y, 'x and y must be valid principal components in pca(mode="scatter", ...)'

            labels = self.__pcalabels
            xdata = self.__transform[:,x-1]
            ydata = self.__transform[:,y-1]
            mode = 'PC'
            perc_weights = self.get_loading_percents()
        elif mode == 'tsne':
            assert self.tsne_trained, 'tSNE model has not been trained, or the model was trained on the wrong chromosome, use tsne()'
            labels = self.__pcalabels
            xdata = self.npos[:, 0]
            ydata = self.npos[:, 1]
            mode = 'tSNE '
            perc_weights = None
        else:
            raise AssertionError('mode "%s" not found' % mode)

        if spot_cols is None: # colour from start of chrom to end
            spot_cols = numpy.arange(0, len(xdata)) # linear colouring
            cmap=cm.inferno
            # I think this is also a system you could use to e.g. put the frequency of something straight on the plot?

        return_data = self.draw.unified_scatter(labels, xdata, ydata, x=x, y=y, filename=filename,
            spot_cols=spot_cols, spots=spots, alpha=alpha, cmap=cmap,
            perc_weights=perc_weights, mode=mode,
            spot_size=spot_size, label_font_size=label_font_size, cut=cut, squish_scales=squish_scales,
            **kargs)

        return(return_data)

    def tsne(self, num_pc, chrom):
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
        assert chrom, 'You must specify a chromosome for tsne()'

        matrix = self.mats[str(chrom).replace('chr', '')][:]
        matrix = numpy.log2(matrix+0.1)
        self.__pcalabels = self.bin_lookup_by_chrom[str(chrom).replace('chr', '')]

        if isinstance(num_pc, int):
            self.__tsne_model = PCA(n_components=num_pc, whiten=True)
            self.__tsne_transform = self.__model.fit_transform(matrix)
            self.__tsne_pcas = self.__transform

        elif isinstance(num_pc, list):
            self.__tsne_model = PCA(n_components=max(num_pc)+1, whiten=True)
            self.__tsne_transform = self.__model.fit_transform(matrix)
            # get only the specific PCs
            self.__tsne_pcas = numpy.array([self.__transform[:,c-1] for c in num_pc]).T

        else:
            raise AssertionError('num_pcs must be either an integer or a list')

        self.__tsne_final_model = TSNE(n_components=2, perplexity=50, init='pca', random_state=0, n_iter=5000,
            n_iter_without_progress=1000)
        self.npos = self.__model.fit_transform(self.__tsne_pcas)

        self.tsne_trained = chrom
        config.log.info("tsne: Trained tSNE")

    def square_plot(self,
        filename:str = None,
        center_anchors = None,
        distal_anchors = None,
        bedpe = None,
        window:int = None,
        num_bins:int = None,
        bracket = None,
        log2 = None,
        colour_map = cm.viridis,
        **kargs
        ):
        """
        **Purpose**
            Draw a square plot for HiC intensity,
            Take the center_anchor points, and then draw out to eith ther distal anchors, or to some <window>
            Take the average of the matrix;

        **Arguemnts**
            filename

            center_anchors

            distal_anchors

            bedpe (Optional, default=None)
                Or, you could pass a genelist with a loc1 and loc2 pair of locations to use for the corners.

        **Returns**
            ?

        """
        assert self.readonly, 'must be readonly to draw a square_plot. set new=False'
        assert filename, 'You need a filename to save the image to'
        assert not (distal_anchors and window), 'Only one of distal_anchors or window can be valid'
        if bedpe: assert num_bins, 'You must specify num_bins for the heatmap if bedpe=True'
        if not window: assert (bedpe or distal_anchors), 'distal_anchors or bedpe not True'
        if not distal_anchors and not bedpe: assert window, 'window not True'
        if bedpe:
            assert 'loc1' in bedpe.keys(), 'need a "loc1" key in bedpe genelist'
            assert 'loc2' in bedpe.keys(), 'need a "loc2" key in bedpe genelist'

        if distal_anchors:
            raise NotImplementedError('distal_anchors not implemented')

        if bedpe:
            __num_skipped = 0
            __intra_chroms = 0
            mat = numpy.zeros((num_bins, num_bins))
            int_range = numpy.arange(num_bins)

            left = [(loc['chr'], loc['left'], loc['right']) for loc in bedpe['loc1']]
            right = [(loc['chr'], loc['left'], loc['right']) for loc in bedpe['loc2']]
            p = progressbar(len(bedpe))
            for aidx, (l, r) in enumerate(zip(left, right)):
                if l[0] != r[0]:
                    __intra_chroms += 1
                    continue # differnent chroms are not supported

                chrom = 'chr{0}'.format(l[0])

                scaled_window = (max([l[2], r[2]]) - min([l[1], r[1]])) // 10

                #print(l, r)
                localLeft1, localRight1 = self.__quick_find_binID_spans(chrom, l[1]-scaled_window, l[2]+scaled_window, do_assert_check=False)
                localLeft2, localRight2 = self.__quick_find_binID_spans(chrom, r[1]-scaled_window, r[2]+scaled_window, do_assert_check=False)

                localLeft = min([localLeft1, localRight1, localLeft2, localRight2])
                localRight = max([localLeft1, localRight1, localLeft2, localRight2])

                #scaled_window = (localRight - localLeft) // 5 # must match the data.shape[0] < 5 below;

                #print(localLeft, localRight)

                localLeft -= scaled_window
                localRight += scaled_window

                #print(localLeft, localRight, scaled_window)

                data = self.mats[chrom][localLeft:localRight, localLeft:localRight]

                if data.shape != mat.shape:
                    if data.shape[0] < 10:
                        # Seems the two loops are too close together for this resolution, skip it;
                        __num_skipped += 1
                        continue
                    #print(data.shape, mat.shape)
                    data = reshap_mats(data, num_bins, num_bins)

                mat += data
                p.update(aidx)

            mat /= len(bedpe)

            m10 = mat.shape[0] // 10

            markers = [m10, mat.shape[0] - m10]

            config.log.info('{0} were on differnet chromsomes and were not used'.format(__intra_chroms))
            config.log.info('{0} loci were too close together for this hiccy resolution and were not used'.format(__num_skipped))

        elif distal_anchors:
            for anchor, distal in zip(center_anchors, distal_anchors):
                pass # Not implemented yet.

        elif window:
            num_bins = (window * 2) // self['bin_size']
            mat = numpy.zeros((num_bins, num_bins))
            int_range = numpy.arange(num_bins)

            p = progressbar(len(center_anchors))
            for aidx, anchor in enumerate(center_anchors):
                chrom = anchor['loc']['chr']
                if 'chr' not in chrom:
                    chrom = 'chr{0}'.format(chrom)
                localLeft, localRight = self.__quick_find_binID_spans(chrom, anchor['loc']['left']-window, anchor['loc']['right']+window)

                1/0 # I feel this is wrong somehow:
                ##data = self.mats[chrom][localLeft:localRight, localLeft:localRight]

                if data.shape != mat.shape:
                    data = reshap_mats(data, num_bins, num_bins)

                mat += data
                p.update(aidx)

            mat /= len(center_anchors)

            markers = [mat.shape[0] // 2, (mat.shape[0] // 2)+1]

        return self.square_plot_heatmap(filename=filename,
            center_anchors=center_anchors, distal_anchors=distal_anchors,
            window=window, bracket=bracket,
            log2=log2,
            mat=mat,
            colour_map=colour_map,
            markers = markers,
            **kargs)

    def square_plot_heatmap(self,
        filename:str = None,
        center_anchors = None,
        distal_anchors = None,
        window:int = None,
        bracket = None,
        log2 = None,
        colour_map = cm.plasma,
        mat = None,
        markers = None,
        **kargs
        ):

        fig = self.draw.getfigure(**kargs)

        # positions of the items in the plot:
        heatmap_location  = [0.05,   0.01,   0.90,   0.90]
        scalebar_location = [0.05,   0.97,   0.90,   0.02]

        if not "colbar_label" in kargs:
            colbar_label = "Density"

        if bracket: # done here so clustering is performed on bracketed data
            #data = self.draw.bracket_data(numpy.log2(self.matrix+0.1), bracket[0], bracket[1])
            if log2:
                with numpy.errstate(divide='ignore'):
                    mat = numpy.log2(mat)
                colbar_label = "Log2(Density)"
            mat = numpy.clip(mat, bracket[0], bracket[1])
            vmin = bracket[0]
            vmax = bracket[1]
        else:
            if log2:
                with numpy.errstate(divide='ignore'):
                    mat = numpy.log2(mat)
                colbar_label = "Log2(Density)"
            mat[mat == -numpy.inf] = 0
            vmin = mat.min()
            vmax = mat.max()

        # ---------------- (heatmap) -----------------------
        ax = fig.add_subplot(121)

        ax.set_position(heatmap_location) # must be done early for imshow
        hm = ax.imshow(mat, cmap=colour_map,
            vmin=vmin, vmax=vmax,
            aspect="auto",
            origin='lower',
            extent=[0, mat.shape[1], 0, mat.shape[0]],
            interpolation=config.get_interpolation_mode(filename))

        #ax3.set_frame_on(True)
        ax.set_position(heatmap_location)
        ax.set_xlim([0,mat.shape[1]])
        ax.set_ylim([0,mat.shape[0]])
        ax.set_yticklabels("")
        ax.set_xticklabels("")

        if markers:
            for m in markers:
                ax.axvline(m, ls=":", lw=0.5, color="grey")
                ax.axhline(m, ls=":", lw=0.5, color="grey")

        ax.tick_params(top=False, bottom=False, left=False, right=False)
        [t.set_fontsize(5) for t in ax.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(5) for t in ax.get_xticklabels()]

        ax0 = fig.add_subplot(122)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)

        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0, cmap=colour_map)
        cb.set_label(colbar_label, fontsize=6)
        cb.ax.tick_params(labelsize=6)
        #[label.set_fontsize(5) for label in ax0.get_xticklabels()]

        actual_filename = self.draw.savefigure(fig, filename, dpi=300)
        config.log.info("heatmap: Saved '%s'" % actual_filename)

        return mat

    def measure_loop_strength(self,
        bedpe=None,
        log=False,
        min_bins:int = 4,
        **kargs):
        '''
        **Purpose**
            Measure the loop strength for the selected bedpe (containig a loc1 and loc2 key)
            that describes a loop of chromatin between two loci. (i.e. take the diagonal)

            Only loops on the same chromsome are supported.

        **Arguments**
            bedpe (Required)
                a genelist containing a loc1 and loc2 key for the start of the loop and the end.

            min_bins (Optional, default=4)
                minimum number of bins distant to measure the loop (i.e. exclude loops that are too close).

        **Returns**
            a new genelist containing all the valid loci with a new key 'loop_strength'

        '''
        assert self.readonly, 'must be readonly to draw a square_plot. set new=False'
        #assert filename, 'You need a filename to save the image to'
        assert bedpe, 'You must specify a bedpe'
        assert 'loc1' in bedpe.keys(), 'need a "loc1" key in bedpe genelist'
        assert 'loc2' in bedpe.keys(), 'need a "loc2" key in bedpe genelist'

        newl = []

        min_pad = 0.01
        __num_skipped = 0
        __intra_chroms = 0

        p = progressbar(len(bedpe))
        for aidx, item in enumerate(bedpe.linearData):
            p.update(aidx)

            l1 = item['loc1'].loc
            l2 = item['loc2'].loc

            if l1['chr'] != l2['chr']:
                __intra_chroms += 1
                continue # differnent chroms are not supported

            chrom = 'chr{}'.format(l1['chr'])

            try:
                localLeft1, localRight1 = self.__quick_find_binID_spans(chrom, l1['left'], l1['right'], do_assert_check=False)
                localLeft2, localRight2 = self.__quick_find_binID_spans(chrom, l2['left'], l2['right'], do_assert_check=False)
            except KeyError:
                # The chrom is in one list, but not in the other, ignore it
                continue

            localLeft = min([localLeft1, localRight1, localLeft2, localRight2])
            localRight = max([localLeft1, localRight1, localLeft2, localRight2])

            #print(localLeft1, localRight1, localLeft2, localRight2, localLeft, localRight, item['loc1'], item['loc2'], self.mats[chrom][localLeft, localRight])

            if localRight - localLeft < min_bins:
                __num_skipped += 1
                continue

            loop_strength = self.mats[chrom][localLeft, localRight]

            item = dict(item)
            if log:
                item['loop_strength'] = math.log2(float(loop_strength)+min_pad)
            else:
                item['loop_strength'] = float(loop_strength)
            newl.append(item)

        gl = genelist()
        gl.load_list(newl)
        gl.name = bedpe.name

        config.log.info('{} were on different chromosomes and were not used'.format(__intra_chroms))
        config.log.info('{} loci were too close together for this hiccy resolution and were not used'.format(__num_skipped))

        return gl

    def contact_probability(self, min_dist, max_dist, anchors=None, filename=None, skip_Y_chromosome=True, **kargs):
        '''
        **Purpose**
            Measure the contact probability from in_dist to max dist,
            draw an image and return the histogram.

        **Arguments**
            min_dist (Required)
                minimum genomic distance to consider.

            max_dist (Required)
                maximum genomic distance to consider.

            anchors (Optional, default=FAlse)
                By default contact_probability considers all genomic locations.
                If this is set to a genelist with a 'loc' key then only those specific loci will be used.

            filename (Optional, default=False)
                filename to save the resulting histogram to.

            skip_Y_chromosome (Optional, default=True)
                HiC data on the Y chromsome is pretty messy, and if your cells are female is not valid. Skip it by default;

        '''
        if anchors: assert 'loc' in anchors.keys(), '"loc" key not found in anchors'
        assert min_dist, 'max_dist must be specified'
        assert max_dist, 'min_dist must be specified'

        # basically you just step through all diagonals, and take the histogram;
        bin_span = math.ceil((max_dist - min_dist) / self['bin_size'])
        min_bin = math.ceil(min_dist / self['bin_size'])

        hist = numpy.zeros(bin_span)

        if anchors: # selected anchors only
            used_bins = set([]) # Don't use same bin twice;

            anchors = anchors['loc']
            p = progressbar(len(anchors))
            for lidx, loc in enumerate(anchors):
                chrom = 'chr{}'.format(loc.loc['chr'])
                cpt = ((loc.loc['left'] + loc.loc['right']) // 2) // self['bin_size']

                if cpt in used_bins:
                    continue
                used_bins.add(cpt)

                left_top_slice = cpt+min_bin+bin_span
                if left_top_slice < self.mats[chrom].shape[1]:
                    left = numpy.ravel(self.mats[chrom][cpt+min_bin:cpt+min_bin+bin_span, cpt:cpt+1])
                    #print('left', left)
                    hist += left

                up_top_slice = cpt-min_bin-bin_span
                if up_top_slice > 1:
                    up = self.mats[chrom][cpt:cpt+1, up_top_slice:cpt-min_bin][0][::-1]
                    #print('up', up)
                    hist += up
                p.update(lidx)

        else: # Whole genome;
            p = progressbar(len(self.mats))
            for cidx, chrom in enumerate(self.mats):
                if chrom == 'chrY' and skip_Y_chromosome:
                    continue # The Y is often a mess skip it;
                for cpt in range(self.mats[chrom].shape[0]):
                    # Don't add the edges, as that includes things like telomeres;

                    ri_top_slice = cpt+min_bin+bin_span
                    if ri_top_slice < self.mats[chrom].shape[1]:
                        ri = numpy.ravel(self.mats[chrom][cpt+min_bin:cpt+min_bin+bin_span, cpt:cpt+1])
                        #print('ri', cpt, cpt+min_bin, cpt+min_bin+bin_span)
                        #print(ri)
                        hist += ri

                    up_top_slice = cpt-min_bin-bin_span
                    if up_top_slice > 1:
                        up = self.mats[chrom][cpt:cpt+1, up_top_slice:cpt-min_bin][0][::-1]
                        #print('up', cpt, up_top_slice, cpt-min_bin)
                        #print(up)
                        hist += up
                p.update(cidx)

        x = range(min_dist, max_dist, self['bin_size'])
        log10x = [math.log10(i) for i in x]

        if filename:
            fig = self.draw.getfigure()
            ax = fig.add_subplot(111)

            hist = numpy.log10(hist)

            ax.plot(log10x, hist)

            self.draw.do_common_args(fig, *kargs)

            actual_filename = self.draw.savefigure(fig, filename, dpi=300)
            config.log.info('Saved {}'.format(actual_filename))

        return hist



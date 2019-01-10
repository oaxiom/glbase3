'''

hic.py

Analysis for HiC data.

TODO:
-----
. merge_hiccys could be adapted to measure variance, and so extract DE?

'''

import pickle, numpy, math
from operator import itemgetter
from shutil import copyfile

import matplotlib.cm as cm
from scipy import ndimage
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from .genelist import genelist
from . import config
from .draw import draw
from .progress import progressbar
from .location import location
from .format import minimal_bed

if config.H5PY_AVAIL:
    import h5py
else:
    raise AssertionError('Asked for a hic, but h5py is not avaialble')

def merge_hiccys(new_hic_filename, name, *hics):
    '''
    **Purpose**
        Take the mean values for two or more hic objects.
        Save the result as a new hiccy File

        This algorithm assumes you know what you are doing, and you are trying to merge
        hiccy objects with the

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
    assert new_hic_filename, 'You must specify a filename in new_hic_filename'
    assert name, 'You must specify a name'
    assert hics, 'a list of hics was not specified'
    assert isinstance(hics, tuple), 'hics must be a list of >=2'
    assert (len(hics) > 1), 'hics must be >=2 length'

    # For the first one to merge, do an OS copy for speed and setup
    copyfile(hics[0], new_hic_filename)

    newhic = hic(filename=new_hic_filename, name=name, new=False)
    newhic.readonly = False # Make it editable

    # I bind all the hics:
    hics = [hic(filename=f, name=f) for f in hics[1:]]

    for chrom in newhic.all_chrom_names:
        # config.log.info('Merging chrom "%s"' % chrom)

        data = numpy.array(newhic.mats[chrom])
        for h in hics:
            data += h.mats[chrom]

    data /= len(hics)
    newhic.mats[chrom] = data

    newhic.close()
    config.log.info('Merged %s matrices' % (len(hics)+1,))

# I'm not using the base_genelist class, so, you need to add in defs as needed.
class hic:
    def __init__(self, filename=None, name='', new=False, inter_chrom_only=True):
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
        self.readonly = True
        self.pca_valid = False
        self.tsne_trained = False
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
            self.all_chrom_names = [] # Made into a set later
            self.draw = draw()
        else: # old
            self.hdf5_handle = h5py.File(filename, 'r')
            # TODO: fetch it back out as a set:
            self.tad_lookup = None

            dat = self.hdf5_handle['all_chrom_names'].value
            dat = [i[0] for i in dat]# flatten
            self.all_chrom_names = [n.decode("ascii", "ignore") for n in dat]

            # Rescue bin_lookup from disk:
            # save the bin data as an emulated dict:
            self.bin_lookup_by_binID = {}
            self.bin_lookup_by_chrom = {}
            for chrom in self.all_chrom_names:
                flat_bin = self.hdf5_handle['bin_lookup/chrom_%s/bins' % chrom].value
                self.bin_lookup_by_chrom[chrom] = [(row[0], row[1], row[2], row[3]) for row in flat_bin]

                for row in flat_bin:
                    self.bin_lookup_by_binID[row[0]] = (chrom, row[1], row[2], row[3]) # row[3] == newBinID

            #print(self.bin_data)

            self.mats = {}
            for chrom in self.all_chrom_names:
                self.mats[chrom] = self.hdf5_handle['matrix_%s/mat' % chrom]

            self.draw = draw()
            config.log.info('Bound "%s" Hic file' % filename)
        return(None)

    def visit(self):
        self.hdf5_handle.visit(print)

    def __len__(self):
        return(self.num_bins)

    def __getitem__(self, index):
        """
        Confers:

        a = expn["condition_name"]

        and inherits normal genelist slicing behaviour
        """
        return(self.hdf5_handle.attrs[index])

    def keys(self):
        return(self.hdf5_handle.keys())

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

        for bin in self.bin_lookup_by_chrom[loc['chr']]:
            # bin = (binID, left, right)
            if bin[2] > loc['left'] and bin[0] < binLeft:
                binLeft = bin[0]
                locLeft = bin[1]
            if loc['right'] > bin[1] and bin[0] > binRight:
                binRight = bin[0]
                locRight = bin[2]

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

        newloc = location(chr=loc['chr'], left=locLeft+1, right=locRight-1) # the +1/-1 stops the bins from overlapping by 1bp, and solves a lot of problems

        binLeft = binLeft - mostLeft
        binRight = binRight - mostLeft
        mostRight = self.bin_lookup_by_binID[mostRight][3] # Convert to newBinID:
        mostLeft = self.bin_lookup_by_binID[mostLeft][3]
        assert binRight-binLeft > 2, 'the genome view (loc) is too small, and contains < 2 bins'

        return(binLeft, binRight, newloc, mostLeft, mostRight)

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

        return(mostLeft, mostRight)

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
        self.all_chrom_names = []
        oh = open(bed_file, 'r')
        for lin in oh:
            # format = chr, left, right, bin#
            lin = lin.strip().split('\t')

            chr = lin[0].replace('chr', '')
            try:
                chr = int(chr)
            except ValueError:
                pass # It's chrX etc.

            bins.append((chr, int(lin[1]), int(lin[2]), int(lin[3])))
            if chr not in self.all_chrom_names:
                self.all_chrom_names.append(chr)

        oh.close()
        self.all_chrom_names = set(self.all_chrom_names)
        self.hdf5_handle.attrs['num_bins'] = len(bins)

        config.log.info('Found %s bins' % self['num_bins'])
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
            grp.create_dataset('bins', (flat_bin.shape), dtype=int, data=flat_bin)

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

        return(True)

    def save_np3_column_matrix(self, filename):
        """
        **Purpose**
            Some other tools ask for a n+3 chromosome matrix,

            In the form
            chrom   left    right   0   0   0   0   0   0    .... #bins

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
            mat = self.mats[chrom].value
            bins = self.bin_lookup_by_chrom[chrom]
            for m, b in zip(mat, bins):
                lin = [chrom_name, b[1], b[2]] + list(m)
                lin = [str(i) for i in lin]
                oh.write('%s\n' % '\t'.join(lin))
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

    def heatmap(self, filename, chr=None, loc=None,
        bracket=None, colour_map=cm.inferno_r, **kargs):
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

            aspect (Optional, default='square')
                image aspect for the heatmap

            colbar_label (Optional, default='log2(Density)')
                Label for the colour bar.

            colour_map (Optional, default=cm.inferno_r)
                matplotlib colour map to use to colour.

            bracket (Optional, default=None)
                clip the data within the ranges [low, high]

        **Returns**
            None, and a file in filename.

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

        if not "aspect" in kargs:
            kargs["aspect"] = "square"
        if not "colbar_label" in kargs:
            kargs["colbar_label"] = "log2(Density)"

        fig = self.draw.getfigure(**kargs)

        # positions of the items in the plot:
        heatmap_location  = [0.05,   0.01,   0.90,   0.90]
        scalebar_location = [0.05,   0.97,   0.90,   0.02]

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
            data[data == -numpy.inf] = 0
            vmin = data.min()
            vmax = data.max()

        # ---------------- (heatmap) -----------------------
        ax3 = fig.add_subplot(121)

        ax3.set_position(heatmap_location) # must be done early for imshow
        hm = ax3.imshow(data, cmap=colour_map, vmin=vmin, vmax=vmax, aspect="auto",
            origin='lower', extent=[0, data.shape[1], 0, data.shape[0]],
            interpolation=config.get_interpolation_mode())

        #ax3.set_frame_on(True)
        ax3.set_position(heatmap_location)
        ax3.set_xlim([0,data.shape[1]])
        ax3.set_ylim([0,data.shape[0]])
        ax3.set_yticklabels("")
        ax3.set_xticklabels("")

        ax3.tick_params(top=False, bottom=False, left=False, right=False)
        [t.set_fontsize(5) for t in ax3.get_yticklabels()] # generally has to go last.
        [t.set_fontsize(5) for t in ax3.get_xticklabels()]

        ax0 = fig.add_subplot(122)
        ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)

        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0, cmap=colour_map)
        cb.set_label(kargs["colbar_label"])
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
            localLeft, localRight, loc, _ = self.__find_binID_spans(loc)

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
            interpolation=config.get_interpolation_mode())

        ax.set_xlim([0,data.shape[1]])
        ax.set_ylim([data.shape[1]/2,(data.shape[1]/2)+data.shape[1]/10]) # get the middle start position, then about 1/10 of the way up
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

        if chr:
            data = self.mats[str(chr).replace('chr', '')]
            if self.tad_lookup:
                tad_calls = self.tad_lookup[loc['chrom']]
        elif loc:
            if not isinstance(loc, location):
                loc = location(loc)

            # Need to get the binIDs and the updated location span
            localLeft, localRight, loc, mostLeft, mostRight = self.__find_binID_spans(loc)
            data = self.mats[loc['chr']][localLeft:localRight, localLeft:localRight]

            # I just assume the bins match
            this_chrom = [0] * (localRight-localLeft+1)
            cindex = expn.getConditionNames().index(expn_cond_name)
            for i in expn:
                if i['loc']['chr'] == loc['chr']:
                    # Take the old BinID and convert it to the new binID:
                    # If inside this part of the chromosome:
                    if loc.qcollide(i['loc']):
                        local_bin_num = (self.bin_lookup_by_binID[i['bin#']][3] - mostLeft) - localLeft
                        this_chrom[local_bin_num] = i['conditions'][cindex]
            plot_y = this_chrom
            plot_x = numpy.arange(0, len(plot_y))

            tad_calls = None
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
            interpolation=config.get_interpolation_mode())

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

        print(max)
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
            perc_weights=self.get_loading_percents()
        elif mode == 'tsne':
            assert self.tsne_trained, 'tSNE model has not been trained, or the model was trained on the wrong chromosome, use tsne()'
            labels = self.__pcalabels
            xdata = self.npos[:, 0]
            ydata = self.npos[:, 1]
            mode = 'tSNE '
            perc_weights=None
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

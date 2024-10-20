"""
flat_heat, part of chipFish, and glbase3

2024 oAxiom

x-axis is the genome, y-axis is the read size. Useful for V plots for ATAC-seq, etc.


"""

import pickle, sys, os, math, zlib
from operator import itemgetter

from .location import location

from .draw import draw
from .progress import progressbar
from .genelist import Genelist
from .flat_track import flat_track, positive_strand_labels, negative_strand_labels
from . import config
from . import utils

import numpy
import matplotlib.pyplot as plot
import matplotlib.cm as cm
import scipy.interpolate as scipyinterp

if config.H5PY_AVAIL:
    import h5py
else:
    raise AssertionError('flatheat requires h5py, but it is not available')

# TODO: Most of this can be inherited. just get() and heatmap() are actually differnet.
class flat_heat:
    def __init__(self,
        filename:str=None,
        name:str=None,
        new:str=False,
        ymax:int=None,
        ybins:int=50, 
        ):
        """
        **Purpose**
            track definition, used for things like sequence reads across the genome

        **Arguments**
            name (string)
                name for the track, if not specified, it takes it from the meta_data
                If a new bd then it defaults to filename

            filename (string)
                directory location of the track file.
                only respected if dir_name is set.

            ybins (int, optional, default=50)
                number of bins for the y-axis
            
            ymax (int, Required if new=True)
                size of the y-axis in base pairs.

        """
        assert filename, 'You must specify a filename'

        self.readonly = True
        self._draw = None

        if new:
            self.readonly = False
            self.ymax = ymax
            self.ybins = ybins
            self.hdf5_handle = h5py.File(filename, 'w')  # This should be x to stop overwriting an exisiting file

            # h5py hierarchy (Made at load_time)
            # Things stored:
            self.hdf5_handle.attrs['name'] = name
            self.hdf5_handle.attrs['bin_format'] = 'float'
            self.hdf5_handle.attrs['num_reads'] = 0
            self.hdf5_handle.attrs['ybins'] = ybins
            self.hdf5_handle.attrs['ymax'] = ymax
            self.hdf5_handle.attrs['version'] = '1.0' # First version of the flatheat format
            self.chrom_names = []
            self.meta_data = self.hdf5_handle.attrs

            self.draw = draw()
            config.log.info(f'Setup new "{filename}" flat file')

        else: # old
            try:
                self.hdf5_handle = h5py.File(filename, 'r')
            except OSError as e:
                print(e)
                if 'file signature not found' in str(e):
                    config.log.error(f'This file "{filename}" is not a flatheat file')
                elif 'Unable to open file' in str(e): # file not found
                    config.log.error('File not found: {}'.format(filename))
                sys.exit()

            self.meta_data = self.hdf5_handle.attrs

            self.ybins = self.meta_data['ybins']
            self.ymax = self.meta_data['ymax']

            self.chrom_names = [i[0] for i in self.hdf5_handle['all_chrom_names'][()]]# flatten
            self.chrom_names = [n.decode("ascii", "ignore") for n in self.chrom_names]

            self.mats = {}
            for chrom in self.chrom_names:
                self.mats[chrom] = self.hdf5_handle[f'matrix_{chrom}/mat']

            self.draw = draw()
            config.log.info(f'Bound "{filename}" flat file')

        # NAme override:
        if name:
            self.name = name
        else:
            self.name = self.hdf5_handle.attrs['name']

    def __repr__(self):
        return "glbase.flatheat"

    def keys(self):
        return self.hdf5_handle.keys()

    def _optimiseData(self):
        pass

    def close(self):
        self.hdf5_handle.flush()
        self.hdf5_handle.close()

    def __getitem__(self, index):
        """
        Confers:

        a = flat["name"] to collect attrs

        """
        if index == 'name':
            return self.name
        return self.hdf5_handle.attrs[index]

    def add_chromosome_array(self, chromosome=None, arr=None):
        '''
        **Purpose**
            Support for loading of complete chromosome arrays.

            This will overwrite any existing data.

        **Arguments**
            chromsome (Required)
                chromsome to insert array into

            arr (Required)
                a Python list of scores for the entire chromosome. Should start at position 0,
                and extend for the complete chromosomes.
        '''
        assert chromosome, 'You must specify a chromosome'
        assert isinstance(arr, numpy.ndarray), 'arr is not a numpy array'

        if 'chr' not in chromosome:
            chromosome = 'chr{}'.format(chromosome)

        grp = self.hdf5_handle.create_group(f'matrix_{chromosome}')
        grp.create_dataset('mat', arr.shape, dtype=float, data=arr, chunks=True, compression='lzf')
        config.log.info(f'Added chrom={chromosome} to table')

        self.chrom_names.append(chromosome)
        return None

    def get_all_chrom_names(self):
        return self.chrom_names

    def set_total_num_reads(self, num_reads):
        """
        **Purpose**
            Get the total number of reads in this library,
            generally for normalization purposes. Number of tags must be set at creation time,
            not always avaialable for all flats
        """
        assert not self.readonly, 'Trying to set_total_num_reads, but read only'
        self.hdf5_handle.attrs['num_reads'] = num_reads
        return None

    def get_total_num_reads(self):
        """
        **Purpose**
            Get the total number of reads in this library,
            generally for normalization purposes. Number of tags must be set at creation time,
            not always available for all flats
        """
        return self.hdf5_handle.attrs['num_reads']

    def get(self, loc, c=None, left=None, rite=None, strand="+", mask_zero=False, **kargs):
        """
        **Purpose**
            get the data between location 'loc'

        **Arguments**
            loc (Required)
                a valid location or string location.

            strand (Optional, default = "+")
                strand, but only valid for stranded tracks

            mask_zero (Optional, default=False)
                return a masked numpy array with zeros masked out.

        **Returns**
            an 'array('i', [0, 1, 2 ... n])' contiginous array
        """
        if loc:
            try:
                pass
            except TypeError: # probably a location string. try to cooerce
                loc = location(loc=loc)
                # Don't catch any exceptions here. Should break.
            c = str(loc.chrom)
            left = int(loc.left) // 10
            rite = int(loc.right) // 10

        if 'chr' not in c:
            c = f'chr{c}'

        ret_array = self.mats[c][left:rite,]

        if mask_zero:
            mask = []
            for dd in ret_array:
                if int(dd*10000000) == 0: # delete if 10 sig figs close to 0
                    mask.append(1)
                else:
                    mask.append(0)
            # This may produce a warning, but apparently I can safely ignore it
            ret_array = numpy.ma.masked_array(ret_array, mask=mask)

        return ret_array

    def get_array_chromosome(self, chrom=None, **kargs): # kargs for compat with trk
        """
        **Purpose**
            Get the enrire array for the chromosome chrom. 
            Note that the x-axis is one column per every 10 bp. So you need to
            multiply coordinates by 10.
            
            This method should probably be taken internal.

        **Arguments**
            chrom (Required)
                chromsome name to collect.

        """
        if 'chr' not in chrom:
            chrom = 'chr{}'.format(chrom)

        assert chrom in self.chrom_names, 'chromosome {0} not in this flat, available chroms are: "{1}"'.format(chrom, ', '.join(self.chrom_names))
        return self.mats[chrom]

    def finalise(self):
        """
        Kept here for placeholder work and API compatability
        """
        assert not self.readonly, 'readonly!'
        # Write the chorm names to the hdf5
        dat = [str(n).encode("ascii", "ignore") for n in self.chrom_names]
        self.hdf5_handle.create_dataset('all_chrom_names', (len(self.chrom_names), 1), 'S10', dat)

    def pileup_heatmap(self,
        genelists=None,
        filename=None,
        scaled=None,
        scaled_view_fraction=0.5,
        window_size=None,
        average=True,
        background=None,
        mask_zero=False,
        respect_strand=True,
        norm_by_read_count=True,
        **kargs):
        """
        **Purpose**
            Draw a set of pileup count scores (averages or cumulative scores)

        **Arguments**
            genelists
                A list of genelists with a "loc" key

            filename
                The filename to save the image to

            scaled (Required)
                If True, then the pileup will be formatted so that the central <scaled_view_fraction>
                will be the coordinates from the genelists, and the flanking regions will be taken
                from  window_size |---|----|---| and fill the remaining space.

                If False, draw a pileup, but center it on the middle of the coordinates in genelist.
                No scaling is performed and the flanking parts are taken from the window_size

            scaled_view_fraction (Optional, default=0.5)
                Only used if scaled == True. This specifies the fraction of the pileup to use
                for the central region (See scaled)

            window_size (Optional, default=None)
                the number of base pairs to use around the centre of the location (if scaled)

                If set to None then it will use the location as specified.

            average (Optional, default=True)
                use the average score if set to True (i.e. divide by the number of items
                in genelist), or use the cumulative score (the total) if set to False

            background (Optional)
                You can supply a list of background coordinates if you want. The lists
                must contain a "loc" key.

            mask_zero (Optional, default=False)
                flat_tracks are continuous and must have a value at all base pairs. I fill
                that value with 0
                However, in some flat_tracks 0 actually means 'no data' and I want to ignore
                those values in the pileup. If that is the case for you, set mask_zero to
                True.

            <other graph-related args>
            pileup also respects:
                xlabel - x-axis label
                ylabel - y-axis label
                title  - title

            respect_strand (Optional, default=False)
                If available, respect the orientation of the strand from the genelist.
                This is useful if you are, say, using the TSS's and want to maintain the
                orientation with respect to the transcription direction.

            norm_by_read_count (Optional, default=True)
                If you are not using a norm_factor for this library then you probably want to set this to True.
                It will divide the resulting number of reads by the total number of reads,
                i.e. it will account for differences in library sizes.

        **Returns**
            (data, background)
            The data from the line graph.

            Retuns a dict in data, with a key for each data[genelist.name]
            background returns an average of the list of background peaks.

        """
        # Common cleanup
        assert genelists, "genelists is None?"

        # flats have lazy setup of draw:
        if not self._draw:
            self._draw = draw(self)

        if not isinstance(genelists, list):
            genelists = [genelists] # make a one item'd list

        if background:
            if not isinstance(background, list):
                background = [background] # make a one item'd list

        read_count = 1.0
        if norm_by_read_count:
            read_count = float(self.get_total_num_reads())
            if read_count <= 0:
                raise AssertionError('norm_by_read_count=True, but this flat_track has no total number of reads')

        all_hists = {}

        fig = self._draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        x = None
        if window_size:
            x = numpy.arange(window_size*2) - (window_size*2)//2

        if window_size:
            loc_span = window_size*2
        else:
            loc_span = len(genelists[0].linearData[0]["loc"]) # I have to assume all locs are identical.

        available_chroms = list(self.mats.keys())
        available_chroms += [c.replace('chr', '') for c in available_chroms] # help with name mangling
        __already_warned = []

        for gl in genelists:
            if window_size:
                hist = numpy.zeros(window_size*2)
                counts = numpy.zeros(window_size*2)
                gl = gl.pointify().expand('loc', window_size)
            else:
                x = numpy.arange(loc_span) # - loc_span//2
                hist = numpy.zeros(loc_span)
                counts = numpy.zeros(loc_span) # used to get the average.

            for i in gl:
                if i['loc']['chr'] not in available_chroms:
                    if i['loc']['chr'] not in __already_warned:
                        config.log.warning('Asked for chromosome {} but not in this flat_track, skipping'.format(i['loc']['chr']))
                        __already_warned.append(i['loc']['chr'])
                    continue

                a = self.get(i["loc"])#[0:window_size*2] # mask_zero is NOT asked of here. because I need to ignore zeros for the average calculation (below)

                if respect_strand:
                    # positive strand is always correct, so I leave as is.
                    # For the reverse strand all I have to do is flip the array.
                    if i["strand"] in negative_strand_labels:
                        a = a[::-1]

                if a.any(): # It's possible that get() will return nothing
                    # For example if you send bad chromosome names or the locations are nonsensical (things like:
                    # chr9_GL000201_RANDOM:-500-1500
                    # Check for a block miss:
                    if len(a) < loc_span: # This should be a very rare case...
                        config.log.warning('Block miss (short)')
                        num_missing = loc_span - len(a)
                        ad = numpy.zeros(num_missing, dtype=float)
                        a = numpy.append(a, ad)

                    hist += a

                if mask_zero: # surely a better way of doing this...
                    t = numpy.zeros(loc_span)
                    for ee, xx in enumerate(a):
                        if xx > 0:
                            t[ee] = 1.0
                    counts += t

            if average and mask_zero:
                hist /= counts
            elif average and not mask_zero:
                hist /= len(gl)

            if norm_by_read_count:
                hist /= read_count

            ax.plot(x, hist, label=gl.name, alpha=0.7)
            all_hists[gl.name] = hist

        bkgd = None
        if background:
            if window_size:
                bkgd = numpy.zeros(window_size*2)
                counts = numpy.zeros(window_size*2)
            else:
                bkgd = numpy.zeros(loc_span)
                counts = numpy.zeros(loc_span)

            bkgd_items = 0
            p = progressbar(len(background))
            for i, back in enumerate(background):
                for b in back:
                    if b['loc']['chr'] not in available_chroms:
                        if b['loc']['chr'] not in __already_warned:
                            config.log.warning('Asked for {} chromosome but not in this flat_track, skipping'.format(b['loc']['chr']))
                            __already_warned.append(b['loc']['chr'])
                        continue

                    if window_size:
                        l = b["loc"].pointify()
                        l = l.expand(window_size)
                        a = self.get(l)[0:window_size*2]
                    else:
                        a = self.get(b["loc"])[0:loc_span]
                    bkgd_items += 1

                    if respect_strand:
                        # positive strand is always correct, so I leave as is.
                        # For the reverse strand all I have to do is flip the array.
                        if b["strand"] in negative_strand_labels:
                            a = a[::-1]

                    bkgd += a
                    if mask_zero:
                        t = numpy.zeros(loc_span)
                        for ee, xx in enumerate(a):
                            if xx > 0:
                                t[ee] = 1.0
                        counts += t

                if average and mask_zero:
                    bkgd /= counts
                elif average and not mask_zero:
                    bkgd /= bkgd_items

                if norm_by_read_count:
                    hist /= read_count

                if i == 0: # add only a single legend.
                    ax.plot(x, bkgd, color="grey", alpha=0.3, label="Random Background")
                else:
                    ax.plot(x, bkgd, color="grey", alpha=0.3)

                # reset arrays
                bkgd = numpy.zeros(len(bkgd))
                counts = numpy.zeros(len(counts))

                p.update(i)

        else:
            bkgd = None

        leg = ax.legend()
        [t.set_fontsize(3) for t in leg.get_texts()]
        ax.set_ylabel("Magnitude")

        if window_size:
            ax.set_xlabel("Base pairs around centre (bp)")
            ax.axvline(0, ls=":", color="grey")
        else:
            ax.set_xlabel('Base pairs (bp)')

        self._draw.do_common_args(ax, **kargs)

        actual_filename = self._draw.savefigure(fig, filename)

        config.log.info("pileup(scaled=False): Saved '{}'".format(actual_filename))

        return (all_hists, bkgd)
    

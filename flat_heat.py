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
            # TODO: self.hdf5_handle.attrs['flat_type'] = 'flat_heat'
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

    def get(self, loc, c=None, left=None, rite=None, strand="+", **kargs):
        """
        **Purpose**
            get the data between location 'loc'

        **Arguments**
            loc (Required)
                a valid location or string location.

            strand (Optional, default = "+")
                strand, but only valid for stranded tracks

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

        ret_array = self.mats[c][left:rite,0:self.ybins]

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
        genelist=None,
        filename=None,
        window_size:int=None,
        average=True,
        respect_strand=True,
        norm_by_read_count=True,
        colour_map = cm.BrBG,
        fast:bool = False,
        **kargs):
        """
        **Purpose**
            Draw a set of pileup count scores (averages or cumulative scores)

        **Arguments**
            genelists
                A genelist with a "loc" key

            filename
                The filename to save the image to

            window_size (Required, default=None)
                the number of base pairs to use around the centre of the location

                If set to None then it will use the location as specified.

            bracket (Optional, default=[min, max])
                set your own bracket for the heatmap.

            average (Optional, default=True)
                use the average score if set to True (i.e. divide by the number of items
                in genelist), or use the cumulative score (the total) if set to False

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

            colour_map (Optional, default=cm.BrBG)
                A matplotlib colormap

        **Returns**
            (data, background)
            The data from the line graph.

            Retuns a dict in data, with a key for each data[genelist.name]
            background returns an average of the list of background peaks.

        """
        # Common cleanup
        assert genelist, "genelist is None?"
        assert window_size, 'You must specify a window size around the center point'
        assert filename, 'You must specify a filename to save the heatmap to'

        # flats have lazy setup of draw:
        if not self._draw:
            self._draw = draw(self)

        read_count = 1.0
        if norm_by_read_count:
            read_count = float(self.get_total_num_reads())
            if read_count <= 0:
                raise AssertionError('norm_by_read_count=True, but this flat_heat has no total number of reads')

        all_hists = {}

        x = numpy.arange(window_size*2) - (window_size*2)//2
        loc_span = window_size*2 // 10
        half_loc = loc_span // 2

        available_chroms = list(self.mats.keys())
        available_chroms += [c.replace('chr', '') for c in available_chroms] # help with name mangling
        __already_warned = []

        gl = genelist

        hist = numpy.zeros((loc_span, self.ybins), dtype=numpy.float64)

        gl = gl.pointify().expand('loc', window_size)

        p = progressbar(len(gl))
        for idx, i in enumerate(gl):
            if i['loc'].chrom not in available_chroms:
                if i['loc'].chrom not in __already_warned:
                    config.log.warning('Asked for chromosome {} but not in this flat_track, skipping'.format(i['loc'].chrom))
                    __already_warned.append(i['loc'].chrom)
                continue

            left = i['loc'].left // 10 # I only use pointfied
            rite = left + half_loc
            left -= half_loc

            #a = self.get(i["loc"]) # No method overhead, this part is too slow already...
            a = self.mats[f"chr{i['loc'].chrom}"][left:rite, 0:self.ybins]

            if respect_strand:
                # positive strand is always correct, so I leave as is.
                # For the reverse strand all I have to do is flip the array.
                if i["strand"] in negative_strand_labels:
                    a = a[::-1,] # flip x-axis

            if a.any(): # It's possible that get() will return an array full of zeros, if so, skip them
                # For example if you send bad chromosome names or the locations are nonsensical (things like:
                # chr9_GL000201_RANDOM:-500-1500
                # Check for a block miss:
                if a.shape[0] < loc_span: # This should be a very rare case...
                    config.log.warning(f'Block miss (short) {i["loc"]}')
                    # TODO: Fill in? Probably better to skip
                    continue

                hist += a

            p.update(idx)

        if average:
            hist /= len(gl)

        if norm_by_read_count:
            hist /= read_count

        if bracket: # done here so clustering is performed on bracketed data
            vmin = bracket[0]
            vmax = bracket[1]
        else:
            vmin = hist.min()
            vmax = hist.max()

        fig = self._draw.getfigure(**kargs)
        ax1 = fig.add_subplot(211)
        hm = ax1.imshow(hist.T, cmap=colour_map, vmin=vmin, vmax=vmax, aspect="auto",
            origin='lower', extent=[0, hist.shape[0], 0, hist.shape[1]],
            interpolation=config.get_interpolation_mode(filename))

        ax0 = fig.add_subplot(212)
        #ax0.set_position(scalebar_location)
        ax0.set_frame_on(False)
        cb = fig.colorbar(hm, orientation="horizontal", cax=ax0)

        ax1.set_xlabel("Base pairs around centre (bp)")
        ax1.axvline(window_size // 10, ls=":", color="grey")

        self._draw.do_common_args(ax1, **kargs)

        actual_filename = self._draw.savefigure(fig, filename)

        config.log.info("pileup_heatmap: Saved '{}'".format(actual_filename))

        return hist


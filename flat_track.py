"""
flat_track2, part of chipFish, and glbase3

2020 oAxiom

New reimplementation using hdf5.

"""

import pickle, sys, os, struct, configparser, math, zlib

from .location import location
from .data import positive_strand_labels, negative_strand_labels

from array import array
from .draw import draw
from .progress import progressbar
from . import config

import numpy
import matplotlib.pyplot as plot
import matplotlib.cm as cm

if config.H5PY_AVAIL:
    import h5py
else:
    raise AssertionError('flat_track now requires h5py, but it is not avaialble')

positive_strand_labels = frozenset(["+", "1", "f", "F", 1])
negative_strand_labels = frozenset(["-", "0", "r", "R", -1, 0, "-1"])

class flat_track():
    def __init__(self, name=None, new=False, filename=None, bin_format='f'):
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

            bin_format (Optional, default=None, required if new=False)
                the format to store the data in,
                This is the same format as the python array, so "i" = integer,
                "f" = float

        """
        self.readonly = True

        if new:
            self.readonly = False
            self.hdf5_handle = h5py.File(filename, 'w')  # This should be x to stop overwriting an exisiting file

            # h5py hierarchy (Made at load_time)
            # Things stored:
            self.hdf5_handle.attrs['name'] = name
            self.hdf5_handle.attrs['bin_format'] = bin_format
            self.hdf5_handle.attrs['num_reads'] = 0
            self.hdf5_handle.attrs['version'] = '5.0' # Fifth major revision in the flat format
            self.chrom_names = []
            self.meta_data = self.hdf5_handle.attrs

            self.draw = draw()
            config.log.info('Setup new "{0}" flat file' .format(filename))

        else: # old
            try:
                self.hdf5_handle = h5py.File(filename, 'r')
            except OSError as e:
                print(e)
                # Check it's the h5py old error
                if 'file signature not found' in str(e):
                    config.log.error('Either this is not a flat_track file, or it is the old version 4 flat tracks')
                    config.log.error('Please regenerate your flat_tracks. The new v5 format is ~4x faster on reads')
                    config.log.error('and about the same to generate the track. File size is about 2x bigger')
                elif 'Unable to open file' in str(e): # file not found
                    config.log.error('File not found: {0}'.format(filename))
                sys.exit()

            self.meta_data = self.hdf5_handle.attrs

            self.chrom_names = [i[0] for i in self.hdf5_handle['all_chrom_names'][()]]# flatten
            self.chrom_names = [n.decode("ascii", "ignore") for n in self.chrom_names]

            self.mats = {}
            for chrom in self.chrom_names:
                self.mats[chrom] = self.hdf5_handle['matrix_%s/mat' % chrom]

            self.draw = draw()
            config.log.info('Bound "%s" flat file' % filename)

    def __repr__(self):
        return "glbase.flat_track"

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

        a = flat["condition_name"]

        """
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

        # TODO: Check chrom does not already exist

        grp = self.hdf5_handle.create_group('matrix_{0}'.format(chromosome))
        grp.create_dataset('mat', arr.shape, dtype=numpy.float32, data=arr, chunks=True, compression='lzf')
        config.log.info('Added chrom={0} to table'.format(chromosome))

        self.chrom_names.append(chromosome)

        return None

    def get_all_chrom_names(self):
        return self.chrom_names

    def set_total_num_reads(self, num_reads):
        """
        **Purpose**
            Get the total number of reads in this library,
            generally for normalization purposes. Number of tags must be set at creation time,
            not always avaialble for all flats
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
                if loc["chr"]: pass
            except TypeError: # probably a location string. try to cooerce
                loc = location(loc=loc)
                # Don't catch any exceptions here. Should break.
            c = str(loc['chr'])
            left = int(loc['left'])
            rite = int(loc['right'])

        if 'chr' not in c:
            c = 'chr{0}'.format(c)

        ret_array = self.mats[c][left:rite]

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

            Shouldn't be needed in this implementation

        **Arguments**
            chrom (Required)
                chromsome name to collect.

        """
        if 'chr' not in chrom:
            chrom = 'chr{0}'.format(chrom)

        assert chrom in self.chrom_names, 'chromosome {0} not in this flat, available chroms are: "{1}"'.format(chrom, ', '.join(self.chrom_names))
        return self.mats[chrom]

    def finalise(self):
        """
        Kept here for placeholder work and API compatability
        """
        assert not self.readonly, 'readonly!'
        # Write the chorm names to the hdf5
        print(self.chrom_names)
        dat = [str(n).encode("ascii", "ignore") for n in self.chrom_names]
        self.hdf5_handle.create_dataset('all_chrom_names', (len(self.chrom_names), 1), 'S10', dat)

    def pileup(self, genelists=None, filename=None, window_size=None, average=True,
        background=None, mask_zero=False, respect_strand=True, norm_by_read_count=False, **kargs):
        """
        **Purpose**
            Draw a set of pileup count scores (averages or cumulative scores)

        **Arguments**
            genelists
                A list of genelist with a "loc" key

            filename
                The filename to save the image to

            window_size (Optional, default=None)
                the number of base pairs to use around the centre of the location
                If set to None then it will use the location as specified.

            average (Optional, default=True)
                use the average score if set to True (i.e. divide by the number of items)
                Or use the cumulative score (the total) if set to False

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

            norm_by_read_count (Optional, default=False)
                If you are not using a norm_factor for this library then you probably want to set this to True.
                It will divide the resulting number of reads by the total number of reads,
                i.e. it will account for differences in library sizes.

        **Returns**
            (data, back)
            The data from the line graph.
            back will be the average of the list of background peaks, but data
            will be the last entry in the peaklists (if a list of peaks) or will correspond to the
            only peaklist provided. e.g.:

            data, back = pileup(genelists=[data1, data2, THISDATAWILLBERETURNED] ...)

            or:

            data, back = pileup(genelists=THISDATAWILLBERETURNED ...)
        """
        assert genelists, "genelists is None?"

        if not isinstance(genelists, list):
            genelists = [genelists] # make a one item'd list

        if background:
            if not isinstance(background, list):
                background = [background] # make a one item'd list

        read_count = 1.0
        if norm_by_read_count:
            read_count = float(self.get_total_num_reads())

        all_hists = {}

        # flats have lazy setup of draw:
        if not self._draw:
            self._draw = draw(self)

        fig = self._draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        x = None
        if window_size:
            x = numpy.arange(window_size*2) - (window_size*2)//2

        if window_size:
            loc_span = window_size*2
        else:
            loc_span = len(genelists[0].linearData[0]["loc"]) # I have to assume all locs are identical.

        for gl in genelists:
            if window_size:
                hist = numpy.zeros(window_size*2)
                counts = numpy.zeros(window_size*2)
                gl = gl.pointify().expand('loc', window_size)
            else:
                x = numpy.arange(loc_span) - loc_span//2
                hist = numpy.zeros(loc_span)
                counts = numpy.zeros(loc_span) # used to get the average.

            for i in gl:
                a = self.get(i["loc"])#[0:window_size*2] # mask_zero is NOT asked of here. because I need to ignore zeros for the average calculation (below)

                if respect_strand:
                    # positive strand is always correct, so I leave as is.
                    # For the reverse strand all I have to do is flip the array.
                    if i["strand"] in negative_strand_labels:
                        a = a[::-1]

                if a: # It's possible that get() will return nothing
                    # For example if you send bad chromosome names or the locations are nonsensical (things like:
                    # chr9_GL000201_RANDOM:-500-1500
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

        ax.axvline(0, ls=":", color="grey")

        leg = ax.legend()
        [t.set_fontsize(7) for t in leg.get_texts()]
        ax.set_ylabel("Magnitude")
        ax.set_xlabel("Base pairs around centre (bp)")

        self._draw.do_common_args(ax, **kargs)

        actual_filename = self._draw.savefigure(fig, filename)

        config.log.info("pileup(): Saved '%s'" % actual_filename)

        return(all_hists, bkgd)

    def heatmap(self, filename=None, genelist=None, distance=1000, read_extend=200, log=2,
        bins=200, sort_by_intensity=True, raw_heatmap_filename=None, bracket=None,
        pointify=True, respect_strand=False, cmap=cm.plasma, norm_by_read_count=False,
        log_pad=None, imshow=True,
        **kargs):
        """
        **Purpose**
            Draw a heatmap of the seq tag density drawn from a genelist with a "loc" key.

        **Arguments**
            genelist (Required)
                a genelist with a 'loc' key.

            filename (Optional)
                filename to save the heatmap into
                Can be set to None if you don't want the png heatmap.

            raw_heatmap_filename (Optional)
                Save a tsv file that contains the heatmap values for each row of the genelist.

            distance (Optional, default=1000)
                Number of base pairs around the location to extend the search

            pointify (Optional, default=True)
                Take the middle point of the 'loc' in the genelist, used in combination with distance

                You can set this to False and distance to 0 and heatmap will use whatever locations are specified in
                the genelist. Note that the locations must all be the same lengths.

            bins (Optional, default=100)
                Number of bins to use. (i.e. the number of columns)

            read_extend (Optional, default=200)
                number of base pairs to extend the sequence tag by.

                NOTE: This value is ignored by flat-tracks

            respect_strand (Optional, default=False)
                Use the strand in the genelist to orientate each genomic location (for example when doing TSS
                or TTS).

            log (Optional, default=2)
                log transform the data, optional, the default is to transform by log base 2.
                Note that this parameter only supports "e", 2, and 10 for bases for log
                if set to None no log transform takes place.

            log_pad (Optional, default=None)
                a Pad value for the log. If None then heatmap() guesses from the min of 0.1 or data.min()

            norm_by_read_count (Optional, default=False)
                Normalise the heatmap by the total number of tags in the seq library.
                Should not be used if the norm_factor system is being used.

            sort_by_intensity (Optional, default=True)
                sort the heatmap so that the most intense is at the top and the least at
                the bottom of the heatmap.

            cmap (Optional, default=?)
                A matplotlib cmap to use for coloring the heatmap

            imshow (Optional, default=True)
                Embed the heatmap as an image inside a vector file.

        **Results**
            file in filename and the heatmap table, and the
            'sorted_original_genelist' (a copy of the genelist) sorted into the same order
            as the heatmap
        """
        assert genelist, "must provide a genelist"
        assert "loc" in list(genelist.keys()), "appears genelist has no 'loc' key"
        assert "left" in list(genelist.linearData[0]["loc"].keys()), "appears the loc key data is malformed"
        assert log in ("e", math.e, 2, 10, None), "this 'log' base not supported"

        table = []

        gl_sorted = genelist.deepcopy()
        #gl_sorted.sort('loc') # No need to sort anymore;
        all_locs = gl_sorted['loc']

        if respect_strand:
            strands = gl_sorted['strand']
        else:
            strands = ['+'] * len(all_locs) # Fake strand labels for code simplicity

        curr_cached_chrom = None
        cached_chrom = None
        bin_size = None

        number_of_tags_in_library = False # For code laziness
        if norm_by_read_count:
            number_of_tags_in_library = self.get_total_num_reads() / float(1e6) # 1e6 for nice numbers

        p = progressbar(len(genelist))
        for idx, read in enumerate(zip(all_locs, strands)):
            l = read[0]
            if pointify:
                l = l.pointify()
            if distance:
                l = l.expand(distance)

            if not bin_size: # I need to sample from an example size.
                # assume all sizes are the same!
                expected_width = len(l)
                bin_size = int(expected_width / float(bins))
                #print 'Binsize:', bin_size
            '''
            # I can dispose and free memory as the locations are now sorted:
            # See if the read_extend is already in the cache:
            if l['chr'] != curr_cached_chrom:
                curr_cached_chrom = l['chr']
                # flat_tracks ignore read_extend, but tracks require it
                cached_chrom = self.get_array_chromosome(l['chr'], read_extend=read_extend) # Will hit the DB if not already in cache
                # if we are a flat_track, we need to put it to a numpy array:
                if isinstance(cached_chrom, list):
                    cached_chrom = numpy.array(cached_chrom, dtype=numpy.float32)

            actual_width = cached_chrom.shape[0] - l['left']

            if actual_width < expected_width:
                #print "!", a.shape
                continue
            else:
                a = cached_chrom[l['left']:l['right']]
            '''
            '''     # This is not working, it pads to the wrong size and eventually leads to obscure errors
            if l['right'] > cached_chrom.shape[0]:
                # Trouble, need to fill in the part of the array with zeros
                # Possible error here it l['left'] is also off the array?

                a = cached_chrom[l['left']:cached_chrom.shape[0]] # stop wrap around
                a = numpy.pad(a, (0,expected_width-actual_width), mode='constant')
                print 'padded'
                #print a, a.shape
                #config.log.error('Asked for part of the chomosome outside of the array')
            '''

            a = self.get(l) # This is much faster than the chrom caching system...

            if respect_strand:
                # positive strand is always correct, so I leave as is.
                # For the reverse strand all I have to do is flip the array.
                if read[1] in negative_strand_labels:
                    a = a[::-1]
            if number_of_tags_in_library:
                #print(a, number_of_tags_in_library)
                a = numpy.array(a, dtype=numpy.float)
                a /= float(number_of_tags_in_library) # This is 1.0 if norm_by_read_count == False

            # bin the data
            ret = utils.bin_data(a, bin_size)
            if ret: # Python list, so testable
                table.append(numpy.array(ret))
            p.update(idx)

        # sort the data by intensity
        # no convenient numpy. So have to do myself.
        mag_tab = []
        for index, row in enumerate(table):
            mag_tab.append({"n": index, "sum": row.max()}) # Shouldn't this be sum?

        data = numpy.array(table)
        if sort_by_intensity:
            mag_tab = sorted(mag_tab, key=itemgetter("sum"))
            data = numpy.array([data[item["n"],] for item in mag_tab])
            data = numpy.delete(data, numpy.s_[-1:], 1) # Nerf the last column.

        # Get the sorted_locs for reporting purposes:
        temp_gl_copy = genelist.deepcopy().linearData # deepcopy for fastness
        sorted_locs = [temp_gl_copy[i['n']] for i in mag_tab] # need to sort them by intensity, if done that way
        sorted_locs.reverse() # Make it into the intuitive order.
        gl = Genelist()
        gl.load_list(sorted_locs)
        sorted_locs = gl

        if log:
            if not log_pad:
                log_pad = min([0.1, max([0.1, numpy.min(data)])])

            if log == "e" or log == math.e:
                data = numpy.log(data+log_pad)
            elif log == 2:
                data = numpy.log2(data+log_pad)
            elif log == 10:
                data = numpy.log10(data+log_pad)

        # draw heatmap

        if not self._draw:
            self._draw = draw()

        if filename:
            if not bracket:
                m = data.mean()
                ma = data.max()
                mi = data.min()
                std = data.std()
                bracket=[m, m+(std*2.0)]
                config.log.info("track.heatmap(): I guessed the bracket ranges as [%.3f, %.3f]" % (bracket[0], bracket[1]))
            elif len(bracket) == 1: # Assume only minimum.
                bracket = [bracket[0], data.max()]

            filename = self._draw.heatmap2(data=data, filename=filename, bracket=bracket, colour_map=cmap,
                imshow=imshow, **kargs)

        if raw_heatmap_filename:
            numpy.savetxt(raw_heatmap_filename, data, delimiter="\t")
            config.log.info("track.heatmap(): Saved raw_heatmap_filename to '%s'" % raw_heatmap_filename)

        config.log.info("track.heatmap(): Saved heatmap tag density to '%s'" % filename)
        return({"data": data, 'sorted_original_genelist': sorted_locs})


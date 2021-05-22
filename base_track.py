"""

base_track

base class for track-like objects (ie. tracks and flats)

"""

import sys, os, sqlite3, time, math, numpy, zlib
from array import array
import matplotlib.cm as cm

from . import config, utils
from .draw import draw
from operator import itemgetter

from . import genelist as Genelist # Namespace mangling :(
from .genelist import genelist
from .location import location
from .errors import FailedToMakeNewDBError
from .progress import progressbar
from .data import negative_strand_labels, positive_strand_labels

class base_track:
    def __init__(self, name=None, new=False, filename=None, norm_factor=1.0, mem_cache=False):
        """
        base track only accepts three arguments,
        the filename, name (this is a legacy thing) and new.
        If new it sets up the meta data, but that's ALL!
        """
        assert filename, "you must specify a filename"

        self.norm_factor = norm_factor

        m = name or filename
        self.meta_data = {"name": m, # Setup a dummy meta_data
            "source_filename": filename,
            "creation_date": time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()),
            "version": "4.0", # forth major revision of the track structures...
            "glbase_version": config.version,
            "bin_format": None}

        # set-up the tables
        if new:
            self.__setup_db(filename) # return an empty track
        else:
            assert os.path.exists(filename), "track '%s' cannot be found" % filename
            self.__load_tables(filename)
            self._load_meta_data()
            # Unpack any meta_data
            if name: # perform override of name metadata (if required)
                self.meta_data["name"] = name

        self._c = None
        self._draw = None # Lazy set-up in track. Don't init draw unless needed.

        self.gl_mem_cache = None
        if mem_cache:
            config.log.info('caching the database "%s" into memory...' % filename)

            from io import StringIO

            # This doens't give a very big speed up in reality:
            self._connection = sqlite3.connect(filename)
            self._connection.text_factory = sqlite3.OptimizedUnicode
            tempfile = StringIO()
            for line in self._connection.iterdump():
                tempfile.write('%s\n' % line)
            self._connection.close()
            tempfile.seek(0)

            # Create a database in memory and import from tempfile
            self._connection = sqlite3.connect(":memory:")
            self._connection.cursor().executescript(tempfile.read())
            self._connection.commit()
            self._connection.row_factory = sqlite3.Row

            self._c = self._connection.cursor()

            # Use an optional genelist-like variant instead
            # This one consumes waaay too much memory. Even a single list murders 10Gb....
            '''
            self.gl_mem_cache = genelist()

            self._connection = sqlite3.connect(filename)
            self._connection.text_factory = sqlite3.OptimizedUnicode
            tempfile = StringIO()
            for line in self._connection.iterdump():
                if 'INSERT' in line:
                    t = line.split(' ')
                    pos = t[3].replace('VALUES', '').strip('();').split(',')
                    chrom = t[2].strip('"').replace('chr_', '')
                    left = int(pos[0])
                    right = int(pos[1])
                    strand = pos[2].strip("'")
                    #print chrom, left, right, strand
                    self.gl_mem_cache.linearData.append({'loc': location(chr=chrom, left=left, right=right), 'strand': strand})

            self.gl_mem_cache._optimise()
            '''

        config.log.info("Bound '%s'" % filename)

    def __getitem__(self, key):
        """
        make meta data accesible like a dict
        """
        if key == "info": # catch this special key
            for k in self.meta_data:
                print("%s\t:\t%s" % (k, self.meta_data[k]))
        else:
            assert key in self.meta_data, "'%s' not found in this track" % key
            return(self.meta_data[key])

    def _load_meta_data(self):
        """
        retrieve the meta data from the
        """
        c = self._connection.cursor()

        # First see if meta data exists.
        # This may be an old db, which has no info table - I'll need to make one then.
        ''' # should be none of these in the wild by now, anyway, just broke the BLOCK_SIZE
        c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='info'")
        if not "info" in [k[0] for k in c.fetchall()]:
            # No meta data at all.
            c.execute("CREATE TABLE info (key TEXT PRIMARY KEY, value TEXT)")
            self._save_meta_data()
            return()
        '''

        # Make certain there are no missing keys.
        # If there are, add them to the table.
        # If the key is there load it into self.meta_data
        c.execute("SELECT key FROM info")
        db_known_keys = [k[0] for k in c.fetchall()]

        # This feels unwise to modify on a load
        for key in self.meta_data:
            if key not in db_known_keys: # This is a new metadata_attribute, save it
                c.execute("INSERT INTO info VALUES (?, ?)", (key, self.meta_data[key]))

        # Okay, get all the other keys as normal
        c.execute("SELECT * FROM info")
        for item in c.fetchall():
            if item[0] == 'pre_build': # Special unpacks
                self.meta_data[item[0]] = list(item[1])
            else:
                self.meta_data[item[0]] = item[1]

        #self._connection.commit()
        c.close()

    def _save_meta_data(self):
        """
        Load a dictionary of meta data into the info table
        """
        c = self._connection.cursor()

        # Work out if there are any new metadata keys to add.
        c.execute("SELECT key FROM info")

        db_known_keys = [k[0] for k in c.fetchall()]

        for key in self.meta_data:
            if key in db_known_keys:
                c.execute("UPDATE info SET value=? WHERE key=?", (self.meta_data[key], key))
            else: # This is a new metadata_attribute, save it
                c.execute("INSERT INTO info VALUES (?, ?)", (key, str(self.meta_data[key])))
        self._connection.commit()
        c.close()

    def __setup_db(self, filename):
        """
        No pre-defined file - I want a new database track.
        """
        # kill any previously exisiting file (Use with care!)
        # make sure the directory is available:
        path = os.path.split(filename)[0]
        if path and not os.path.exists(path):
            os.makedirs(path)

        try:
            if os.path.exists(filename): # overwrite old file.
                os.remove(filename)
                # This could potentially fail - I should report and fail
                # nicely... At the moment it just throws an exception.
        except Exception:
            raise FailedToMakeNewDBError(filename,)

        self.__load_tables(filename)

        c = self._connection.cursor()

        c.execute("CREATE TABLE info (key TEXT PRIMARY KEY, value TEXT)")
        for key in self.meta_data: # Load in meta data
            c.execute("INSERT INTO info VALUES (?, ?)", (key, self.meta_data[key]))

        self._connection.commit()
        c.close()

    def __repr__(self):
        return("glbase.base_track")

    def __load_tables(self, filename):
        """
        just load in the tables.
        (basically, fill __connection)
        """
        self._connection = sqlite3.connect(filename)
        self._connection.text_factory = sqlite3.OptimizedUnicode

    def _format_data(self, data):
        """
        array('i', []) --> whatever it's stored as in db
        """
        return(sqlite3.Binary(zlib.compress(data.tostring())))

    def _unformat_data(self, data):
        """
        whatever stored as in db --> array('i', [])
        """
        #print "ret:",[d for d in data], ":"
        a = array(self.bin_format)
        a.fromstring(zlib.decompress(data))
        return(a)

    def finalise(self):
        """
        finalise the database (shrink unused edit space)
        dump useless bits etc.
        You must call this! to finalise the db.
        Or get() will not work!
        This copies the cache onto disk and closes the db.
        """
        # do a commit
        self._save_meta_data()
        self._connection.commit()

    def old_heatmap(self, filename=None, genelist=None, distance=1000, read_extend=200, log=2,
        bins=20, sort_by_intensity=True, **kargs):
        """
        TO DEPRECATE!


        **Purpose**
            Draw a heatmap of the seq tag density drawn from a genelist with a "loc" key.

        **Arguments**
            filename (Required)
                filename to save the heatmap into

            genelist (Required)
                a genelist with a 'loc' key.

            distance (Optional, default=1000)
                Number of base pairs around the location to extend the search

            bins (Optional, default=20)
                Number of bins to use. (i.e. the number of columns)

            read_extend (Optional, default=200)
                number of base pairs to extend the sequence tag by.

            log (Optional, default=2)
                log transform the data, optional, the default is to transform by log base 2.
                Note that this parameter only supports "e", 2, and 10 for bases for log

            sort_by_intensity (Optional, default=True)
                sort the heatmap so that the most intense is at the top and the least at
                the bottom of the heatmap.

        **Results**
            file in filename and the heatmap table
        """
        assert filename, "must specify a filename"
        assert genelist, "must provide a genelist"
        assert "loc" in list(genelist.keys()), "appears genelist has no 'loc' key"
        assert "left" in list(genelist.linearData[0]["loc"].keys()), "appears the loc key data is malformed"
        assert log in (False, "e", math.e, 2, 10), "this 'log' base not supported"

        table = []
        bin_size = int((distance*2) / bins)

        for item in genelist.linearData:
            l = item["loc"].pointify().expand(distance)

            row = self.get(l, read_extend=read_extend)

            # bin the data
            row = numpy.array(utils.bin_data(row, bin_size))

            table.append(row)

        # sort the data by intensity
        # no convenient numpy. So have to do myself.
        mag_tab = [{"n": index, "sum": row.max()} for index, row in enumerate(table)]
        if sort_by_intensity:
            mag_tab = sorted(mag_tab, key=itemgetter("sum"))

        data = numpy.array(table)+1

        newt = [data[item["n"],] for item in mag_tab]
        data = numpy.array(newt)

        if log:
            if log in ["e", math.e]:
                data = numpy.log(data)-1
            elif log == 2:
                data = numpy.log2(data)-1
            elif log == 10:
                data = numpy.log10(data)-1

        # draw heatmap

        if not self._draw:
            self._draw = draw()

        filename = self._draw.heatmap2(data=data, filename=filename, bracket=[data.min(), data.max()], **kargs)

        config.log.info("Saved pileup tag density to '%s'" % filename)
        return({"data": data})


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

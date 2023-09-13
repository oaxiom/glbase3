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
            assert os.path.isfile(filename), f"track '{filename}' cannot be found"
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
        if path and not os.path.isdir(path):
            os.makedirs(path)

        try:
            if os.path.isfile(filename): # overwrite old file.
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

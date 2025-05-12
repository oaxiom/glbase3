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

from .genelist import genelist
from .location import location
from .errors import FailedToMakeNewDBError
from .progress import progressbar
from .data import negative_strand_labels, positive_strand_labels

class base_sql:
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
            return self.meta_data[key]

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
        return "glbase.base_track"

    def __load_tables(self, filename):
        """
        just load in the tables.
        (basically, fill __connection)
        """
        self._connection = sqlite3.connect(filename)

    def _format_data(self, data):
        """
        array('i', []) --> whatever it's stored as in db
        """
        return sqlite3.Binary(zlib.compress(data.tostring()))

    def _unformat_data(self, data):
        """
        whatever stored as in db --> array('i', [])
        """
        #print "ret:",[d for d in data], ":"
        a = array(self.bin_format)
        a.fromstring(zlib.decompress(data))
        return a

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

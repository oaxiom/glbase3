"""
A special ordered list of genes.
behaves like a normal list, but each element contains a heterogenous set of data.

"""

import sys, os, csv, copy, random, pickle, re, numpy, scipy, gzip

from operator import itemgetter

from . import config
from . import utils
from .flags import *
from .helpers import *
from .location import location
from .draw import draw
from .history import historyContainer
from .errors import AssertionError, UnRecognisedCSVFormatError, UnrecognisedFileFormatError, ArgumentError
from .progress import progressbar
from .base_genelist import _base_genelist
from .format import sniffer, sniffer_tsv, _load_hmmer_tbl

class Genelist(_base_genelist): # gets a special uppercase for some dodgy code in map() I don't dare refactor.
    """
    **Purpose**
        This is a class container for any arrangement of heterogenous data.

        it is good for dealing with csv/tsv files with arbitrary columns - arranging
        them by keys and then allowing cross-matching and data extraction.
        genelist is the generic implementation. Other derived classes are available for
        e.g. peaklists, genome annotations and microarray (or expression) based data.

    **Arguments**
        name (Optional)
            Will default to the name of the file if a file is loaded,
            otherwise name will be set to "Generic List" by default.
            us the name argument to give it a custom nam,e

        filename (Optional)
            the name of a file to attempt to load.

        force_tsv (Optional)
            if you specify a filename to load then
            setting this argument to True will force the file to be treated
            as a tsv (tab-separated) rather than the default csv (comma
            separated).

        loadable_list (Optional, default=None)
            If you supply a list of dicts then glbase will attempt to build a genelist out of it.

        gzip (Optional, default=False)
            The input file is gzipped

        format (Required, or use format.sniffer to let glbase guess)
            Format specifiers are a mini language that can describe any valid TSV file.

            They should be a dictionary, in the form:

            {'column_name1': 1,
            'column_name2': 2}

            where the key of the dictionary will be the key name, and the value will be the
            column number.

            location data is handled a little differently, as the values may be split across
            several columns. In this case glbase is looking for a special 'location' tag with a specific format.

            Suppose the chromosome was in column 3, the left coordinate was in column 4 and right in column 5,
            to add a genome location to our format, we would add a 'loc' key, containing this info:

            {'column_name1': 1,
            'column_name2': 2,
            "loc": "location(chr=column[3], left=column[3], right=column[5])"}

            To help deal with unusal syntax in TSV or CSV files there are a list of reserved
            key names that perform some special function:


            duplicates_key
                ?

            skiplines
                do not start loading from the file until you get to line number 'skiplines': value

            debug
                print out a debug load of the file, stopping at 'debug': X line numbers

            special
                ?

            skiptill
                do not start loading the file until you see a line starting with 'skiptill' and
                start loading the file from the next line.

            force_tsv
                forse the loader to assume the file is a TSV, rather than the defaul CSV

            gtf_decorators
                This specifies the column number that contains GTF decorators, which will be split into key:value and added to the genelist

            endwith
                Stop loading the file if you see a line that contains the value specified in endwith

            __description__
                ?

            commentlines
                Ignore lines that start with this string (e.g. 'commentlines': '#' is quite commont)

            keepifxin
                ?

            __column_must_be_used
                ?

            __ignore_empty_columns
                Ignore a column if there is no value in the column, this is for when TSVs/CSVs
                are incomplete and are missing columns on specific lines, but you don't want to have
                to sanitise the TSV/CSV, and would prefer to just fill in the blank with nothing.

            As an example, here is the full format for a complete BED file:

            {"loc": "location(chr=column[0], left=column[1], right=column[2])",
            "name": 3, "score": 4, "strand": 5, "thickStart": 6, "thickEnd": 7,
            "itemRgb": 8, "blockCount": 9, "blockSizes": 10, "blockStarts": 11,
            "force_tsv": True, "skiplines": -1, "commentlines": "#"}

            see also glbase/format.py for a list of already defined format specifiers
            that you can call using:

            gl = genelist(..., format=format.full_bed)

    """
    def __init__(self, filename=None, loadable_list=None, gzip=False, **kargs):
        # This call signature is used in a few places, so modify with care
        valig_args = ["name", "format", "force_tsv",]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError(self.__init__, k)

        self.linearData = []
        self.dataByChr = None # this is private, use get by loc.
        self.debug = False
        self.draw = draw(self)
        self.name = "Generic List"
        self.metadata = {} # container for various metadata's to extract figures from.
        self.__deathline = None # Error reporting in load_CSV()
        self.__deathindx = None

        format = sniffer
        if "format" in kargs:
            format = kargs["format"] # I expect a filename = is coming.

        if "force_tsv" in kargs and kargs["force_tsv"]:
            format["force_tsv"] = True

        if filename:
            if "format" in kargs:
                self.load(filename=filename, format=format, gzip=gzip)
            else:
                raise AssertionError('Due to excessive ambiguity the sniffing function of genelists has been removed and you now MUST provide a format argument, you can reenable this feature by specifying the sniffer: format=format.sniffer')

            config.log.info("genelist(): loaded '%s' found %s items" % (filename, len(self.linearData)))
        elif loadable_list:
            self.load_list(loadable_list)

        if "name" in kargs: # Here so it overrides anything above.
            self.name = kargs["name"]

    def load(self, filename=None, format=None, gzip=False, **kargs):
        """
        **Purpose**

        load a file into the genelist. load will attempt to load the file
        based on the filename, unless a format is specified.

        **Arguments**

        filename
            absolute filename (including path) to the actual file.
            can include path short cuts (e.g. "./", "../" etc)

        format (Optional, default = "sniffer" (ie. guess))
            format specifer, see format.py, flags.py and helpers.py and the
            documentation on how to write a valid format specifier

        **Result**

        fills the genelist with the data from the file as specified by
        the format specifier.
        """
        assert filename, "No filename specified"
        assert os.path.exists(os.path.realpath(filename)), "File %s not found" % filename

        self.path = os.path.split(os.path.realpath(filename))[0]
        self.filename = os.path.split(os.path.realpath(filename))[1]
        self.fullfilename = filename
        if self.filename.find(".") != -1:
            self.name = "".join(self.filename.split(".")[:-1])
        else:
            self.name = self.filename

        if format:
            if "special" in format: # special loads
                if format["special"] == "fasta":
                    self.linearData = utils.convertFASTAtoDict(filename=filename, gzip_input=gzip)
                    # See if I can parse names into a location?
                    try:
                        for item in self.linearData:
                            item["loc"] = location(loc=item["name"])
                    except Exception:
                        pass
                    self._optimiseData()
                    return(True)
                if format["special"] == "hmmer_tbl":
                    self.linearData = _load_hmmer_tbl(filename)
                    self._optimiseData()
                    return(True)
        else:
            raise AssertionError('Due to excessive ambiguity the sniffing function of genelists has been removed and you now MUST provide a format argument')

        csv_headers = frozenset(["csv", "xls", "tsv", "txt", "bed"])
        if filename.split(".")[-1].lower() in csv_headers: # check the last one for a csv-like header
            self.loadCSV(filename=filename, format=format, gzip=gzip, **kargs)
        elif filename.split(".")[-1] in ["glb"]:
            self = glload(filename) # will this work?
        else:
            self.loadCSV(filename=filename, format=format, gzip=gzip, **kargs)

        if "force_tsv" not in kargs and "force_tsv" not in format and len(list(self.keys())) == 1:
            config.log.warning("List contains only a single key, are you sure this is not a tsv?")

        return(True) # must have made it to one - if it fails it should trigger

    def loadCSV(self, filename=None, format=None, **kargs):
        """
        **Purpose**

        load a CSV file into the genelist

        **Arguments**

        filename
            absolute filename (including path) to the actual file.
            can include path short cuts (e.g. "./", "../" etc)

        format (Optional, default = "sniffer" (ie. guess))
            format specifer, see flags.py and helpers.py and the
            documentation on how to write a valid format specifier

        force_tsv (Optional, default=False)
            if you don't send a format specifier, but want to force the
            genelist to load as a tsv, you can set this flag to True.
            NOTE: If you send a format argument then this argument is ignored!

            As another option, you can use the sniffer_tsv format specifier

        name (Optional, Default = based on the filename)
            name of the genelist, used intermittently as labels in
            figures and the like.

        **Result**

        fills the genelist with the CSV table.
        """
        assert os.path.exists(os.path.realpath(filename)), "File %s not found" % filename

        self.path = os.path.split(os.path.realpath(filename))[0]
        self.filename = os.path.split(os.path.realpath(filename))[1]
        self.fullfilename = filename
        self.name = self.filename.split(".")[0]

        # decide whether to respect the force_tsv arg.
        if not format:
            if "force_tsv" in kargs and kargs["force_tsv"]:
                format = sniffer_tsv
            else:
                format = sniffer

        if "debug" in format and format["debug"]:
            print("--------")
            print("DEBUG load:")
            self._loadCSV(filename=self.fullfilename, format=format, **kargs)
        else:
            try:
                self._loadCSV(filename=self.fullfilename, format=format, **kargs)
            except Exception:
                try: # try again, guessing it might be a tsv
                    config.log.warning("Failed to load as a 'csv'. Try to load as a 'tsv' file")
                    format["dialect"] = csv.excel_tab
                    self._loadCSV(filename=self.fullfilename, format=format, **kargs)
                except Exception: # oh dear. Die.
                    if self.__deathline:
                        config.log.error("Died on line: '%s'" % self.__deathline)
                    raise UnRecognisedCSVFormatError("'%s' appears mangled, the file does not fit the format specifier" % self.fullfilename, self.fullfilename, format)

    def _loadCSV(self, **kargs):
        """
        (Internal)

        Actual loadCSV()

        This is mainly so the try/except block above doesn't get completely out of control
        and allows debug code.
        """
        assert "filename" in kargs, "No filename specified"
        assert "format" in kargs, "_loadCSV requres a format specifier"
        assert os.path.exists(kargs["filename"]), "%s file not found" % kargs["filename"]

        filename = kargs["filename"]
        format = kargs["format"]

        temp_data = []
        if 'gzip' in kargs and kargs['gzip']:
            oh = gzip.open(filename, "rt")
        else:
            oh = open(filename, "rt")

        if "force_tsv" in kargs and kargs["force_tsv"]: # force_tsv takes priority
            reader = csv.reader(oh, dialect=csv.excel_tab)
        elif "force_tsv" in format and format["force_tsv"]:
            reader = csv.reader(oh, dialect=csv.excel_tab)
        elif "dialect" in format and format["dialect"]:
            reader = csv.reader(oh, dialect=format["dialect"])
        else:
            reader = csv.reader(oh)

        if "skiplines" in format:
            skiplines = format["skiplines"]
        else:
            skiplines = 0 # skip any header row by default.

        if "skiptill" in format and format["skiptill"]:
            skiptill = format["skiptill"]
        else:
            skiptill = "Done" # This is done as truth testing fails as format["skiptill"] != None

        if "sniffer" in format:
            # I need to construct my own format
            format = {}
            for top_line in reader:
                for index, key in enumerate(top_line): # get all the potential keys.
                    format[key] = index
                skiplines = -1 # if the sniffer is used then when it gets to the below
                # index will = 0 = skiplines causing the first line to be missed.
                break

        debug_line = 0

        # sanitise format[key] data for security.
        newf = {}
        for k in format:
            if isinstance(format[k], dict) and "code" in format[k]:
                newf[k] = {"code": format[k]["code"].replace("__", "").replace(" sys,", "")} # Some security suggestions:, {"__builtins__":None, column, location}, {})
            elif isinstance(format[k], str) and "location" in format[k]:
                newf[k] = format[k].replace("__", "").replace(" sys,", "")
            else:
                newf[k] = format[k]
        format = newf

        for index, column in enumerate(reader):
            # This is cryptically called column, when it is actually row.\
            # there is a reason for that, it is so that in the formats it appears:
            # "refseq": column[1] # see :)
            #print index, column # debug for when all else fails!
            self.__deathline = column # For debugging purposes
            self.__deathindx = index

            if not column: # if row is completely empty, so just omit.
                continue

            if index <= skiplines or skiptill != "Done":
                if column and True in [skiptill in item for item in column]:
                    skiptill = "Done"
                continue

            if "__column_must_be_used" in format and format["__column_must_be_used"]:
                if not column[format["__column_must_be_used"]]:
                    continue # Omit data when this particular column is blank

            if "endwith" in format and format["endwith"]:
                if True in [format["endwith"] in item for item in column]:
                    break

            if "debug" in format and format["debug"]:
                debug_line += 1
                print("%s:'%s'" % (index, column))
                if isinstance(format["debug"], int) and debug_line > format["debug"]:
                    break # If an integer, collect that many items.

            if column[0] not in typical_headers:
                if "commentlines" in format and format["commentlines"]:
                    if column[0][0] == format["commentlines"]:
                        continue

                # passed all the tests
                temp_data.append(self._processKey(format, column))

            #print temp_data[-1] # tells you what got loaded onto the list.
        oh.close()

        self.linearData = temp_data

        self._optimiseData()
        return(True)

    def _optimiseData(self):
        """
        (Internal)
        Call me after modifying the data to bin and build the internal tables.
        """
        self.dataByChr = None
        if not self.linearData: # list is empty, not possible to optimise anything...
            return(False)

        # Guess a loc key
        loc_key = None
        if "tss_loc" in self.linearData[0]: # always use tss_loc in preference of loc, if available
            loc_key = "tss_loc"
        elif "loc" in self.linearData[0]:
            loc_key = "loc" # Don't change this though. annotate() relies on the bucket system using tss_loc

        if "tss_loc" in self.linearData[0] and "loc" in self.linearData[0]:
            config.log.warning("List contains both 'tss_loc' and 'loc'. By default glbase will use 'tss_loc' for overlaps/collisions/annotations")

        if loc_key in self.linearData[0]: # just checking the first entry.
            self.dataByChr = {}
            self.dataByChrIndexLookBack = {}
            self.buckets = {}
            for n, item in enumerate(self.linearData): # build the chromosome quick search maps.
                chr = item[loc_key]["chr"]
                if not chr in self.dataByChr:
                    self.dataByChr[chr] = []
                    self.dataByChrIndexLookBack[chr] = []
                self.dataByChr[chr].append(item)
                self.dataByChrIndexLookBack[chr].append(n) # oh sweet, sweet dirty hack...
                # I can't remember what this look-back is for, but you
                # can use it to get the linearData index even though looking at the
                # dataByChr data It is not documented for a reason!
                # New bucket system to go in here.

                if not chr in self.buckets:
                    self.buckets[chr] = {}
                # work out the bucket(s) for the location.
                # which bucket is left and right in?
                left_buck = (item[loc_key]["left"]//config.bucket_size)*config.bucket_size
                right_buck = ((item[loc_key]["right"]+config.bucket_size)//config.bucket_size)*config.bucket_size
                buckets_reqd = list(range(left_buck, right_buck, config.bucket_size))

                #print n, item[loc_key], buckets_reqd, left_buck, right_buck, len(buckets_reqd)

                for b in buckets_reqd:
                    if not b in self.buckets[chr]:
                        self.buckets[chr][b] = []
                    self.buckets[chr][b].append(n) # use index to maintain uniqueness.

        # Then to do a collision =
        """
        left_buck = (item[loc_key]["left"]//config.bucket_size)*config.bucket_size
        right_buck = (item[loc_key]["right"]//config.bucket_size)*config.bucket_size
        buckets_reqd = range(left_buck, right_buck, config.bucket_size)

        for buck in buckets_reqd:
            cbuck = set(self.buckets[chr][buck]) # unique ids
        """

        """
        # Build a quickfinder:

        This is dict representation of the list, and allows you to do things like:

        indeces = self.qkeyfind[key][value]
        """
        self.qkeyfind = {}
        for index, item in enumerate(self.linearData):
            for key in item:
                if not key in self.qkeyfind:
                    self.qkeyfind[key] = {}

                try:
                    if not item[key] in self.qkeyfind[key]:
                        self.qkeyfind[key][item[key]] = []
                    self.qkeyfind[key][item[key]].append(index)
                except TypeError:
                    # The item in unhashable and cannot be added to the qkeyfind
                    # This should be pretty rare if not impossible.
                    #print '!Unhashable key: %s for qkeyfind system' % key
                    pass

            # Now to do a find you just go:
            # item_indeces = self.qkeyfind["name"]["Stat3"]

        return(True)

    def isChromosomeAvailable(self, chromosome):
        """
        you must check me before trying to access dataByChr[]
        """
        if chromosome in self.dataByChr:
            return(True)
        else:
            return(False)
        return(False)

    def getAllUnorderedData(self):
        """
        (Obselete?)
        return the list of unordered data
        """
        return(self.linearData)

    def _findDataByKeyLazy(self, key, value): # override????? surely find?
        """
        (Internal)

        find the first matching entry of key value pair

        This version is lazy, so I take the min() and return that item
        """
        if key in self.qkeyfind:
            if value in self.qkeyfind[key]:
                return self.linearData[min(self.qkeyfind[key][value])]
        return None # not found;

    def _findDataByKeyGreedy(self, key, value): # override????? surely finditer?
        """
        finds all - returns a list
        """
        ret = []
        item_indeces = None
        if key in self.qkeyfind:
            if value in self.qkeyfind[key]:
                item_indeces = self.qkeyfind[key][value]

        if item_indeces:
            return([self.linearData[i] for i in item_indeces])
        return(None)

    def get(self, key, value, mode="greedy"):
        """
        **Purpose**
            get all items that match "value" in "key"
            you can also use find() to get the first matching item, or set mode to "lazy" to
            achieve the same effect

            get will always return a genelist, even if there is only 1 item

        **Arguments**
            key
                the key of the genelist to search in

            value
                the value to look for in key

            mode (Optional, default="greedy")
                the mode of search. "greedy" searches will
                find all examples, "lazy" searches will only find
                the first in the list. Most of the time you mean
                "greedy"

        **Returns**
            A new genelist or None
        """
        assert key in self.keys(), '"{0}" key not found in this list'.format(key)

        if mode == "greedy":
            r = self._findDataByKeyGreedy(key, value)
        elif mode == "lazy":
            r = self._findDataByKeyLazy(key, value) # lazy just returns a single dict
            if r:
                r = [r]
        else:
            raise AssertionError("mode '%s' for get() is not known" % mode)

        # The internal methods return vanilla lists for compatability with er... internal stuffs
        # repackage to a new gl.
        if r:
            newl = self.shallowcopy() # shallowcopy for once as we will load in our own list.
            newl.load_list(r)
            return(newl)
        return(None)

    def index(self, key, value):
        """
        **Purpose**
            A bit like the Python index() but for key:value pairs

            NOTE: The method is lazy and finds the first index and returns that.

        **Arguments**
            Key (Required)
                the key to search in

            Value (Required)
                the value to find

        **Returns**
            The index of the list the item is contsained at.
        """
        assert key, "must send key"
        assert value, "must send value"

        return(min(self.qkeyfind[key][value]))

    def _findAllLabelsByKey(self, key):
        """
        (Internal)
        Returns a 1D list of all the values under Key.
        Most useful for things like geneList["Symbol"]
        geneList["entrez"]
        """
        return([x[key] for x in self.linearData])

    def _findByLabel(self, key, value):
        """
        (Internal)
        key is the key to search with
        toFind is some sort value to compare

        There is a subtle problem here if the user tries to find a list or other non-hashable
        object. But then It would be kind of wierd to do that...

        (This is used in at least draw._heatmap_and_plot())

        This version is deprectaed. The official internal methods are:
        _findDataByKeyLazy|Greedy()
        """
        return(self._findDataByKeyLazy(key, value)) # not found;

    def _findByLoc(self, key, loc):
        """
        (internal)
        coords should be in formal format, chrX:int(left)-int(right)
        key is unused

        lazy.
        """
        ret = []

        if self.isChromosomeAvailable(loc["chr"]):
            for item in self.dataByChr[loc["chr"]]:
                if utils.qcollide(loc["left"], loc["right"], item[key]["left"], item[key]["right"]):
                    ret.append(item)
        return(ret)

    def saveTSV(self, filename=None, **kargs):
        """
        **Purpose**
            save the geneList or similar object as a tsv
            Note: This is not always available.
            As the geneList becomes more complex it loses the ability to be
            expressed simply as a csv-file. In those cases you must use
            the save() method to save a binary representation.

            saveTSV is *generally* 'consistent' if you can succesfully save
            as a tsv then you can reload the same list as that particular
            type. Be careful though. Behaviour like this will work fine::

                microarray.saveTSV(filename="file.csv")
                anewlist = genelist(filename="file.csv")
                anewlist # this list is now a genelist and not an expression.
                         # You must load as an expression:
                rnaseq = expression(filename="file.csv")

        **Arguments**
            filename
                filename to save, including the full path.

            key_order (Optional, default=None)
                send a list, specifying the order you'd like to write the keys
                by default saveTSV() will write to the file in an essentially
                random order. But if you want to specify the order then
                send a list of key names and it will save them in that order.

                Also, you need only specify the left side of the column.
                Any unspecified columns will be appended to the right
                hand side in a random order.

                Note that saveTSV() always saves all columns of data.

            tsv (True|False)
                save as a tsv file see also saveCSV()

            no_header (Optional, default=False)
                Don't write the first line header for this file. Usually it's the list
                of keys, then the second line will be the data. If no_header is set to False
                then the first line of the file will begin with the data.

        **Result**
            returns None.
            saves a TSV representation of the genelist.
        """
        self.saveCSV(filename, tsv=True, **kargs)

    def saveCSV(self, filename=None, no_header=False, **kargs):
        """
        **Purpose**

            save the geneList or similar object as a csv
            Note: This is not always available.
            As the geneList becomes more complex it loses the ability to be
            expressed simply as a csv-file. In those cases you must use
            the save() method to save a binary representation.

            saveCSV is guaranted to be 'consistent' if you can succesfully save
            as a csv then you can reload the same list as that particular
            type. Be careful though. Behaviour like this will work fine::

                microarray.saveCSV(filename="file.csv")
                anewlist = genelist(filename="file.csv")
                anewlist # this list is now a genelist and not a microarray.
                         # You must load as a microarry:
                amicroarry = microarray(filename="file.csv")

            saving to a csv will will blank the history, and any other meta-
            data generated about the list.

        **Arguments**

            filename
                filename to save, including the full path.

            key_order (List)
                send a list, specifying the order you'd like to write the keys
                by default saveTSV() will write to the file in an essentially
                random order. But if you want to specify the order then
                send a list of key names and it will save them in that order.

                Also, you need only specify the left side of the column.
                Any unspecified columns will be appended to the right
                hand side in a random order.

                Note that saveTSV() always saves all columns of data.

            tsv (True|False)
                save as a tsv file see also saveTSV()

            no_header (Optional, default=False)
                Don't write the first line header for this file. Usually it's the list
                of keys, then the second line will be the data. If no_header is set to False
                then the first line of the file will begin with the data.

        **Result**
            returns None.
            saves a CSV representation of the geneList.
        """
        assert filename, "No filename specified"

        valig_args = ["filename", "key_order", "tsv", "no_header"]
        for k in kargs:
            if k not in valig_args:
                raise ArgumentError(self.saveCSV, k)

        oh = open(filename, "w")
        if not self.linearData: # data is empty, fail graciously.
            config.log.error("csv file '%s' is empty, no file written" % filename)
            oh.close()
            return(None)

        if "tsv" in kargs and kargs["tsv"]:
            writer = csv.writer(oh, dialect=csv.excel_tab)
        else:
            writer = csv.writer(oh)

        # work out key order and the like:
        write_keys = []
        if "key_order" in kargs:
            write_keys = kargs["key_order"]
            # now add in any missing keys to the right side of the list:
            for item in list(self.keys()):
                if item not in write_keys:
                    write_keys.append(item)
        else:
            # just selece them all:
            write_keys = list(self.keys())

        if not no_header:
            writer.writerow(write_keys) # write the header row.

        for data in self.linearData:
            line = []
            for key in write_keys: # this preserves the order of the dict.
                if key in data:
                    line.append(data[key])
                else:
                    line.append("") # a blank key, fail gracefully.
            writer.writerow(line)
        oh.close()
        config.log.info("Saved '%s'" % filename)
        return(None)

    def saveFASTA(self, filename=None, seq_key="seq", **kargs):
        """
        **Purpose**

            Save data as a FASTA file. You can only do this if the
            data has a "seq" tag. The 'name' of the fasta is derived from
            the name="key" value you specify.

        **Arguments**

            filename
                absolute path to the file (including path)

            name (Optional, tries to use "seq_loc" key otherwise defaults to "null_n")
                a key in the list to use as a name for each fasta entry.
                defaults to "null_n", where n is the position in the list

            seq_key (Optional, default="seq")
                The key name that contains the sequence data. defaults to "seq" if not set (this
                is the usual place glbase stores sequence data). Note that there is no sanity checking on
                the key you use, and glbase will save whatever is contained in seq_key.
                This way you can make some really strange looking FASTA files if you like.

        **Result**

            returns True if complete.
            Saves a fasta file of the sequence data to filename.

        """
        valid_args = ["filename", "name"]
        for key in kargs:
            assert key in valid_args, "saveFASTA() - Argument '%s' is not recognised" % key

        assert filename, "No filename specified"
        assert self.linearData[0]["seq"], "No sequence data available in this list"

        oh = open(filename, "w")

        if "name" in kargs:
            name = kargs["name"]
            if not isinstance(name, list):
                name = [name]
        else:
            name = "null_"
            if "seq_loc" in self: # Default to seq_loc if available
                name = ["seq_loc"]

        for index, item in enumerate(self.linearData):
            if name == "null_":
                save_name = "null_%s" % index
            else:
                save_name = "_".join([str(item[n]) for n in name])

            oh.write(">%s\n" % save_name)
            oh.write("%s\n" % item[seq_key])
        oh.close()
        config.log.info("Saved FASTA file: %s" % filename)
        return(True)

    def saveBED(self, filename=None, extra_keys=None, id=None, score=None, uniqueID=False, loc_only=False,
        **kargs):
        """
        **Purpose**
            save the genelist in bed format
            list must have a valid loc key

            BED format files are quite diverse in their organisation.

            The official definition is here:
            http://genome.ucsc.edu/FAQ/FAQformat

            glbase saveBED will save location in columns 0, 1 and 2.

            glbase will also enforce values in columns 3, 4 and 5. A lot of downstream
            programs require at least something in columns 3, 4 and 5. glbase will make spoof
            values for these columns. Using "+" for all strand columns and 0 in columns 4 and 5.

            You can modify the data in column 4 by sending a key name with the id= argument.
            Similarly you can modify column 5 (score) by sending a key-name to score.

            Finally if your data has a strand key that will be used for column 6.

            extra_keys will then be placed in separate columns after column 6.

        **Arguments**
            filename
                Name of the file to save the bed to

            id (Optional)
                A key in the genelist to use for the id column (column 4)

            uniqueID (Optional, default=False)
                Sometime the id column must contain a unique value. If you set this to True, glbase
                will spoof a unique id but adding <name>-n to the id column, where <name> is the
                name of the genelist and n is the n'th entry in the genelist.

                if id is set to a key value then instead of the name of the genelist the key will
                be used instead and -n appended to whatever is in that key.

            score (Optional)
                A key in the genelist to use for the score column (column 5)

            extra_keys = []
                A list of extra keys to save into the bed. These will appear after column 6 (strand).

            loc_only (Optional, default=False)
                If set to True then output a BED file containing only the loc (i.e. the frist three columns:
                'chrom, left, right')

                Note that anything in extra_keys will still be respected.

            gzip (Optional, default=False)
                gzip the output file

        **Returns**
            A saved bed file and None
        """
        if 'gzip' in kargs and kargs['gzip']:
            oh = gzip.open(filename, "wt")
        else:
            oh = open(filename, "w")

        for index, item in enumerate(self.linearData):
            todo = ["chr%s" % str(item["loc"]["chr"]), str(item["loc"]["left"]), str(item["loc"]["right"])]
            if not loc_only:
                if uniqueID:
                    if id:
                        todo += ["%s-%s" % (str(item[id]), index)]
                    else:
                        todo += ["%s-%s" % (self.name, index)]
                else:
                    if id:
                        todo += [str(item[id])]
                    else:
                        todo += ["0"]

                if score:
                    todo += [str(item[score])]
                else:
                    todo += ["0"]

                if "strand" in item:
                    todo += [str(item["strand"])]
                else:
                    todo += ["+"]

            if extra_keys:
                todo += [str(item[k]) for k in extra_keys]

            oh.write("%s\n" % ("\t".join(todo)))

        oh.close()
        config.log.info("Saved '%s' BED file" % filename)
        return(filename)

    def saveGTF(self, filename=None, strand=None, source=None, feature=None, score=None, frame=None,
        loc=None, **kargs):
        """
        **Purpose**
            save the data as a gtf.

            NOTE: descriptive options are NOT supported at hte start of the gtf.

        **Arguments**
            filename (Required)
                the filename to save to.

            The next arguments correspond to the columns of the gtf file:
            1: loc["chrom"]
            2: source
            3: feature
            4: loc["left"]
            5: loc["right"]
            6: score
            7: strand
            8: frame (+/-/.)
            9: group (decorators) - all the other keys in the genelist.

            loc (Required)
                the key name to extract the location data from

            strand (Optional, default="+")
                if specified this key will be used to extract strand information

            source (Optional, default="glbase")
                if specified, take the values from this key.
                Usually, this key is used for the programme that generated this feature.

            feature (Optional, default="exon")
                if specified, take the feature description from this key

            score (Optional, default=".")
                a key name to extract 'score' data from. Must be an integer between 0..1000
                or "."

            NOTES:
                All other keys will be added to the group field, in the form:
                    key_name "value"; other_key "other_value";

        **Return**
            None and a succesfully saved file in filename

        """
        assert filename, "you must specify a filename"
        assert loc, "you must specify a loc key"
        assert loc in self, "'%s' not found in this list" % loc

        # work out the decorator keys:
        keys = list(self.keys())
        keys.remove(loc)
        if strand:
            keys.remove(strand)
        if source:
            keys.remove(source)
        if feature:
            keys.remove(feature)
        if score:
            keys.remove(score)

        oh = open(filename, "w")
        for item in self:
            # structured part of line:
            out = [None, "glbase", "exon", None, None, ".", "+", ".", []]
            out[0] = "chr%s" % item[loc]["chr"]
            out[3] = str(item[loc]["left"])
            out[4] = str(item[loc]["right"])
            if strand:
                out[6] = str(item[strand])
            if source:
                out[1] = str(item[source])
            if feature:
                out[2] = str(item[feature])
            if score:
                out[5] = str(item[score])

            # gtf formats can tolerate missing decorators:
            for k in keys:
                if k in item:
                    out[8].append('%s "%s"' % (k, item[k]))

            # Join all the bits together
            out[8] = "%s;" % ("; ".join(out[8]))
            out = "\t".join(out)

            #print out
            oh.write("%s\n" % out)
            #break
        oh.close()
        config.log.info("Saved '%s' GTF file" % filename)
        return(None)

    def sort(self, key=None, reverse=False):
        """
        Sort the data into a particular order based on key.
        This sorts the list in-place in the same style as Python.
        ie.

        ret = list.sort() - is True and not the sorted list

        list.sort() - now list contains the sorted list

        **Arguments**

        key
            the key to use to sort by.
            must be some sort of sortable item

        reverse (Optional, default=False)
            By default the list is sorted smallest to largest.
            reverse = True sorts largest to smallest.

        **Result**

        returns True if it completes.
        sorts the list IN PLACE.
        """
        assert key, "No such key '%s'" % key
        assert key in self.linearData[0], "Data does not have key '%s'" % key

        self.linearData = sorted(self.linearData, key=itemgetter(key))
        if reverse:
            self.linearData.reverse()
        self._optimiseData()
        return(True)

    def shuffle(self, key=None, reverse=False):
        """
        Randomly shuffle the order of the genelist

        **Arguments**
            None

        **Result**
        Returns None, and shuffles the list IN PLACE
        """
        new_order = list(range(0, len(self)))
        random.shuffle(new_order)

        newl = []
        for p in new_order:
            newl.append(self.linearData[p])
        self.linearData = newl

        self._optimiseData()
        return(True)

    def multi_sort(self, keys):
        """
        **Purpose**
            Sort a genelist using multiple keys.

        **Arguments**
            keys (Required)
                A list of key names to sort. Sorting will first sort keys[0] then key[1] through key[n]

        **Returns**
            returns True if it completes.
            sorts the list IN PLACE.
        """
        #assert key, "No such key '%s'" % key
        #assert key in self.linearData[0], "Data does not have key '%s'" % key

        comparers = [ ((itemgetter(k[1:].strip()), -1) if k.startswith('-') else (itemgetter(k.strip()), 1)) for k in keys]
        def comparer(left, right):
            for fn, mult in comparers:
                result = cmp(fn(left), fn(right))
                if result:
                    return mult * result
            else:
                return 0
        self.linearData = sorted(self.linearData, cmp=comparer)
        self._optimiseData()
        return(True)

    def reverse(self):
        """
        reverse the order of the list, in place.

        **Arguments**

        None

        **Result**

        returns True if okay or false.
        """
        self.linearData.reverse()
        self._optimiseData() # just in case.
        return(True)

    def getValuesInRange(self, key=None, low=None, high=None):
        """
        **Purpose**
            make a new list with all values in key that are between low and high

        **Arguments**
            key (Required)
                The key to search

            low, high (Required)
                low and high values to select items.

                Items are selected based on this logic:

                low >= x <= high

        **Returns**
            A new genelist containing only the items that pass the criteria.
        """
        assert key, "must specify a key from which to extract values"
        assert low, "'low' number not valid"
        assert high, "'high' number not valid"

        newl = self.shallowcopy()
        newl.linearData = []
        for item in self.linearData:
            if item[key] >= low and item[key] <= high:
                newl.linearData.append(item)
        newl._optimiseData()
        config.log.info("getValuesInRange: '%s' items passed (%s >= x <= %s) criteria" % (len(newl), low, high))
        return(newl)

    #------------------------------ Overrides --------------------------

    def __contains__(self, item):
        """
        (Override)
        There may be some obscure failures with this item - to do with
        returned lists. IF you send a [ {} ... {} ] like object
        derived from a genelist then it will fail (but then it should?)
        but if you use slices it should be okay:
        a = genelist[0] # returns a single dict
        a = genelist[0:10] # returns a new genelist
        """
        if not self.linearData:
            return(False)

        if item in self.linearData[0]:
            return(True)
        return(False)

    def __repr__(self):
        """
        (Override)
        report the underlying representation
        """
        return("glbase.genelist")

    def __str__(self):
        """
        (Override)
        give a sensible print out.
        """
        if len(self.linearData) > config.NUM_ITEMS_TO_PRINT:
            out = []
            # welcome to perl
            for index in range(config.NUM_ITEMS_TO_PRINT):
                out.append("%s: %s" % (index, ", ".join(["%s: %s" % (key, self.linearData[index][key]) for key in self.linearData[index]])))
            out = "%s\n... truncated, showing %s/%s" % ("\n".join(out), config.NUM_ITEMS_TO_PRINT, len(self.linearData))

            if config.PRINT_LAST_ITEM:
                out = "%s\n%s" % (out, "%s: %s" % (len(self.linearData), ", ".join(["%s: %s" % (key, self.linearData[-1][key]) for key in self.linearData[-1]])))

        elif len(self.linearData) == 0:
            out = "This list is Empty"

        else: # just print first entry.
            out = []
            for index in range(len(self.linearData)):
                out.append("%s: %s" % (index, ", ".join(["%s: %s" % (key, self.linearData[index][key]) for key in self.linearData[index]])))
            out = "%s\nShowing %s/%s" % ("\n".join(out), len(self.linearData), len(self.linearData))

        return(out)

    def all(self):
        """
        So you can do things like print gl.all()
        """
        out = []
        # welcome to perl
        for index in range(len(self.linearData)): # Yeah, funky way to do it, but maintains compatability with __str__()
            out.append("%s: %s" % (index, ", ".join(["%s: %s" % (key, self.linearData[index][key]) for key in self.linearData[index]])))

        return("\n".join(out))

    def _collectIdenticalKeys(self, gene_list):
        """
        (Internal)
        What it says, returns a list of valid keys in common between this list and gene_list
        """
        return(list(set(self.keys()) & set(gene_list.keys())))

    def getColumns(self, return_keys=None):
        """
        **Purpose**
            return a new genelist only containing the columns specified in return _keys (a list)
        """
        assert isinstance(return_keys, list), "return_keys must have a list"

        newl = self.shallowcopy()
        newl.linearData = []

        for item in self.linearData:
            newl.linearData.append(dict((k, item[k]) for k in return_keys)) # hack for lack of dict comprehension
            # assumes all keys are in the dict

        newl._optimiseData()
        config.log.info("getColumns: got only the columns: %s" % ", ".join(return_keys))
        return(newl)

    def getGeneSet(self, key=None, list_of_items=None, use_re=True, **kargs):
        """
        Deprecated, see getRowsByKey
        """
        return(self.getRowsByKey(key=key, list_of_values=list_of_items, use_re=use_re, **kargs))

    def getRowsByKey(self, key=None, values=None, use_re=True, case_sensitive=True, **kargs):
        """
        **Purpose**
            extract all rows from a genelist for which the values in key are in the
            list_of_items

            You can send regular expressions and they will be
            interpreted correctly.

            NOTE that getRowsByKey() is a SLOW look up.

            If you need speed use get(), which is basically a free lookup of the list
            (even for greedy searches) but does not support regular expressions

        **Arguments**
            values (Required)
                a list of items to collect

            key (Optional, default=None)
                the key to search in.
                If None, all keys are searched.

            case_sensitive (Optional, default=True)
                Set to True to make the search case sensitive. Only works if use_re=True

            use_re (Optional, default=True)
                Unset this if you want exact matches or are getting unexpected behaviour
                from the regular expressions.

        **Returns**
            A new genelist containing only those items.
        """
        assert values, "getRowsByKey: 'values' argument cannot be None"
        if not case_sensitive:
            assert use_re, 'use_re must be True if case_sensitive is False'

        if not isinstance(values, list):
            values = [values]

        # This should be made super fast with qkeyfind

        newl = self.shallowcopy()
        newl.linearData = []

        if use_re: # re-ise the items.
            if case_sensitive:
                list_of_items = [re.compile(i) for i in values]
            else:
                list_of_items = [re.compile(i, re.IGNORECASE) for i in values]

            if not key: # split here for clarity.
                for item in self.linearData:
                    for k in item: # no key specified, check all.
                        for r in list_of_items:
                            if r.search(str(item[k])):
                                newl.linearData.append(utils.qdeepcopy(item))
            else:
                for item in self.linearData:
                    for r in list_of_items:
                        if r.search(str(item[key])): # sometimes gene names accidentally get stored as numbers
                            newl.linearData.append(utils.qdeepcopy(item))
        else: # no re.
            if not key: # split here for clarity.
                for item in self.linearData:
                    for k in item: # no key specified, check all.
                        for r in values:
                            if r == item[k]:
                                newl.linearData.append(item.copy())
            else:
                for item in self.linearData:
                    for r in values:
                        if r == item[key]:
                            newl.linearData.append(item.copy())

        if newl:
            newl._optimiseData()
        else:
            config.log.info("getRowsByKey: Found 0 items")
            return(None)

        config.log.info("getRowsByKey: Found %s items" % len(newl))
        return(newl)

    def filter_by_value(self, key=None, evaluator=None, value=None, **kargs):
        """
        **Purpose**
            Filter data based on a key with some numeric data.

            If you skip the keyword arguments you can write nice things like this:

            newdata = expn.filter_by_value("q-value", "<", 0.05)

        **Arguments**
            key (Required)
                The name of the key to use to filter the data on.
                or a name of a condition in the expression object to filter on.

            evaluator (Required, values=["gt", "lt", "gte", "lte", "equal"])
                The comparator.
                    gt = '>' greater than value
                    lt = '<' less than value
                    gte = '>=' greater than or equal to value
                    lte = '<=' less than or equal to value
                    equal = "==" equal to value

                    You can also send the > < >= <= or == as a string symbol as well.

            value (Required)
                The value of change required to pass the test.

        **Returns**
            A new expression-object containing only the items that pass.
        """
        assert key, "filter_by_value(): 'key' argument is required"
        assert evaluator in ("gt", "lt", "gte", "lte", "equal", ">", "<", ">=", "<=", "=="), "filter_by_value(): evaluator argument '%s' not recognised" % evaluator

        if self.__repr__() == "glbase.expression" and key in self._conditions:
            assert key in self._conditions, "filter_by_value():'%s' not found in this expression object" % key
            its_a_condition = True
        else:
            assert key in list(self.keys()), "filter_by_value(): no key named '%s' found in this genelist object" % key
            its_a_condition = False

        new_expn = []

        conv_dict = {"gt": ">", "lt": "<", "gte": ">=", "lte": " <=", "equal": "=="}
        if evaluator in ("gt", "lt", "gte", "lte", "equal"):
            evaluator = conv_dict[evaluator]

        if its_a_condition:
            cond_index = self._conditions.index(key)
            for item in self.linearData:
                if eval("%s %s %s" % (item["conditions"][cond_index], evaluator, value)):
                    new_expn.append(item)
        else: # filter on a normal key.
            for item in self.linearData:
                if eval("%s %s %s" % (item[key], evaluator, value)):
                    new_expn.append(item)

        ret = self.shallowcopy()
        ret.load_list(new_expn) # In case I make optimisations to load_list()

        config.log.info("filter_by_value: Filtered expression for ['%s' %s %s], found: %s" % (key, evaluator, value, len(ret)))
        return(ret)

    def map(self, genelist=None, peaklist=None, microarray=None, genome=None, key=None,
        greedy=True, logic="and", silent=False, **kargs):
        """
        **Purpose**
            map() merges two genelist-like objects and outputs a new genelist.

            It matches, by the key, each item that overlap in the two genelist and
            returns a new genelist only containing the matching items between the two lists.

            The new genelist will inherit from 'the right', for
            example if you have a expression-object you should perform the map in this
            order::

                result = gene_list.map(genelist=expn, key="refseq") # expn is an expresion object

            'result' will now be a expression object with all the appropriate
            methods.

            If however, you did this by mistake::

                result = expn.map(genelist=gene_list, key="refseq") # expn is an expression object

            It will still work fine, but now, trying::

                result.heatmap(filename="foo.png")

            will fail, because the result is a vanilla genelist rather than
            an expression-object as you might intend.

            Also note, this method is 'greedy' by default and and will take all matching
            entries it finds. This can be changed by setting greedy=False.

        **Arguments**
            genelist
                some sort of genelist-like object,
                examples inglude genelist, expression, genome, etc

            key
                a key in common between the two lists you can use to map
                them against each other.

            image_filename (Optional)
                save a venn diagram

            venn_proportional (Optional)
                enable experimental proportional venn diagrams.

                Note that for a proper venn, both lists should be unique for
                the key you are using to match. glbase does not check that this is
                the case. This can be useful to estimate the Venn overlap so glbase remains
                silent for non-unique lists, but can occasionally give bizarre results,
                such as negative numbers in a particular condition.

            title (Optional, default = None)
                title for the venn diagram.

            greedy (Optional, default=True)
                set to True to collect all matching entries, including duplicates. (The default
                behaviour)
                If set to False then the search finds the first matching entry only

            logic (Optional, default="and")
                a logic operator to apply to the map
                Accepted operators are:

                "and" = only keep the item if it appears in both lists
                "notright" = for each item in the right hand list, only keep it if it is NOT in the left hand list

                Be aware of the strange syntax of 'notright'. This tests the item in the right list
                and only keeps it if it is NOT in the left list.

        **Result**

            returns a new genelist-like object containing the overlapping
            objects, inheriting methods from the
            right hand side of the function.

        """
        valid_args = ("genelist", "peaklist", "microarray", "key",
            "image_filename", "title", "venn_proportional", "greedy")
        for k in kargs:
            if not k in valid_args:
                raise ArgumentError(self.map, k)

        assert logic in ("and", "notright"), "logic '%s' not supported" % logic

        if repr(genelist) == "glbase.delayedlist": # delayedlists will fail an assertion
            gene_list = genelist
        else:
            assert genome or genelist or peaklist or microarray, "map(): No valid genelist specified"

        if genelist:
            gene_list = genelist
        if peaklist:
            gene_list = peaklist
        if microarray:
            gene_list = microarray
        if genome:
            gene_list = genome

        __warning_assymetric_errs = False

        assert key, "Must specify a 'key' to map the two lists"
        #assert key in gene_list.linearData[0], "'%s' key not found in provided genelist '%s'" % (key, self.name)
        #assert key in self.linearData[0], "'%s' key not found in self '%s'" % (key, self.name)
        map_key = key

        p = progressbar(len(gene_list)) # leave as len() for compatability with delayedlists
        # speed up with a triangular search?
        newl = gene_list.shallowcopy()
        if repr(genelist) == "glbase.delayedlist": # Special exception for delayedlists, need to convert to vanilla genelist:
            newl = Genelist()
            newl.name = gene_list.name

        newl.linearData = []
        for index, item in enumerate(gene_list):
            if greedy:
                results = self._findDataByKeyGreedy(map_key, item[map_key])
            else:
                results = self._findDataByKeyLazy(map_key, item[map_key])
                if results:
                    results = [results] # coerce to a single member list to simplify code below

            if results:
                if logic == "and":
                    for r in results:
                        new_entry = utils.qdeepcopy(r) # inherit from the right
                        new_entry.update(item) # Key items inherit from the right hand side

                        # add a special case for expression objects:
                        # test that both objects are actually expression objects with _conditions and ["conditions"]:
                        # Funky syntax in case I ever derive a descendent of expression:
                        if "conditions" in item and "conditions" in r: # See if conditions in both genelists:
                            # The below line will escape the rare occasions a genelist is sent that has "conditions" but no _conditions
                            if hasattr(self, '_conditions') and hasattr(gene_list, '_conditions'): # I think we can safely assume both are bonafide expression
                                if self._conditions != gene_list._conditions: # DONT do this if the conditions are identical.
                                    new_entry["conditions"] = item["conditions"] + r["conditions"]
                                    newl._conditions = gene_list._conditions + self._conditions # will update multiple times, whoops.

                                    # only look at the err keys if I am merging the conditions
                                    if "err" in item and "err" in r:
                                        if self._conditions != gene_list._conditions: # DONT do this if the conditions are identical.
                                            new_entry["err"] = item["err"] + r["err"]
                                    elif "err" in new_entry: # Only one list has an err key, output a warning and kill it.
                                        if not __warning_assymetric_errs:
                                            __warning_assymetric_errs = True
                                            config.log.warning("map(): Only one of the two lists has an 'err' key, deleting it")
                                        del new_entry["err"]

                        newl.linearData.append(new_entry)
            elif logic == "notright":
                newl.linearData.append(utils.qdeepcopy(item)) # only inherit from the right, can't inheret from the left, as no matching map

            p.update(index)

        if "image_filename" in kargs and kargs["image_filename"]:
            if not greedy:
                config.log.warning("map(): greedy=False, this can lead to inaccurate results in the Venn diagram")

            venn_proportional = False
            if "venn_proportional" in kargs and kargs["venn_proportional"]:
                venn_proportional = True

            title = "".join(kargs["image_filename"].split(".")[:-1]) # guess a title
            if "title" in kargs:
                title = kargs["title"]

            if "." in kargs["image_filename"]:
                filename_root = "".join(kargs["image_filename"].split(".")[:-1])
            else:
                filename_root = kargs["image_filename"]

            labels = {"left": self.name, "right": gene_list.name, "title": title}
            # and modify the output and draw the venn
            self.draw._vennDiagram2(len(self)-len(newl), len(gene_list)-len(newl), len(newl),
                filename="%s_venn.png" % filename_root, proportional=venn_proportional,
                labels=labels)

        if not silent:
            if logic == "notright":
                config.log.info("map: '%s' vs '%s', using '%s' via '%s', kept: %s items" % (self.name, gene_list.name, map_key, logic, len(newl)))
            else:
                config.log.info("map: '%s' vs '%s', using '%s', found: %s items" % (self.name, gene_list.name, map_key, len(newl)))

        if len(newl.linearData):
            newl._optimiseData()
            return(newl)
        return(None)

    def _findNearbyGenes(self, coords, distance=10000):
        """
        expects:
        coords = chr1:10000-10002
        distance = distance from the coords;
        # relies on you list having a valid tss_loc

        reports dist_to_tss relative to the strand of the gene.
        """

        #assert coords, "Cannot annotate: %s" % coords
        assert "tss_loc" in self.linearData[0], "no available tss_loc key"

        if not self.isChromosomeAvailable(str(coords["chr"])):
            return(False)

        peakCentre = (coords["left"] + coords["right"]) / 2

        ret = []
        for line in self.dataByChr[str(coords["chr"])]:
            line_loc = line["tss_loc"]
            tss_start = line_loc["left"]

            if (peakCentre >= (tss_start-distance)) and (peakCentre <= (tss_start+distance)):
                # add a new _dist_to_tss tag;
                try:
                    strand = line["strand"]
                except Exception:
                    strand = None
                if strand:
                    if strand in positive_strand_labels:
                        line["dist_to_tss"] = peakCentre-tss_start
                    elif strand in negative_strand_labels:
                        line["dist_to_tss"] = tss_start-peakCentre
                else:
                    line["dist_to_tss"] = peakCentre-tss_start
                #print line, peakCentre, tss_start
                ret.append(line)
        return(ret)

    def annotate(self, genelist=None, key_to_match="loc", distance=10000, window=2000,
        image_filename=None, closest_only=False, **kargs):
        """
        **Purpose**
            Annotate some other genelist with this list (usually a genome or expression data set).

            By default annotate will look for a 'tss_loc' key in the genome object. If it can't find
            one then it will use loc, and pointify() it before use

        **Arguments**
            genelist
                A genelist-like object to annotate with the information in this list.

            key_to_match
                must be some sort of location tag::

                    e.g.
                    * loc "chr1:1110100-1110200"
                    * location "chr1:1110100-1110200"
                    * tss_loc "chr1:1110100-1110200"
                    * etc.

            distance
               The distance (in base pairs) to look for a match.

            image_filename (Optional)
                If set to a filename it will cause annotate to save a moving average window of
                the distances annotated.

            window (Optional)
                annotate() draws a histogram of the distance to the
                nearest tss_loc creates an image and saves it to image_filename.
                The window argument specifies the size of the moving window
                used in the calculation of the graph.

            closest_only (Optional, default=False)
                keep the closest annotation only, by default annotate() will collect all
                genes within distance of the peak.

                Note that as many transcripts share the same TSS, only the first
                annotation will be collected. This can be semi-random in which
                actual annotation gets collected. So be careful that is actually what you
                want.

        **Result**
            returns a new genelist-like object inherited from the right-side.
            i.e.::

                result = genelist.map(expression, "tss_loc")

            result will be a microarray.

            however::

                result = expression.map(chip_list, "loc")

            result will be a chip_list

            If image_filename is True:
                Saves the frequency of hits relative to the location tag as a PNG image

        """
        if "peaklist" in kargs: # override if sent a peaklist
            genelist = kargs["peaklist"]

        assert genelist, "you must specify a genelist"
        assert key_to_match in genelist.linearData[0], "The '%s' key is not found in the genelist" % key_to_match
        genome_loc_key = "tss_loc"
        #if genome_loc_key not in self.keys():
        #    genome_loc_key = "loc" # default to loc if genome object
        assert genome_loc_key in list(self.keys()), "The annotation list does not have a valid transcription start-site key, are the genelist and genome the wrong way around?"

        t = [] # blank the data though

        headerLabels = set(["chr", "loc", "name"])
        dist_hist = []

        total_hits = 0

        # self is the genome, genelist has buckets, genome does not
        newl = []
        p = progressbar(len(genelist))
        for index, item in enumerate(genelist.linearData):

            loc = location(loc=item[key_to_match]) # copy
            loc = loc.pointify()
            loc = loc.expand(distance)

            # work out which of the buckets required:
            left_buck = ((loc["left"]-1)//config.bucket_size)*config.bucket_size
            right_buck = ((loc["right"])//config.bucket_size)*config.bucket_size
            buckets_reqd = list(range(left_buck, right_buck+config.bucket_size, config.bucket_size)) # make sure to get the right spanning and left spanning sites

            # get the ids reqd.
            loc_ids = set()
            if buckets_reqd:
                for buck in buckets_reqd:
                    if loc["chr"] in self.buckets:
                        if buck in self.buckets[loc["chr"]]:
                            loc_ids.update(self.buckets[loc["chr"]][buck]) # set = unique ids
            # loc_ids now contains all of the indeces of the items in linearData that need checking

            # Now I check through the loc_ids and collect all peaks within distance
            peaks = []

            anns = []
            for index in loc_ids:
                annotation = self.linearData[index]
                if loc.qcollide(annotation[genome_loc_key]):
                    new_entry = {}
                    for key in annotation:
                        new_entry[key] = annotation[key]

                    # dist_to_tss must be corrected for strand:
                    if "strand" in annotation:
                        if annotation["strand"] in positive_strand_labels:
                            new_entry["dist_to_tss"] = loc.qdistance(annotation[genome_loc_key])
                        elif annotation["strand"] in negative_strand_labels:
                            new_entry["dist_to_tss"] = -loc.qdistance(annotation[genome_loc_key])
                    else:
                        new_entry["dist_to_tss"] = loc.qdistance(annotation[genome_loc_key])
                    anns.append(new_entry)

            if anns:
                if closest_only:
                    closest = anns[0] # This could be done faster by decorate and sort?
                    for i in anns:
                        if abs(i["dist_to_tss"]) < abs(closest["dist_to_tss"]):
                            closest = i

                    new_entry = {} # merge the dictionaries
                    for key in closest:
                        new_entry[key] = closest[key]

                    for key in item:
                        new_entry[key] = item[key] # from the genelist, genelists overwrite genome

                    newl.append(new_entry)

                    if image_filename: # Although result can be biased?
                        dist_hist.append(new_entry["dist_to_tss"])
                else:
                    for annotation in anns:
                        # I want to merge the two lists in a new dict;
                        new_entry = {}

                        for key in annotation:
                            new_entry[key] = annotation[key] # from the genome

                        for key in item:
                            new_entry[key] = item[key]

                        newl.append(new_entry)

                        if image_filename:
                            dist_hist.append(new_entry["dist_to_tss"])

        if not newl:
            config.log.warning("Nothing nearby to annotate for list '%s'" % genelist.name)
            return(None)

        newgl = genelist.shallowcopy() # make a new copy, I need to copy to preserve the attributes.
        # or if it is a derived class.
        newgl.load_list(newl)

        if image_filename:
            # first I need to bin the data (by 1) - then do the moving average.
            linData = numpy.zeros(distance*2) # set-up a blank array
            for item in dist_hist:
                linData[int(item - distance)] += 1
            x, m_average = utils.movingAverage(linData, window)

            self.draw._plot(image_filename,
                m_average, x=x,
                title="%s - loc of site around tss, moving window (%s)" % (newgl.name, window),
                xlabel="Distance to transcription start site (base pairs)",
                ylabel="Frequency (Raw)", **kargs)

        config.log.info("Annotated %s, found: %s within: %s bp" % (newgl.name, len(newgl), distance))
        return(newgl)

    def pointify(self, key="loc"):
        """
        Convert all of the loc coordinates to a single point, centred
        around the middle of the coordinates

        Uses a 'loc' key as the default location to pointify

        **Arguments**

            key (default = "loc")
                specify a location key to pointify.
                defaults to loc


        **Result**

            Returns a new list with 'pointified' coordinates - the coordinates
            are now a single base pair centred around the middle of the
            coordinates.
        """
        assert key in self.linearData[0], "'%s' not in this list" % key

        newl = self.deepcopy()
        for item in newl:
            item[key] = item[key].pointify()

        newl._optimiseData()

        config.log.info("Pointified peaklist '%s'" % self.name)
        return(newl)

    def addFakeKey(self, key=None, value=None):
        """
        **Purpose**
            You need to add a fake key ot the list so that it becomes compatible with some downstream function.
            But you don't care what is in the key, or just want to set the key to a specific value for the entire list

            This is the one to use::

                gl = genelist("a_bed.bed", format=format.minimal_bed)

                print gl
                0: loc: chr1:1000-2000

                # Argh! I need a strand key for the downstream, but I don't actually care what's in strand!

                gl = gl.addFakeKey("strand", "+")

                # Phew! That's better!

                print gl
                0: loc: chr1:1000-2000, strand: +

        **Arguments**
            key (Required)
                the key to name to add

            value (Optional, default=None)
                A value to fill into each new key.

        **Returns**
            A new genelist with the added key.
        """
        assert key , "You must specify a new key name"

        newl = self.deepcopy()
        for item in newl:
            item[key] = value

        newl._optimiseData()
        config.log.info("addFakeKey(): Added a new key '%s'" % key)
        return(newl)

    def expand(self, key="loc", base_pairs=None, side="both"):
        """
        Add base_pairs to the left and right of the location specified in 'key'

        Uses a 'loc' key as the default location to pointify

        **Arguments**

            key (default = "loc")
                specify a location key to pointify.
                defaults to loc

            base_pairs (Required)
                Number of base pairs

            side (Optional, default="both")
                The side to use to expand the location by.

                "both"
                loc = chromosome, left+base_pairs, right+base_pairs

                "left"
                loc = chromosome, left+base_pairs, right

                "right"
                loc = chromosome, left, right+base_pairs

        **Result**

            Returns a new list with 'pointified' coordinates - the coordinates
            are now a single base pair centred around the middle of the
            coordinates.
        """
        assert key in self.linearData[0], "'%s' not in this list" % key

        newl = self.deepcopy()
        for item in newl:
            if side == "both":
                item[key] = item[key].expand(base_pairs)
            elif side == "left":
                item[key] = item[key].expandLeft(base_pairs)
            elif side == "right":
                item[key] = item[key].expandRight(base_pairs)

        newl._optimiseData()

        config.log.info("Expanded '%s' in genelist '%s' by %s base pairs" % (key, self.name, base_pairs))
        return(newl)

    def pointLeft(self, key="loc"):
        """
        pointify the location so that it is set to the left-most base pair.

        i.e.
        loc = (chromosome, left, left)

        **Arguments**

            key (default = "loc")
                specify a location key to pointify.
                defaults to loc

        **Result**
            Returns a new list
        """
        assert key in self.linearData[0], "'%s' not in this list" % key

        newl = self.deepcopy()
        for item in newl:
            item[key] = item[key].pointLeft()

        newl._optimiseData()

        config.log.info("pointLeft genelist %s" % (self.name))
        return(newl)

    def pointRight(self, key="loc"):
        """
        pointify the location so that it is set to the left-most base pair.

        i.e.
        loc = (chromosome, right, right)

        **Arguments**
            key (default = "loc")
                specify a location key to pointify.
                defaults to loc

        **Result**
            Returns a new list
        """
        assert key in self.linearData[0], "'%s' not in this list" % key

        newl = self.deepcopy()
        for item in newl:
            item[key] = item[key].pointRight()

        newl._optimiseData()

        config.log.info("pointRight genelist %s" % (self.name))
        return(newl)

    def collide(self, compare_mode="Collide", loc_key="loc", delta=200, title=None, bins=20,
        add_tags=False, image_filename=None, keep_rank=False, genelist=None, **kargs):
        """
        **Purpose**

            collide two lists using some sort of location tag or key.
            takes two lists, the parent list and one of microarray, peaklist
            or other geneList-like object with a "location"-like key.

            Be careful in your usage of collide() and overlap(), they both merge
            coordinates, but collide() uses the centre of the peak and expands
            it by delta, whereas
            overlap() only expands the left and right border by delta

        **Arguments**
            genelist
                must be some sort of geneList-like object

            loc_key (Optional, Default: "loc")
                the key to use as a location tag.

            delta (Optional, Default: 200)
                the collision fuzzy window +- around the middle of the loc_key tag

            logic (Optional, Default: "and")
                can be one of the below:
                "and" = keep only collisions in both lists.
                "not" = perform a 'not collision', keeping elements that do not collide.
                "notinleft" = Only keep items in the right list if the are 'NOT in the left list'
                "notinright" = Only keep items in the left list if they are 'NOT in the right list'

            image_filename (Optional)
                save a dotplot of the ranks of the colliding peaks. ({filename}_dots.png)
                and save a venn diagram of the overlap. ({filename}_venn.png)
                and a frequency histogram showing the collision distance ({filename}_distance.png)

            bins (Optional, default=20)
                Sets the number of bins in the frequency histogram showing the collision
                distances.

            merge (Optional, Default: False)
                merge the keys of the two lists keys together. This is useful for
                heterogenous lists. For homogenous lists, don't set this
                or set it to False and collide() will keep the resulting list
                homogenous too.
                This will only merge 'extra' keys, it does not overwrite the
                original lists keys where two values exist.

            title (Optional, default = None)
                title for the image figures.

            keep_rank (Optional, default=False)
                if set to true then I will add a new key "ranks" that contains two ranks, the
                ranks of the peak in the (old, new) lists.

            add_tags (Optional, default=False)
                If I find a key with the name specified in 'add_tags' I will add the numbers together,
                If set to False then I will take the tag score from self.

        **Result**

            returns a new genelist derived from the intersect of the keys in
            the original list and the new gene list
            the resulting "loc" key will be a single base pair mid-way
            between the centres of the two colliding coordinates.
        """
        assert genelist, "You must specify a genelist"

        newl = self._unified_collide_overlap(compare_mode, loc_key, delta, title, bins, add_tags, image_filename, keep_rank, genelist=genelist, **kargs)

        len_res = 0 # because newl can be None.
        if newl:
            len_res = len(newl)

        if not config.SILENT:
            if "logic" in kargs and kargs["logic"]:
                config.log.info( "%s (%s) two lists using '%s', found: %s non-overlapping sites" % (compare_mode, kargs["logic"], loc_key, len_res))
            else:
                try:
                    #print len(self), len(genelist), len_res
                    perc = (len_res/ float(len(self)+len(genelist)-len_res) ) # This is actually the Jaccard index
                except ZeroDivisionError: # This can occur if you collide two identical lists
                    perc = 1.0
                config.log.info("%s(): found: %s (Jaccard=%.3f) overlaps in [%s&%s] with '%s' key" % (compare_mode.lower(), len_res, perc, self.name, genelist.name, loc_key))

        return(newl)

    def overlap(self, compare_mode="Overlap", loc_key="loc", delta=200, title=None, bins=20, add_tags=False, image_filename=None,
        keep_rank=False, **kargs):
        """
        **Purpose**

            overlap two lists using some sort of location tag or key.
            takes two lists, the parent list and one of microarray, peaklist
            or other geneList-like object with a "location"-like key.

            Be careful in your usage of collide() and overlap(), they both merge
            coordinates, but collide() uses the centre of the peak, whereas
            overlap() uses the peak span.

        **Arguments**
            genelist
                must be some sort of geneList-like object

            loc_key (Optional, Default: "loc")
                the key to use as a location tag. Only used in this list, in the other list it
                will use whatever key is built by the bucket system. In almost all cases this
                will be tss_loc first, followed by loc

            delta (Optional, Default: 200)
                the collision fuzzy window +- around the left and right
                points of the loc_key tag

            logic (Optional, Default: "and")
                can be one of the below:
                "and" = keep only collisions in both lists.
                "not" = perform a 'not collision', keeping elements that do not collide.
                "notinleft" = Only keep items in the right list if the are 'NOT in the left list'
                "notinright" = Only keep items in the left list if they are 'NOT in the right list'

            image_filename (Optional)
                save a dot-plot map of the colliding lists.

            bins (Optional, default=20)
                Sets the number of bins in the frequency histogram showing the collision
                distances.

            merge (Optional, Default=False)
                merge the two lists keys together. This is useful for
                heterogenouse lists. For homogenous lists, don't set this
                or set it to False and collide() will keep the resulting list
                homogenous too.

            title (Optional, default = None)
                title for the image figures.

            keep_rank (Optional, default=False)
                if set to true then I will add a new key "ranks" that contains two ranks, the
                ranks of the peak in the (old, new) lists.

            add_tags (Optional, default=False)
                If I find a key with the name specified in 'add_tags' I will add the numbers together,
                If set to False then I will take the tag score from self.

        **Result**

            returns a new genelist derived from the intersect of the keys in
            the original list and the new gene list
            the resulting "loc" key will be the left and right most extreme points
            of the two overlapping coordinates.
        """
        newl = self._unified_collide_overlap(compare_mode, loc_key, delta, title, bins, add_tags, image_filename, keep_rank, **kargs)

        len_res = 0 # because newl can be None.
        if newl:
            len_res = len(newl)

        if not config.SILENT:
            if "logic" in kargs and kargs["logic"]:
                config.log.info( "%s (%s) two lists using '%s', found: %s non-overlapping sites" % (compare_mode, kargs["logic"], loc_key, len_res))
            else:
                config.log.info("%s two lists using '%s' key, found: %s overlaps" % (compare_mode, loc_key, len_res))

        return(newl)

    def _unified_collide_overlap(self, compare_mode=None, loc_key="loc", delta=200, title=None, bins=20, add_tags=False, image_filename=None,
        keep_rank=False, **kargs):
        """
        A new unified overlap() collide() code.

        Hopefully less bugs, more features. Or at the very least easier to fix and test.
        """
        assert compare_mode, "compare_mode cannot be False"

        if image_filename:
            config.log.warning("_unified_collide_overlap(): use of image_filename to draw a Venn is not recommended. The values in the Venn")
            config.log.warning("                   may not match the actual overlap values (due to multiple overlaps for a single location)")
            if "logic" in kargs and kargs["logic"] != "and":
                config.log.warning("_unified_collide_overlap(): The output to the Venn (image_filename) when logic != 'and' is not correct")

        __add_tag_keys_warning = False
        __add_tag_keys_float_warning = False
        __add_tag_keys_int_warning = False

        if "genelist" in kargs:
            gene_list = kargs["genelist"]
        else:
            raise AssertionError("_unified_collide_overlap(): No valid genelist-like object") # bodge the error as I do loading at the same time.
            # yes it could be done in a one line assert, but how then do I load gene_list?

        assert loc_key in self[0], "_unified_collide_overlap(): The 'loc_key' '%s' not found in this list" % (loc_key, )
        assert loc_key in gene_list[0], "_unified_collide_overlap(): The 'loc_key' '%s' not found in the other list" % (loc_key, )

        merge = False
        if "merge" in kargs:
            merge = kargs["merge"]

        mode = "and"
        if "logic" in kargs and kargs["logic"] != "and":
            assert kargs["logic"] in ("and", "not", "notinleft", "notinright"), "'%s' logic not found/supported" % kargs["logic"]
            mode = kargs["logic"]
            foundA = [False for x in range(len(self))]
            foundB = [False for x in range(len(gene_list))]

        newl = gene_list.shallowcopy()
        newl.linearData = []

        if image_filename:
            result_array = [0 for x in range(max(len(gene_list), len(self)))]
            x_data = []
            y_data = []
            scatter_data = {"found": [], "not_found": []}
            dist_array = []

        # get a progress bar
        p = progressbar(len(self))

        for indexA, item in enumerate(self):
            if gene_list.isChromosomeAvailable(item[loc_key]["chr"]):
                locA = item[loc_key]
                if compare_mode == "Collide":
                    locA = locA.pointify().expand(delta)
                elif compare_mode == "Overlap" and delta != 0:
                    locA = locA.expand(delta)

                # work out which of the buckets required:
                left_buck = ((locA["left"]-1-delta)//config.bucket_size)*config.bucket_size
                right_buck = ((locA["right"]+delta)//config.bucket_size)*config.bucket_size
                buckets_reqd = list(range(left_buck, right_buck+config.bucket_size, config.bucket_size)) # make sure to get the right spanning and left spanning sites

                # get the ids reqd.
                loc_ids = set()
                if buckets_reqd:
                    for buck in buckets_reqd:
                        if locA["chr"] in gene_list.buckets:
                            if buck in gene_list.buckets[locA["chr"]]:
                                loc_ids.update(gene_list.buckets[locA["chr"]][buck]) # set = unique ids
                # loc_ids now contains all of the indeces of the items in linearData that need checking

                for indexB in loc_ids:
                    other = gene_list.linearData[indexB]

                    if compare_mode == "Collide":
                        locB = gene_list.linearData[indexB][loc_key].pointify().expand(delta)
                    elif compare_mode == "Overlap":
                        locB = gene_list.linearData[indexB][loc_key].expand(delta)

                    #print item["name"],"vs", other["name"], locA, locB, locA.qcollide(locB)

                    if locA.qcollide(locB):
                        if mode != "and":  # only matters if logic="not*"
                            foundB[indexB] = True # gene_list
                            foundA[indexA] = True # self

                        if merge:
                            a = utils.qdeepcopy(other)
                            for k in item:
                                if not k in a:
                                    a.update({k: item[k]})
                        else:
                            a = utils.qdeepcopy(other)

                        a["dist"] = locA.qdistance(locB) # add a new key.

                        # Take the mid-point between the two location centres.
                        a["loc"] = location(chr=locA["chr"], left=int((locA["left"] + locB["left"])/2), right=int((locA["right"] + locB["right"])/2))
                        a["loc"] = a["loc"].pointify().expand(delta)

                        if keep_rank:
                            a["ranks"] = (indexA, indexB)

                        if add_tags and add_tags in other and add_tags in item:
                            if add_tags in item and add_tags in other:
                                try:
                                    a[add_tags] = float(item[add_tags]) + float(other[add_tags])
                                except TypeError:

                                    if not __add_tag_keys_float_warning:
                                        config.log.warning("add_tag key float addition failed")
                                        __add_tag_keys_float_warning = True
                                    try:
                                        a[add_tags] = int(item[add_tags]) + int(other[add_tags])
                                    except TypeError as ValueError:

                                        if not __add_tag_keys_int_warning:
                                            __add_tag_keys_int_warning = True
                                            config.log.warning("add_tag key int addition failed")

                                        a[add_tags] = item[add_tags]
                                        if not __add_tag_keys_warning:
                                            config.log.warning("add_tag key could not be added together")
                                            __add_tag_keys_warning = True

                        newl.linearData.append(a)

                        if image_filename: # collect the data for the images.
                            scatter_data["found"].append((indexA, indexB))
                            result_array[indexA] += 1
                            result_array[indexB] += 1
                            x_data.append(indexA)
                            y_data.append(indexB)
                            dist_array.append(a["dist"])

                    # only gets here if not found.
                    if image_filename: # indexB are not labelled as not Found. Bug or feature?
                        scatter_data["not_found"] = indexA
            p.update(indexA)

        # if logic = "not", I want the opposite of all the items found
        if mode in ("not", "notinleft", "notinright"):
            newl.linearData = []
            if mode == "not" or mode == "notinleft":
                for index, item in enumerate(gene_list):
                    if not foundB[index]:
                        a = utils.qdeepcopy(item)
                        a["loc"] = a["loc"].pointify() # pointify result
                        newl.linearData.append(a)

            if mode == "not" or mode == "notinright":
                for index, item in enumerate(self):
                    if not foundA[index]:
                        a = utils.qdeepcopy(item)
                        a["loc"] = a["loc"].pointify() # pointify result
                        newl.linearData.append(a)

        if image_filename and scatter_data["found"]:
            filename_root = image_filename

            if not title:
                title = filename_root

            x_array, y_array = list(zip(*scatter_data["found"])) # neat trick to unzip list

            # I need to add indexB not founds...

            self.draw.nice_scatter(x=x_array, y=y_array, filename="%s_dot.png" % image_filename,
                do_best_fit_line=True, print_correlation="pearsonr2",
                xlabel=self.name, ylabel=gene_list.name)

            labels = {"left": self.name, "right": gene_list.name, "title": title}

            # and modify the output and draw the venn
            self.draw._vennDiagram2(len(self)-len(newl), len(gene_list)-len(newl), len(newl),
                filename="%s_venn.png" % filename_root, proportional=False, labels=labels)

            #self.draw._qhist(data=dist_array, filename="%s_distance.png" % filename_root, bins=bins, title=title)

        if len(newl) == 0:
            return(None)

        newl._optimiseData()
        return(newl)

    def removeDuplicatesByLoc(self, mode, key="loc", delta=200, use_strand=False):
        """
        **Purpose**
            Remove duplicates in the list based on a location tag.
            Scans all of the loc tags and removes items that overlap.

            Performs a pointify() and expand(delta) on the locations.

            This method will keep the first encountered item of a particular duplicate.
            (i.e. it is a lazy search) As such this method is not expected to give identical
            answers:

            The results you get may be influenced by the order of the list.

            Also note that any pair of peaks can only be merged a single time.
            This stops a condition where peaks can just keep merging into larger and larger peaks
            this is particularly a problem in the 'overlap' mode.

        **Arguments**
            mode (Required)
                You must set one of two modes either:
                    'pointify_expand':
                        first pointify() (i.e. take the mid point between the left and right coords)
                        and then expand(delta) (i.e. symmetrically expand the left and right by 'delta')
                        This is the old default mode
                    'overlap':
                        use the left and right span of the coords to perform an overlap.
                        Note that the new peak will take the min(l1, l2) and the max(r1, r2) where
                        l = left, r = right and 1 = the first genelist and 2 = the second genelist.

            key (Optional, default="loc")
                the loc key to use

            delta (Optional, default=200)
                The number of base pairs to search across.

                Note that delta is NOT used when the mode = 'overlap'

            use_strand (Optional, default=False)
                use the 'strand' key in determining uniqueness

                Note: that it assumes the format of the two strand keys are the same in each
                list (ie. one of '+/-', '-1/+1', etc.)

        **Returns**
            A new genelist containing only unique sites
        """
        if key == 'loc':
            all_keys = list(self.keys())
            assert 'tss_loc' not in all_keys, 'removeDuplicatesByLoc will give incorrect results if your genelist contains both a loc and tss_loc key. Please delete one or the other key'

        assert mode in ('pointify_expand', 'overlap'), "mode is not one of 'pointify_expand', or 'overlap'"

        if use_strand and 'strand' not in self.keys():
            raise AssertionError('use_strand=True, but no "strand" key in this list')

        if not use_strand and 'strand' in self.keys():
            config.log.warning('use_strand=False, but this list has a "strand" key')

        if mode == 'pointify_expand':
            mask = [0] * len(self.linearData) # keep track of masked entries
            newl = []
            for index, item in enumerate(self.linearData):
                locA = item[key].pointify().expand(delta)
                if mask[index] == 0: # only do if not already masked/searched
                    # Do a collision check

                    # work out which of the buckets required:
                    left_buck = ((locA["left"]-1-delta)//config.bucket_size)*config.bucket_size # Add an extra delta for accurate bucket spanning overlaps
                    right_buck = ((locA["right"]+1+delta)//config.bucket_size)*config.bucket_size
                    buckets_reqd = list(range(left_buck, right_buck+config.bucket_size, config.bucket_size)) # make sure to get the right spanning and left spanning sites

                    # get the ids reqd.
                    loc_ids = set()
                    if buckets_reqd:
                        for buck in buckets_reqd:
                            if locA["chr"] in self.buckets:
                                if buck in self.buckets[locA["chr"]]:
                                    loc_ids.update(self.buckets[locA["chr"]][buck]) # set = unique ids
                    # loc_ids now contains all of the indeces of the items in linearData that need checking

                    for indexB in loc_ids:
                        if indexB != index:
                            #other = self.linearData[indexB]
                            locB = self.linearData[indexB][key].pointify().expand(delta)
                            if locA.qcollide(locB):
                                if use_strand:
                                    strandA = self.linearData[index]['strand']
                                    strandB = self.linearData[indexB]['strand']
                                    if strandA == strandB: # Assumes both lists use the same strand format...
                                        mask[indexB] = 1 # found an overlap
                                else:
                                    mask[indexB] = 1

                    mask[index] = 1
                    newl.append(item) # add in the item

            ov = self.shallowcopy()
            ov.load_list(newl)
        elif mode == 'overlap':
            # get the maximum peak size for a decent estimate of delta:
            delta = max([len(i['loc']) for i in self.linearData])
            config.log.info('removeDuplicatesByLoc: delta used from bucket searching {0}'.format(delta))

            if delta > 2000:
                config.log.warning("removeDuplicatesByLoc: The maximum location size is >2000 bp, performance may be poor for removeDuplicatesByLoc()")

            mask = [0] * len(self.linearData) # keep track of masked entries
            newl = []
            for index, item in enumerate(self.linearData):
                locA = item[key]
                if mask[index] == 0: # only do if not already masked/searched
                    # Do a collision check

                    # work out which of the buckets required:
                    left_buck = ((locA["left"]-1-delta)//config.bucket_size)*config.bucket_size # Add an extra delta for accurate bucket spanning overlaps
                    right_buck = ((locA["right"]+1+delta)//config.bucket_size)*config.bucket_size
                    buckets_reqd = list(range(left_buck, right_buck+config.bucket_size, config.bucket_size)) # make sure to get the right spanning and left spanning sites

                    # get the ids reqd.
                    loc_ids = set()
                    if buckets_reqd:
                        for buck in buckets_reqd:
                            if locA["chr"] in self.buckets:
                                if buck in self.buckets[locA["chr"]]:
                                    loc_ids.update(self.buckets[locA["chr"]][buck]) # set = unique ids
                    # loc_ids now contains all of the indeces of the items in linearData that need checking

                    for indexB in loc_ids:
                        if indexB != index:
                            #other = self.linearData[indexB]
                            locB = self.linearData[indexB][key]
                            if locA.qcollide(locB):
                                if use_strand:
                                    strandA = self.linearData[index]['strand']
                                    strandB = self.linearData[indexB]['strand']
                                    if strandA == strandB: # Assumes both lists use the same strand format...
                                        mask[indexB] = 1 # found an overlap
                                else:
                                    mask[indexB] = 1

                    mask[index] = 1
                    newl.append(item) # add in the item

            ov = self.shallowcopy()
            ov.load_list(newl)

        config.log.info("Removed {0} duplicates, {1} remain".format((len(self) - len(ov)), len(ov)))
        return ov

    def removeDuplicates(self, key=None, **kargs):
        """
        **Purpose**
            remove the duplicates in the list and returns a new list;
            keeps the first example it finds

            This will only delete duplicates within the 'key'. For example,
            these three entries in a genelist:

            1: name: Stat3, score: 20, splicing: canonical
            2: name: Stat3, score: 30, splicing: alternate
            3: name: Nanog, score: 40, splicing: alternate

            gl = gl.removeDuplicates("name")

            will give:

            1: name: Stat3, score: 20, splicing: canonical
            3: name: Nanog, score: 40, splicing: alternate

            whilst

            gl = gl.removeDuplicates("splicing")

            will result in:

            1: name: Stat3, score: 20, splicing: canonical
            2: name: Stat3, score: 30, splicing: alternate

        **Arguments**
            key
                The key in which to make search for duplicates.

        **Returns**
            The new list with the duplicates removed.
        """
        assert key, "No key specified"
        assert key in list(self.keys()), "the key '%s' was not found in this genelist" % key

        newl = self.shallowcopy()
        newl.linearData = []
        count = 0
        p = progressbar(len(self))

        for item in self.qkeyfind[key]:
            newl.linearData.append(self.linearData[min(self.qkeyfind[key][item])]) # grab first
            # Will only apply a single item (the earliest) even if there
            # IS only one of these items.

        newl._optimiseData()

        config.log.info("removeDuplicates(): %s duplicates, list now %s items long" % (len(self) - len(newl), len(newl)))
        return(newl)

    def removeExactDuplicates(self):
        """
        **Purpose**
            removes exact duplicates where all of the keys match. Keeping the first
            found copy

        **Returns**
            The new list with the duplicates removed.
        """
        newl = self.shallowcopy()
        newl.linearData = []
        count = 0

        # TODO: Could be done with: list(map(dict, frozenset(frozenset(i.items()) for i in marked_for_deletion)))
        # Lazy at the moment, this function is very rarely used. Better to speed up *ByKey()
        unq = set()
        kord = list(self.linearData[0].keys())# fix the key order

        for item in self.linearData:
            valstr = "".join([str(item[k]) for k in kord])
            if valstr not in unq:
                unq.add(valstr)
                newl.linearData.append(item) # add first item found

        newl._optimiseData()

        config.log.info("removeExactDuplicates(): %s exact duplicates" % (len(self) - len(newl)))
        return(newl)

    def removeEmptyDataByKey(self, key=None):
        """
        **Purpose**
            remove any entry that has empty data within 'key'

            For example, consider this data::

                1: name: Nanog, annot: ,      score: 20
                2: name: Sox2,  annot: yep,   score: 30
                3: name: Stat3, annot: kinda, score:

            You can see some columns with empty data::

                data.removeEmptyDataByKey("annot")

            will result in::

                1: name: Nanog, annot: ,      score: 20
                2: name: Stat3, annot: kinda, score:

            Notice that although score is also empty, this function only considers
            data in the annot: key

            Similarly, using "score" will delete the Stat3 entry only::

                data.removeEmptyDataByKey("score")
                print data

                1: name: Nanog, annot: ,      score: 20
                2: name: Sox2,  annot: yep,   score: 30

        **Arguments**
            key
                The key in which to delete empty data

        **Returns**
            The new list with empty data in key removed
        """
        assert key, "No key specified"

        newl = self.deepcopy()
        oldl = newl.linearData # preserve the copies
        newl.linearData = []
        count = 0
        p = progressbar(len(self))

        for item in oldl:
            if item[key]:
                newl.linearData.append(item) # grab first
            # Will only apply a single item (the earliest) even if there
            # IS only one of these items.

        newl._optimiseData()

        config.log.info("Removed empty data in %s key: %s entries" % (key, len(self) - len(newl)))
        return(newl)

    def act(self, operation, key1, key2, result_key=None):
        """
        **Purpose**
            perfrom some sort of action (or operation) on two keys
            this is an inter-list operation that performs action on key1
            and key2 and then stores the result in result_key

        **Arguments**
            operation
                They type of operation to perform.
                valid operations are:
                ["add", "sub", "div", "pow", "mul", "cat"]
                Most of them are self explanatory and act on numbers.
                The exception is cat, which concatenates two strings into
                one string. The order is always 'key1key2'. cat will
                also happily concatente two numbers into a string
                make sure this is what you actually want.

                add: key1 + key2 # make sure values are numbers
                sub: key1 - key2 # make sire values are numbers
                div: key1 / key2 # fills Division By Zero Errors with 0.0 and returns a warning
                pow: key1**key2
                mul: key1 * key2
                cat: concatenate two strings 'key1key2'

            key1, key2
                the two key names to act on.

            result_key (Optional)
                the name of the key to store the result in.
                (defualts to '<key1><operation><key2>')
                Will overwrite a key if it is already present.

        **Returns**
            This is an IN PLACE operation and returns None.
        """
        valid_operations = ["add", "sub", "div", "pow", "mul", "cat"]

        op_str_map = {"add": "+", "sub": "-", "div": "/", "pow": "^", "mul": "*",
            "cat": "<cat>"} # used in history and for a default storage key

        if operation not in valid_operations:
            raise BadOperationError("%s.act()" % self.__repr__(), operation)

        if not result_key:
            result_key = "%s%s%s" % (key1, op_str_map[operation], key2)

        for item in self.linearData:
            #try:
            if operation == "add":
                item[result_key] = item[key1] + item[key2]
            elif operation == "sub":
                item[result_key] = item[key1] - item[key2]
            elif operation == "div":
                try:
                    item[result_key] = item[key1] / item[key2]
                except ZeroDivisionError:
                    item[result_key] = 0.0
                    config.log.warning("Divide by zero error in the data, value set to '0.0'!")
            elif operation == "mul":
                item[result_key] = item[key1] * item[key2]
            elif operation == "pow":
                item[result_key] = item[key1] ** item[key2]
            elif operation == "cat":
                item[result_key] = "%s%s" % (item[key1], item[key2])

        return(None)

    def load_list(self, list_to_load, name=False):
        """
        **Purpose**
            You've generated your own [{ ... }, { ...}] like list
            (A list of dicts) and you want to either reload it into
            a genelist-like object or load it into an empty genelist.
            This is the method to do that officially.

            This method should be used with great care. Some sanity
            checking is done. But not very much.

        **Arguments**
            list_to_load
                must be a list of dicts.

            name
                Allows you to change the name of the list. By default it will keep
                the previous name.

        **Returns**
            None. This is one of the few IN PLACE methods and returns
            None.
        """
        try:
            list_to_load[0]
            i = list_to_load.__iter__()
        except TypeError:
            raise AssertionError("Type Error, the list appears not to be actually a list")

        try:
            item = list_to_load[0]
            i = [item[k] for k in item]
        except Exception:
            raise AssertionError("Type Error, the list of items appears not to contain a dictionary item")

        #try:
        self.linearData = utils.qdeepcopy(list_to_load)
        #except Exception:
        #    self.linearData = copy.deepcopy(list_to_load)

        # See if we have a name:
        if name:
            self.name = name

        self._optimiseData()
        return(None)

    def raw_data(self):
        """
        **Purpose**
            The internal representation of the data will/may be changing in the
            future - to get access to the underlying list you can get it here.
            In future the actual data may be a lot more messy.
            This will return whatever representation is an iterable()
            containing dictionaries. Like, say, a list of dictionaries
            (the current implementation)

        **Returns**
            The underlying data object. At the moment
            it returns a list of dictionaries, but in the future this may
            change.
        """
        return(self.linearData)

    def unfold(self, key=None):
        """
        **Purpose**
            This method 'unfolds' the list based on a key that also contains
            a list. It's easiest just to illustrate how it works::

                a = [{"type": "one", "data": [1,2]}, {"type": "two", "data": [60,70]}]

                b = a.unfold(key="data")

                b = [{"type": "one", "data": 1},
                    {"type": "one", "data": 2},
                    {"type": "two", "data": 60},
                    {"type": "two", "data": 70}]

            WARNING: This method can be very dangerous in the wrong hands.
            Be certain that the output is actually what you are after!

        **Arguments**
            key
                The key name of a list that will unfold the list

        **Returns**
            The unfolded list

        """
        assert isinstance(self.linearData[0][key], list), "the unfold key '%s' is not a list"

        newl = self.shallowcopy__()
        newl.linearData = []

        for item in self:
            for u in item[key]:
                newitem = utils.qdeepcopy(item)
                newitem[key] = u
                newl.linearData.append(newitem)

        newl._optimiseData()
        return(newl)

    def sample(self, number_to_get=None, seed=None):
        """
        **Purpose**
            Sample n random samples from this genelist

            NOTE: the number of samples to get cannot be larger than
            the size of the list.

        **Arguments**
            number_to_get
                The number of samples to get from the list

            seed
                Rand seed generator, set to a specific seed for code reproducibility

        **Returns**
            A new genelist containing the sampled list.
        """
        assert number_to_get, "Invalid number of samples"
        assert number_to_get < len(self), "the number to get is larger than the list, not possible to sample"

        if seed:
            random.seed(seed)

        shuf = list(range(len(self)))
        random.shuffle(shuf)
        to_get = shuf[0:number_to_get]
        samples = []

        for i in to_get:
            samples.append(utils.qdeepcopy(self.linearData[i]))

        newgl = self.shallowcopy()
        newgl.linearData = samples
        newgl._optimiseData()

        return(newgl)

    def find(self, value):
        """
        **Purpose**
            find value in any key in the genelist.
            Note that find is lazy and returns only the first matching item it finds.
            Use map() for more accurate mapping, or get to collect all

        **Arguments**
            value
                find this 'value' anywhere in the genelist (can be in any key). Can also be any value, number
                string etc.

        **Returns**
            either the first item it finds or False
        """
        for i in self:
            for k in i:
                if i[k] == value:
                    return(i)

        return(False)

    def renameKey(self, old_key_name, new_key_name, keep_old_key=False, replace_existing_key=False):
        """
        **Purpose**
            rename a key with a new name

        **Arguments**
            old_key_name, new_key_name
                the old and new key names

            keep_old_key (Optional, default = False)
                set to True if you want to keep the old key name and value

            replace_existing_key (Optional, default=False)
                Force renameKey() to replace a key with new_key_name.
                For example if your list already has a 'loc' key, and you do this::

                    newgl = nemgl.renameKey("other_loc", "loc")

                glbase will break with the critical error::

                    CRITICAL: new_key_name 'loc' is present in the list already!

                set replace_exisiting_key=True and glbase will no longer complain and will
                silently overwrite the old key

        **Returns**
            A new list with a renamed key
        """
        assert old_key_name, "you must specify an old key name"
        assert new_key_name, "you must specify a new key name"
        assert old_key_name in self.linearData[0], "old_key_name '%s' not found in the list" % old_key_name
        assert old_key_name != new_key_name, 'Cannot replace the old key with the same new key'
        if not replace_existing_key:
            assert new_key_name not in self.linearData[0], "new_key_name '%s' is present in the list already! You might want to set replace_existing_key=True if you are sure" % new_key_name

        newl = []

        newl = utils.qdeepcopy(self.linearData) # copy
        [newitem.update({new_key_name: newitem[old_key_name]}) for newitem in newl] # in-place modify
        if not keep_old_key:
            [newitem.pop(old_key_name) for newitem in newl]

        """
        for item in self.linearData:
            newitem = item.copy()
            newitem[new_key_name] = newitem[old_key_name]
            if not keep_old_key:
                del newitem[old_key_name]
            newl.append(newitem)
        """

        newgl = self.shallowcopy()
        newgl.linearData = newl
        newgl._optimiseData()
        config.log.info("Renamed key '%s' to '%s'" % (old_key_name, new_key_name))
        return(newgl)

    def repairKey(self, key_to_repair, fill_in_key, **kargs):
        '''
        **Purpose**
            genelists will tolerate 'holes' (missing key:values) in individual entries.

            A specific example is loading things like a gtf file, which in that format will tolerate missing attirbute tags

            glbase is quite happy with this, but it may cause problems downstream if you try to grab a single key from a genelist

            This method will fill in the holes in 'key_to_repair' by dragging data from 'fill_in key'
        '''
        assert key_to_repair in list(self.keys()), 'key_to_repair: "%s" not found' % key_to_repair
        assert fill_in_key in list(self.keys()), 'fill_in_key: "%s" not found' % fill_in_key

        replaced = 0
        newl = self.deepcopy()
        for item in newl:
            if key_to_repair not in item or not item[key_to_repair]:
                item[key_to_repair] = item[fill_in_key]
                replaced += 1
        newl._optimiseData()

        config.log.info('repairKey: Repaired %s keys' % replaced)
        return(newl)


    def splitKeyValue(self, key, key_sep=" ", val_sep=":"):
        """
        **Purpose**
            split() the values of key into key:value pairs and add them back into the genelist

            A good example is this:

            After loading a fasta  the entries are like this:

            0: name: ENSP00000451042 pep:known chromosome:GRCh37:14:22907539:22907546:1 gene:ENSG00000223997 transcript:ENST00000415118 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene, seq: EIV

            Note that the name value contains the long string:

            "ENSP00000451042 pep:known chromosome:GRCh37:14:22907539:22907546:1 gene:ENSG00000223997 transcript:ENST00000415118 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene"

            It would be more useful to split this up into key:value pairs so glbase can get at the values:

            fasta = fasta.splitKeyValue("name", " ", ":")

            Results in the more useful:

            0: seq: EIV, transcript_biotype: TR_D_gene, pep: known, gene_biotype: TR_D_gene, gene: ENSG00000223997, transcript: ENST00000415118, chromosome: GRCh37:14:22907539:22907546:1

            NOTE: This will work for a relatively simple example, but will fail rapidly in more complex versions.
            Notice how the first ENSP00000451042 is lost as it does not have a val_sep (:).

        **Arguments**
            key (Required)
                The key to split

            key_sep (Required, default=" ")
                the separator to discriminate key:value pairs

            val_sep (Required, default=",")
                the separator inbetween the key and value.

        **Returns**
            A new genelist
        """

        newl = self.shallowcopy()
        newl.linearData = []

        for row in self:
            newk = dict(row)
            kk = newk.pop(key)

            kvs = kk.split(key_sep)

            for i in kvs:
                if val_sep in i: # no val_sep so omit it;
                    t = i.split(val_sep) # potential error if more than one sep.
                    k = t[0]
                    v = val_sep.join(t[1:])

                    newk[k] = v

            newl.linearData.append(newk)

        newl._optimiseData()
        config.log.info("splitKeyValue(): split '%s' into ~'%s'%s'%s' key value pairs" % (key, len(k), val_sep, len(v)))
        return(newl)

    def splitKey(self, key, key_names, keep_original=False):
        """
        **Purpose**
            split a key that is a list into a new set of keys using key_names.

            For example:
                a = [{"data": [0, 5, 20]}]

                b = a.splitKey("data", ["data_one", "data_two", "data_three"])

                b = [{"data_one": 0, "data_two": 5, "data_three": 20}]

            Some restrictions:
            1. The length of key_names must equal the length of the list stored in Key
            2. All lists stored in 'key' must be the same length.
            3. The data stored in key is assumed to always be in the same order

        **Arguments**
            key
                the key to split

            key_names
                A list of key names to use to split the key up by

            keep_original (Optional, default=False)
                keep the original key as well. Default behaviour deletes the key and its data

        **Returns**
            The new list
        """
        assert key in self.linearData[0], "'%s' not found in this genelist" % key
        assert len(key_names) == len(self.linearData[0][key]), "the key_names list is not the same length as the data in the genelist"

        newl = self.deepcopy()

        for item in newl:
            for i, newk in enumerate(key_names):
                item[newk] = item[key][i]

            if not keep_original:
                del item[key]

        newl._optimiseData()
        return(newl)

    def joinKey(self, new_key_name, formatter, keyA, keyB, keep_originals=False):
        """
        **Purpose**
            Perform a string formatting operation on two keys to jon them together.

            Does this:

            item[new_key_name] = formatter % (keyA, keyB)

        **Arguments**
            new_key_name (Required)
                The new key name.

            formatter (REquired)
                string format operation to perform.

            keyA, keyB (Required)
                the two keys to format.

            keep_originals (Optional, default=False)
                keep the original keys.

        **Returns**
            The new genelist
        """
        assert keyA in self.linearData[0], "keyA '%s' not found in this genelist" % keyA
        assert keyB in self.linearData[0], "keyB '%s' not found in this genelist" % keyB

        newl = self.deepcopy()

        for item in newl:
            item[new_key_name] = formatter % (item[keyA], item[keyB])

            if not keep_originals:
                if not new_key_name == keyA: # Don't delete if it is also the new key.
                    del item[keyA]
                if not new_key_name == keyB:
                    del item[keyB]

        newl._optimiseData()
        return(newl)

    def pie(self, key=None, filename=None, font_size=12, **kargs):
        """
        **Purpose**
            Draw a pie chart based on the values found in 'key'.

            Counts all of the duplicates in pie and then draws a pie-chart and then returns the raw
            data.

        **Arguments**
            key (Required)
                key name to use to construct the pie.

            filename (Required)
                filename to save the image to.

            draw_percents (Optional, default=False)
                draw percent in each category within the pie chart.

            colours (Optional)
                a list of colours, either valid matplotlib names or a list of
                RGB colours (e.g. "#ff0000" for red, "#00ff00" for green etc).
                Note that the number of colours must match the number of categories.

            cmap (Optional)
                An optional matplotlib colourmap to colour the segments with
                Overrides entries in colours

            font_size (Optional, default=12)
                The font size of the pie-chart labels

        **Returns**
            None and an image in filename and the raw data in the form of a dictionary::

                {"class": <count>, "class2": <count>, ...}
        """
        data = {}
        for item in self.qkeyfind[key]:
            data[item] = len(self.qkeyfind[key][item])

        labels = list(data.keys())
        fracs = [data[k] for k in labels]

        newfilename = self.draw.pie(fracs, labels, filename, font_size=font_size, **kargs)

        config.log.info("Saved pie to '%s'" % newfilename)
        return(data)

    def frequencyAgainstArray(self, filename=None, match_key=None, expression=None,
        spline_interpolate=False,
        step_style=False, window=None, **kargs):
        """
        Draw a peaklist and compare against an array.
        Draws a three panel figure showing an array heatmap, the binding events
        and a moving average of the binding events.

        **Arguments**
            expression (Required)
                must be an expression-like object containing numerical data.

            filename (Required)
                a file name for the image including the path.

            window (Optional)
                size of the sliding window to use.
                Defaults to 10% of the length of the microarray.
                If you specify a number here it is in items on the expression array
                and not a percentage.

            bracket (Optional)
                bracket the expression data within a range (e.g. (0,2))

            match_key (Required)
                the key to match between the expression data and the peaks.

            tag_key (Optional)
                use a 'tag score' or some other key:value pair to use
                as a metric for the peaklist.

                By default the list is treated as a binary list, either
                it is or is not a binding site. More sophisticated analysis
                may use the number of seq reads, or fold change over
                background or some other metric to give a more analogue-like
                score. You can specify the score to use using this argument

            spline_interpolate (Optional, default=False)
                spline_interpolate the line graph, valid interpolations are:

                 'slinear', 'quadratic', 'cubic'

                 This will override the default moving average window if set.

            draw_frames (Optional, default=False)
                draw a frame around each of the elements in the figure.

        **Result**
            returns a new peaklist containing the overlapping sites only. The microarray condition value
            will be added to the data and you can get a new microarray out using something like:

            o = p.frequencyAgainstArray()

            newm = expression(loadable_list=o, expn="conditions", ...)

            plots a multi-panel image of the array and a moving average plot of the density of
            chip-seq peaks. Saves the image to filename.
        """
        assert filename, "must specify a filename to save as"
        assert expression, "must provide some expression data"
        assert match_key, "'match_key' is required"
        assert match_key in list(self.linearData[0].keys()), "match_key '%' not found in this list"
        assert match_key in list(expression.linearData[0].keys()), "match_key '%' not found in expression object"
        if spline_interpolate:
            assert spline_interpolate in ('slinear', 'quadratic', 'cubic' ), "'%s' spline_interpolate method not found" % spline_interpolate

        tag_key = None
        if "tag_key" in kargs and kargs["tag_key"]:
            tag_key = kargs["tag_key"]

        # get arrange and sort the data.
        arraydata = expression.getExpressionTable() # for compatability
        #arraydata = numpy.array([array_dict_data[key] for key in array_dict_data]) # needs to be sorted already.
        # There is a potential bug here, with the column names if multiple data is sent.
        if "bracket" in kargs:
            arraydata = self.draw.bracket_data(arraydata, kargs["bracket"][0], kargs["bracket"][1])

        bin = [0 for x in arraydata]

        newgl = self.shallowcopy()
        newl = []
        for index, item in enumerate(expression.linearData):
            found = self._findDataByKeyGreedy(match_key, item[match_key])
            if found:
                newf = []
                for i in found:
                    newi = utils.qdeepcopy(i)
                    newi.update({"conditions": item["conditions"]})
                    newf.append(newi)

                    if "tag_key" in kargs and kargs["tag_key"]:
                        bin[index] += i[kargs["tag_key"]]
                    else:
                        bin[index] = 1

                newl += newf

        if not newl:
            raise AssertionError("frequencyAgainstArray: no matches were found, it's possible this is correct, but it is highly unlikely. I suspect you have not specified match_key correctly")

        newgl.load_list(newl)

        if spline_interpolate:
            f = scipy.interpolate.interp1d(list(range(len(bin))), bin, kind=spline_interpolate)
            xnew = numpy.linspace(0, 40, 40) # A space to interpolate into
            peak_data = list(f(xnew)) # dump out the falues from formula
        else:
            if not window: # get 10% of the list
                window = int(len(bin) * 0.1) # bin is the same length as the data
            peak_data = utils.movingAverage(bin, window, normalise=True)[1]

        # reload/override not really a good way to do this...
        # I should reload a new dict... As I may inadvertantly override another argument?
        kargs["filename"] = filename
        kargs["arraydata"] = arraydata.T
        kargs["row_names"] = expression["name"]
        kargs["col_names"] = expression.getConditionNames()
        if not "bracket" in kargs:
            kargs["bracket"] = [0, 1]

        # Sorry this is a bit of a mess. I was still getting the hang of kargs and
        # Now I can't be arsed to clean it up.
        # Actually, I just spent quite a while getting it into
        # A more sensible shape. Now my enthusiasm has evaporated...
        actual_filename = self.draw._heatmap_and_plot(peakdata=peak_data, bin=bin, row_label_key=match_key, **kargs)

        config.log.info("frequencyAgainstArray: Saved '%s'" % actual_filename)
        return(newgl)

    def islocinlist(self, loc, key="loc", mode="collide", delta=200):
        """
        **Purpose**
            check if the specified loc is in this peaklist.

        **Arguments**
            loc
                the location to check for

            key (Optional, default=False)
                the location key to search in

            mode (Optional, default = False)
                the collision mode to use ("collide" or "overlap")

            delta (Optional, default = 200 bp)
                the delta around the peak to check

        **Returns**
            True or False
        """
        assert loc, "loc must be valid"
        if not isinstance(loc, location):
            loc = location(loc=loc)

        # Easy rejection if the chr bucket is not available:
        if loc["chr"] not in self.buckets:
            return(False)

        if mode == "collide":
            loc = loc.pointify().expand(delta)
        elif mode == "overlap":
            if delta:
                loc = loc.expand(delta)

        # work out which of the buckets is required:
        left_buck = ((loc["left"]-1-delta)//config.bucket_size)*config.bucket_size
        right_buck = ((loc["right"]+delta)//config.bucket_size)*config.bucket_size
        buckets_reqd = list(range(left_buck, right_buck+config.bucket_size, config.bucket_size)) # make sure to get the right spanning and left spanning sites

        # get the ids reqd.
        loc_ids = set()
        if buckets_reqd:
            for buck in buckets_reqd:
                if buck in self.buckets[loc["chr"]]:
                    loc_ids.update(self.buckets[loc["chr"]][buck]) # set = unique ids

        # loc_ids is a set, and has no order.
        #print loc_ids
        for index in loc_ids:
            #print loc.qcollide(self.linearData[index]["loc"]), loc, self.linearData[index]["loc"]
            if loc.qcollide(self.linearData[index]["loc"]):
                return(True)
        return(False)

    def islocinlist_dist(self, loc, key="loc", mode="collide", delta=200):
        """
        **Purpose**
            This will return the distance (in base pairs)

            Note that if the locations are exactly on top of each other the
            distance = 0. Which you may test as no overlap. So, this function returns
            a dictionary with two values: "found" and "distance"

            You should test for location in the following manner::

                res = genelist.islocinlist_dist(loc)
                if res["found"]:
                    # distance is valid
                    # do something with distance here
                    distances += res["distance"]

        **Arguments**
            loc
                the location to check for

            key (Optional, default=False)
                the location key to search in

            mode (Optional, default = False)
                the collision mode to use ("collide" or "overlap")

            delta (Optional, default = 200 bp)
                the delta around the peak to check

        **Returns**
            returns a dictionary, containing
            {"found": <True|False>, "distance": <in bp>}
            returns the distance to the closest
        """
        assert loc, "loc must be valid"
        if not isinstance(loc, location):
            loc = location(loc=loc)

        if mode == "collide":
            loc = loc.pointify().expand(delta)
        elif mode == "overlap":
            if delta:
                loc = loc.expand(delta)

        # work out which of the buckets is required:
        left_buck = ((loc["left"]-1-delta)//config.bucket_size)*config.bucket_size
        right_buck = ((loc["right"]+delta)//config.bucket_size)*config.bucket_size
        buckets_reqd = list(range(left_buck, right_buck+config.bucket_size, config.bucket_size)) # make sure to get the right spanning and left spanning sites

        # get the ids reqd.
        loc_ids = set()
        if buckets_reqd:
            for buck in buckets_reqd:
                if loc["chr"] in self.buckets:
                    if buck in self.buckets[loc["chr"]]:
                        loc_ids.update(self.buckets[loc["chr"]][buck]) # set = unique ids

        # loc_ids is a set, and has no order.
        for index in loc_ids:
            if loc.qcollide(self.linearData[index]["loc"]):
                return({"found": True, "distance": loc.qdistance(self.linearData[index]["loc"])})
        return({"found": False, "distance": None})

    def distribution(self, genelist=None, random_lists=None, filename=None, **kargs):
        """
        **Purpose**
            draw a distribution plot, work out the frequency
            based on the frequency of the nearest annotated gene
            found in genome

        **Arguments**
            genome (Required)
                Some sort of reference set or collection of items you want to
                annotate against. Must contain a valid "loc" key

            random_lists (Required)
                A set of genelist-like objects with a valid "loc" key

            filename (Required)
                the filename to use to save the resulting figure to.

        **Returns**
            The values:
                hist, back, back_err, labels

            in a tuple

            hist = the numbers of genes in each bin
            back = the numbers of genes in each bin, from the background
            back_err = the standard error of the background data.
            labels = labels for the items, in order with hist and back

            saves and image to filename

            prints a table of distributions to the console (stdout)
        """

        assert filename, "no filename specified"

        if not isinstance(random_lists, list):
            background = [random_lists] # make a one item'd list

        hist = {"desert": 0, (-200, -100): 0, (-100, -50): 0, (-50, -10): 0, (-10, 0): 0,
            (0, 10): 0, (10, 50): 0, (50, 100): 0, (100, 200): 0}

        back = {"desert": [], (-200, -100): [], (-100, -50): [], (-50, -10): [], (-10, 0): [],
            (0, 10): [], (10, 50): [], (50, 100): [], (100, 200): []}

        for i, glist in enumerate([self] + random_lists):
            ann = genome.annotate(genelist=glist, distance=int(0.2e6), closest_only=True)
            #
            tback = {(-200, -100): 0, (-100, -50): 0, (-50, -10): 0, (-10, 0): 0,
                (0, 10): 0, (10, 50): 0, (50, 100): 0, (100, 200): 0}

            for d in ann:
                if i == 0:
                    for k in hist:
                        c = d["dist_to_tss"] / 1000
                        if c >= k[0] and c < k[1]:
                            hist[k] += 1
                else:
                    for k in back:
                        c = d["dist_to_tss"] / 1000
                        if c >= k[0] and c < k[1]:
                            tback[k] += 1

            # append the back
            # and work out desert
            if i == 0:
                hist["desert"] = len(self) - len(ann)
                s = sum([hist[k] for k in hist])
            else:
                back["desert"].append(float(len(glist) - len(ann))/len(glist))
                for k in tback:
                    back[k].append(float(tback[k])/len(glist))
                s = sum([tback[k] for k in tback] + [len(glist) - len(ann)])

        kord = ["desert",
            (-200, -100),
            (-100, -50),
            (-50, -10),
            (-10, 0),
            (0, 10),
            (10, 50),
            (50, 100),
            (100, 200)]

        # arrange a table:
        print("class\tann\tmean\tstdev\tstderr")
        for k in kord:
            t = [k, hist[k], numpy.mean(back[k]), numpy.std(back[k]), numpy.std(back[k])/math.sqrt(len(back[k]))]
            print("\t".join([str(i) for i in t]))

        labels = ["Gene desert", "-200 kb to -100 kb", "-100 kb to -50 kb", "-50 kb to -10 kb", "-10 kb to 0 kb",
            "0 kb to 10 kb", "10 kb to 50 kb", "50 kb to 100 kb", "100 kb to 200 kb"]

        data = numpy.array([hist[k] for k in kord], dtype=numpy.float64)
        rand = numpy.array([numpy.mean(back[k]) for k in kord])

        total = sum(data)

        data = (data / total)*100
        rand = rand*100
        err = [(numpy.array(back[k]))*100 for k in kord]
        err = [numpy.std(i) for i in err]

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        ax.set_position([0.1, 0.2, 0.8, 0.70])

        x_bar = numpy.arange(len(labels))
        width = 0.35
        ax.bar(x_bar, data, width, color="orange", label=self.name, ec="none")
        ax.bar(x_bar + width, rand, width, color="black", yerr=err, ecolor="black", label="Randomly selected background sites")
        ax.set_xticklabels(labels)
        ax.set_ylim([0,40])
        ax.legend()

        ax.set_ylabel("Percent in category", size=20)

        ax.set_xticks(x_bar+width)
        fig.autofmt_xdate()
        [t.set_fontsize(16) for t in ax.get_xticklabels()]
        [t.set_fontsize(16) for t in ax.get_yticklabels()]

        actual_filename = self.draw.savefigure(fig, filename)
        return(data, rand, err, labels)

    def genome_distribution(self, genome, random_lists, filename=None, **kargs):
        """
        **Purpose**
            draw a genome distribution plot, based on the frequency of the nearest annotated gene
            found in genome

        **Arguments**
            genome (Required)
                the genome to use as the annotation, must contain a valid "tss_loc" key and must
                be a "genome" object with an "annotate()" method

            random_lists (Required)
                A set of genelist-like objects with a valid "loc" key

            filename (Required)
                the filename to use to save the resulting figure to.

        **Returns**
            The values:
                hist, back, back_err, labels

            in a tuple

            hist = the numbers of genes in each bin
            back = the numbers of genes in each bin, from the background
            back_err = the standard error of the background data.
            labels = labels for the items, in order with hist and back

            saves and image to filename

            prints a table of distributions to the console (stdout)
        """
        if random_lists:
            if not isinstance(random_lists, list):
                background = [random_lists] # make a one item'd list

        hist = {"desert": 0, (-200, -100): 0, (-100, -50): 0, (-50, -10): 0, (-10, 0): 0,
            (0, 10): 0, (10, 50): 0, (50, 100): 0, (100, 200): 0}

        back = {"desert": [], (-200, -100): [], (-100, -50): [], (-50, -10): [], (-10, 0): [],
            (0, 10): [], (10, 50): [], (50, 100): [], (100, 200): []}

        if random_lists:
            todos = [self] + random_lists
        else:
            todos = [self]

        for i, glist in enumerate(todos):
            ann = genome.annotate(genelist=glist, distance=200000, closest_only=True)

            tback = {(-200, -100): 0, (-100, -50): 0, (-50, -10): 0, (-10, 0): 0,
                (0, 10): 0, (10, 50): 0, (50, 100): 0, (100, 200): 0}

            for d in ann:
                if i == 0:
                    for k in hist:
                        if k == 'desert':
                            continue

                        c = d["dist_to_tss"] / 1000
                        if c >= k[0] and c < k[1]:
                            hist[k] += 1
                else:
                    for k in back:
                        if k == 'desert':
                            continue

                        c = d["dist_to_tss"] / 1000
                        if c >= k[0] and c < k[1]:
                            tback[k] += 1

            # append the back
            # and work out desert
            if i == 0:
                hist["desert"] = len(self) - len(ann)
                s = sum([hist[k] for k in hist])
            else:
                back["desert"].append(float(len(glist) - len(ann))/len(glist))
                for k in tback:
                    back[k].append(float(tback[k])/len(glist))
                s = sum([tback[k] for k in tback] + [len(glist) - len(ann)])

        kord = ["desert",
            (-200, -100),
            (-100, -50),
            (-50, -10),
            (-10, 0),
            (0, 10),
            (10, 50),
            (50, 100),
            (100, 200)]

        labels = ["Gene desert", "-200 kb to -100 kb", "-100 kb to -50 kb", "-50 kb to -10 kb", "-10 kb to 0 kb",
            "0 kb to 10 kb", "10 kb to 50 kb", "50 kb to 100 kb", "100 kb to 200 kb"]

        dlabels = dict(list(zip(kord, labels)))

        # arrange a table:
        #print "class\tann\tmean\tstdev\tstderr"
        for k in kord:
            t = [dlabels[k], hist[k], numpy.mean(back[k]), numpy.std(back[k]), numpy.std(back[k])/math.sqrt(len(back[k]))]
            #print "\t".join([str(i) for i in t])

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        ax.set_position([0.1, 0.2, 0.8, 0.70])

        data = numpy.array([hist[k] for k in kord], dtype=numpy.float64)

        total = sum(data)

        data = (data / total)*100

        x_bar = numpy.arange(len(labels))
        width = 0.35
        ax.bar(x_bar, data, width, color="orange", label=self.name, ec="none")

        if random_lists:
            rand = numpy.array([numpy.mean(back[k]) for k in kord])
            rand = rand*100
            err = [(numpy.array(back[k]))for k in kord]
            err = [numpy.std(i)*100 for i in err]
            ax.bar(x_bar + width, rand, width, color="black", yerr=err, ecolor="black", label="Background")
        else:
            rand = None # spoof entries for the return()
            err = None

        ax.set_xticklabels(labels)
        ax.set_ylim([0,max(max(rand), max(data))])
        ax.legend()

        ax.set_ylabel("Percent in category", size=20)

        ax.set_xticks(x_bar+width)
        fig.autofmt_xdate()
        [t.set_fontsize(7) for t in ax.get_xticklabels()]
        [t.set_fontsize(7) for t in ax.get_yticklabels()]

        self.draw.do_common_args(ax, **kargs)
        if filename:
            actual_filename = self.draw.savefigure(fig, filename)
        return(data, rand, err, labels)

    def hist(self, key=None, filename=None, range=None, suppress_zeros=False, log=None,
        kde=False, covariance=0.2, **kargs):
        """
        **Purpose**
            Draw a normal histogram of the expression values

        **Arguments**
            filename (Required)
                the filename to save the resulting image to.

            key (Required)
                key to use for numeric data.

            covariance (Optional, default=0.2)
                undocumented wierdness (Sorry)

            range (Optional, default=(0, max(expression)))
                This is the maximum value to build the KDE over. (ie. 0 ... mmax).
                You probably really want to chage this value, as small outliers will distort the
                distribution massively.
                Play around with the value until it gets you a nice normal-like distribution.

            suppress_zeros (Optional, default=False)
                If this is set to True expression values of zero are not included in the
                histogram.

            log (Optional, default=False)
                log transform the data.
                At the moment only log2 is supported.

            kde (Optional, default=False)
                use kernel density estimation to smooth the data

            covariance (Optional, default=0.2)
                Value for KDE smoothing.

        **Returns**
            The values used to build the histogram
            and an image in filename
        """
        assert filename, "Need filename, leh!"
        assert key, "You must specify a key that contains numerical data"

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        values = self[key]
        if not range: # sample a range if none specified
            range = (0, max(expn_values))

        if suppress_zeros:
            values = [v for v in expn_values if int(v*10000) > 0]

        values = numpy.array(values)

        if log:
            values = numpy.log2(values)

        if kde:
            values = utils.kde(values, range=range, covariance=covariance, bins=100)

        ax.hist(values, bins=200, range=range, normed=True, histtype='step', label=k)

        ax.legend(ncol=1)

        self.draw.do_common_args(ax, **kargs)
        real_filename = self.draw.savefigure(fig, filename)

        config.log.info("Saved '%s'" % real_filename)
        return(values)

    def bar_chart(self, filename=None, labels=None, data=None, percents=False,
        **kargs):
        """
        **Purpose**
            Draw a barchart of all of the values in the gene list.

        **Arguments**
            data (Required)
                The key name to use to count the different classes.

            labels (Required)
                The key in the genelist to use for labels.

            filename (Required)
                The name to save the image to.

            percents (Optional, default=True)
                present as percentage of the appropriate list (i.e. normalise for list length)

        **Returns**
            None
        """

        data = self[data]
        labels = self[labels]
        if percents:
            data = numpy.array(data)
            data = (data * 100.0) / data.sum()

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)

        x_bar = numpy.arange(len(labels))
        width = 0.35

        ax.bar(x_bar, data, width, color="black", label=self.name, ec="black")

        ax.set_xticklabels(labels)
        ax.set_ylim([0,max(data)+2])
        ax.set_xlim([x_bar[0]-(1.0/len(x_bar)), x_bar[-1]+1])

        if percents:
            ax.set_ylabel("Percent", size=20)
        else:
            ax.set_ylabel("Number", size=20)
        #ax.set_xlabel("Category name", size=20)

        ax.set_xticks(x_bar+width)

        self.draw.do_common_args(ax, **kargs)

        fig.autofmt_xdate()
        fig.tight_layout()

        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("bar_chart(): Saved '%s'" % actual_filename)

        return(None)

    def frequency_bar_chart(self, filename=None, random_backgrounds=None, key=None, percents=True,
        horizontal=False, **kargs):
        """
        **Purpose**
            Draw a bar_chart of all the different items found in key, the x-axis are all
            of the different keys and the y-axis is the frequency.

            An optional random_background can be sent for comparison.

        **Arguments**
            key (Required)
                The key name to use to count the different classes.

            filename (Required)
                The name to save the image to.

            random_backgrounds (Optional, default=None)
                A list of random_backgrounds to plot a comparison frequency against.

            percents (Optional, default=True)
                present as percentage of the appropriate list (i.e. normalise for list length)

            horizontal (Optional, default=False)
                rotate the bar chart so that frequency is on the x-axis and category is
                on the y-axis.

        **Returns**
            A dictionary containing the frequencies of the data and the background (or None if not
            background data used).
        """
        assert filename, "need a filename!"
        assert key, "need a key!"
        assert key in list(self.keys()), "key not found in this data"

        if random_backgrounds:
            if not isinstance(random_backgrounds, list):
                random_backgrounds = [random_backgrounds]
            for item in random_backgrounds:
                assert key in item.linearData[0], "key '%s' not found in the random background" % key # check the first entry only.

        # get the counts:
        # All stored in qkeyfind

        counts = {}
        for k in self.qkeyfind[key]:
            if percents:
                counts[k] = len(self.qkeyfind[key][k]) * 100.0 / len(self)
            else:
                counts[k] = len(self.qkeyfind[key][k])

        rand_counts = None
        rand_stderr = None
        if random_backgrounds:
            rand_counts = {}
            for rand in random_backgrounds:
                for k in rand.qkeyfind[key]:
                    if k not in rand_counts:
                        rand_counts[k] = []

                    if percents: # convert or not:
                        rand_counts[k].append((len(rand.qkeyfind[key][k]) * 100.0) / len(rand))
                    else:
                        rand_counts[k].append(len(rand.qkeyfind[key][k]))

            # flatten rand_counts and get stderr
            rand_stderr = {}
            for k in rand_counts:
                rand_stderr[k] = numpy.std(rand_counts[k])
                rand_counts[k] = numpy.mean(rand_counts[k])

        if not random_backgrounds:
            labels = list(counts.keys())
        else: # collect eh rand keys too.
            labels = [] # make sure all possible labes are collected:

            for k in counts:
                if k not in labels:
                    labels.append(k)

            for k in rand_counts:
                if k not in labels:
                    labels.append(k)

        if horizontal and not "aspect" in kargs: # set aspect if horizontal set and no aspect specified by user
            kargs["aspect"] = "long"

        fig = self.draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        ax.set_position([0.1, 0.2, 0.8, 0.70])

        data = []
        rand_data = []
        rand_err = []
        for k in labels: # unfold data:
            if k in counts:
                data.append(counts[k])
            else:
                data.append(0) # make sure rands with zero counts still get added.

            if random_backgrounds:
                if k in rand_counts:
                    rand_data.append(rand_counts[k])
                    rand_err.append(rand_stderr[k])
                else:
                    rand_data.append(0) # block error for missing keys
                    rand_err.append(0)

        x_bar = numpy.arange(len(labels))
        width = 0.35

        if horizontal:
            # matplotlib buggers up the location, adjust it;
            ax.set_position([0.3, 0.08, 0.65, 0.88])

            ax.barh(x_bar, data, width, color="none", label=self.name, ec="black") # barh auto-rotates
            ax.set_ylim([x_bar[0]-(1.0/len(x_bar)), x_bar[-1]+1])

            if rand_counts:
                ax.barh(x_bar+width, rand_data, width, color="black", xerr=rand_err, ecolor="black", label="Random background")

            ax.set_yticklabels(labels)

            if rand_data:
                ax.set_xlim([0,max([max(rand_data)+1, max(data)+1])])
            else:
                ax.set_xlim([0,max(data)+1])

            if percents:
                ax.set_xlabel("Percent in category", size=20)
            else:
                ax.set_xlabel("Number in category", size=20)
            ax.set_ylabel("Category name", size=16)

            ax.set_yticks(x_bar+width)

        else:
            ax.bar(x_bar, data, width, color="none", label=self.name, ec="black")

            if rand_counts:
                ax.bar(x_bar + width, rand_data, width, color="black", yerr=rand_err, ecolor="black", label="Random background")
            ax.set_xticklabels(labels)

            if rand_data:
                ax.set_ylim([0,max([max(rand_data)+2, max(data)+2])])
            else:
                ax.set_ylim([0,max(data)+2])

            ax.set_xlim([x_bar[0]-(1.0/len(x_bar)), x_bar[-1]+1])

            if percents:
                ax.set_ylabel("Percent in category", size=20)
            else:
                ax.set_ylabel("Number in category", size=20)
            ax.set_xlabel("Category name", size=20)

            ax.set_xticks(x_bar+width)

        ax.legend()

        if not horizontal:
            fig.autofmt_xdate()
        [t.set_fontsize(16) for t in ax.get_xticklabels()]
        [t.set_fontsize(16) for t in ax.get_yticklabels()]

        self.draw.do_common_args(ax, **kargs)

        actual_filename = self.draw.savefigure(fig, filename)
        config.log.info("frequency_bar_chart(): Saved '%s'" % actual_filename)

        # turn the result into a genelist:
        gl = genelist(name="bar_chart result")
        loadlist = []

        if random_backgrounds: # This could probably be all done much nicer if I thought about it better.
            all_keys = set(counts.keys()) | set(rand_counts.keys())
        else:
            all_keys = set(counts.keys())

        if percents:
            ppkey = "percent"
        else:
            ppkey = "count"

        for k in all_keys: # build a new genelist summarising the results
            if k in counts: # sometimes back has, but not counts
                item = {"name": k, ppkey: counts[k]}
            else:
                item = {"name": k, ppkey: 0}
            if random_backgrounds:
                if k in rand_counts:
                    item["background_%s" % ppkey] = rand_counts[k]
                    item["background_stderr"] = rand_stderr[k]
                else: # no key.
                    item["background_%s" % ppkey] = 0
                    item["background_stderr"] = 0
            loadlist.append(item)
        gl.load_list(loadlist)

        return(gl)

    def remove(self, key=None, value=None):
        """
        remove a particular item (or items) by key and value

        **Arguments**
            key (Required)
                The key to look for value in

            value (Required)
                if this value is found in key. Then remove the item
        **Retuns**
            A new genelist, with the item(s) removed
        """
        assert key, "remove(): You must specify a key"
        assert value, "remove(): You must specify a value"
        assert key in list(self.keys()), "remove(): key '%s was not found in this genelist" % key

        newl = self.shallowcopy() # Just use a view as we are not really modifying the data.
        newl.linearData = [] # get a new view
        removed = 0

        for item in self.linearData:
            if not item[key] == value:
                newl.linearData.append(item)
            else:
                removed += 1

        newl._optimiseData()
        config.log.info("remove(): removed %s items" % removed)
        return(newl)

    def cumulative_annotation_plot(self, filename, peaklists, randoms=None, annotation_range=(100, 50000, 500), **kargs):
        """
        **Purpose**
            Plot the cumulative number of your set of genes with a peak nearby.

            This list should be a list of genes, with a valid "tss_loc" key.

        **Arguments**
            filename (Required)
                The filename to save to.

            peaklist (Required)
                a list of genelists containing lists of peaks.

            randoms (Optional, default=None)
                A list of randoms for comparison.

            annotation_range (Optional, default=(100, 50000, 500))
                The genomic range to annotate over.

        **Returns**
            A file saved in filename.
        """
        assert filename, "cumulative_annotation_plot(): You must provide a filename"

        xs = list(range(100, 50000, 500))
        result = []
        annotated_genes = self.annotate(genelist=peaklists, distance=annotation_range[1]).removeDuplicates("enst")

        p = progressbar(len(xs))
        for n, d in enumerate(xs):
            this_amount = 0
            for gene in annotated_genes:
                #print gene["dist_to_tss"], abs(gene["dist_to_tss"]) < d
                if abs(gene["dist_to_tss"]) < d:
                    this_amount += 1
            result.append((this_amount/float(len(self))) * 100.0) # get within 100%

            if randoms:
                for r in random_stat1:
                    r_ann = enst.annotate(genelist=r, distance=d)
                    if r_ann: # can sometimes be zero.
                        mapped = r_ann.map(genelist=wtko, key="ensg")
                        if mapped:
                            mapped = mapped.removeDuplicates("ensg") # ... bad map key

                    if not r.name in rands_st1:
                        rands_st1[r.name] = []
                    if r_ann and mapped:
                        rands_st1[r.name].append(len(mapped)/len(wtko) * 100.0)
                    else:
                        rands_st1[r.name].append(0.0)
            p.update(n)

        fig = self.draw.getfigure(**kargs)

        ax = fig.add_subplot(111)
        ax.plot(xs, result, color="red")
        if randoms:
            for k in randoms:
                ax.plot(xs, rands_st3[k], color="grey", alpha=0.5)
        ax.set_title(self.name)
        ax.legend(loc=2)

        self.draw.savefigure(fig, filename, **kargs)
        return(None)

genelist = Genelist # Basically used only im map() for some dodgy old code I do not want to refactor.

"""
Think of me as a delayed version of geneList

* this is the only special case of geneList (delayed)
"""
import sys, os, time, copy, csv, gzip

from . import config
from . import utils
from .data import *
from array import array as qarray
from .draw import draw
from .genelist import genelist
from .history import historyContainer
from .errors import AssertionError, NotSupportedError, DelayedListError
from .location import location

# try to load non-standard libs.
try:
    import matplotlib.pyplot as plot
    MATPLOTLIB_AVAIL = True
except:
    print("Warning: matplotlib not available or not installed")
    MATPLOTLIB_AVAIL = False

class delayedlist(genelist):
    """
    **Purpose**
    
        delayedlist is very similar to a genelist - except the data never
        makes it into memory. This allows you to iterate over huge files from disk.
    
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

        gzip (Optional, default=False)
            The input file is a gzip file. 

        format
            format specifier. (see docs... complex)

    """
    
    def __init__(self, filename=None, format=None, force_tsv=False, gzip=False, **kargs):
        genelist.__init__(self) # no kargs though. I want the mpty version.
        
        self.__len_estimate = None

        assert filename, "No Filename"
        assert os.path.exists(filename), "%s not found" % (filename)
        assert format, "You must provide a format for delayedlist. I cannot guess its format."

        self.path = os.path.split(os.path.realpath(filename))[0]
        self.filename = os.path.split(os.path.realpath(filename))[1]
        self.fullpath = filename
        self.filehandle = None
        self.format = format # override default
        self.gzip = gzip

        if force_tsv: 
            self.format["force_tsv"] = True

        config.log.info("Bound '%s' as a delayedlist" % filename)
        self._optimiseData()

    def __repr__(self):
        return("glbase.delayedlist")

    def collide(self, **kargs):
        """
        Note: the 'logic' command is not supported for delayedlists
        only "and" can be performed.
        """
        self._optimiseData()
        if "logic" in kargs: raise NotSupportedError("'logic' commands not supported for delayedlist.collide()")

        assert "peaklist" in kargs or "genelist" in kargs or "microarray" in kargs, "You must provide a genelist-like object"
        assert "loc_key" in kargs, "You must provide a 'loc_key' name"
        
        if "peaklist" in kargs: 
            gene_list = kargs["peaklist"]
        elif "genelist" in kargs: 
            gene_list = kargs["genelist"]
        elif "microarray" in kargs:
            gene_list = kargs["microarray"]
        
        assert kargs["loc_key"] in gene_list[0]
        assert kargs["loc_key"] in next(self.__iter__()) # get an item and test it
        self._optimiseData()

        delta = 200
        if "delta" in kargs: delta = kargs["delta"]

        return(genelist.collide(self, genelist=gene_list, loc_key=kargs["loc_key"], delta=delta, merge=True))

    def overlap(self, delta=200, **kargs):
        """
        Note: the 'logic' command is not supported for delayedlists
        only "and" can be performed.
        """
        self._optimiseData()

        assert "peaklist" in kargs or "genelist" in kargs or "microarray" in kargs, "You must provide a genelist-like object"
        assert "loc_key" in kargs, "You must provide a 'loc_key' name"
        assert 'logic' not in kargs, "the 'logic' system is not supported if one of the genelists is a delayedlist"

        # get the genelist object:
        if "peaklist" in kargs: gene_list = kargs["peaklist"]
        elif "genelist" in kargs: gene_list = kargs["genelist"]
        elif "microarray" in kargs: gene_list = kargs["microarray"]

        assert kargs["loc_key"] in gene_list[0]
        assert kargs["loc_key"] in next(self.__iter__()) # get an item and test it

        self._optimiseData() # reset the __iter__

        return(genelist.overlap(self, genelist=gene_list, loc_key=kargs["loc_key"], delta=delta, merge=True))

    def map(self):
        raise AssertionError('delayedlists cannot be mapped in this direction, try the other way: genelist.map(genelist=delayedlist, key="...", ...)')

    def __len__(self):
        # I need to collect an estimate
        if not self.__len_estimate:
            if not self.gzip:
                f = open(self.fullpath, 'rb')           
                lines = 0
                for line in f: lines += 1

                self.__len_estimate = lines-1 # start from 0
                
            else: # gzipped file variant
                f = gzip.open(self.fullpath, 'rb') # must be rb :(  
                lines = 0
                for line in f.readlines(): lines += 1
                self.__len_estimate = lines-1 # start from 0
                
        return(self.__len_estimate) # 

    def __getitem__(self, index):
        """
        (Override)
        confers a = geneList[0] behaviour
        This is broken. It returns only the very first entry,
        whatever 'index' is sent to it, it disregards it.
        This only continues to exist for compatability with some internal
        routines.
        """
        self._optimiseData()
        return(next(self.__iter__()))
        self._optimiseData()

    def __iter__(self):
        """
        (Override)
        make the geneList behave like a normal iterator (list)
        """
        try:
            column = next(self.__reader) # get started
            self.cindex += 1
            
            while column:
                d = None
                while not d:
                    if column: # list is empty, so omit.  
                        if "commentlines" in self.format:
                            if column[0][0] == self.format["commentlines"]: # csv will split the table and returns a list
                                column = None # force a skip of this row, don't use continue it will just hang
                                d = None
                            elif column[0].startswith(self.format["commentlines"]):
                                column = None
                                d = None
                    
                    if column:
                        if "keepifxin" in self.format:
                            if True in [self.format["keepifxin"] in i for i in column]:
                                if (not (column[0] in typical_headers)):
                                    d = self._processKey(self.format, column)
                            else:
                                d = None # not present, skip this line

                        else: # just do normally
                            if (not (column[0] in typical_headers)):
                                d = self._processKey(self.format, column)
                    
                    if not d: # d is bad, grab another
                        column = next(self.__reader)
                        self.cindex += 1
                
                # I do quoting = NONE now, so I need to manually deal with containing quotes.  
                for k in d: # d must be valid to get here
                    if isinstance(d[k], str): # only for strings 
                        if d[k][0] == "'" and d[k][-1] == "'":
                            d[k] = d[k].strip("'")
                        if d[k][0] == '"' and d[k][-1] == '"': 
                            d[k] = d[k].strip('"')
                
                yield d # d should be valid
                column = next(self.__reader) # get the next item
                self.cindex += 1
        except StopIteration:
            self._optimiseData()
            raise StopIteration

    def _optimiseData(self):
        """
        (Override)
        Impossible to optimise the data.
        so we just reset the entry point.
        This makes the iterator work like you would expect:
        a new iterator will go back to the beginning of the list.
        """
        if self.filehandle: 
            self.filehandle.close()

        if not self.gzip:
            self.filehandle = open(self.fullpath, "rt")
        else:
            self.filehandle = gzip.open(self.fullpath, 'rt') # must be rb :(
        
        
        if "force_tsv" in self.format and self.format["force_tsv"]:
            self.__reader = csv.reader(self.filehandle, dialect=csv.excel_tab, quoting=csv.QUOTE_NONE)
        elif "dialect" in self.format:
            self.__reader = csv.reader(self.filehandle, dialect=self.format["dialect"], quoting=csv.QUOTE_NONE)
        else:
            self.__reader = csv.reader(self.filehandle, quoting=csv.QUOTE_NONE)

        if "skiplines" in self.format:
            if self.format["skiplines"] != -1: # no skipped lines, good to go.
                for i, x in enumerate(self.__reader):
                    if i == self.format["skiplines"]:
                        break
        else: # default behaviour of genelist is to always skip the first line
            next(self.__reader)

        if "skiptill" in self.format:
            done = False
            while not done:
                for line in self.__reader:
                    if self.format["skiptill"] in line:
                        done = True
                        break

        self.linearData = self.__iter__()
        self.cindex = 0
        return(True)

    def __str__(self):
        self._optimiseData()

        # load in a bunch of data and dump it into self.linearData
        temp_data = []
        for index, item in enumerate(self):
            temp_data.append(item)
            if index > config.NUM_ITEMS_TO_PRINT-2:
                break # only get the first n data.s
        self.linearData = temp_data
        ret = genelist.__str__(self)
        self._optimiseData()
        return("%s\nThis is a delayedlist - only the first %s entries are shown" %(ret, config.NUM_ITEMS_TO_PRINT))

    def save(self):
        raise NotSupportedError("Cannot save a binary representation of a delayedlist")

    def saveCSV(self, **kargs):
        raise NotSupportedError("delayedlists do not support saveCSV()")

    def getChIPSeqTags(self, gene_list, bins=None, bSaveMergedImages=True):
        """
        **Purpose**

        Count the number of chip seq tags that lie under some arrangement of the 'bins'

        NOTE: Although this functionality remains in delayedlist it has
        been deprecated in favour of tracks. This method is astonishingly
        slow and the user is strongly encouraged to look up tracks
        and their usage.

        **Arguments**

        gene_list
            your genelist or genelist-like object

        bins
            an array of bins spannign the range of base pairs you want to bin.
            e.g. [x for x in range(-5000, 5000, 100)] (This is the default)

        bSaveMergedImages
            undocumented.

        **Result**

        returns a taglist list (a descendent of geneList)
        that behaves very much like a genelist.
        However, any sepcial methods from the input gene_list will be lost
        in preference to taglist's methods.
        the chip_seq tag frequency is stored in the "chip_tags" key
        """
        print("Warning: getChIPSeqTags() is slow (depending upon the size of your tag file)")
        print("         press Ctrl+C to interupt.")

        # test gene_list has a valid sorted
        assert gene_list.dataByChr["1"], "List does not have a valid location tag"
        assert gene_list.dataByChr["1"][0]["tss_loc"] , "List does not have a valid tss_loc tag"

        if not bins:
            bins = [x for x in range(-5000, 5000, 100)]

        newl = taglist(bins, gene_list.name) # loses any non-standard methods... :(
        newl.linearData = []

        #def mapAgainstSeqFile(path, seq_filename, listofgenes, bins=[-5000, -4000, -3000, -2000, -1000, 0, 1000, 2000, 3000, 4000, 5000], locCol=0, symmetric=False):
        oh = open(self.fullpath, "rU")

        results = [] # fill a blank array to store the scores.
        for item in range(len(gene_list)+1):
            results.append(qarray("L", [0 for item in bins]))
        flatBin = qarray("L", [0 for item in bins])

        binLeft = bins[0]
        binRight = bins[-1]
        n = 0 # counters
        t = 0
        f = 0

        start = time.time()

        for line in oh:
            s = line.split("\t") # avoid the overhead of csv.reader.
            n += 1
            chr = s[self.format["chr"]].replace("chr", "")
            left = int(s[self.format["left"]])
            right = int(s[self.format["right"]])
            if gene_list.isChromosomeAvailable(chr):
                for data in gene_list.dataByChr[chr]:
                    loc = utils.getLocation(data["tss_loc"])
                    #print data

                    tss = loc["left"] # always for tss_loc

                    mid_tag = left + ((right - left)/2)
                    #print tss, binLeft, mid_tag, binRight
                    if ((tss+binLeft) < mid_tag) and ((tss+binRight) > mid_tag):
                        if data["strand"] == "+":
                            local_tag_loc = (mid_tag + binLeft) - (tss + binLeft)
                        elif data["strand"] == "-":
                            local_tag_loc = (tss + binLeft) - (mid_tag + binLeft)
                        #print local_tag_loc, mid_tag

                        l = bins[0]
                        for i, b in enumerate(bins[1:]):
                            #print "l:%s b:%s" % (l,b)
                            if (local_tag_loc > l and local_tag_loc < b):
                                binSet = results[data["n"]] # get the location in the original list;
                                binSet[i] += 1
                                flatBin[i] += 1
                                f += 1
                                #print "g:", tss, item[1], "t:", mid_tag, "b:", bins[i], "-", bins[i+1]
                                #break
                            l = b

            if n > 1000000: # counter
                t += 1
                print("Info: Done: %s,000,000 - Found: %s tags" % (t, f))
                n = 0
                #break

        oh.close() # finished with the file.

        # load my newl
        for index, item in enumerate(gene_list.linearData):
            c = copy.deepcopy(item) # if you don't copy, it will mangle the original list...
            if "chip_tags" not in c:
                c["chip_tags"] = {}
            c["chip_tags"][self.name] = results[index]
            newl.linearData.append(c)

        # draw some pictures of flatBin.

        end = time.time()
        config.log.info("Time taken: %i mins" % ((end - start) / 60))

        config.log.info("All done, found:", f)
        if MATPLOTLIB_AVAIL:
            self._drawMerged(flatBin, gene_list.path)

        newl._optimiseData()
        return(newl)

    def _drawMerged(self, flatBin, path, window=1):
        """
        (Internal)
        """
        plot.cla()
        if window > 1:
            n = utils.movingAverage(flatBin, window)
        else:
            n = flatBin
        plot.plot(n)
        plot.savefig(os.path.join(path, "ChIP_merge_%s.png" % self.name))
        config.log.info("Saved a merge of the ChIP-peak to: %s" % os.path.join(path, "ChIP_merge_%s.png" % self.name))

    def reset(self):
        """
        **Purpose**
            reset the delayedlist to the 0th element.
            This is a bit of a hack. for most purposes
            delayedlist will be correctly reset. The exception is this case:

            for item in delayedlist:
                ... some code
                break

            for item in delayedlist:
                !!! Error! The list
                continues from where it left off, not from the zeroth
                element as expected.

            Actually, (not tested) I think iterating over a delayedlist
            twice is generally broken, and you should reset.
            However, there is no way for delayedlist to know if the next
            iteration is actually the first iteration of the list
            and not a continuing iteration.

        **Arguments**
            None

        **Results**
            resets the list to the zeroth entry.
        """
        self._optimiseData()

    def getSequences(self, FASTAfilename, genome, loc_key="loc", delta=False):
        """
        **Purpose**
            write out a FASTA file on the fly. 
            This is so that delayed lists can cope with 
            sequence retrieval.
            
            The current arrangement genome.getSequences() will
            copy the delayedlist into memory as it adds sequence to
            each key. That kind of defeats the object of the 
            delayedlists. 
            
            This function will write out a FASTA file without loading
            the delayedlist into memory.
            
            However, this function will always return 'None'
            
        **Arguments**
            FASTAfilename
                name of the FASTA file to save the sequence to.
                
            genome
                a genome object with a bound sequence attached.
                
            loc_key (Optional, default="loc")
                the location key to use for the genome spans
            
            name_key (Optional, default=None)
                a key to use as the name of the fasta entry.
                default is None. This will use an index number for
                each fasta entry, of the form: "<genome.name>_<index>"
            
            delta (Optional, default=False)
                an integer describing the symmetric span around the centre
                of the region to take
                set to False to use only the coords specified in the loc_key
                
        **Returns**
            None
        """
        assert loc_key in self[0], "no '%' loc_key found in this list" % loc_key
        
        oh = open(FASTAfilename, "w")
        
        for index, item in enumerate(self):
            loc = item[loc_key]

            if delta:
                loc = loc.pointify()
                loc = loc.expand(delta)
            if loc["left"] > 0:
                seq = genome.getSequence(loc=loc)
            
            if seq:
                oh.write(">%s_%s\n" % (genome.name, index))
                oh.write("%s\n" % seq)
        oh.close()
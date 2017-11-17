"""
track, part of glbase

Brand new bucket reworking of tracks.

This approach is currently about 30% faster than the old track (10 s versus 15 s). 
And could potentially be expanded for PE reads.

Additionally, this system will be more efficient for larger libraries.

This approach should also win in organisms with only a few chromsomes and large numbers of reads (e.g. Drosophila)

Finally, this approach still has a lot of room for optimisation, whereas track.py is mined out.
For track.py 1/2 the time is spent in SQL and the other half in a very tight get() loop.

Here a lot more time is spent in Python which could potentially be optimised.
Then as seq depths get larger this approach will win.

Unfortunately this approach has one fatal flaw at the moment, in that it takes about 3-10x as long
to initially build the database. In fact on real world data and not simulated data it even seems to fail
eventually (presumably as the CACHE_SIZE is filled, but it's not clear exactly why). 

"""

from __future__ import division

import cPickle, sys, os, struct, math, sqlite3, zlib, time, csv

from progress import progressbar
from errors import AssertionError
from location import location
import utils, config
from data import positive_strand_labels, negative_strand_labels
from draw import draw
import matplotlib.cm as cm
import matplotlib.pyplot as plot
import scipy.stats as stats
from base_track import base_track
import array as pyarray

import numpy
from numpy import array, zeros, set_printoptions, int32, append, linspace, argmax, amax, delete, integer

TRACK_CACHE_SIZE = 10 # number of track segments to cache. (Deprecated?)

"""
This number below may become a significant problem later. 

As the number of reads in each block can potentially be unlimited we can consume all RAM
filling up a single block. In cases like that user would need to reduce the block sizes...
"""
CACHE_SIZE = 100000 # number of blocks per chrom before we start trimming

class track(base_track):
    """
    track definition, used for things like sequence reads across the genome
    
    **Arguments**
        name (string)
            name for the track (defaults to filename)

        filename (string)
            directory location of the track file.
            only respected if dir_name is set.
        
        new (Optional, default=False)
            Use seqToTrk() in preference of this. But if you know what you are 
            doing then this will generate a new (empty) db.
            
        paired_end (Optional, default=False)
            This is a paired-end library.
            
            NOTE: PE support is experimental and incomplete. 
            
        bucket_size (Optional, default=1000)
            The size of buckets to partition the data into
            It is strongly recommended that you do not change this value.

    """
    def __init__(self, name=None, new=False, filename=None, bucket_size=1000, paired_end=False, **kargs):
        base_track.__init__(self, name, new, filename)
        
        self.__bucket_size=bucket_size # should be meta data
        self.paired_end = paired_end # a PE library or not
        self.cache = {}
        
        if new:
            self.__setup_tables(filename)
            self.__next_id = 0 # a running unique id for each read

    def __repr__(self):
        return("glbase.track3 (experimental)")

    def __setup_tables(self, filename):
        """
        No pre-defined file - I want a new database track.
        """
        # If i make it here, then a lot of grunt work already done in base_track
        c = self._connection.cursor()

        c.execute("CREATE TABLE main (chromosome TEXT PRIMARY KEY, seq_reads INT)")
        c.execute("CREATE TABLE bucket_ids (chromosome TEXT, buck_id INT)") # neither can be indexed as neither is unique...
        
        self._connection.commit()
        c.close()

    def __add_chromosome(self, chromosome):
        """
        add a chromosome to the main table.
        add a chromosome table.

        returns True if succesfully created/already present.
        """
        c = self._connection.cursor()
        # check chromosome is not already present.
        if self.__has_chromosome(chromosome):
            return(True)

        c.execute("INSERT INTO main VALUES (?, ?)", (chromosome, 0)) # add chr to master table.

        table_name = "buck_%s" % (str(chromosome),)
        c.execute("CREATE TABLE %s (buck_id INT, set_data TEXT)" % table_name) # and corresponding bucket table
        #c.execute("CREATE INDEX %s_idx ON %s(buck_id)" % (table_name, table_name))
        
        # Add a cache entry:
        self.cache[chromosome] = {}
        
        c.close()
        return(True)

    def __add_bucket(self, chromosome, buck_id):
        """
        add a chromosome to the main table.
        add a chromosome table.

        returns True if succesfully created/already present.
        """
        c = self._connection.cursor()
        # check bucket is not already present.
        if self.__has_bucket(chromosome, buck_id):
            return(True)

        c.execute("INSERT INTO bucket_ids VALUES (?, ?)", (chromosome, buck_id)) # add chr to master table.

        # make the new chromsome table:
        table_name = "buck_%s" % (str(chromosome),)
        c.execute("INSERT INTO %s VALUES (?, ?)" % table_name, (buck_id, sqlite3.Binary(cPickle.dumps(set([]))))) # add in spoof data
        
        #print "added bucket", table_name, buck_id
        # stick it on the cache
        self.cache[chromosome][buck_id] = set([])
        
        c.close()
        return(True)

    def __has_chromosome(self, chromosome):
        """
        do we have that chromosome?
        """
        c = self._connection.cursor()

        c.execute("SELECT chromosome FROM main WHERE chromosome=?", (chromosome, ))

        result = c.fetchone()

        c.close()

        if result:
            return(True)
        return(False)

    def __has_bucket(self, chromosome, bucket_id):
        """
        Test if the bucket exits
        """
        c = self._connection.cursor()

        # Check if we have the master chrom table?
        c.execute("SELECT * FROM bucket_ids WHERE chromosome=? AND buck_id=?", (chromosome, bucket_id))
        result = c.fetchone()
        if result:
            return(True)
        return(False)  

    def add_location(self, loc, pe_loc=None, strand="+", pe_strand=None):
        """
        **Purpose**
            Add a location to the track.
            Increments the score by 'increment' from loc["left"] to
            loc["right"]

        **Arguments**
            loc

            pe_loc
                Support for paired-end reads. Insert the location of the paired end
                here.

            strand
                strand of the read
                
            pe_strand 
                strand of the paired-end.

            increment

        **Returns**
            True, if completes succesfully, or exception.
        """
        if not self.__has_chromosome(loc["chr"]):
            self.__add_chromosome(loc["chr"])

        # convert strand to +, -
        if strand in positive_strand_labels:
            strand = "+"
        elif strand in negative_strand_labels:
            strand = "-"

        c = self._connection.cursor()

        # insert location into new array.
        table_name = "chr_%s" % str(loc["chr"])

        read = (self.__next_id, loc["chr"], loc["left"], loc["right"], strand)
        buckets_reqd = self.__get_buck_ids(loc["left"], loc["right"])

        if buckets_reqd: # This could easily be done faster... patience, problems don't lie here ATM.
            for buck in buckets_reqd:
                if not self.__has_bucket(loc["chr"], buck):
                    self.__add_bucket(loc["chr"], buck)
                
                buck_set = self.__get_buck(loc["chr"], buck)
                buck_set.add(read)
                
                self.__update_buck(loc["chr"], buck, buck_set)

        self.__next_id += 1
        self.__flush_cache()
        c.close()
        return(True)

    def __get_buck(self, chrom, bucket_id):
        """
        get the data within chrom/bucket_id
        
        Does not test if the bucket exists...
        """
        # see if I'm on the cache first
        if bucket_id in self.cache[chrom]:
            return(self.cache[chrom][bucket_id])
        
        # Okay, get from db:
        c = self._connection.cursor()
        table_name = "buck_%s" % str(chrom)
        c.execute("SELECT set_data FROM %s WHERE buck_id=?" % table_name, (bucket_id, ))
        result = c.fetchone()
        c.close()

        #print "got:", chrom, bucket_id

        # stick it on the cache
        self.cache[chrom][bucket_id] = cPickle.loads(str(result[0]))

        if result:
            return(self.cache[chrom][bucket_id])
        raise Exception, "No Block!"

    def __update_buck(self, chrom, buck_id, buck_set):
        """
        Update the bucket with a new set
        """
        #if buck_id in self.cache[chrom]: # Must be on cache as get_buck was run
        self.cache[chrom][buck_id] = buck_set

    def __commit_buck(self, chrom, buck_id, buck_set):
        """
        Update the bucket with a new set
        """
        data_to_put = sqlite3.Binary(cPickle.dumps(buck_set, 1))
        
        c = self._connection.cursor()
        table_name = "buck_%s" % str(chrom)
        c.execute("UPDATE %s SET set_data=? WHERE buck_id=?" % table_name, 
                (data_to_put, buck_id))
        c.close()

    def __flush_cache(self, cache_size=CACHE_SIZE):
        """
        Clean out the cache
        
        cache_size is items per chrom.
        """
        for chrom in self.cache:
            while len(self.cache[chrom]) > cache_size:
                item = self.cache[chrom].popitem() # get a random block
                self.__commit_buck(chrom, item[0], item[1])
    
    def __get_buck_ids(self, left, right):
        """
        return a list of all the bucket_ids this read spans
        """
        left_buck = int((left-1)/self.__bucket_size)*self.__bucket_size
        right_buck = int((right)/self.__bucket_size)*self.__bucket_size
        buckets_reqd = range(left_buck, right_buck+self.__bucket_size, self.__bucket_size) # make sure to get the right spanning and left spanning sites

        return(buckets_reqd)

    def get(self, loc, resolution=None, read_extend=0, kde_smooth=False, 
        view_wid=0, **kargs):
        """
        **Purpose**
            get the data between location 'loc' and return it formatted as
            a nbp resolution array

        **Arguments**
            loc (Required)
                a valid location or string location.

            resolution (Optional, default = 1bp)
                nbp resolution required (you should probably send a float for accurate rendering)

            read_extend (Optional, default = 0)
                extend the read length to 'fill in the peak'
                if the original reads are 36bp, then add ~70bp to give an
                estimated size of the peak.
                If the reads are end-based, then set this to the estimated
                size of the DNA shear.
            
            kde_smooth (Experimental)
                perform kde smoothng on the data, using the integer specified as an option.
                In this case the read_extend acts as a tag shift instead of a read_extend
                Hence set that to half of the expected shear size.

        **Returns**
            an 'numpy.array([0, 1, 2 ... n])' contiginous array
            or a tuple containing two arrays, one for each strand.
        """
        if not isinstance(loc, location):
            loc = location(loc=loc)
        extended_loc = loc.expand(read_extend)

        result = self.get_reads(extended_loc)
        
        if kde_smooth:
            return(self.__kde_smooth(loc, reads, resolution, 0, view_wid, read_extend))

        # make a single array
        #a = zeros(int( (loc["right"]-loc["left"]+resolution)/resolution ), dtype=integer)
        if resolution:
            a = [0] * int( (loc["right"]-loc["left"]+resolution)/resolution ) # Fast list allocation
        else: 
            a = [0] * (loc["right"]-loc["left"]+1) # Fast list allocation
            #a = pyarray.array('l', a)
        # Python lists are much faster for this

        loc_left = loc["left"]
        loc_right = loc["right"]
        len_a = len(a)
        
        for r in result:
            id, chrom, read_left, read_right, strand = r
            
            if strand == "+":
                read_right += (read_extend + 1) # coords are open
            elif strand == "-" :
                read_left -= read_extend
                read_right += 1 # coords are open 
            
            if resolution:
                rel_array_left = (read_left - loc_left) // resolution
                rel_array_right = (read_right - loc_left) // resolution  
            else:
                rel_array_left = (read_left - loc_left)
                rel_array_right = (read_right - loc_left)          
            
            if rel_array_left <= 0:
                rel_array_left = 0
            if rel_array_right > len_a:
                rel_array_right = len_a
            
            #a[rel_array_left:rel_array_right] += 1 # Numpy only
            # The below is a very tiny amount faster
            
            for array_relative_location in xrange(rel_array_left, rel_array_right, 1):
                a[array_relative_location] += 1
            
        return(numpy.array(a))


    def __kde_smooth(self, loc, reads, resolution, bandwidth, view_wid, read_shift=100):
        """
        Internal abstraction for kde smoothing
        
        Returns a new array 
        """
        # Shift the reads
        newr = []
        for r in reads:
            if r[2] in positive_strand_labels:
                newr.append(float((r[0] + read_shift) - loc["left"])) # Must be floats for kde to work
            elif r[2] in negative_strand_labels:
                newr.append(float((r[1] - read_shift) - loc["left"]))
        
        a = linspace(0, loc["right"] - loc["left"], view_wid)
        
        # Hack gaussian_kde()
        def covariance_factor(self):
            return 0.02
        
        kde = stats.gaussian_kde(newr)
        setattr(kde, 'covariance_factor', covariance_factor.__get__(kde, type(kde)))
        kde._compute_covariance()

        kk = kde.evaluate(a) * 1000000 # resacle to get in integer range.
        res = array(kk, dtype=integer)
        
        return(res)
                
    def get_array_chromosome(self, chromosome, strand=None, resolution=1, read_extend=0, **kargs):
        """
        **Purpose**
            get the entire array data for the chromosome

        **Arguments**
            chromosome (Required)
                a number '1', string "1", or X, Y

            strand (Optional, default = None, ie. collect and merge both strands)
                strand, but only valid for stranded tracks
                if "+" return only that strand, if '-' return only the negative
                strand (will recognise several forms of strand, e.g. F/R, +/-

            resolution (default = 1bp)
                nbp resolution required (you should probably send a float for accurate rendering)

            read_extend (Optional, default = 0)
                extend the read length to 'fill in the peak'
                if the original reads are 36bp, then add ~70bp to give an
                estimated size of the peak.
                If the reads are end-based, then set this to the estimated
                size of the DNA shear.

        **Returns**
            an 'numpy.array([0, 1, 2 ... n], dtype=integer)' contiginous array of integers
            or a tuple containing two arrays, one for each strand.
        """
        if strand: 
            raise NotImplementedError, "Eh... strand not supported yet..."

        if not self._c:
            self._c = self._connection.cursor()
            
        table_name = "chr_%s" % chromosome
        self._c.execute("SELECT * FROM %s" % table_name)
        reads = self._c.fetchall()

        # I need to find the right most read to estimate the size of the track array.
        right_most = 0
        for i in reads:
            if right_most < i[1]+read_extend:
                right_most = i[1]+read_extend

        # make an array.
        a = zeros(int(right_most+int(resolution), int(resolution)), dtype=integer)

        for r in reads:
            # read_extend
            if r[2] in positive_strand_labels:
                left = r[0]
                right = r[1] + read_extend + 1
            elif r[2] in negative_strand_labels:
                left = r[0] - read_extend
                right = r[1] + 1

            # check for array wrap arounds:
            if left < 0: left = 0

            for loc in xrange(left, right, 1):
                if resolution <= 1: # force 1bp read # this may be incorrect as it may add the read several times until it increments?
                    # this is the source of the artivacts when resolution < 1.0?
                    a[loc] += 1
                else:
                    a[int(loc/resolution)] += 1

        #print "array_len", len(a)

        return(a)

    def get_reads(self, loc, strand=None):
        """
        **Purpose**
            get all of the sequence reads between location 'loc' and return
            it formatted as a list of tuples: (left, right, strand), seq reads.

        **Arguments**
            loc (Required)
                a valid location or string location.

        **Returns**
            a list containing all of the reads between loc.
        """
        if not isinstance(loc, location):
            loc = location(loc=loc)
            
        if not self._c:
            self._c = self._connection.cursor()

        # Should test for chrom here?

        if len(loc["chr"]) < 30: # small security measure.
            table_name = "chr_%s" % loc["chr"]

        reads = set() 

        # bucket_ids I need to interrogate
        bucks_reqd = self.__get_buck_ids(loc["left"], loc["right"])
        for bucket_id in bucks_reqd:
            table_name = "buck_%s" % loc["chr"]
            #result = self._connection.execute("SELECT set_data FROM %s WHERE buck_id=?" % table_name, (buck,)) 
            
            if self.__has_bucket(loc["chr"], bucket_id):
                self._c.execute("SELECT set_data FROM %s WHERE buck_id=? LIMIT 1" % table_name, (bucket_id, ))
                result = cPickle.loads(str(self._c.fetchone()[0]))
                if result:
                    reads |= result

        result = []

        for r in reads: # Now they are unique, with no overlapping reads.
            # retrieve only the overlapping reads
            # test here to avoid method-call overhead:
            if loc["right"] >= r[2] and loc["left"] <= r[3]:
            #if loc.qcollide(location(chr=loc["chr"], left=r[2], right=r[3])):
                result.append(r)
                #print loc, r, loc.qdistance(location(chr=loc["chr"], left=r[0], right=r[1]))
                
        return(result)

    def get_read_count(self, loc):
        """
        **Purpose**
            get the number of reads within the location specified

        **Arguments**
            loc (Required)
                a valid location or string location.

        **Returns**
            an integer (or 0) containing the number of reads falling within
            the location string.
        """
        if not self._c:
            self._c = self._connection.cursor()

        table_name = "chr_%s" % loc["chr"]

        self._c.execute("SELECT * FROM %s WHERE (?>=left AND ?<=right) OR (?>=left AND ?<=right) OR (left<=? AND right>=?) OR (?<=left AND ?>=right)" % table_name,
            (loc["left"], loc["left"], loc["right"], loc["right"], loc["left"], loc["right"], loc["left"], loc["right"]))

        # this is here in the hope that in the future cursor.rowcount
        # will be correctly supported...
        # at the moment probably provides little or no benefit.

        return(len(self._c.fetchall()))

    def get_chromosome_names(self):
        """
        **Purpose**
            Return a list of all the valid chromosome names in the database
            
        **Arguments**
            None
            
        **Returns**
            A list of strings containing all of the chromosome names in the track
        """
        if not self._c:
            self._c = self._connection.cursor()
        
        self._c.execute("SELECT chromosome FROM main")
        r = [i[0] for i in self._c.fetchall()]
        return(r)

    def finalise(self):
        """
        finalise the database (shrink unused edit space)
        dump useless bits etc.
        You must call this! to finalise the db.
        Or get() will not work!
        This copies the cache onto disk and closes the db.
        """
        # do a commit
        self.__flush_cache(0) # empty
        self.meta_data["num_tags"] = self.__next_id
        self._save_meta_data()
        self._connection.commit()
        self._connection.execute("VACUUM")
        self._connection.execute("ANALYZE")
        self._connection.commit()

    def _debug__print_all_tables(self):
        c = self._connection.cursor()
        c.execute("SELECT * FROM main")
        result = c.fetchall()

        print "Main:"
        for item in result:
            print item

        print "Chr_Tables:"
        for item in result:
            table_name = "chr_%s" % str(item[0])[0] # stop injections.
            print " Table", table_name
            c.execute("SELECT * FROM %s" % table_name)
            chr_table_res = c.fetchall()
            for i in chr_table_res:
                print " ", i
        c.close()

    def pileup(self, genelist=None, key="loc", filename=None, heatmap_filename=None, 
        bin_size=500, window_size=5000, read_extend=200, use_kde=False, simple_cleanup=False,
        only_do=False, stranded=False, respect_strand=True, raw_tag_filename=None,
        **kargs): 
        """
        **Purpose**
            draw cumulative 'pileups' of the tag density in the track based on a genelist
            containing a "loc" tag
        
        **Arguments**
            genelist (Required)
                A genelist-like object containing a "loc"-like key tag
            
            key (Optional, default="loc")
                the key to use for the locations. Must contain valid location data:
                chr1:1000-1011 etc. draw_pileup() will use the centre of the location if it is a 
                span.
                
            filename (Required)
                the name of the image file to save the pileup graph to.
                
            heatmap_filename (Optional)
                draw a heatmap, where each row is a single gene. And red is the tag density.
                Please see heatmap() for a better, more customizable means of drawing 
                heatmaps of tags
            
            raw_tag_filename (Optional)
                Save a tsv file that contains the heatmap values for each row of the genelist.
            
            bin_size (Optional, default=500)
                bin_size to use for the heatmap
                
            window_size (Optional, default=5000)
                The window size +- around the centre of the peak to build the tag density map
                from
                
            read_extend (Optional, default=200)
                extend the read x bp either 5' or 3' depending upon the strand of the read.
                If use_kde is true then this will be the 'tag shift' value instead.
                
            use_kde (Optional)
                Use KDE versions of the tracks instead (Experimental)
                
            simple_cleanup (False)
                remove rows from the pileup where sum(tag_count) < simple_cleanup
                
            stranded (Optional, default=False)
                build a stranded pileup, with + reads in blue and - reads in red
                
            respect_strand (Optional, default=True)
                If available, respect the orientation of the strand from the genelist.
                This is useful if you are, say, using the TSS's and want to maintain the
                orientation with respect to the transcription direction.
                            
        **Returns**
            If succesful returns a list of lists containing the a single entry for each
            entry in the original genelist (in the same order as the genelist).
        """
        assert filename, "you must specify a filename"
        assert genelist, "you must specify a genelist"
        assert key in genelist.linearData[0], "the genelist appears to be lacking a loc key"
        
        if stranded:
            return(self.__draw_pileup_stranded(genelist, filename, window_size, **kargs))
        else:
            return(self.__draw_pileup_normal(genelist, key, filename, heatmap_filename, 
                bin_size, window_size, read_extend, use_kde, simple_cleanup, 
                only_do=only_do, raw_tag_filename=raw_tag_filename, respect_strand=respect_strand, **kargs))
                
    def __draw_pileup_stranded(self, genelist=None, filename=None, bandwidth=300, **kargs):
        """
        **Purpose**
            Build a histogram density plot of the paired reads. 
            This will estimate the approximate observed shear size in your chip-seq data.
            It pairs the data together for all pairs of all reads within a specified bandwidth
            then outputs a plot of the resulting histogram.
            
        **Arguments**
            filename
                the filename to save the image(s) to.
                
            genelist (Required)
                some sort of genelistlike object with a 'loc' key containing genomic coordinates
            
            bandwidth (Optional, default=300)
                area around the centre of the peak to build the cumulative distributions.
        
        **Returns**
            None
            and some image files in <base_filename>.<draw_mode>
        """      
        # step along each chromosome. Quit if there are no reads in the window
        
        if not self._c:
            self._c = self._connection.cursor()
        
        hist = {"+": zeros(bandwidth+1), "-": zeros(bandwidth+1)}
        
        p = progressbar(len(genelist))
        for n, read in enumerate(genelist):
            loc = read["loc"].pointify().expand(bandwidth//2)
            for s in self.get_reads(loc): # get reads returns all overlapping reads. I need to trim 
                # the edges to stop negative array positions
                
                loc_left = s[0] - loc["left"]
                loc_right = s[1] - loc["left"]

                if s[2] == "+" and loc_left > 0:
                    loc_left = s[0] - loc["left"]
                    hist["+"][loc_left] += 1
                elif s[2] == "-" and loc_right < bandwidth-1:
                    loc_right = s[1] - loc["left"]
                    hist["-"][loc_right] += 1

            p.update(n)
            
        if not self._draw:
            self._draw = draw(self)
            
        # now plot:
        fig = self._draw.getfigure()
        ax = fig.add_subplot(111)
        ax.plot(hist["+"], color="blue")
        ax.plot(hist["-"], color="red")
        
        max_left = argmax(hist["+"])
        max_right = argmax(hist["-"])
        
        realfilename = self._draw._saveFigure(fig, filename)
        config.log.info("Saved shear_size_pileup '%s'" % realfilename)
        return(hist)
    
    def __draw_pileup_normal(self, genelist=None, key="loc", filename=None, heatmap_filename=None, 
        bin_size=500, window_size=5000, read_extend=200, use_kde=False, simple_cleanup=False,
        only_do=False, respect_strand=True, raw_tag_filename=None, mask_zero=False,
        **kargs): 
        """
        The normal pileup views
        """
        # See if there is a proper stransd key in there somewhere:
        if respect_strand and (not "strand" in genelist.linearData[0]):
            config.log.warning("I could not find the 'strand' key, setting respect_strand to False")
            respect_strand = False
            
        n = 0
        h = 0
        pileup = None
        sampled = False
        binned_data = None
        setup_bins = False
        
        if not simple_cleanup:
            simple_cleanup = 0

        p = progressbar(len(genelist))
        for i, g in enumerate(genelist):
            l = g[key].pointify().expand(window_size)
            if not use_kde:
                a = self.get(l, read_extend=read_extend) # read_extend=shear size
            if use_kde:
                a = self.get(l, read_extend=read_extend, kde_smooth=True, view_wid=window_size) # read_extend=tag shift
                
            if respect_strand:
                # positive strand is always correct, so I leave as is.
                # For the reverse strand all I have to do is flip the array.
                if g["strand"] in negative_strand_labels:
                    a = a[::-1]

            if sum(a) > simple_cleanup: # Only keep if tag count is > simple_cleanup
                if not sampled: # numpy __nonzero__ retardedness
                    pileup = a
                    sampled = True
                else:
                    pileup = pileup + a

                if heatmap_filename or raw_tag_filename:
                    if not setup_bins:
                        binned_data = array([utils.bin_data(a, bin_size)])
                        setup_bins = True
                    else:
                        #print binned_data, [utils.bin_data(a, bin_size)]
                        binned_data = append(binned_data, [utils.bin_data(a, bin_size)], axis=0)
                
            if only_do and n > only_do: 
                #print only_do, n
                break
            p.update(i)

        if not self._draw:
            self._draw = draw()

        # matplotlib pileup graph
        fig = self._draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        ax.plot(pileup)
        real_filename = self._draw.savefigure(fig, filename)
        config.log.info("Saved pileup tag density to '%s'" % filename)
               
        other_args = {}
        if "vmax" in kargs:
            other_args["vmax"] = kargs["vmax"]
        if "vmin" in kargs:
            other_args["vmin"] = kargs["vmin"]

        # spin this off into a .heatmap() method?
        if not setup_bins:
            config.log.info("No tags found, heatmap_filename or raw_tag_filename not used")
        else:
            if heatmap_filename or raw_tag_filename:
                print binned_data
                binned_data = numpy.delete(binned_data, numpy.s_[-1:], 1) # kill the rightmost empty col.
                
                if raw_tag_filename:
                    print binned_data
                    binned_data = numpy.delete(binned_data, numpy.s_[-1:], 1)
                    oh = open(raw_tag_filename, "w")
                    writer = csv.writer(oh, dialect=csv.excel_tab)
                    print binned_data
                    writer.writerows(binned_data)
                    oh.close()
                    config.log.info("saved raw_tag_file to '%s'" % raw_tag_filename)
                    
                if heatmap_filename:
                    if other_args: 
                        real_filename = self._draw._simple_heatmap(data=binned_data, filename=heatmap_filename, 
                            dpi=150, figsize=(6, 24), aspect="long", **other_args)
                    else:
                        real_filename = self._draw._simple_heatmap(data=binned_data, filename=heatmap_filename, 
                            dpi=150, figsize=(6, 24), aspect="long")
                    config.log.info("saved heatmap to '%s'" % real_filename)

        return(pileup)

    def heatmap(self, filename=None, genelist=None, distance=1000, read_extend=200, log=2, 
        bins=20, sort_by_intensity=True, **kargs):
        """
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
        assert filename, "must specifty a filename"
        assert genelist, "must provide a genelist"
        assert "loc" in genelist.keys(), "appears genelist has no 'loc' key"
        assert "left" in genelist.linearData[0]["loc"].keys(), "appears the loc key data is malformed"
        assert log in ("e", math.e, 2, 10), "this 'log' base not supported"
        
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
        mag_tab = []
        for index, row in enumerate(table):
            mag_tab.append({"n": index, "sum": row.max()})
        
        if sort_by_intensity:
            mag_tab = sorted(mag_tab, key=itemgetter("sum"))
        
        data = numpy.array(table)+1
        
        newt = []
        for item in mag_tab:
            newt.append(data[item["n"],])
        data = numpy.array(newt)
        
        if log:
            if log == "e" or log == math.e:
                data = numpy.log(data)-1
            elif log == 2:
                data = numpy.log2(data)-1
            elif log == 10:
                data = numpy.log10(data)-1
               
        # draw heatmap
        
        if not self._draw:
            self._draw = draw()
            
        filename = self._draw.heatmap2(data=data, filename=filename, bracket=[1, data.max()], **kargs)
        
        config.log.info("Saved pileup tag density to '%s'" % filename)
        return({"data": data})

def build():
    """
    temp func for profiling
    """
    t = track(filename="testnew.trk3", name="test", new=True)
    for n in xrange(0, 10000): #test = 10000
        l = random.randint(0, 100000)
        t.add_location(location(chr="1", left=l, right=l+35), strand="+")
    t.finalise()

if __name__ == "__main__":
    """
    track .py:   1.0 s and 14 s
    tracknew.py: 2.5s and 36 s
    track3.py:   2.9s and 9.0s !
    """

    import cProfile, pstats    
    import random, time
    from location import location
    from genelist import genelist
    
    s = time.time()
    print "Building...",
    cProfile.run("build()", "profile.pro")
    p = pstats.Stats("profile.pro")
    p.strip_dirs().sort_stats("time").print_stats()

    e = time.time()
    print e-s, "s"
    
    s = time.time()
    print "Fake bed..."
    # fake a bed
    newb = []
    for n in xrange(0, 1000):
        l = random.randint(1000, 100000) # 1000 is for the window size. -ve locs are real bad.
        newb.append({"loc": location(chr="1", left=l, right=l+200), "strand": "+"})
    bed = genelist()
    bed.load_list(newb)
    e = time.time()
    print e-s, "s"
    
    t = track(filename="testnew.trk3")
    
    print "Pileup..."
    #cProfile.run("t.pileup(genelist=bed, filename='test.png', bin_size=10, window_size=1000)", "profile.pro")
    p = pstats.Stats("profile.pro")
    p.strip_dirs().sort_stats("time").print_stats()

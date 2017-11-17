"""
track, part of glbase

"""



import pickle, sys, os, struct, math, sqlite3, zlib, time, csv, zlib

from operator import itemgetter

from .progress import progressbar
from .errors import AssertionError
from .location import location
from . import genelist as Genelist
from . import utils, config
from .data import positive_strand_labels, negative_strand_labels
from .draw import draw
import matplotlib.cm as cm
import matplotlib.pyplot as plot
import scipy.stats as stats
from scipy.stats import pearsonr
from .base_track import base_track

import numpy
from numpy import array, zeros, set_printoptions, int32, append, linspace, argmax, amax, delete

TRACK_CACHE_SIZE = 10 # number of track segments to cache.

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

        norm_factor (Optional, default = 1.0)
            An optional normalization factor. Data is multiplied by this number before display.
            Can only be specified at creation time, cannot be modified later.
            
        mem_cache (Optional, default=False)
            Instead of doing the whole thing from disk, first cache the DB in memory for
            super fast access. Just make sure you have enough memory!
        
        pre_build (Optional, default=[0, 100, 200])
            prebuild genome arrays for various read_extend parameters for fast lookups.
            By default glbase builds indexes for 100 and 200 bp read_extends and unextended 
            (i.e. the complete frag size, useful for e.g. paired-end data). Raw frag sizes are returned
            when read_extend and pre_build are 0
            
    """
    def __init__(self, name=None, new=False, filename=None, norm_factor=None, mem_cache=False, pre_build=[0, 100, 200], **kargs):
        base_track.__init__(self, name, new, filename, norm_factor, mem_cache)
                                
        if new:
            if norm_factor is None:
                norm_factor = 1.0
            self.meta_data['norm_factor'] = str(norm_factor)
            self.pre_build = pre_build
            self.meta_data['pre_build'] = pre_build
            self.__setup_tables(filename)
            config.log.info('Generating new track')
        else:
            assert not norm_factor, 'norm_factor can only be modified at creation time'
            # if not new, get pre_build from the metadatum
        
        try:
            self.pre_build =  self.meta_data['pre_build']
        except KeyError:
            raise AssertionError('meta data not found in trk file, this suggests the trk file is incomplete, please check your trk file and regenerate if required')
            
        self.norm_factor = float(self.meta_data['norm_factor'])
        if self.norm_factor != 1.0:
            config.log.info('track: Using norm_factor=%.3f' % self.norm_factor)
            #print norm_factor, self.norm_factor, str(norm_factor) != str(self.norm_factor)
            #if str(norm_factor) != str(self.norm_factor):
            #    config.log.error('the norm_factor supplied here does not match the norm_factor used during the creation of the track!')
            #    raise AssertionError, 'norm_factor != norm_factor (%.2f != %.2f) stored in the trk file.' % (norm_factor, self.norm_factor)
        self.__warned_about_zlib = False # To deprecate

    def __repr__(self):
        return("glbase.track")

    # I use these ones here as tracks prefer numpy arrays
    # LATER: port the flat_tracks to use numpy arrays
    def _format_data(self, data):
        """
        array('i', []) --> whatever it's stored as in db
        """
        return(sqlite3.Binary(zlib.compress(data.dumps(), 1)))

    def _unformat_data(self, data):
        """
        whatever stored as in db --> array('i', [])
        """
        #print "ret:",[d for d in data], ":"
        try:
            a = numpy.loads(zlib.decompress(data))
        except pickle.UnpicklingError:
            a = numpy.loads(zlib.decompress(data))
            if not self.__warned_about_zlib:
                config.log.warning('Tracks are no no longer Zlib compressed by default. Tracks will need to be regenerated')
                config.log.warning('Benefit is that they are now faster!')
                config.log.warning('In future versions of glbase this warning will be removed and an error will be produced')
                self.__warned_about_zlib = True
        return(a)

    def __setup_tables(self, filename):
        """
        No pre-defined file - I want a new database track.
        """
        # If i make it here, then a lot of grunt work already done in base_track
        c = self._connection.cursor()

        c.execute("CREATE TABLE main (chromosome TEXT PRIMARY KEY, seq_reads INT)")
        
        c.execute("CREATE TABLE pre_build (chromosome TEXT, read_extend TEXT, array_blob BLOB)")

        self._connection.commit()
        c.close()

    def finalise(self):
        c = self._connection.cursor()
    
        # Do the prebuilds:
        list_of_all_chroms_in_db = self.get_chromosome_names()
        for read_extend in self.pre_build:
            for chrom in list_of_all_chroms_in_db:
                a = self.get_array_chromosome(chrom, read_extend=read_extend, _silent=True) # Force a cache miss
                c.execute('INSERT INTO pre_build VALUES (?, ?, ?)', (chrom, read_extend, self._format_data(a)))
                config.log.info('Cached chrom=%s, read_extend=%s' % (chrom, str(read_extend)))
        #print self.meta_data
        
        base_track.finalise(self)
        
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

        # make the new chromsome table:
        table_name = "chr_%s" % str(chromosome)
        c.execute("CREATE TABLE %s (left INT, right INT, strand TEXT)" % (table_name, ))
        #c.execute("CREATE INDEX %s_com_idx ON %s(left, right)" % (table_name, table_name))
        #c.execute("CREATE INDEX %s_lef_idx ON %s(left)" % (table_name, table_name))
        #c.execute("CREATE INDEX %s_rig_idx ON %s(right)" % (table_name, table_name))

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

    def add_location(self, loc, strand="+", increment=1):
        """
        **Purpose**
            Add a location to the track.
            Increments the score by 'increment' from loc["left"] to
            loc["right"]

        **Arguments**
            loc

            strand

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

        if not self._c:
            self._c = self._connection.cursor()

        # insert location into new array.
        table_name = "chr_%s" % str(loc["chr"])

        # get the old number of seq_reads
        self._c.execute("SELECT seq_reads FROM main WHERE chromosome=?", (loc["chr"], ))
        current_seq_reads = self._c.fetchone()[0] # always returns a tuple

        self._c.execute("UPDATE main SET seq_reads=? WHERE chromosome=?", (current_seq_reads+1, loc["chr"]))

        # add the location to the seq table:
        self._c.execute("INSERT INTO %s VALUES (?, ?, ?)" % table_name, (loc["left"], loc["right"], strand))

        #c.close()

    def get(self, loc, resolution=1, read_extend=0, kde_smooth=False, 
        view_wid=0, strand=False, **kargs):
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

            strand (Optional, default=False)
                collect only reads on the specified strand. (track will use read strand 
                information intelligently, if present).
    
        **Returns**
            an 'numpy.array([0, 1, 2 ... n])' contiginous array
            or a tuple containing two arrays, one for each strand.
        """
        if not isinstance(loc, location):
            loc = location(loc=loc)
        extended_loc = loc.expand(read_extend)
        
        result = self.get_reads(extended_loc, strand=strand)
        
        if kde_smooth:
            return(self.__kde_smooth(loc, result, resolution, 0, view_wid, read_extend))

        loc_left = loc["left"]
        loc_right = loc["right"]

        # make a single array
        a = [0] * int( (loc_right-loc_left+resolution)/resolution ) # Fast list allocation
        # Python lists are much faster for this than numpy or array

        len_a = len(a)
        
        for r in result:
            read_left, read_right, strand = r
            if strand == "+":
                read_right += (read_extend + 1) # coords are open
            elif strand == "-" :
                read_left -= read_extend
                read_right += 1 # coords are open 
            
            rel_array_left = int((read_left - loc_left) // resolution)
            rel_array_right = int((read_right - loc_left) // resolution)        
            
            if rel_array_left <= 0:
                rel_array_left = 0
            if rel_array_right > len_a:
                rel_array_right = len_a
            
            for array_relative_location in range(rel_array_left, rel_array_right, 1):
                a[array_relative_location] += 1
            
            #a[rel_array_left:rel_array_right] += 1 # Why is this slower than the for loop? # could be done with num_expr?
            
            #[a[array_relative_location].__add__(1) for array_relative_location in xrange(rel_array_left, rel_array_right, 1)] # just returns the exact item, a is unaffected?
        return(numpy.array(a)*self.norm_factor)

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
        res = array(kk)
        
        return(res)
    
    def get_all_reads_on_chromosome(self, chrom, strand=None):
        """
        **Purpose**
            Get all of the reads on chromosomes.
            
        **Arguments**
            chromosome (Required)  
                the name of the chromosome to pull from.
                
            strand (Optional)
                selct + or - strands and only collect those.
                
        **Returns**
            An iterator to collect the data from.
            
            You must process the data in some sort of for loop:
            for read in trk.get_all_reads_on_chromosome("1")
        """
        assert chrom, "You must provide a chromosome"
        assert chrom in self.get_chromosome_names(), "chromsome '%s' not found in this track" % chromosome
        
        if not self._c:
            self._c = self._connection.cursor()
            
        if len(chrom) < 30: # small security check
            table_name = "chr_%s" % chrom
        
        if strand:
            result = self._c.execute("SELECT * FROM %s WHERE strand=?" % table_name, strand)
        else:
            result = self._c.execute("SELECT * FROM %s" % table_name)
            
        #reads = self._c.fetchall()
        return(result)                

    def get_array_chromosome(self, chrom, read_extend=0, strand=None, resolution=1, _silent=False, **kargs):
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
                
                NOT SUPPORTED AT THIS TIME

            resolution (default = 1bp)
                nbp resolution required (you should probably send a float for accurate rendering)

            read_extend (Optional, default = 0)
                extend the read length to 'fill in the peak'
                if the original reads are 36bp, then add ~70bp to give an
                estimated size of the peak.
                If the reads are end-based, then set this to the estimated
                size of the DNA shear.
                
                Use a read_extend of 0 to return the actual frags.

        **Returns**
            an 'numpy.array([0, 1, 2 ... n], dtype=integer)' contiginous array of integers
            or a tuple containing two arrays, one for each strand.
        """
        if strand: 
            raise NotImplementedError("Eh... strand not supported yet...")
        
        c = self._connection.cursor()

        # Find out if we already have this array:
        c.execute("SELECT chromosome FROM pre_build WHERE (chromosome=? AND read_extend=?)", (chrom, read_extend))
        result = c.fetchone()

        if not result:
            if not _silent: # It purposely misses the cache when building the track
                config.log.warning('Cache miss on chromosome=%s, read_extend=%s' % (chrom, read_extend))
            return(self.__cache_miss_get_array_chromosome(chromosome=chrom, read_extend=read_extend)) # Don't have... :(
            # The above is already * self.norm_factor

        # have a changed copy:
        c.execute("SELECT array_blob FROM pre_build WHERE (chromosome=? AND read_extend=?)", (chrom, read_extend))
        return(self._unformat_data(c.fetchone()[0])) # DO NOT multiply the below result by norm_factor
        # The prebuilt causes a cache _miss and a return above. norm_factor is applied at the end of 
        # __cache_miss_get_array_chromosome() and the result is stored in the array_blob
                
    def __cache_miss_get_array_chromosome(self, chromosome, strand=None, resolution=1, read_extend=0, **kargs):
        # Generate the chromosome array for a cache miss
        if not self._c:
            self._c = self._connection.cursor()
            
        table_name = "chr_%s" % chromosome
        self._c.execute("SELECT * FROM %s" % table_name)
        reads = sorted(self._c.fetchall(), key=itemgetter(0)) # shouldn't this be 1?
        # I need to find the right most read to estimate the size of the track array.
        right_most = reads[-1][1]+read_extend+1000 # Add a large enough pad onto the end, particularly for weird data with ultra long reads

        # make an array.
        a = [0] * int( (right_most+resolution)/resolution ) # Fast list allocation
        # Python lists are much faster for this than numpy or array
        len_a = len(a) # == right_most+resolution

        for r in reads:
            read_left, read_right, strand = r
            if read_extend > 0: # if == 0 then use the total frag size
                if strand == "+":
                    read_right += (read_extend + 1) # coords are open
                elif strand == "-" :
                    read_left -= read_extend
                    read_right += 1 # coords are open 
            
            rel_array_left = read_left
            rel_array_right = read_right
            if resolution != 1:
                rel_array_left = int(read_left // resolution)
                rel_array_right = int(read_right // resolution)        
            
            if rel_array_left <= 0:
                rel_array_left = 0
            if rel_array_right > len_a: # This should never happen?
                rel_array_right = len_a

            # fold up to 1 liner           
            # This one liner does not work for some reason.
            #[a[array_relative_location] + 1 for array_relative_location in xrange(rel_array_left, rel_array_right, 1)]
            for array_relative_location in range(rel_array_left, rel_array_right, 1):
                a[array_relative_location] += 1

        #print "array_len", len(a)

        return(numpy.array(a)*self.norm_factor)

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
        
        if self.gl_mem_cache: # Use the mem cache if available
            # work out which of the buckets is required:
            left_buck = int((loc["left"]-1-delta)/config.bucket_size)*config.bucket_size
            right_buck = int((loc["right"]+delta)/config.bucket_size)*config.bucket_size
            buckets_reqd = list(range(left_buck, right_buck+config.bucket_size, config.bucket_size)) # make sure to get the right spanning and left spanning sites

            # get the ids reqd.                
            loc_ids = set()
            if buckets_reqd:
                for buck in buckets_reqd:
                    if buck in self.gl_mem_cache.buckets[loc["chr"]]:
                        loc_ids.update(self.gl_mem_cache.buckets[loc["chr"]][buck]) # set = unique ids
       
            # loc_ids is a set, and has no order. 
            #print loc_ids
            for index in loc_ids:
                #print loc.qcollide(self.linearData[index]["loc"]), loc, self.linearData[index]["loc"]
                if loc.qcollide(self.linearData[index]["loc"]):
                    # result expected in form :read_left, read_right, strand
                    result.append((self.linearData[index]["loc"]['left'], self.linearData[index]["loc"]['right'], self.linearData[index]["strand"])) 

        #if len(loc["chr"]) < 30: # small security measure.
        table_name = "chr_%s" % loc["chr"]

        #result = self._connection.execute("SELECT * FROM %s WHERE (?>=left AND ?<=right) OR (?>=left AND ?<=right) OR (left<=? AND right>=?) OR (?<=left AND ?>=right)" % table_name,
        #    (loc["left"], loc["left"], loc["right"], loc["right"], loc["left"], loc["right"], loc["left"], loc["right"]))
        
        # This is the code used in location.collide():
        #self["right"] >= loc["left"] and self["left"] <= loc["right"]
        result = self._connection.execute("SELECT left, right, strand FROM %s WHERE (right >= ? AND left <= ?)" % table_name,
            (loc["left"], loc["right"]))
    
        #result = None       
        result = result.fetchall() # safer for empty lists and reusing the cursor
        
        if result and strand: # sort out only this strand
            if strand in positive_strand_labels:
                strand_to_get = positive_strand_labels
            elif strand in negative_strand_labels:
                strand_to_get = negative_strand_labels
            
            newl = []
            for r in result:
                if r[2] in strand_to_get:
                    newl.append(r)
            result = newl                    

        return(result)

    def get_read_count(self, loc):
        """
        **Purpose**
            get the number of reads within the location specified

        **Arguments**
            loc (Required)
                a valid location or string location.

        **Returns**
            an float (or 0.0) containing the number of reads falling within
            the location string.
        """
        if not self._c:
            self._c = self._connection.cursor()
            
        if not isinstance(loc, location):
            loc = location(loc=loc)

        table_name = "chr_%s" % loc["chr"]

        self._c.execute("SELECT left, right, strand FROM %s WHERE (right >= ? AND left <= ?)" % table_name,
            (loc["left"], loc["right"]))
            
        return(len(self._c.fetchall())*self.norm_factor)

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
        return(set(r))

    def get_numreads_on_chromosome(self, name):
        """
        **Purpose**
            Return the number of reads on chromosme name
            
        **Arguments**
            name (Required)
                get the number of reads on chromsome 'name'
            
        **Returns**
            An integer containing the number of reads
        """
        if not self._c:
            self._c = self._connection.cursor()

        self._c.execute("SELECT chromosome, seq_reads FROM main WHERE chromosome=?", (str(name), ))
        r = self._c.fetchone()
        return(r[1])

    def get_total_num_reads(self):
        """
        **Purpose**
            Return the number total number of reads for this track.
            
        **Arguments**
            None
            
        **Returns**
            An integer containing the number of reads
        """
        if not self._c:
            self._c = self._connection.cursor()

        self._c.execute("SELECT chromosome, seq_reads FROM main")
        r = [int(i[1]) for i in self._c.fetchall()]
        
        return(sum(r))

    def _debug__print_all_tables(self):
        c = self._connection.cursor()
        c.execute("SELECT * FROM main")
        result = c.fetchall()

        print("Main:")
        for item in result:
            print(item)

        print("Chr_Tables:")
        for item in result:
            table_name = "chr_%s" % str(item[0])[0] # stop injections.
            print(" Table", table_name)
            c.execute("SELECT * FROM %s" % table_name)
            chr_table_res = c.fetchall()
            for i in chr_table_res:
                print(" ", i)
        c.close()

    def saveBedGraph(self, filename, bin_size=100, read_extend=0):
        '''
        **Purpose**
            Save the track as a BedGraph
            Will take into account the norm_factor if available
        
        **Arguments**
            filename (Required)
                filename to save BG to
                
            bin_size (Optional, default=100)
                the size for each bin (resolution) of the BedGraph
                
            read_extend (Optional, default=0)
                extend the reads on the 3' end by this many base pairs.
                set to 0 if your reads are the total fragments (e.g. paired-end data).
        
        **Returns**
            None
        '''
        assert filename, 'You must provide a filename'
        
        oh = open(filename, 'w')
        for chrom in sorted(self.get_chromosome_names()):
            this_chrom = self.get_array_chromosome(chrom, read_extend=read_extend)
            config.log.info("Doing Chromosome '%s'" % chrom)
            min_position = 0 # Could guess, but the below code will trim the padding zeros
            max_position = len(this_chrom)
            for l in range(min_position, max_position, bin_size):
                value = numpy.mean(this_chrom[l:l+bin_size])
                if value > 0.0: # If zero then it is okay to skip the loc.
                    oh.write('chr%s\t%s\t%s\t%s\n' % (chrom, l, l+bin_size, numpy.mean(this_chrom[l:l+bin_size]))) # data is already norm_factor corrected  
        oh.close()
        config.log.info("saveBedGraph: Saved '%s'" % filename)
        return(None)
    
    def pileup(self, genelist=None, key="loc", filename=None, heatmap_filename=None, 
        bin_size=500, window_size=5000, read_extend=200, use_kde=False, simple_cleanup=False,
        only_do=False, stranded=False, respect_strand=True, raw_tag_filename=None,
        norm_by_read_count=False, pointify=True,
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
            
            normalize (Optional, default=True)
                IMPORTANT
                
                If you are using the norm_factor system, you MUST set this to False!
            
            bin_size (Optional, default=500)
                bin_size to use for the heatmap
                
            pointify (Optional, default=True)
                convert ythe genomic locations in 'genelist' to a single point
                (Usually used in combination with window_size).
            
            window_size (Optional, default=5000)
                The window size +- around the centre of the peak to build the tag density map
                from
                
            read_extend (Optional, default=200)
                extend the read x bp either 5' or 3' depending upon the strand of the read.
                If use_kde is true then this will be the 'tag shift' value instead.
                
            use_kde (Optional)
                Use KDE versions of the tracks instead (Experimental)
                
            simple_cleanup (False)
                remove rows from the pileup that have < simple_cleanup tag counts
                
            stranded (Optional, default=False)
                build a stranded pileup, with + reads in blue and - reads in red
                
            respect_strand (Optional, default=True)
                If available, respect the orientation of the strand from the genelist.
                This is useful if you are, say, using the TSS's and want to maintain the
                orientation with respect to the transcription direction.
                
            norm_by_read_count (Optional, default=False)
                If you are not using a norm_factor for this library then you probably want to set this to True. 
                It will divide the resulting number of reads by the total number of reads, 
                i.e. it will account for differences in library sizes. 
                            
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
            return(self.__draw_pileup_normal(genelist=genelist, key=key, filename=filename, 
                heatmap_filename=heatmap_filename, 
                bin_size=bin_size, window_size=window_size, read_extend=read_extend, use_kde=use_kde, 
                simple_cleanup=simple_cleanup, pointify=pointify,
                norm_by_read_count=norm_by_read_count, 
                only_do=only_do, raw_tag_filename=raw_tag_filename, respect_strand=respect_strand, 
                **kargs))
    '''
    def __draw_pileup_normal(self, genelist=None, key="loc", filename=None, heatmap_filename=None, 
        bin_size=500, window_size=5000, read_extend=200, use_kde=False, simple_cleanup=False,
        only_do=False, respect_strand=True, raw_tag_filename=None, mask_zero=False, 
        norm_by_read_count=False, normalize=False,
        **kargs): 
    '''         
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
        
        # get a sorted list of all the locs I am going to use.
        all_locs = genelist['loc']
        all_locs.sort()        
                
        p = progressbar(len(genelist))
        for n, read in enumerate(all_locs):
            loc = read.pointify().expand(bandwidth//2)
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
        norm_by_read_count=False, pointify=True,
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
        binned_data = None
        setup_bins = False
        number_of_tags_in_library = 1.0 # For code laziness
        if norm_by_read_count:
            number_of_tags_in_library = self.get_total_num_reads()/float(1e6) # div 1e6 for number nicness.

        # get a sorted list of all the locs I am going to use.
        # I need to do this:
        gl_sorted = genelist.deepcopy()
        gl_sorted.sort(key)
        all_locs = gl_sorted[key]      
        strands = ['+'] * len(all_locs)
        if respect_strand:
            strands = gl_sorted['strand']
        
        curr_cached_chrom = None
        cached_chrom = None
        
        p = progressbar(len(genelist))
        for i, read in enumerate(zip(all_locs, strands)):
            l = read[0]
            if pointify:
                l = l.pointify()
            if window_size:
                l = l.expand(window_size)
            
            # I can dispose and free memory as the locations are now sorted:
            # See if the read_extend is already in the cache:
            if l['chr'] != curr_cached_chrom:
                curr_cached_chrom = l['chr']
                cached_chrom = self.get_array_chromosome(l['chr'], read_extend) # auto-deals with cahce issues
            
            # UseKDE is now unimplemented:
            '''
            if not use_kde:
                a = self.get(l, read_extend=read_extend) # read_extend=shear size
            else:
                a = self.get(l, read_extend=read_extend, kde_smooth=True, view_wid=window_size) # read_extend=tag shift
            '''
            assert l['left'] >= 0, 'left genome coordinate is less than zero "%s"' % l
                        
            if l['right'] > cached_chrom.shape[0]:
                # Trouble, need to fill in the part of the array with zeros
                # Possible error here it l['left'] is also off the array?
                expected_width = l['right'] - l['left']
                actual_width = cached_chrom.shape[0] - l['left']
                a = cached_chrom[l['left']:cached_chrom.shape[0]] # stop wrap around
                a = numpy.pad(a, (0,expected_width-actual_width), mode='constant')
                #print a, a.shape
                config.log.warning('Asked for part of the chomosome outside of the array "%s", skipping this loc' % str(l))
                continue
            else:
                a = cached_chrom[l['left']:l['right']]
            
            #print read
            
            if respect_strand:
                # positive strand is always correct, so I leave as is.
                # For the reverse strand all I have to do is flip the array.
                if read[1] in negative_strand_labels:
                    a = a[::-1]
            
            # It is possible to ask for an array longer than the length of the array
            # NOFIX?
            
            #print l, a.shape
            
            if pileup is None: # numpy __nonzero__ retardedness
                pileup = a
                binned_data = array([utils.bin_data(a, bin_size)])
                setup_bins = True
            else:        
                if sum(a) > simple_cleanup: # Only keep if tag count is > simple_cleanup     
                    pileup = pileup + a

                if heatmap_filename or raw_tag_filename: 
                    if setup_bins:
                        #print binned_data, [utils.bin_data(a, bin_size)]
                        binned_data = append(binned_data, [utils.bin_data(a, bin_size)], axis=0)
                
            if only_do and n > only_do: 
                #print only_do, n
                break
            p.update(i)
       
        if not self._draw:
            self._draw = draw()
                
        if pileup is None: # numpy iszero testing:
            raise AssertionError('no data found, either the bed is empty, has no regions or the trk file is empty')
        
        if norm_by_read_count:
            config.log.info('Normalized by read count')
            pileup /= float(number_of_tags_in_library) # This one should not be used if norm_factor is also used
        # This one SHOULD be used, even if norm_factor is True
        pileup /= float(len(genelist)) # convert it back to a relative tag density.

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
        if heatmap_filename or raw_tag_filename:
            binned_data = numpy.delete(binned_data, numpy.s_[-1:], 1) # kill the rightmost empty col.
            
            if raw_tag_filename:
                binned_data = numpy.delete(binned_data, numpy.s_[-1:], 1)
                oh = open(raw_tag_filename, "w")
                writer = csv.writer(oh, dialect=csv.excel_tab)
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

        ret = {"pileup": pileup}
        if heatmap_filename or raw_tag_filename: #if binned_data:
            # __nonzero__ not set in numpy arrays, so assume binned_data is valid
            # if doing heatmap
            ret["binned_data"] = binned_data

        return(ret)

    def heatmap(self, filename=None, genelist=None, distance=1000, read_extend=200, log=2, 
        bins=200, sort_by_intensity=True, raw_heatmap_filename=None, bracket=None, 
        pointify=True, respect_strand=False, cmap=cm.BuPu, norm_by_read_count=False,
        imshow=True,
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
                
            respect_strand (Optional, default=False)
                Use the strand in the genelist to orientate each genomic location (for example when doing TSS 
                or TTS).
                
            log (Optional, default=2)
                log transform the data, optional, the default is to transform by log base 2.
                Note that this parameter only supports "e", 2, and 10 for bases for log
                if set to None no log transform takes place.
            
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
        
        # get a sorted list of all the locs I am going to use.
        gl_sorted = genelist.deepcopy()
        gl_sorted.sort('loc')
        all_locs = gl_sorted['loc']   
           
        if respect_strand:
            strands = gl_sorted['strand']
        else:
            strands = ['+'] * len(all_locs) # Fake strand labels for code simplicity
        
        curr_cached_chrom = None
        cached_chrom = None
        bin_size = None

        number_of_tags_in_library = 1.0 # For code laziness
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
            
            # I can dispose and free memory as the locations are now sorted:
            # See if the read_extend is already in the cache:
            if l['chr'] != curr_cached_chrom:
                curr_cached_chrom = l['chr']
                cached_chrom = self.get_array_chromosome(l['chr'], read_extend) # Will hit the DB if not already in cache
        
            actual_width = cached_chrom.shape[0] - l['left']              

            if actual_width < expected_width:
                #print "!", a.shape
                continue
            else:
                a = cached_chrom[l['left']:l['right']]
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

            if respect_strand:
                # positive strand is always correct, so I leave as is.
                # For the reverse strand all I have to do is flip the array.
                if read[1] in negative_strand_labels:
                    a = a[::-1]
            
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
        
        if sort_by_intensity:
            mag_tab = sorted(mag_tab, key=itemgetter("sum"))
        
        if log:
            data = numpy.array(table)+1 # Should be in tag count, so a reasonable pad
        else:
            data = numpy.array(table)
        
        newt = []
        for item in mag_tab:
            newt.append(data[item["n"],])
        data = numpy.array(newt)
        #print data.shape
        data = numpy.delete(data, numpy.s_[-1:], 1) # Nerf the last column.
        
        # Get the sorted_locs for reporting purposes:
        temp_gl_copy = genelist.deepcopy().linearData # deepcopy for fastness
        sorted_locs = [temp_gl_copy[i['n']] for i in mag_tab] # need to sort them by intensity, if done that way
        sorted_locs.reverse() # Make it into the intuitive order.
        gl = Genelist.genelist()
        gl.load_list(sorted_locs)
        sorted_locs = gl
        
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

    def measure_frip(self, genelist=None, sample=None, delta=None, pointify=False):
        """
        **Purpose**
            Measure the FRiP 'Fraction of Reads in Peaks' as defined by Landt et al., 2012; 
            Gen Res, 22:1813-1831.
            
            Essentially the fraction of reads inside a list of peaks (from the genelist).
            
        **Arguments**
            genelist (Required)
                A list of peaks, must have a "loc" key containing the location data.
                Ideally this is the peak spans reported by your peak discovery tool,
                but you can provide delta=xxx bp argument to expand the locations.
                
            sample (Optional)
                sample only the first n peaks. By default, all peaks are used.
                
            delta (Optional, default=None)
                a value to expand the locations by + and - delta.
                
            pointify (Optional, default=False)
                'pointify' (Convert a span of base pairs to the centre point).
                Executed before 'delta'
            
        **Returns**
            The FRiP score, the total number of reads and the number of reads inside the peaks
        """
        assert genelist, "You must provide a genelist"
        assert "loc" in genelist.linearData[0], "The genelist does not appear to have a 'loc' key"
        
        if sample:
            gl = gl[sample]
        else:
            gl = genelist.deepcopy() # get a copy as I may mess with it.
        
        if pointify:
            gl = gl.pointify("loc")
            
        if delta:
            gl = gl.expand("loc", delta)
        
        # work out the total number of reads in this library
        chr_names = self.get_chromosome_names()
        total_reads = 0
        for chrom in chr_names:
            total_reads += self.get_numreads_on_chromosome(chrom)
                    
        num_reads = 0
        p = progressbar(len(gl))
        for idx, item in enumerate(gl):
            num_reads += self.get_read_count(item["loc"])
            p.update(idx)
        
        return({"FRiP": num_reads/float(total_reads), "Total number of reads": total_reads, "Number of reads in peaks": num_reads})

    def qc_encode_idr(self, chrom_sizes=None, filename=None, max_shift=400, **kargs):
        """
        **Purpose**
            Perform QC for ChIP-seq libraries, as explained 
            
            https://sites.google.com/site/anshulkundaje/projects/idr
            
            and in Landt et al., 2012, Gen Res, 22:1813-1831.
        
        **Arguments**
            chromosome_sizes (Required)
                You must provide a dict, containing the chromosome names and the 
                number of base pairs.
                
                For mm9 this data is available as part of glbase and can be specified:
                
                trk.qc_encode_idr(chromosome_sizes=gldata.chromsizes["mm9"])
                
                Potentially, hg18, hg19, mm8 and mm9 will be available too. maybe.
                
            filename (Required)
                filename to save the plot to
                
        **Returns**
            NSC and RSC values. (See Landt et al., 2012; Gen Res, 22:1813-1831) for 
            details.
        """ 
        assert chrom_sizes, "You must provide chromosome sizes"
        assert filename, "You must provide a filename"

        if not self._draw:
            self._draw = draw()
        
        # I only need to generate the + strand once.
        plus_strands = {}
        minu_strands = {}
        
        # constructing a numpy array is excessively large. I only need to store pairs of reads
        
        all_chroms = set(self.get_chromosome_names()) & set([i.replace("chr", "") for i in list(chrom_sizes.keys())]) # only iterate ones in both list
        all_p = numpy.array([])
        all_m = numpy.array([])
        res = []
        pears = numpy.zeros([max_shift, len(all_chroms)])
        
        for idx, chrom in enumerate(all_chroms):
            this_p = numpy.array([r[0] for r in self.get_all_reads_on_chromosome(chrom, "+")])
            this_m = numpy.array([r[1] for r in self.get_all_reads_on_chromosome(chrom, "-")])
                        
            p = progressbar(max_shift)
            for n in range(max_shift):                   
                this_m = this_m - 1    
                                 
                union = numpy.union1d(this_p, this_m) # only ones I will have to look at
                #union = numpy.intersect1d(this_p, this_m) 
                #xor_union = numpy.setxor1d(this_p, this_m)
                #union = numpy.append(union, xor_union)
    
                pair_p = numpy.bincount(this_p, minlength=max(this_p.max(), this_m.max())+1)[union]
                pair_m = numpy.bincount(this_m, minlength=max(this_p.max(), this_m.max())+1)[union]
                
                pears[n, idx] = pearsonr(pair_p, pair_m)[0]
                p.update(n)
                """
                fig = self._draw.getfigure()
                ax = fig.add_subplot(111)
                ax.plot(pair_p, pair_m, 'o', mec="none", alpha=0.2)
                ax.set_title("Pearson: %.3f" % pears[n, idx])
                fig.savefig("plots/test_%s_%s.png"% (chrom, n))
                """
            print("Done chromosome '%s'" % chrom)

        print(pears)            
        for row in pears:
            res.append(numpy.average(row))
        print(res)
        
        fig = self._draw.getfigure(**kargs)
        ax = fig.add_subplot(111)
        ax.plot(numpy.arange(len(res)), res)
        self._draw.savefigure(fig, filename)            
        
        ret = {"NSC": 0.0, "RSC": 0.0}
        return(ret)

if __name__ == "__main__":
    """

    Current 15 s

    """
    import random, time
    from .location import location
    from .genelist import genelist
    
    s = time.time()
    print("Building...")
    t = track(filename="testold.trk2", name="test", new=True)
    for n in range(0, 10000):
        l = random.randint(0, 100000)
        t.add_location(location(chr="1", left=l, right=l+35), strand="+")
    t.finalise()
    e = time.time()
    print(e-s, "s")
    #t.finalise()
    
    print(t.get_reads('chr1:100-200'))
    
    s = time.time()
    print("Fake bed...")
    # fake a bed
    newb = []
    for n in range(0, 1000):
        l = random.randint(1000, 100000) # 1000 is for the window size. -ve locs are real bad.
        newb.append({"loc": location(chr="1", left=l, right=l+200), "strand": "+"})
    bed = genelist()
    bed.load_list(newb)
    e = time.time()
    print(e-s, "s")
    
    t = track(filename="testold.trk2")
    print("Pileup...")
    import cProfile, pstats
    cProfile.run("t.pileup(genelist=bed, filename='test.png', bin_size=10, window_size=1000)", "profile.pro")
    p = pstats.Stats("profile.pro")
    p.strip_dirs().sort_stats("time").print_stats()

    print(t.pileup(genelist=bed, filename='/tmp/test2.png', respect_strand=True))
    print(t.pileup(genelist=bed, filename='/tmp/test2.png', pointify=False, respect_strand=True))
    
    print(bed.all())
    print(t.pileup(genelist=bed, filename='/tmp/test2.png', pointify=False, window_size=0, respect_strand=True))

    print(t.heatmap(genelist=bed, raw_heatmap_filename="/tmp/test.tsv", filename='/tmp/test.png', bin_size=10, window_size=1000))
    print(t.heatmap(genelist=bed, raw_heatmap_filename="/tmp/test.tsv", filename='/tmp/test.png', bin_size=10, window_size=1000, log=None))
    print(t.heatmap(genelist=bed, raw_heatmap_filename="/tmp/test.tsv", filename='/tmp/test.png', bin_size=10, window_size=1000, log=None, respect_strand=True))
    

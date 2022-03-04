"""

An SQL replacement for genome

This loses a lot of the flexibility of genelists, but could potentially be abstracted into
a genelist-like interface.

TODO:
-----
. This should be abstracted down to a base_sql class which mimics the base functionality of base_genelist
. This can be specified as putting the sql into memory, so not requireing an in disk store.

"""



import pickle, sys, os, struct, math, sqlite3, zlib, time, csv

from operator import itemgetter

from .progress import progressbar
from .errors import AssertionError
from .location import location
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

class genome_sql(base_track):
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

    """
    def __init__(self, name=None, new=False, filename=None, norm_factor=1.0, **kargs):
        base_track.__init__(self, name, new, filename, norm_factor)

        if new:
            self.__setup_tables(filename)

    def __repr__(self):
        return("glbase.genome_sql")

    def __setup_tables(self, filename):
        """
        No pre-defined file - I want a new database track.
        """
        # If i make it here, then a lot of grunt work already done in base_track
        c = self._connection.cursor()

        c.execute("CREATE TABLE main (chromosome TEXT PRIMARY KEY, num_features INT)")

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

        # make the new chromsome table:
        table_name = "chr_%s" % str(chromosome)
        c.execute("CREATE TABLE %s (transcript_left INT, transcript_right INT, cds_left INT, cds_right INT, exonStarts TEXT, exonEnds TEXT, name TEXT, strand TEXT, feature_type TEXT)" % (table_name, ))

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

        return bool(result)

    def add_feature(self, loc, cds_loc, exonCounts, exonStarts, exonEnds, name, strand, feature_type='gene'):
        """
        **Purpose**
            Add a location to the track.
            Increments the score by 'increment' from loc["left"] to
            loc["right"]

        **Arguments**
            loc

            strand

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

        # get the old number of seq_reads
        c.execute("SELECT num_features FROM main WHERE chromosome=?", (loc["chr"], ))
        current_seq_reads = c.fetchone()[0] # always returns a tuple

        c.execute("UPDATE main SET num_features=? WHERE chromosome=?", (current_seq_reads+1, loc["chr"]))

        # add the location to the seq table:
        insert_data = (loc['left'], loc['right'], cds_loc['left'], cds_loc['right'], str(exonStarts), str(exonEnds), str(name), str(strand), feature_type)
        c.execute("INSERT INTO %s VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)" % table_name, insert_data)

        c.close()

    def find(self, name=None, greedy=True):
        """
        **Purpose**
            Find an element in the 'name' key that matches the provided name or return None

        **Arguments**
            name (Required)
                a name (e.g. ENST name, gene symbol, etc. )

        **Results**
            A lits of matches or None
        """
        assert name, 'You must provide a name'

        all_results = []
        for chrom in self.get_chromosome_names():
            table_name = "chr_%s" % chrom
            result = self._connection.execute("SELECT * FROM %s WHERE (name == ?)" % table_name, (name,))
            result = result.fetchall()
            if result:
                all_results += (self.__format_results(result, chrom))

        if not all_results:
            all_results = None
        return all_results

    def __format_results(self, results, chrom):
        '''
        Internal function to pack the results into a nice dict
        '''
        newres = []
        for r in results:
            newr = {'loc': location(chr=chrom, left=r[0], right=r[1]),
                'cds_loc': location(chr=chrom, left=r[2], right=r[3]),
                'exonStarts': eval(r[4]),
                'exonEnds': eval(r[5]),
                'name': r[6], 'type': r[8], 'strand': r[7]}
            newres.append(newr)
        return newres

    def getFeatures(self, loc=None, **kargs):
        """
        **Purpose**
            get all of the genomic features (probably genes) within a
            certain location span.

        **Arguments**
            loc
                either a location() or a cooercable location().

        **Returns**
            A vanilla list containing a bunch of keys describing any features
            sitting in the genomic span specified by location.
        """
        assert loc, "no location provided"

        try:
            loc = location(loc=loc)
        except (TypeError, IndexError):
            raise AssertionError("cannot cooerce location into correct form. Location is mangled?")

        table_name = "chr_%s" % loc["chr"]

        result = self._connection.execute("SELECT * FROM %s WHERE (transcript_right >= ? AND transcript_left <= ?)" % table_name,
            (loc["left"], loc["right"]))

        result = result.fetchall() # safer for empty lists and reusing the cursor

        if result:
            result = self.__format_results(result, loc['chr'])
        if not result: # Compatability with chipFish
            result = []

        return(result)

    def get_chromosome_names(self):
        """
        **Purpose**
            Return a list of all the valid chromosome names in the database

        **Arguments**
            None

        **Returns**
            A list of strings containing all of the chromosome names in the genome_sql
        """
        if not self._c:
            self._c = self._connection.cursor()

        self._c.execute("SELECT chromosome FROM main")
        r = [i[0] for i in self._c.fetchall()]
        return(set(r))

    def get_feature_count(self):
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

        self._c.execute("SELECT chromosome, num_features FROM main")
        r = [int(i[1]) for i in self._c.fetchall()]

        return(sum(r))

    def __iter__(self):
        """
        (Override)
        make the geneList behave like a normal iterator (list)
        """
        all_chrom_names = self.get_chromosome_names()

        for c in all_chrom_names:
            table_name = "chr_%s" % c

            result = self._connection.execute("SELECT * FROM %s" % table_name)

            r = True # Survive first while

            while r:
                r = result.fetchone() # safer for empty lists and reusing the cursor

                if r:
                    # This needs to be abstracted away
                    # Repack item into a nice format:
                    # (57049987, 57050281, 57049987, 57050281, '[1]', '[1]', 'SINE-AluJb', '-', 'SINE')
                    r = {'loc': location(chr=c, left=r[0], right=r[1]),
                        'cds_loc': location(chr=c, left=r[2], right=r[3]),
                        'exonStarts': eval(r[4]),
                        'exonEnds': eval(r[4]),
                        'name': r[6], 'type': r[8], 'strand': r[7]}
                    yield r

    def __len__(self):
        return(self.get_feature_count())

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

    def bindSequence(self, path=None):
        """
        **Purpose**

        Bind genome fasta files so that this genome object will recognise
        the sequence. This step is required if you want to use fastalists
        and genome.getSequence()

        **Arguments**

        path
            path specifying the locations of the FASTA files that make
            up the sequence data. They usually come in the form "chr1.fa"
            for human and mouse genomes.

            bindSequence will only work with multi-fasta files, i.e. the fasta genome should
            be in the form:

            FILE: chr1.fa::

                >chr1
                NNNNNNNNNNNNNNNNNNNNNNNN

            FILE: chr2.fa::

                >chr2
                NNNNNNNNNNNNNNNNNNNNNNNN

            etc.

            The names of the chromosomes will come from the names of the fasta files before the
            period and with 'chr' removed (if present), so for example:

            chr1.fa

            will result in '1' entries in the db.

            And similarly '1.fa' will result in chr names of '1'.

        **Result**

        returns True if complete.
        genome.getSequence(loc="chrN:left-right") will now correctly return
        the sequence specified by the location.
        """
        raise NotImplementedError

    def getSequence(self, loc=None, **kargs):
        """
        **Purpose**

        get the sequence under coords...

        **Arguments**

        loc or coords (one is Required)
            genomic coordinates of the form "chr1:100100-100200"

        strand
            Use, +, f, top, forward, 0, for the top strand
            Use, -, r, bottom, reverse, 1, for the reverse complement strand
            If the strand is not specified then the + strand will be returned.

        mask (Optional, default=False)
            'repeat mask' the returned sequence (i.e. convert lower-case
            acgt to NNNN). DOES NOT perform repeat masking. Only converts lower-cased
            bases to NNN.

        **Result**

        returns a string containing the sequence at 'coords'
        """
        raise NotImplementedError

    def getSequences(self, genelist=None, loc_key='loc', replace_loc_key=True, strand_key=False,
        mask=False, pointify=False, delta=False, **kargs):
        """
        **Purpose**

        Get all of the sequences from a gene_list-like object (e.g. a peaklist,
        microarray, etc) with some sort of valid location key (e.g. "chr1:10000-20000")

        I've checked this extensively - if you provide the correct location
        it will send back the correct sequence. Any further errors are
        generally from wrong locations. So Heh.

        **Arguments**

        genelist (Required)
            some sort of genelist-like object

        loc_key (Required)
            the name of the location key (e.g. "loc", "tss_loc", "coords", etc..)

        Optional Arguments

        replace_loc_key (Optional, default=True)
            If you specify a delta then the sequence will cover a different region than that
            specified in the loc key. hence getSequences() replaces the "loc" key with
            the new loc for which the sequence was taken from. If you set this to False
            then the old loc key is NOT overwritten.

        strand_key
            If you want the list to respect the strand you must tell it the name of
            a 'strand' key.

        deltaleft=n
            expand the coordinates left by n base pairs
            Will respect the orientation of the strand.
            (You must specify a 'strand_key')

        deltaright=n
            expand the coordinates rightwards by n base pairs.
            Will respect the orientation
            of the strand.
            (You must specify a 'strand_key')

        delta=n
            expand the coordinates by n, added onto the left and right

        pointify (True|False, default=False)
            turn the location into a single base pair based on the centre
            of the coordinates (best used in combination with delta to
            expand reads symmetrically, pointify will be performed before
            the expansion of the coordinates)

        mask (default=False)
            use the upper and lower case of the fasta files to 'mask'
            the sequence. This will turn acgt to NNNN.

            This is not a proper repeat masker and relies on your genome
            being repeat masked. For human and mouse this is usually true,
            but for other genomes has not been tested.

        **Result**

        returns a copy of the original genelist-like object
        with the new key "seq" containing the sequence.
        Will add a new key "seq_loc" that contains the new location that
        the seq spans across.
        """
        raise NotImplementedError

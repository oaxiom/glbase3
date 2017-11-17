#! /usr/bin/python
"""

Part of glbase,
converts a sequence file to a track (graph)-like display

"""

import sys, os, csv, time, cProfile, pstats

from .. import track
from .. import delayedlist # get delayedlist from glbase
from .. import genelist
from .. import format
from .. import config

def seqToTrk(infilename=None, outfilename=None, name=None, stranded=True, format=format.minimal_bed, 
    gzip=False, **kargs):
    """
    **Purpose**
        Convert a list of genomic coordinates to a trk database (Actually an SQL
        db)

    **Arguments**
        infilename
            the name of the filename to read from. This should be a csv.
            The data is loaded using a delayedlist(), so memory
            is not a big problem.
            
            You can also send a list of filenames, 

        outfilename
            The filename of the trk file to write to.

        name (Required)
            A name describing the track.

        format (Optional, default=format.minimal_bed)
            a format specifier in the glbase format describing the csv
            file.

        stranded (Optional, default=True)
            expects a "strand" key and will store the strand.

        norm_factor (Optional, default = 1.0)
            An optional normalization factor. Data is multiplied by this number before display

        gzip (Optional, default=False)
            The input file(s) is a gzip file.

    **Returns**
        True on completion
        and a trk file in outfilename

    """
    assert infilename, "seqToTrk(): You must specify infilename"
    assert outfilename, "seqToTrk(): You must specify outfilename"
    assert name, "seqToTrk(): You must specify a name for the track"
    
    __strand_checked = False
    
    if not isinstance(infilename, list):
        infilename = [infilename]

    n = 0
    m = 0
    total = 0
        
    t = track(filename=outfilename, stranded=stranded, new=True, name=name, **kargs)

    s = time.time()
    for file in infilename:
        #config.log.info("Started %s -> %s" % (file, outfilename)) # Don't need as the delayedlist binding will put some output
        seqfile = delayedlist(filename=os.path.realpath(file), format=format, gzip=gzip)
        for item in seqfile:
            if stranded: 
                if not __strand_checked: # purposely throw an error
                    assert 'strand' in item, 'stranded=True, but no "strand" key found'
                    __strand_checked = True
                    
                t.add_location(item["loc"], strand=item["strand"])
            else:
                t.add_location(item["loc"])

            n += 1
            total += 1
            if n > 1000000:
                m += 1
                n = 0
                config.log.info("%s,000,000 tags read" % m) # How to correct for 1000 M tags?
                #break
     
    # 1000 = averages 8-9 s
    # 3.65 s cache.
    # 0.61 s better cacheing, less commits
    # 10,000 used to check now.
    # 5.98 s, 22.8Mb (Binary(array.tostring())
    # 6.0 s 17.1 Mb (str/eval)
    # 6.0 s 259kb (zlib/Binary()/tostring())
    # track_new
    # 2.5 s # overhead is almost all delayedlist now...

    config.log.info("Library contains '%s' tags" % total)
    config.log.info('Building cache, finalizing...')
    t.finalise()
    e = time.time()
    config.log.info("Took: %.1f seconds" % (e-s))
    return(True)

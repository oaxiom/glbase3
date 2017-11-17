#! /usr/bin/python
"""

Part of glbase,
converts a sequence file to a track (graph)-like display



"""

import sys, os, csv, time
from glob import glob

from .. import flat_track
from .. import config
from .. import location

def gerp_to_flat(path, outfilenameA, outfilenameB, name, bin_format="f", **kargs):
    """
    **Purpose**
        Convert a list of genomic coordinates to a trk database (Actually an SQL
        db)

    **Arguments**
        path
            gerp files come in the from chr1.maf.rates
            This should point to the path

        outfilenameA (Required)
        outfilenameB (Required)
            GERP data comes in two columns, the 'neutral rate' and the 'RS score'. 
            flats can only store a single value for each bp, so I have to generate two 
            flats. A will store the neutral value and B will store the RS score.
            
            This does lead to limitations later as the flats must be treated as separate tracks.

            Don't forget though that the data returned by pileup() is a numpy array.
            So you could do neat things like::
            
                rs_score = f.pileup(...)
                neutral = f.pileup(...)
                
                average = rs_score - neutral
                
                then use maplotlib::
                pyplot.plot(average)
                pyplot.savefig(...)
            
            or:: 
            
                pyplot.plot(rs_score, label=average)
                pyplot.plot(neutral, label="neutral")
                pylot.legend()
                pyplot.savefig(...)

        name
            A name describing the track

    **Returns**
        True on completion
        and a trk file in outfilename

    """
    assert os.path.realpath(path), "no filename specified"
    assert os.path.realpath(outfilenameA), "no save filename specified"
    assert os.path.realpath(outfilenameB), "no save filename specified"

    n = 0
    m = 0
    total = 0

    #fa = flat_track(filename=outfilenameA, new=True, name=name, bin_format=bin_format)
    fb = flat_track(filename=outfilenameB, new=True, name=name, bin_format=bin_format)

    config.log.info("Started %s -> %s and %s" % (path, outfilenameA, outfilenameB))
    s = time.time()
    step = 1
    
    for filename in glob(os.path.join(path, "*.maf.rates")):
        config.log.info("Doing %s" % filename)
        oh = open(filename, "rU")
        chrom = filename.split(".")[0].replace("chr", "")
        cleft = 0

        for line in oh:
            if line: # Just in case there are some empty lines
                d = line.split()
                #print d[0], cleft
                #fa.add_score(chromosome=chrom, left=cleft, right=cleft+step, score=float(d[0]))
                fb.add_score(chromosome=chrom, left=cleft, right=cleft+step, score=float(d[1]))
                cleft += step 
                
                if n>1e6:
                    m += 1
                    print "%s,000,000 bp" % m
                    n = 0
                n += step
        oh.close()

    e = time.time()

    config.log.info("Finalise library...")
    #fa.finalise()
    fb.finalise()
    config.log.info("Took: %s seconds" % (e-s))
    return(True)
    

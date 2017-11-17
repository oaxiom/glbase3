#! /usr/bin/python
"""

Part of glbase,
converts a sequence file to a track (graph)-like display



"""

import sys, os, csv, time

from .. import flat_track
from .. import config
from .. import location

def wigstep_to_flat(infilename, outfilename, name, bin_format=None, **kargs):
    """
    **Purpose**
        Convert a list of genomic coordinates to a flat file (Actually an SQL
        db)

    **Arguments**
        infilename
            the name of the wiggle step filename to read from. 

        outfilename
            The filename of the flat file to write to.

        name
            A name describing the track
            
        bin_format
            the format to use to store the data. Valid values are:
            
            i = integers
            f = floats

    **Returns**
        True on completion
        and a flat file in outfilename

    """
    assert os.path.realpath(infilename), "no filename specified"
    assert os.path.realpath(outfilename), "no save filename specified"

    n = 0
    m = 0
    total = 0

    f = flat_track(filename=outfilename, new=True, name=name, bin_format=bin_format)

    config.log.info("Started %s -> %s" % (infilename, outfilename))

    s = time.time()
    oh = open(infilename, "rU")
    cleft = 0
    for line in oh:
        if line: # Just in case there are some empty lines
            if "fixedStep" in line:
                t = line.split()
                chrom = t[1].split("=")[1]
                cleft = int(t[2].split("=")[1])
                step = int(t[3].split("=")[1])
            else:
                #f.add_score(loc=location("%s:%s-%s" % (chrom, cleft, cleft+step)), score=float(line))
                f.add_score(chromosome=chrom.replace("chr", ""), left=cleft, right=cleft+step, score=float(line))
                cleft += step    
            
            if n>1e6:
                m += 1
                print("%s,000,000 bp" % m)
                n = 0
            n += step

    e = time.time()

    config.log.info("Finalise library. Contains '%.1e' bps of data" % (int(m*1e6)))
    f.finalise()
    config.log.info("Took: %s seconds" % (e-s))
    return(True)
    

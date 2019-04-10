#! /usr/bin/python
"""

Part of glbase,
converts a sequence file to a track (graph)-like display

This module is incorrectly named.

"""

import sys, os, csv, time

import gzip as opengzip
from .. import flat_track
from .. import config
from .. import location

def bedgraph_to_flat(infilename, outfilename, name, gzip=None, all_in_mem=False, **kargs):
    """
    **Purpose**
        Convert a bedGraph file to a flat file (Actually an SQL
        db)

        A bedGraph file should look like this:

        chr1	3000881	3000986	1
        chr1	3000986	3001081	2
        chr1	3001081	3001157	1
        chr1	3001157	3001186	2
        chr1	3001186	3001357	1
        chr1	3001535	3001735	1
        chr1	3002009	3002209	1
        chr1	3003894	3003917	1
        chr1	3003917	3004094	2


    **Arguments**
        infilename
            the name of the bedGraph file to read from.

        outfilename
            The filename of the flat file to write to.

        name
            A name describing the track

        gzip (Optional, default=False)
            The input file(s) is a gzip file.

        all_in_mem (Optional, default=False)
            If you have a lot of memory, set this to True and the flat creation is done
            all in memory and only committed to disk right at the end.

    **Returns**
        True on completion
        and a flat file in outfilename

    """
    assert os.path.realpath(infilename), "no filename specified"
    assert os.path.realpath(outfilename), "no save filename specified"

    bin_format = 'f' # Only one supported
    n = 0
    m = 0
    total = 0

    f = flat_track(filename=outfilename, new=True, name=name, bin_format=bin_format)

    config.log.info("Started %s -> %s" % (infilename, outfilename))

    s = time.time()
    if not gzip:
        oh = open(infilename, "rt")
    else:
        oh = opengzip.open(infilename, 'rt')

    cleft = 0
    for line in oh:
        #print(line)
        if not "#" in line:
            line = line.split()
            f.add_score(chromosome=line[0].replace("chr", ""),
                left=int(line[1]),
                right=int(line[2]),
                score=float(line[3].strip()),
                all_in_mem=all_in_mem) # Do it all in memory. I hope you have enough!

            if n % 1000000 == 0:
                print('Processed: {:,} bp'.format(n))
                #break
            n += int(line[2]) - int(line[1])

    e = time.time()

    config.log.info("Finalise library. Contains '%.1e' bps of data" % (int(m*1e6)))
    f.finalise()
    config.log.info("Took: %s seconds" % (e-s))
    return(True)


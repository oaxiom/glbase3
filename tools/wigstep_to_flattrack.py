#! /usr/bin/python
"""

Part of glbase,
converts a sequence file to a track (graph)-like display



"""

import sys, os, csv, time, numpy
import gzip as opengzip
from .. import flat_track
from .. import config
from .. import location

def wigstep_to_flat(infilename, outfilename, name, bin_format=None, gzip=False, skip_non_standard_chroms=False, **kargs):
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

            f = floats

        skip_non_standard_chroms (Optional, default=False)
            only use canonical chromsome names (i.e. chr1, chrX, chrM)
            and not scaffolds and unplaced (chr14_KIA... etc).

    **Returns**
        True on completion
        and a flat file in outfilename

    """
    assert os.path.realpath(infilename), "no filename specified"
    assert os.path.realpath(outfilename), "no save filename specified"

    if not gzip:
        open_mode = open
    else:
        open_mode = opengzip.open

    n = 0
    m = 0
    total = 0

    f = flat_track(filename=outfilename, new=True, name=name, bin_format='f')

    config.log.info("Started {} -> {}".format(infilename, outfilename))

    s = time.time()
    oh = open_mode(infilename, "rt")
    cleft = 0
    newchrom = None
    lastchrom = None

    for line in oh:
        if 'track' in line: # Ignore headers'
            continue

        if line: # Just in case there are some empty lines
            if "fixedStep" in line:
                t = line.split()
                chrom = t[1].split("=")[1]

                # if change of chrom, save it to the flat;
                if lastchrom and lastchrom != chrom and newchrom:
                    if skip_non_standard_chroms and '_' in lastchrom:
                        continue
                    else:
                        f.add_chromosome_array(lastchrom, numpy.array(newchrom))
                        config.log.info('Finished {} with {:,} bp of data'.format(lastchrom, len(newchrom)))
                    newchrom = []

                lastchrom = chrom
                cleft = int(t[2].split("=")[1])
                step = int(t[3].split("=")[1])

                if not newchrom:# setup new array;
                    newchrom = [0] * cleft
                else: # pad to the new location cleft
                    newchrom += [0] * (cleft-len(newchrom))
            else:
                #f.add_score(loc=location("%s:%s-%s" % (chrom, cleft, cleft+step)), score=float(line))
                newchrom.append(float(line))
                cleft += step

            if n > 1e7: # because n+= step;
                m += 1
                print("{:,} bp".format(int(m*1e7)))
                n = 0
            n += step

    e = time.time()

    config.log.info("Finalise library. Contains approx. {:,} bps of data".format(int(m*1e6)))
    f.finalise()
    config.log.info("Took: %s seconds" % (e-s))
    return True

